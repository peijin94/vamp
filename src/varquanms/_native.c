#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <math.h>
#include <fenv.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <zlib.h>
#include <time.h>

static double thread_cpu(void) {struct timespec t;clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t);return t.tv_sec+t.tv_nsec*1e-9;}

static void put32(uint8_t *p,uint32_t value) {for(int k=0;k<4;k++) p[k]=(uint8_t)(value>>(8*k));}
static void header(uint8_t *p,int bits,int compressed,uint32_t n,uint32_t exceptions,float step,int mapping) {
    p[0]=bits;p[1]=compressed;p[2]=mapping;p[3]=0;put32(p+4,n);put32(p+8,exceptions);memcpy(p+12,&step,4);
}

/* Evaluate the complete block search inside one GIL-free region. A generous
 * bound margin prevents reduction-order false rejections; Python checks the
 * winner strictly and retries its reference search on a near-boundary miss. */
static PyObject *search(PyObject *self,PyObject *args) {
    Py_buffer a={0},f={0},widths={0},scales={0},bounds={0};int channels,level,mapping=0;
    if(!PyArg_ParseTuple(args,"y*y*y*y*y*ii|i",&a,&f,&widths,&scales,&bounds,&channels,&level,&mapping)) return NULL;
    Py_ssize_t n=a.len/8, candidates=widths.len;
    int valid=a.len%8==0 && n>0 && n<=262144 && f.len==n && channels>0 && n%channels==0 &&
        candidates>0 && candidates<=4096 && scales.len==4*candidates && bounds.len==32 && level>=-1 && level<=9 && (mapping==0||mapping==1) && fegetround()==FE_TONEAREST;
    const uint8_t *width=(const uint8_t *)widths.buf;
    for(Py_ssize_t k=0;k<candidates && valid;k++) valid=width[k]==2||width[k]==4||width[k]==6||width[k]==8||width[k]==12||width[k]==16;
    if(!valid) {PyErr_SetString(PyExc_ValueError,"Invalid native search buffers/parameters");goto cleanup_search;}
    uint8_t *raw=malloc(16*n), *compressed=malloc(compressBound(16*n)), *best=malloc(8*n+16);
    float *decoded=malloc(8*n),*best_decoded=malloc(8*n);
    if(!raw||!compressed||!best||!decoded||!best_decoded) {
        free(raw);free(compressed);free(best);free(decoded);free(best_decoded);PyErr_NoMemory();goto cleanup_search;
    }
    const float *x=(const float *)a.buf,*steps=(const float *)scales.buf;
    const double *bound=(const double *)bounds.buf;
    const uint8_t *flags=(const uint8_t *)f.buf;
    size_t best_size=8*n+16; int failed=0;double compression_cpu=0,quantization_cpu=0,cpu_stamp;
    Py_BEGIN_ALLOW_THREADS
    header(best,0,0,(uint32_t)n,0,0,0);memcpy(best+16,a.buf,8*n);memcpy(best_decoded,a.buf,8*n);
    if(level>=0) {
        cpu_stamp=thread_cpu();
        uLongf size=compressBound(8*n);
        if(compress2(compressed,&size,(const Bytef *)a.buf,8*n,level)!=Z_OK) failed=1;
        else if(size<(size_t)8*n) {header(best,0,1,(uint32_t)n,0,0,0);memcpy(best+16,compressed,size);best_size=size+16;}
        compression_cpu+=thread_cpu()-cpu_stamp;
    }
    for(Py_ssize_t k=0;k<candidates && !failed;k++) {
        int bits=width[k],maximum=(1<<(bits-1))-1;float step=steps[k];
        double center=maximum+(mapping?0.5:0.0),shift=mapping?center:0;
        int lower=mapping?0:-maximum,upper=mapping?(1<<bits)-1:maximum;
        if(!isfinite(step)||step<0) continue;
        cpu_stamp=thread_cpu();
        size_t packed=(2*n*bits+7)/8,ne=0,byte=0;
        double limit=(double)step*center,sum=0,sum_r=0,sum_i=0,row_sum=0,max_square=0,max_row=0;
        Py_ssize_t count=0,row_count=0;uint64_t accumulator=0;int available=0,overflow=0;
        /* Exception positions are staged after maximum code capacity; values
         * are copied into the final body after the exception count is known. */
        uint32_t *positions=(uint32_t *)(raw+4*n);
        for(Py_ssize_t j=0;j<n;j++) {
            float r=x[2*j],im=x[2*j+1];int finite=isfinite(r)&&isfinite(im);
            int exact=flags[j]||!finite||(r==0&&signbit(r))||(im==0&&signbit(im))||fabs((double)r)>limit||fabs((double)im)>limit;
            int qr=mapping?maximum+1:0,qi=qr;
            if(!exact&&step!=0) {qr=(int)nearbyint((double)r/step+shift);qi=(int)nearbyint((double)im/step+shift);}
            if(qr<lower||qr>upper||qi<lower||qi>upper) {overflow=1;break;}
            float cr=(float)qr-(float)shift,ci=(float)qi-(float)shift;
            decoded[2*j]=cr*step;decoded[2*j+1]=ci*step;
            if(!isfinite(decoded[2*j])||!isfinite(decoded[2*j+1])) {overflow=1;break;}
            uint32_t pair[2]={(uint32_t)(qr+(mapping?0:maximum)),(uint32_t)(qi+(mapping?0:maximum))};
            for(int h=0;h<2;h++) {accumulator|=(uint64_t)pair[h]<<available;available+=bits;
                while(available>=8) {raw[byte++]=(uint8_t)accumulator;accumulator>>=8;available-=8;}}
            if(exact) {positions[ne++]=(uint32_t)j;memcpy(decoded+2*j,x+2*j,8);}
            if(!flags[j]&&finite) {
                double er=(double)decoded[2*j]-r,ei=(double)decoded[2*j+1]-im,square=er*er+ei*ei;
                sum+=square;sum_r+=er;sum_i+=ei;row_sum+=square;count++;row_count++;
                if(square>max_square) max_square=square;
            }
            if((j+1)%channels==0) {double mse=row_count?row_sum/(2*row_count):0;if(mse>max_row) max_row=mse;row_sum=0;row_count=0;}
        }
        if(overflow) {quantization_cpu+=thread_cpu()-cpu_stamp;continue;}
        if(available) raw[byte++]=(uint8_t)accumulator;
        double max_error=sqrt(max_square),mse=count?sum/(2*count):0,bias=count?hypot(sum_r/count,sum_i/count):0;
        if(mse>bound[0]*(1+1e-9)||max_error>bound[1]*(1+1e-9)||max_row>bound[3]*(1+1e-9)||bias>bound[2]+max_error*1e-9) {quantization_cpu+=thread_cpu()-cpu_stamp;continue;}
        /* Move the position list before copying original values, so neither
         * operation overwrites the staged indexes before they are consumed. */
        memmove(raw+packed,positions,4*ne);
        for(size_t j=0;j<ne;j++) {uint32_t index;memcpy(&index,raw+packed+4*j,4);memcpy(raw+packed+4*ne+8*j,x+2*index,8);}
        size_t body_size=packed+12*ne;uint8_t *body=raw;int compressed_id=0;
        quantization_cpu+=thread_cpu()-cpu_stamp;
        cpu_stamp=thread_cpu();
        if(level>=0) {uLongf size=compressBound(body_size);
            if(compress2(compressed,&size,raw,body_size,level)!=Z_OK) {failed=1;break;}
            if(size<body_size) {body=compressed;body_size=size;compressed_id=1;}}
        compression_cpu+=thread_cpu()-cpu_stamp;
        if(body_size+16<best_size) {best_size=body_size+16;header(best,bits,compressed_id,(uint32_t)n,(uint32_t)ne,step,mapping);memcpy(best+16,body,body_size);memcpy(best_decoded,decoded,8*n);}
    }
    Py_END_ALLOW_THREADS
    PyObject *result=NULL;
    if(failed) PyErr_SetString(PyExc_RuntimeError,"Native zlib compression failed");
    else result=Py_BuildValue("y#y#dd",best,(Py_ssize_t)best_size,best_decoded,(Py_ssize_t)(8*n),quantization_cpu,compression_cpu);
    free(raw);free(compressed);free(best);free(decoded);free(best_decoded);
    PyBuffer_Release(&a);PyBuffer_Release(&f);PyBuffer_Release(&widths);PyBuffer_Release(&scales);PyBuffer_Release(&bounds);return result;
cleanup_search:
    PyBuffer_Release(&a);PyBuffer_Release(&f);PyBuffer_Release(&widths);PyBuffer_Release(&scales);PyBuffer_Release(&bounds);return NULL;
}

/* Optional CPU kernels. Buffer inputs are checked before releasing the GIL.
 * No fast-math or fused multiply/add: reconstruction is float32 multiply. */
static PyObject *candidate(PyObject *self, PyObject *args) {
    Py_buffer a={0}, f={0}; int bits, channels,mapping=0; double step_arg;
    if (!PyArg_ParseTuple(args,"y*y*idi|i",&a,&f,&bits,&step_arg,&channels,&mapping)) return NULL;
    Py_ssize_t n=a.len/8;
    if (a.len%8 || n<1 || n>262144 || f.len!=n || channels<1 || n%channels ||
        !(bits==2 || bits==4 || bits==6 || bits==8 || bits==12 || bits==16) ||
        !isfinite(step_arg) || step_arg<0 || !(mapping==0||mapping==1) || fegetround()!=FE_TONEAREST) {
        PyBuffer_Release(&a); PyBuffer_Release(&f);
        return PyErr_Format(PyExc_ValueError,"Invalid native candidate buffers/parameters/rounding mode");
    }
    float step=(float)step_arg;
    if (!isfinite(step)) { PyBuffer_Release(&a); PyBuffer_Release(&f); return PyErr_Format(PyExc_ValueError,"Invalid native scale"); }
    Py_ssize_t size=(2*n*bits+7)/8;
    PyObject *decoded=PyBytes_FromStringAndSize(NULL,n*8), *packed=PyBytes_FromStringAndSize(NULL,size);
    PyObject *positions=PyBytes_FromStringAndSize(NULL,n*4), *exceptions=PyBytes_FromStringAndSize(NULL,n*8);
    if (!decoded || !packed || !positions || !exceptions) goto fail;
    const float *x=(const float *)a.buf; const uint8_t *flag=(const uint8_t *)f.buf;
    float *out=(float *)PyBytes_AS_STRING(decoded);
    uint8_t *codes=(uint8_t *)PyBytes_AS_STRING(packed), *pos=(uint8_t *)PyBytes_AS_STRING(positions), *exc=(uint8_t *)PyBytes_AS_STRING(exceptions);
    int maximum=(1<<(bits-1))-1, failed=0;
    double center=maximum+(mapping?0.5:0.0),shift=mapping?center:0;
    int lower=mapping?0:-maximum,upper=mapping?(1<<bits)-1:maximum;
    double limit=(double)step*center, sum=0, row_sum=0, max_square=0, max_row=0;
    Py_ssize_t count=0, row_count=0, ne=0, byte=0;
    uint64_t accumulator=0; int available=0;
    Py_BEGIN_ALLOW_THREADS
    for (Py_ssize_t j=0;j<n;j++) {
        float r=x[2*j], im=x[2*j+1];
        int finite=isfinite(r)&&isfinite(im);
        int exact=flag[j] || !finite || (r==0 && signbit(r)) || (im==0 && signbit(im)) || fabs((double)r)>limit || fabs((double)im)>limit;
        int qr=mapping?maximum+1:0, qi=qr;
        if (!exact && step!=0) {
            qr=(int)nearbyint((double)r/(double)step+shift);
            qi=(int)nearbyint((double)im/(double)step+shift);
        }
        if (qr<lower || qr>upper || qi<lower || qi>upper) { failed=1; break; }
        float cr=(float)qr-(float)shift,ci=(float)qi-(float)shift;
        out[2*j]=cr*step; out[2*j+1]=ci*step;
        if (!isfinite(out[2*j]) || !isfinite(out[2*j+1])) { failed=1; break; }
        uint32_t pair[2]={(uint32_t)(qr+(mapping?0:maximum)),(uint32_t)(qi+(mapping?0:maximum))};
        for(int k=0;k<2;k++) {
            accumulator|=(uint64_t)pair[k]<<available; available+=bits;
            while(available>=8) { codes[byte++]=(uint8_t)accumulator; accumulator>>=8; available-=8; }
        }
        if (exact) {
            uint32_t index=(uint32_t)j;
            for(int k=0;k<4;k++) pos[4*ne+k]=(uint8_t)(index>>(8*k));
            memcpy(exc+8*ne,(const char *)a.buf+8*j,8);
            memcpy(out+2*j,x+2*j,8); ne++;
        }
        if (!flag[j] && finite) {
            double er=(double)out[2*j]-r, ei=(double)out[2*j+1]-im;
            double squared=er*er+ei*ei;
            sum+=squared; row_sum+=squared; count++; row_count++;
            if(squared>max_square) max_square=squared;
        }
        if ((j+1)%channels==0) {
            double mse=row_count ? row_sum/(2*row_count) : 0;
            if(mse>max_row) max_row=mse;
            row_sum=0; row_count=0;
        }
    }
    if(available) codes[byte++]=(uint8_t)accumulator;
    Py_END_ALLOW_THREADS
    PyBuffer_Release(&a); PyBuffer_Release(&f);
    if(failed) { Py_DECREF(decoded);Py_DECREF(packed);Py_DECREF(positions);Py_DECREF(exceptions);Py_RETURN_NONE; }
    if(_PyBytes_Resize(&positions,ne*4)<0 || _PyBytes_Resize(&exceptions,ne*8)<0) goto alloc_fail;
    return Py_BuildValue("NNNNddd",decoded,packed,positions,exceptions,count?sum/(2*count):0,sqrt(max_square),max_row);
fail:
    PyBuffer_Release(&a); PyBuffer_Release(&f);
alloc_fail:
    Py_XDECREF(decoded); Py_XDECREF(packed); Py_XDECREF(positions); Py_XDECREF(exceptions); return NULL;
}

static PyObject *quality_kernel(PyObject *self, PyObject *args) {
    Py_buffer a={0}, b={0}, f={0}; Py_ssize_t rows, channels, products;
    if(!PyArg_ParseTuple(args,"y*y*y*nnn",&a,&b,&f,&rows,&channels,&products)) return NULL;
    Py_ssize_t n=a.len/8;
    if(a.len%8 || b.len!=a.len || f.len!=n || rows<1 || channels<1 || channels>65536 || products<1 || products>32 || n%(channels*products) || rows!=n/(channels*products)) {
        PyBuffer_Release(&a);PyBuffer_Release(&b);PyBuffer_Release(&f);
        return PyErr_Format(PyExc_ValueError,"Invalid native quality buffers/dimensions");
    }
    Py_ssize_t sizes[10]={products*7*8,rows*products*8,rows*products*8,channels*products*8,
        channels*products*8,(channels-1)*products*6*8,(channels-1)*products*8,products*products*6*8,products*products*8,n*8};
    PyObject *result=PyTuple_New(11); if(!result) goto fail_quality;
    for(int k=0;k<10;k++) {
        PyObject *value=PyBytes_FromStringAndSize(NULL,sizes[k]);
        if(!value) {Py_DECREF(result);goto fail_quality;}
        if(k!=9) memset(PyBytes_AS_STRING(value),0,sizes[k]);
        PyTuple_SET_ITEM(result,k,value);
    }
    double *global=(double *)PyBytes_AS_STRING(PyTuple_GET_ITEM(result,0));
    double *bmse=(double *)PyBytes_AS_STRING(PyTuple_GET_ITEM(result,1)), *bmax=(double *)PyBytes_AS_STRING(PyTuple_GET_ITEM(result,2));
    double *csum=(double *)PyBytes_AS_STRING(PyTuple_GET_ITEM(result,3));
    uint64_t *ccount=(uint64_t *)PyBytes_AS_STRING(PyTuple_GET_ITEM(result,4));
    double *adj=(double *)PyBytes_AS_STRING(PyTuple_GET_ITEM(result,5));
    uint64_t *acount=(uint64_t *)PyBytes_AS_STRING(PyTuple_GET_ITEM(result,6));
    double *cross=(double *)PyBytes_AS_STRING(PyTuple_GET_ITEM(result,7));
    uint64_t *xcount=(uint64_t *)PyBytes_AS_STRING(PyTuple_GET_ITEM(result,8));
    double *absolute=(double *)PyBytes_AS_STRING(PyTuple_GET_ITEM(result,9));
    const float *x=(const float *)a.buf, *y=(const float *)b.buf;
    const uint8_t *flag=(const uint8_t *)f.buf;
    Py_ssize_t protected_count=0; int failed=0;
    Py_BEGIN_ALLOW_THREADS
    for(Py_ssize_t r=0;r<rows && !failed;r++) {
        double prev_r[32]={0},prev_i[32]={0},row_sum[32]={0};
        uint64_t row_count[32]={0}; int prev_valid[32]={0};
        for(Py_ssize_t c=0;c<channels;c++) {
            double er[32]={0},ei[32]={0}; int valid[32]={0};
            for(Py_ssize_t p=0;p<products;p++) {
                Py_ssize_t j=(r*channels+c)*products+p, cp=c*products+p, rp=r*products+p;
                float xr=x[2*j],xi=x[2*j+1];
                int finite=isfinite(xr)&&isfinite(xi);
                int protect=flag[j] || !finite || (xr==0 && signbit(xr)) || (xi==0 && signbit(xi));
                if(protect) {protected_count++; if(memcmp(x+2*j,y+2*j,8)!=0) failed=1;}
                valid[p]=!flag[j] && finite;
                if(valid[p]) {
                    if(!isfinite(y[2*j]) || !isfinite(y[2*j+1])) {failed=1;break;}
                    er[p]=(double)y[2*j]-xr; ei[p]=(double)y[2*j+1]-xi;
                    double rr=er[p]*er[p],ii=ei[p]*ei[p], squared=rr+ii, magnitude=hypot(er[p],ei[p]);
                    double *g=global+7*p;
                    absolute[p*rows*channels+(Py_ssize_t)g[6]]=magnitude;
                    g[0]+=er[p];g[1]+=ei[p];g[2]+=rr;g[3]+=ii;g[6]++;
                    if(magnitude>g[4]) g[4]=magnitude;
                    row_sum[p]+=squared;row_count[p]++;csum[cp]+=squared;ccount[cp]++;
                    if(magnitude>bmax[rp]) bmax[rp]=magnitude;
                }
                if(c && valid[p] && prev_valid[p]) {
                    Py_ssize_t k=(c-1)*products+p; double *s=adj+6*k;
                    s[0]+=prev_r[p]*er[p]+prev_i[p]*ei[p];s[1]+=prev_i[p]*er[p]-prev_r[p]*ei[p];
                    s[2]+=prev_r[p];s[3]+=prev_i[p];s[4]+=er[p];s[5]+=ei[p];acount[k]++;
                }
                prev_r[p]=er[p];prev_i[p]=ei[p];prev_valid[p]=valid[p];
            }
            for(Py_ssize_t p=0;p<products;p++) if(valid[p]) {
                for(Py_ssize_t q=p;q<products;q++) if(valid[q]) {
                    Py_ssize_t k=p*products+q;double *s=cross+6*k;
                    s[0]+=er[p]*er[q]+ei[p]*ei[q];s[1]+=ei[p]*er[q]-er[p]*ei[q];
                    s[2]+=er[p];s[3]+=ei[p];s[4]+=er[q];s[5]+=ei[q];xcount[k]++;
                }
            }
        }
        for(Py_ssize_t p=0;p<products;p++) {
            double mse=row_count[p]?row_sum[p]/(2*row_count[p]):0;
            bmse[r*products+p]=mse;
            if(mse>global[7*p+5]) global[7*p+5]=mse;
        }
    }
    for(Py_ssize_t p=0;p<products;p++) for(Py_ssize_t q=0;q<p;q++) {
        double *s=cross+6*(p*products+q), *t=cross+6*(q*products+p);
        s[0]=t[0];s[1]=-t[1];s[2]=t[4];s[3]=t[5];s[4]=t[2];s[5]=t[3];
        xcount[p*products+q]=xcount[q*products+p];
    }
    Py_END_ALLOW_THREADS
    PyBuffer_Release(&a);PyBuffer_Release(&b);PyBuffer_Release(&f);
    if(failed) {Py_DECREF(result);return PyErr_Format(PyExc_ValueError,"Protected sample changed or nonfinite reconstruction");}
    PyObject *count=PyLong_FromSsize_t(protected_count);
    if(!count) {Py_DECREF(result);return NULL;}
    PyTuple_SET_ITEM(result,10,count);return result;
fail_quality:
    PyBuffer_Release(&a);PyBuffer_Release(&b);PyBuffer_Release(&f);return NULL;
}

/* Power levels are evaluated with the same ordered binary64 operations as the
 * reference decoder. No pow(), cached external codebook, or fused arithmetic. */
static float power_value(int q,int maximum,float radius,float center,float gamma) {
    double u=(double)(2*q-maximum)/maximum,magnitude=fabs(u);
    if(gamma!=1) {magnitude=magnitude*magnitude;if(gamma==.25f) magnitude=magnitude*magnitude;}
    double residual=copysign(magnitude,u)*(double)radius;
    return (float)(residual+(double)center);
}

static int power_code(double residual,double target,int maximum,float radius,float center,float gamma,float *value) {
    if(radius==0) {*value=center;return (maximum+1)/2;}
    double u=residual/(double)radius,magnitude=fabs(u);
    if(gamma!=1) {magnitude=sqrt(magnitude);if(gamma==.25f) magnitude=sqrt(magnitude);}
    double coordinate=(copysign(magnitude,u)+1)*(maximum/2.0);
    int lo=(int)floor(coordinate);if(lo<0) lo=0;if(lo>maximum) lo=maximum;
    int hi=lo<maximum?lo+1:maximum;
    float a=power_value(lo,maximum,radius,center,gamma),b=power_value(hi,maximum,radius,center,gamma);
    double da=fabs((double)a-target),db=fabs((double)b-target);
    if(db<da||(db==da&&hi%2==0)) {*value=b;return hi;}
    *value=a;return lo;
}

static PyObject *power_search(PyObject *self,PyObject *args) {
    Py_buffer a={0},f={0},widths={0},parameters={0},bounds={0};int channels,level;
    if(!PyArg_ParseTuple(args,"y*y*y*y*y*ii",&a,&f,&widths,&parameters,&bounds,&channels,&level)) return NULL;
    Py_ssize_t n=a.len/8,candidates=widths.len;
    int valid=a.len%8==0&&n>0&&n<=262144&&f.len==n&&channels>0&&n%channels==0&&candidates>0&&candidates<=4096&&
        parameters.len==16*candidates&&bounds.len==32&&level>=-1&&level<=9&&fegetround()==FE_TONEAREST;
    const uint8_t *width=(const uint8_t *)widths.buf;
    const float *params=(const float *)parameters.buf;
    for(Py_ssize_t k=0;k<candidates&&valid;k++) {
        const float *p=params+4*k;
        valid=(width[k]==2||width[k]==4||width[k]==6||width[k]==8||width[k]==12||width[k]==16)&&
              isfinite(p[0])&&p[0]>=0&&isfinite(p[1])&&isfinite(p[2])&&(p[3]==1||p[3]==.5f||p[3]==.25f);
    }
    if(!valid) {PyErr_SetString(PyExc_ValueError,"Invalid native power search buffers/parameters");goto cleanup_power;}
    uint8_t *raw=malloc(16*n),*compressed=malloc(compressBound(16*n)),*best=malloc(8*n+28);
    float *decoded=malloc(8*n),*best_decoded=malloc(8*n);
    if(!raw||!compressed||!best||!decoded||!best_decoded) {
        free(raw);free(compressed);free(best);free(decoded);free(best_decoded);PyErr_NoMemory();goto cleanup_power;
    }
    const float *x=(const float *)a.buf;const uint8_t *flags=(const uint8_t *)f.buf;const double *bound=(const double *)bounds.buf;
    size_t best_size=8*n+16;int failed=0;double quant_cpu=0,compress_cpu=0,stamp;
    Py_BEGIN_ALLOW_THREADS
    header(best,0,0,(uint32_t)n,0,0,0);memcpy(best+16,a.buf,8*n);memcpy(best_decoded,a.buf,8*n);
    if(level>=0) {stamp=thread_cpu();uLongf size=compressBound(8*n);
        if(compress2(compressed,&size,(const Bytef *)a.buf,8*n,level)!=Z_OK) failed=1;
        else if(size<(size_t)8*n) {header(best,0,1,(uint32_t)n,0,0,0);memcpy(best+16,compressed,size);best_size=size+16;}
        compress_cpu+=thread_cpu()-stamp;}
    for(Py_ssize_t k=0;k<candidates&&!failed;k++) {
        int bits=width[k],maximum=(1<<bits)-1;const float *p=params+4*k;
        float radius=p[0],cr=p[1],ci=p[2],gamma=p[3];
        size_t packed=(2*n*bits+7)/8,ne=0,byte=0;uint32_t *positions=(uint32_t *)(raw+4*n);
        double sum=0,sum_r=0,sum_i=0,row_sum=0,max_square=0,max_row=0;
        Py_ssize_t count=0,row_count=0;uint64_t accumulator=0;int available=0,rejected=0;
        stamp=thread_cpu();
        for(Py_ssize_t j=0;j<n;j++) {
            float r=x[2*j],im=x[2*j+1];int finite=isfinite(r)&&isfinite(im);
            double rr=(double)r-cr,ri=(double)im-ci;
            int exact=flags[j]||!finite||(r==0&&signbit(r))||(im==0&&signbit(im))||fabs(rr)>radius||fabs(ri)>radius;
            uint32_t pair[2];
            pair[0]=power_code(exact?0:rr,exact?cr:r,maximum,radius,cr,gamma,decoded+2*j);
            pair[1]=power_code(exact?0:ri,exact?ci:im,maximum,radius,ci,gamma,decoded+2*j+1);
            if(!isfinite(decoded[2*j])||!isfinite(decoded[2*j+1])) {rejected=1;break;}
            for(int h=0;h<2;h++) {accumulator|=(uint64_t)pair[h]<<available;available+=bits;
                while(available>=8) {raw[byte++]=(uint8_t)accumulator;accumulator>>=8;available-=8;}}
            if(exact) {positions[ne++]=(uint32_t)j;memcpy(decoded+2*j,x+2*j,8);}
            if(!flags[j]&&finite) {
                double er=(double)decoded[2*j]-r,ei=(double)decoded[2*j+1]-im,square=er*er+ei*ei;
                sum+=square;sum_r+=er;sum_i+=ei;row_sum+=square;count++;row_count++;
                if(square>max_square) max_square=square;
                /* Prefix sums have only nonnegative terms. Full cell/block
                 * counts overestimate valid denominators: safe early rejection. */
                if(square>bound[1]*bound[1]*(1+1e-9)||sum>bound[0]*(2*n)*(1+1e-9)||row_sum>bound[3]*(2*channels)*(1+1e-9)) {rejected=1;break;}
            }
            if((j+1)%channels==0) {double mse=row_count?row_sum/(2*row_count):0;if(mse>max_row) max_row=mse;row_sum=0;row_count=0;}
        }
        if(rejected) {quant_cpu+=thread_cpu()-stamp;continue;}
        if(available) raw[byte++]=(uint8_t)accumulator;
        double max_error=sqrt(max_square),mse=count?sum/(2*count):0,bias=count?hypot(sum_r/count,sum_i/count):0;
        if(mse>bound[0]*(1+1e-9)||max_error>bound[1]*(1+1e-9)||max_row>bound[3]*(1+1e-9)||bias>bound[2]+max_error*1e-9) {quant_cpu+=thread_cpu()-stamp;continue;}
        memmove(raw+packed,positions,4*ne);
        for(size_t j=0;j<ne;j++) {uint32_t index;memcpy(&index,raw+packed+4*j,4);memcpy(raw+packed+4*ne+8*j,x+2*index,8);}
        size_t body_size=packed+12*ne;uint8_t *body=raw;int compressed_id=0;
        quant_cpu+=thread_cpu()-stamp;stamp=thread_cpu();
        if(level>=0 && best_size>29 && body_size>1) {
            /* A compressed candidate must beat BOTH its raw representation and
             * the current winner including the 28-byte header. A bounded
             * output buffer stops deflate once it cannot win; Z_BUF_ERROR is
             * an ordinary losing candidate, not a damaged archive. */
            uLongf size=body_size-1;
            if(size>best_size-29) size=best_size-29;
            int rc=compress2(compressed,&size,raw,body_size,level);
            if(rc==Z_OK) {body=compressed;body_size=size;compressed_id=1;}
            else if(rc!=Z_BUF_ERROR) {failed=1;break;}
        }
        compress_cpu+=thread_cpu()-stamp;
        if(body_size+28<best_size) {best_size=body_size+28;header(best,bits,compressed_id,(uint32_t)n,(uint32_t)ne,radius,2);
            memcpy(best+16,p+1,12);memcpy(best+28,body,body_size);memcpy(best_decoded,decoded,8*n);}
    }
    Py_END_ALLOW_THREADS
    PyObject *result=NULL;
    if(failed) PyErr_SetString(PyExc_RuntimeError,"Native power compression failed");
    else result=Py_BuildValue("y#y#dd",best,(Py_ssize_t)best_size,best_decoded,(Py_ssize_t)(8*n),quant_cpu,compress_cpu);
    free(raw);free(compressed);free(best);free(decoded);free(best_decoded);
    PyBuffer_Release(&a);PyBuffer_Release(&f);PyBuffer_Release(&widths);PyBuffer_Release(&parameters);PyBuffer_Release(&bounds);return result;
cleanup_power:
    PyBuffer_Release(&a);PyBuffer_Release(&f);PyBuffer_Release(&widths);PyBuffer_Release(&parameters);PyBuffer_Release(&bounds);return NULL;
}

static PyMethodDef methods[]={
    {"candidate",candidate,METH_VARARGS,"Fused deterministic quantization, packing and conservative error precheck."},
    {"search",search,METH_VARARGS,"Complete adaptive candidate search with zlib, without holding the GIL."},
    {"quality",quality_kernel,METH_VARARGS,"Streaming error moments, paired covariances and protected-value checks."},
    {"power_search",power_search,METH_VARARGS,"Median-centered signed-power quantization with strict error selection."},
    {NULL,NULL,0,NULL}
};
static struct PyModuleDef module={PyModuleDef_HEAD_INIT,"_native",NULL,-1,methods};
PyMODINIT_FUNC PyInit__native(void) {
    PyObject *m=PyModule_Create(&module);
    if(m && PyModule_AddIntConstant(m,"MIDRISE_V1",1)<0) {Py_DECREF(m);return NULL;}
    if(m && PyModule_AddIntConstant(m,"POWER_V1",1)<0) {Py_DECREF(m);return NULL;}
    return m;
}
