"""Copy and safely extract a supplied MS tar, checking its manifest bitwise."""
import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import tarfile
import time


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(8*1024**2),b''):h.update(block)
    return h.hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('source',type=Path);p.add_argument('--raw',type=Path,required=True)
    p.add_argument('--working',type=Path,required=True);p.add_argument('--report',type=Path,required=True)
    args=p.parse_args();started=time.perf_counter()
    source=args.source.resolve();manifest_path=source.with_name(source.name+'.manifest.json')
    manifest=json.loads(manifest_path.read_text());root=manifest['root_directory']
    if PurePosixPath(root).parts!=(root,) or root in ('.','..'):
        raise ValueError('Manifest root must be one safe directory name')
    checksum=source.with_name(source.name+'.sha256').read_text().split()[0]
    if checksum!=manifest['archive_sha256'] or source.stat().st_size!=manifest['archive_bytes']:
        raise ValueError('Provided archive length/checksum records disagree')
    target=args.raw/source.name;output=args.working/root
    if target.exists() or output.exists() or args.report.exists():raise ValueError('Outputs must be new')
    args.raw.mkdir(parents=True,exist_ok=True);args.working.mkdir(parents=True,exist_ok=True)
    partial=target.with_name(target.name+'.partial');h=hashlib.sha256()
    with source.open('rb') as src,partial.open('xb') as dst:
        for block in iter(lambda:src.read(8*1024**2),b''):dst.write(block);h.update(block)
        dst.flush();os.fsync(dst.fileno())
    if h.hexdigest()!=checksum or digest(partial)!=checksum:
        raise ValueError('Archive copy SHA-256 mismatch; partial file retained')
    os.link(partial,target);partial.unlink()
    for src in (manifest_path,source.with_name(source.name+'.sha256'),source.parent/'config.json'):
        if src.exists():
            with src.open('rb') as f,(args.raw/src.name).open('xb') as out:shutil.copyfileobj(f,out)
    staging=args.working/(root+'.partial');staging.mkdir()
    with tarfile.open(target,'r:') as archive:
        members=archive.getmembers();seen=set();files={}
        for member in members:
            path=PurePosixPath(member.name)
            if path.is_absolute() or '..' in path.parts or not path.parts or path.parts[0]!=root or path in seen:
                raise ValueError('Unsafe or duplicate tar member: '+member.name)
            seen.add(path)
            if not (member.isfile() or member.isdir()):raise ValueError('Unsupported tar entry: '+member.name)
            if member.isfile():files[str(PurePosixPath(*path.parts[1:]))]=member
        if set(files)!=set(manifest['files']) or len(files)!=manifest['file_count']:
            raise ValueError('Tar inventory differs from manifest')
        for relative,member in files.items():
            destination=staging/relative;destination.parent.mkdir(parents=True,exist_ok=True)
            record=manifest['files'][relative]
            if member.size!=record['bytes']:raise ValueError('Tar member size mismatch: '+relative)
            with archive.extractfile(member) as src,destination.open('xb') as dst:shutil.copyfileobj(src,dst,8*1024**2)
            if digest(destination)!=record['sha256']:raise ValueError('Extracted SHA-256 mismatch: '+relative)
            destination.chmod(member.mode & 0o777)
    staging.rename(output)
    result={'schema':'varquanms.manifest-import/1','status':'verified','source_tar':str(source),
            'local_tar':str(target),'extracted_ms':str(output),'tar_sha256':checksum,
            'archive_bytes':target.stat().st_size,'file_count':len(files),'payload_bytes':sum(v['bytes'] for v in manifest['files'].values()),
            'files':{k:{'bytes':v['bytes'],'sha256':v['sha256']} for k,v in manifest['files'].items()},
            'elapsed_seconds':time.perf_counter()-started,'note':'Source was read only. Copied tar and every extracted regular file match supplied SHA-256 records; metadata/layout/calibration inspection is separate.'}
    with args.report.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:result[k] for k in ('status','archive_bytes','file_count','payload_bytes','elapsed_seconds')},indent=2))


if __name__=='__main__':main()
