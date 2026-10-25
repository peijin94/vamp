"""Matching single-Stokes reference runs for the standard I/V validation."""
import fcntl
import json
from pathlib import Path
import subprocess
import sys

from run_wsclean import snapshot

STANDARD_STOKES=('I','V')


def reference_path(root,stokes):
    if stokes not in STANDARD_STOKES:raise ValueError('Expected Stokes I or V')
    return Path(root)/('reference' if stokes=='I' else 'reference-v')


def image_paths(root,name,stokes,repeat=False):
    if stokes not in STANDARD_STOKES:raise ValueError('Expected Stokes I or V')
    stem=name+('' if stokes=='I' else '-v')
    return (Path(root)/(stem+('-repeat-images' if repeat else '-images')),
            Path(root)/(stem+('-repeat-comparison.json' if repeat else '-comparison.json')))


def restored_image(directory):
    paths=sorted(Path(directory).glob('image*-image.fits'))
    if len(paths)!=1:raise ValueError('Expected one restored single-Stokes image')
    return paths[0]


def ensure_references(source,root,stokes,env):
    """Create missing references once, or verify their source and settings.

    An advisory lock coordinates independent first/repeat experiment processes;
    incomplete reference directories fail explicitly and are never overwritten.
    """
    root=Path(root);root.mkdir(parents=True,exist_ok=True)
    source_hashes=snapshot(Path(source))
    result={}
    for pol in stokes:
        directory=reference_path(root,pol)
        with (root/f'.reference-{pol.lower()}.lock').open('a') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            if not (directory/'complete.json').exists():
                if directory.exists():raise ValueError('Incomplete reference requires a new reference-root: '+str(directory))
                command=[sys.executable,str(Path(__file__).with_name('run_wsclean.py')),str(source),
                         '--output',str(directory),'--data-column','CORRECTED_DATA',
                         '--clean-profile','calibrated','--threads','8','--polarizations',pol]
                with (root/f'reference-{pol.lower()}-driver.log').open('a') as log:
                    subprocess.run(command,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
            record=json.loads((directory/'complete.json').read_text())
            expected=dict(size=2048,scale='4amin',polarizations=pol,data_column='CORRECTED_DATA',
                          clean_profile='calibrated',weighting='briggs 0',multiscale=True,threads=8,
                          niter=150000,mgain=.95,auto_mask=3,auto_threshold=.5,
                          multiscale_scale_bias=.8,local_rms=True,horizon_mask='10deg',
                          reorder=True,memory_percent=40,update_model_required=False)
            if (record.get('status')!='complete' or record.get('source_unchanged') is not True
                    or record.get('source_sha256')!=source_hashes
                    or any(record['parameters'].get(k)!=v for k,v in expected.items())):
                raise ValueError('Reference source or calibrated imaging settings differ: '+str(directory))
            restored_image(directory)
            result[pol]=directory
    return result
