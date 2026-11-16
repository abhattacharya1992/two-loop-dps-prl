#!/usr/bin/env python3
"""Fetch checksum-pinned published VVamp input and rebuild Python-local translation.

Only generated inputs/build products under this directory are created.
No third-party source is committed. An existing verified local archive can
be supplied with --archive for offline use. Rebuild on changing Python versions:
algebra.marshal is an implementation-specific generated cache, not source code.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile
import urllib.request

HERE=Path(__file__).resolve().parent
URL="https://lcgpackages.web.cern.ch/lcgpackages/tarFiles/sources/ggvvamp-1.0.tar.gz"
SHA256="1c69a6beb3be85b07225af867586f242e999337a1ccf1fb940016a5341dfc091"


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):
            h.update(block)
    return h.hexdigest()


def safe_unpack(path,target):
    with tarfile.open(path,'r:gz') as arc:
        # Reject special files and links as well as path traversal.
        members=arc.getmembers()
        for entry in members:
            p=Path(entry.name)
            if p.is_absolute() or '..' in p.parts or not (entry.isfile() or entry.isdir()):
                raise ValueError(f"Unsafe archive entry: {entry.name}")
        if hasattr(tarfile,'data_filter'):
            arc.extractall(target,members=members,filter='data')
        else:
            arc.extractall(target,members=members)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--archive',type=Path,help='previously downloaded ggvvamp-1.0.tar.gz')
    args=p.parse_args()
    inputs=HERE/'inputs'; build=HERE/'build'
    inputs.mkdir(exist_ok=True);build.mkdir(exist_ok=True)
    archive=args.archive.resolve() if args.archive else inputs/'ggvvamp-1.0.tar.gz'
    if not archive.exists():
        if args.archive:
            raise FileNotFoundError(archive)
        with tempfile.NamedTemporaryFile(dir=inputs,delete=False) as f:
            tmp=Path(f.name)
        try:
            request=urllib.request.Request(URL,headers={'User-Agent':'DPS-unequal-virtuality-reproduction/1'})
            with urllib.request.urlopen(request,timeout=120) as response,tmp.open('wb') as out:
                shutil.copyfileobj(response,out)
            if digest(tmp)!=SHA256:
                raise ValueError('Downloaded archive checksum mismatch; file not used')
            tmp.replace(archive)
        finally:
            tmp.unlink(missing_ok=True)
    if digest(archive)!=SHA256:
        raise ValueError('Archive checksum mismatch; file not used')
    safe_unpack(archive,build)
    translator=HERE.parent/'zz/scripts/translate_ggvvamp.py'
    subprocess.run([sys.executable,str(translator),str(build/'ggvvamp-1.0'),str(build/'ggv')],check=True)
    info={'url':URL,'archive_sha256':SHA256,'python':sys.version,'translator_sha256':digest(translator),
          'crossing_source_sha256':digest(build/'ggvvamp-1.0/ggvvampprec.cpp')}
    (build/'input_provenance.json').write_text(json.dumps(info,indent=2)+'\n')
    print('Verified input and translated evaluator:',build/'ggv')


if __name__=='__main__':
    main()
