"""Build the original VVamp special-function stage with native GiNaC/CLN.

Requires g++, GiNaC and CLN development headers/libraries. The optional --prefix
points to an unpacked installation's usr directory. The build directory must
be new; no externally generated object or Python special-function value is used.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import re
import subprocess
import tarfile

ARCHIVE_SHA='1c69a6beb3be85b07225af867586f242e999337a1ccf1fb940016a5341dfc091'
HERE=Path(__file__).resolve().parent

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def replace_once(path,old,new):
    s=path.read_text()
    if s.count(old)!=1:raise ValueError('Unexpected original source: '+str(path))
    path.write_text(s.replace(old,new))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('archive',type=Path)
    p.add_argument('build',type=Path)
    p.add_argument('--prefix',type=Path)
    p.add_argument('--workers',type=int,default=2)
    args=p.parse_args();build=args.build.resolve()
    if build.exists():raise ValueError('Use a new native build directory')
    if sha(args.archive)!=ARCHIVE_SHA:raise ValueError('VVamp archive checksum mismatch')
    if args.workers<1:raise ValueError('Positive worker count required')
    build.mkdir(parents=True)
    with tarfile.open(args.archive) as tar:tar.extractall(build,filter='data')
    source=build/'ggvvamp-1.0';header=source/'ggvvamp.h'
    # These are author-provided compile-time precision switches.
    replace_once(header,'#define QUAD_PRECISION','// #define QUAD_PRECISION')
    replace_once(header,'#define DOUBLE_PRECISION','// #define DOUBLE_PRECISION')
    funcs=source/'p0002.cpp';s=funcs.read_text()
    if not re.search(r'\blst\(',s):raise ValueError('Expected legacy GiNaC calls absent')
    funcs.write_text('#include "native_ginac_compat.h"\n'+re.sub(r'\blst\(','legacy_lst(',s))
    units=[source/'ggvvamp.cpp',funcs]
    # Real GiNaC headers must precede the arithmetic-only checker's thin stub.
    compat=HERE/'cxx/native_ginac_compat.h'
    (source/compat.name).write_bytes(compat.read_bytes())
    include=['-I'+str(source)];link=['-Wl,--gc-sections','-lginac','-lcln']
    if args.prefix:
        prefix=args.prefix.resolve();lib=prefix/'lib/x86_64-linux-gnu'
        include=['-I'+str(prefix/'include'),*include]
        link=['-L'+str(lib),'-Wl,-rpath,'+str(lib),*link]
    flags=['g++','-O0','-std=gnu++11','-fPIC','-ffunction-sections','-fdata-sections',*include]
    def compile_unit(unit):
        obj=unit.with_suffix('.o')
        subprocess.run([*flags,'-c',str(unit),'-o',str(obj)],check=True)
        return obj
    with ThreadPoolExecutor(args.workers) as pool:
        futures=[pool.submit(compile_unit,u) for u in units]
        for i,future in enumerate(as_completed(futures),1):
            print(i,'/',len(units),future.result().name,flush=True)
    executable=build/'native_functions'
    subprocess.run([*flags,str(HERE/'cxx/native_functions.cpp'),
                    *map(str,(u.with_suffix('.o') for u in units)),*link,'-o',str(executable)],check=True)
    manifest={'archive_sha256':ARCHIVE_SHA,'executable_sha256':sha(executable),
              'compiler':subprocess.run(['g++','--version'],check=True,capture_output=True,text=True).stdout,
              'scope':'original special-function stage only','compile_flags':flags,'link_flags':link,
              'adaptations':['CLN precision switch','modern GiNaC list construction',
                             'unused amplitude sections discarded by linker'],
              'adapter_sources':{str(f.relative_to(HERE)):sha(f) for f in
                                (Path(__file__),HERE/'cxx/native_functions.cpp',HERE/'cxx/native_ginac_compat.h')}}
    (build/'native_build.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('Built',executable,flush=True)

if __name__=='__main__':main()
