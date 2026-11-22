"""Build the optional C17 SCIP/QSopt_ex libraries from two locked source archives."""
import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_LOCK = ROOT / 'config/scip-integer-v15.lock.json'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def unpack(archive, digest, target):
    if sha(archive) != digest:
        raise ValueError(f'Unexpected archive digest: {archive}')
    target.mkdir()
    with tarfile.open(archive) as bundle:
        if sys.version_info >= (3, 12):
            bundle.extractall(target, filter='data')
        else:
            bundle.extractall(target)
    source, = target.iterdir()
    return source


def cmake(source, build, options, target, jobs):
    subprocess.run(['cmake', '-S', source, '-B', build, '-DCMAKE_BUILD_TYPE=Release',
                    '-DCMAKE_EXPORT_COMPILE_COMMANDS=ON', *options], check=True)
    subprocess.run(['cmake', '--build', build, '--target', target, '-j', str(jobs)], check=True)


def qs_sources(source, generated):
    generated.mkdir()
    make = (source/'Makefile.am').read_text().replace('\\\n', ' ')
    def names(key):
        return re.search('^'+key+r'\s*=\s*(.+)$', make, re.M)[1].split()
    for p in (source/'qsopt_ex').iterdir():
        if p.suffix in ('.c', '.h'):
            shutil.copy2(p, generated/p.name)
    units = [Path(f).name for f in names('MAIN_SOURCE_FILES')]
    for mode, typename in [('dbl','double'), ('mpq','mpq_t'), ('mpf','mpf_t')]:
        for key in ['TEMPLATE_SOURCE_FILES','TEMPLATE_PUBLIC_HEADER_FILES','TEMPLATE_PRIVATE_HEADER_FILES']:
            for name in names(key):
                original = source/name
                dst = generated/(original.stem+'_'+mode+original.suffix)
                dst.write_text(original.read_text().replace('EGLPNUM_TYPENAME', mode)
                               .replace('EGLPNUM_TYPE', typename))
                if dst.suffix == '.c':
                    units.append(dst.name)
    (generated/'config.h').write_text('''#define DEBUG 0
#define VERBOSE_LEVEL 0
#define HAVE_SYS_RESOURCE_H 1
#define HAVE_SYS_TIMES_H 1
#define HAVE_SYS_PARAM_H 1
#define HAVE_GETRUSAGE 1
#define HAVE_TIMES 1
#define HAVE_CLOCK 1
#define PACKAGE_NAME "QSopt_ex"
#define PACKAGE_VERSION "2.5.10.3"
#define PACKAGE_STRING "QSopt_ex 2.5.10.3"
''')
    (generated/'CMakeLists.txt').write_text('''cmake_minimum_required(VERSION 3.20)
project(qsopt C)
set(CMAKE_C_STANDARD 17)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)
find_path(GMP_INCLUDE_DIR gmp.h REQUIRED)
add_library(qsopt STATIC '''+' '.join(units)+''')
target_include_directories(qsopt PRIVATE . ${GMP_INCLUDE_DIR})
target_compile_definitions(qsopt PRIVATE HAVE_CONFIG_H=1 _DARWIN_C_SOURCE _DEFAULT_SOURCE)
set_source_files_properties(eg_lpnum.c PROPERTIES COMPILE_DEFINITIONS EG_LPNUM_MEMSLAB=0)
''')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scip-archive', type=Path, required=True)
    parser.add_argument('--qsopt-archive', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--gmp-prefix', type=Path)
    parser.add_argument('--lock', type=Path, default=DEFAULT_LOCK)
    parser.add_argument('--jobs', type=int, default=4)
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    lock_path = args.lock.resolve()
    config = json.loads(lock_path.read_text())
    for name, digest in config['patches'].items():
        if sha(ROOT/name) != digest:
            raise ValueError(f'Unexpected patch digest: {name}')
    scip_spec, qs_spec = config['sources']
    scip = unpack(args.scip_archive, scip_spec['sha256'], out/'upstream-scip')
    qs = unpack(args.qsopt_archive, qs_spec['sha256'], out/'upstream-qsopt')
    overlays = config.get('overlays', {})
    overlay_paths = {p for paths in overlays.values() for p in paths}
    base_patches = [p for p in config['patches'] if p.endswith('.patch') and p not in overlay_paths]
    scip_patch, = [p for p in base_patches if 'scip-9.2.4' in p]
    qs_patch, = [p for p in base_patches if 'qsopt-ex' in p]
    for tree, spec, base in [(scip, scip_spec, scip_patch), (qs, qs_spec, qs_patch)]:
        for patch in [base, *overlays.get(spec['name'], [])]:
            if patch not in config['patches']:
                raise ValueError(f'Overlay missing digest: {patch}')
            subprocess.run(['patch', '-p1', '-i', ROOT/patch], cwd=tree, check=True)
    # The disabled PaPILO and symmetry implementations consist entirely of C.
    for base in ['scip/presol_milp','symmetry/compute_symmetry_none']:
        shutil.copy2(scip/'src'/(base+'.cpp'), scip/'src'/(base+'.c'))
    common = ['-DCMAKE_PREFIX_PATH='+str(args.gmp_prefix.resolve())] if args.gmp_prefix else []
    generated = out/'qs-build'
    qs_sources(qs, generated)
    cmake(generated, out/'qs-obj', common, 'qsopt', args.jobs)
    compat = out/'qso-compat'
    compat.mkdir()
    for header in generated.glob('*.h'):
        shutil.copy2(header, compat/header.name)
    shutil.copy2(ROOT/'config/patches/scip-qsopt-ex-compat.h', compat/'qsopt.h')
    qslib = out/'qs-obj/libqsopt.a'
    build = out/'scip-obj'
    cflags = os.getenv('CFLAGS', '')+' -D_DARWIN_C_SOURCE -D_DEFAULT_SOURCE'
    options = ['-DAUTOBUILD=OFF','-DLPS=qso','-DSYM=none','-DTPI=none','-DEXPRINT=none',
               '-DIPOPT=OFF','-DPAPILO=OFF','-DREADLINE=OFF','-DAMPL=OFF','-DWORHP=OFF',
               '-DZIMPL=OFF','-DSHARED=OFF','-DBUILD_TESTING=OFF','-DCMAKE_C_FLAGS='+cflags,
               '-DQSO_INCLUDE_DIR='+str(compat),'-DQSO_LIBRARY='+str(qslib)]
    cmake(scip, build, common+options, 'libscip', args.jobs)
    compiled = [r for r in json.loads((build/'compile_commands.json').read_text())
                if 'libscip.dir/' in r['command']]
    if not compiled or any(not r['file'].endswith('.c') or '-std=c17' not in r['command']
                           for r in compiled):
        raise ValueError('Expected a C17-only SCIP library')
    lib, = (build/'lib').glob('libscip.a')
    sciplib = out/'libscip-c17.a'
    shutil.copy2(lib, sciplib)
    for source, spec in [(scip,scip_spec), (qs,qs_spec)]:
        shutil.copy2(source/spec['license_file'], out/(spec['name']+'-LICENSE'))
    cache = out/'lab-cache.cmake'
    cache.write_text('\n'.join(f'set({k} [=[{val}]=] CACHE {kind} "" FORCE)'
        for k,val,kind in [('LAB_SCIP_INTEGER','ON','BOOL'),('SCIP_INCLUDE_DIR',scip/'src','PATH'),
                           ('SCIP_CONFIG_INCLUDE_DIR',build,'PATH'),('SCIP_C17_LIBRARY',sciplib,'FILEPATH'),
                           ('QSOPT_C17_LIBRARY',qslib,'FILEPATH')])+'\n')
    (out/'build.json').write_text(json.dumps({'lock_sha256':sha(lock_path),
            'scip_c17_units':len(compiled),'libraries':{str(p):sha(p) for p in [sciplib,qslib]}},indent=2)+'\n')
    print(f'Configure MILP Lab with cmake -C "{cache}" -S . -B build-scip')


if __name__ == '__main__':
    main()
