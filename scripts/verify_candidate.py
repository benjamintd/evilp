"""Verify the pinned feature build and record the executable identities tested by CTest."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
from xml.etree import ElementTree
from build_scip_integer import locked_patch


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ('milpsolve', 'milpbench', 'milpcheck')
FEATURES = ('LAB_SCIP_INTEGER', 'LAB_GLPK_BASIS', 'LAB_KISSAT')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def identity(build):
    cache = {}
    for line in (build / 'CMakeCache.txt').read_text().splitlines():
        if line and not line.startswith(('#', '//')) and ':' in line and '=' in line:
            key, value = line.split('=', 1)
            cache[key.split(':', 1)[0]] = value
    lock = ROOT / (ROOT / 'config/candidate-scip-lock.txt').read_text().strip()
    if not all(cache.get(key) == 'ON' for key in FEATURES):
        raise ValueError('candidate requires SCIP, GLPK and Kissat')
    backend = Path(cache['SCIP_C17_LIBRARY']).parent / 'build.json'
    record = json.loads(backend.read_text())
    if record['lock_sha256'] != sha(lock):
        raise ValueError('SCIP build does not match the designated candidate lock')
    source = Path(cache['CMAKE_HOME_DIRECTORY'])
    for patch, expected in json.loads(lock.read_text())['patches'].items():
        locked_patch(source, patch, expected)
    libraries = {}
    for key in ('SCIP_C17_LIBRARY', 'QSOPT_C17_LIBRARY', 'GLPK_LIBRARY', 'KISSAT_LIBRARY'):
        path = Path(cache[key]).resolve()
        libraries[key] = sha(path)
        if key in ('SCIP_C17_LIBRARY', 'QSOPT_C17_LIBRARY') and record['libraries'].get(str(path)) != libraries[key]:
            raise ValueError(f'{key} differs from dependency build record')
    return {'lock': lock.name, 'lock_sha256': sha(lock), 'features': list(FEATURES),
            'build_type': cache['CMAKE_BUILD_TYPE'], 'project_sanitizers': cache['LAB_SANITIZE'] == 'ON',
            'binaries': {tool: sha(build / tool) for tool in TOOLS}, 'libraries': libraries}


def verified(build, receipt):
    data = json.loads(receipt.read_text())
    if data['identity'] != identity(build) or not data['passed']:
        raise ValueError('campaign binary or configuration differs from the passing verification')
    junit = receipt.with_suffix('.xml')
    if sha(junit) != data['contracts_sha256']:
        raise ValueError('candidate verification evidence changed')
    return data


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--build', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--jobs', type=int, default=2)
    args = p.parse_args()
    build, output = args.build.resolve(), args.out.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    xml = output.with_suffix('.xml')
    before = identity(build)
    # Reserve both artifacts. Interrupted verification cannot authorize a campaign.
    with output.open('x') as stream:
        xml.touch(exist_ok=False)
        test = subprocess.run(['ctest', '--test-dir', build, '--output-on-failure', '--no-tests=error',
                               '-j', str(args.jobs), '--output-junit', xml])
        cases = list(ElementTree.parse(xml).iter('testcase'))
        required = {'scip_mip_contract', 'contracted_integer_mip_contract', 'scip_bound_evidence_contract',
                    'binary_sat_contract', 'qso_row_bounds_contract', 'v16_engine_contract'}
        passed = (test.returncode == 0 and required <= {case.get('name') for case in cases}
                  and all(case.get('status') == 'run' and not any(case.find(tag) is not None
                          for tag in ('failure', 'error', 'skipped')) for case in cases)
                  and identity(build) == before)
        json.dump({'identity': before, 'contracts': len(cases), 'contracts_sha256': sha(xml),
                   'passed': passed}, stream, indent=2)
        stream.write('\n')
    return 0 if passed else 1


if __name__ == '__main__':
    raise SystemExit(main())
