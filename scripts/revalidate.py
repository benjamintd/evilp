#!/usr/bin/env python3
"""Rerun independent validation without changing historical solver runs."""
import argparse
import collections
import hashlib
import json
import pathlib
import re
import subprocess
import time

from run_evaluation import EVALUATION_SHA256, assess_evidence, strict_json, read_witness_ledgers, witness_index

ROOT = pathlib.Path(__file__).resolve().parents[1]
EVIDENCE_FILES = {'--solution': 'solution.sol', '--dual': 'dual.sol', '--farkas': 'infeasible.farkas'}
VALID_STATUSES = {0: {'valid'}, 1: {'invalid'}, 2: {'error'}, 3: {'unsupported'},
                  4: {'no_solution', 'unavailable'}}


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def checker_command(row, prefix, model, checker):
    solution = prefix / 'solution.sol'
    dual = prefix / 'dual.sol'
    recorded = row.get('checker_argv')
    if recorded:
        # The campaign directory may have moved. Rebind each input by its role
        # and keep the original argument order, format and numerical options.
        args = list(recorded)
        args[0] = str(checker)
        inputs = {'--model': model, **{flag: prefix / name for flag, name in EVIDENCE_FILES.items()}}
        for index, flag in enumerate(recorded[:-1]):
            if flag in inputs:
                args[index + 1] = str(inputs[flag])
        return args

    # Older runner versions recorded dual_file and the checker outcome, but
    # only signalled primal presence by writing solution.sol.
    if 'checker_argv' in row and row.get('checker_exit_code') is None:
        return None
    has_dual = bool(row.get('dual_file'))
    has_primal = (solution.is_file() or row['checker'].get('objective') is not None
                  or row['checker'].get('status') == 'no_solution'
                  or row['checker'].get('checked_feasible') is True
                  or (not has_dual and row.get('checker_exit_code') is not None))
    if not has_primal and not has_dual:
        return None
    args = [str(checker), '--model', str(model)]
    if has_primal:
        fmt = 'simple' if row['solver'] == 'lab' else row['solver']
        args.extend(['--solution', str(solution), '--format', fmt])
    if row.get('relaxation', False):
        args.append('--relax')
    if has_dual:
        args.extend(['--dual', str(dual)])
    return args


def validate(args, prefix, hashes):
    if args is None and hashes:
        return {'status': 'error', 'checked_feasible': False,
                'message': 'Recorded evidence has no checker invocation'}, None, 0, None
    if args is None:
        return {'status': 'unavailable', 'checked_feasible': False}, None, 0, None
    code = None
    invoked = None
    started = None
    try:
        if hashes is not None:
            for flag, name in EVIDENCE_FILES.items():
                if (flag in args) != (name in hashes):
                    raise ValueError(f'checker arguments disagree with recorded evidence for {name}')
                if name in hashes and sha(prefix / name) != hashes[name]:
                    raise ValueError(f'retained evidence SHA-256 mismatch for {name}')
        started = time.monotonic()
        invoked = args
        process = subprocess.run(args, capture_output=True, timeout=120)
        code = process.returncode
        if code not in VALID_STATUSES:
            raise ValueError(f'checker exited with code {code}')
        result = strict_json(process.stdout)
        status = result.get('status') if isinstance(result, dict) else None
        if status == 'valid' and result.get('relaxation') is not ('--relax' in args):
            raise ValueError('checker report disagrees with requested relaxation scope')
        if not isinstance(status, str) or status not in VALID_STATUSES[code]:
            raise ValueError(f'checker report disagrees with exit code {code}')
    except (ValueError, subprocess.TimeoutExpired, OSError) as error:
        result = {'status': 'error', 'checked_feasible': False, 'message': str(error)}
    return result, code, time.monotonic() - started if started is not None else 0, invoked


def write_new(path, record):
    with path.open('x') as stream:
        json.dump(record, stream, indent=2)
        stream.write('\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=pathlib.Path)
    parser.add_argument('--checker', type=pathlib.Path, default=pathlib.Path('build/milpcheck'))
    parser.add_argument('--witness-runs', action='append', type=pathlib.Path, default=[])
    parser.add_argument('--name', default='revalidated')
    args = parser.parse_args()
    if not re.fullmatch(r'[A-Za-z0-9_-]+', args.name):
        parser.error('Validation name must use letters, digits, underscore or hyphen')

    directory = args.directory.resolve()
    checker = args.checker.resolve()
    identity = {'checker_sha256': sha(checker), 'driver_sha256': sha(pathlib.Path(__file__))}
    identity['evaluation_sha256'] = EVALUATION_SHA256
    record = 'checker-' + args.name + '.json'
    archived, identity['witness_ledgers'] = read_witness_ledgers(args.witness_runs)
    checked_rows = []
    counts = collections.Counter()
    with (directory / 'runs.jsonl').open() as rows, (directory / (args.name + '.jsonl')).open('x') as out:
        for line in rows:
            if not line.strip():
                continue
            row = strict_json(line)
            prefix = directory / row['instance'] / row['solver'] / ('seed-' + str(row['seed']))
            model = pathlib.Path(row['model_path'])
            if not model.exists():
                model = ROOT / 'data/miplib2017/instances' / (row['instance'] + '.mps.gz')
            if sha(model) != row['model_sha256']:
                raise SystemExit('Model changed: ' + row['instance'])
            command = checker_command(row, prefix, model, checker)
            result, code, wall, invoked = validate(command, prefix, row.get('evidence_sha256'))
            original = {key: row[key] for key in ('checker_sha256', 'checker', 'checker_exit_code',
                                                   'checker_wall_s', 'checker_argv') if key in row}
            row.update(original_validation=original, checker_sha256=identity['checker_sha256'],
                       checker=result, checker_argv=invoked, checker_exit_code=code, checker_wall_s=wall,
                       revalidation={**identity, 'solver_times_unchanged': True, 'record': record})
            checked_rows.append((prefix, row))
        witnesses = witness_index([*(row for _, row in checked_rows), *archived])
        for prefix, row in checked_rows:
            row['evaluation'] = assess_evidence(row, row, relaxation=row.get('relaxation', False),
                                                witnesses=witnesses)
            write_new(prefix / record, {**identity, 'checker_argv': row['checker_argv'],
                                        'exit_code': row['checker_exit_code'], 'checker': row['checker'],
                                        'evaluation': row['evaluation']})
            out.write(json.dumps(row) + '\n')
            counts[row['evaluation']['checker_outcome']] += 1

    if sha(checker) != identity['checker_sha256']:
        raise SystemExit('Checker changed during validation')
    metadata = 'revalidation.json' if args.name == 'revalidated' else args.name + '.meta.json'
    write_new(directory / metadata, {**identity, 'runs': sum(counts.values()), 'outcomes': counts,
                                     'solver_times_unchanged': True})
    print(json.dumps(counts, indent=2))
    return int(counts['error'] > 0)


if __name__ == '__main__':
    raise SystemExit(main())
