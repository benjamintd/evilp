#!/usr/bin/env python3
"""Independent test oracle for exported algebra, not a cut-validity checker."""
import argparse
from fractions import Fraction as Q
import hashlib
import json
import pathlib


def matrix(raw):
    result = [[Q(0) for _ in range(raw['cols'])] for _ in range(raw['rows'])]
    for j in range(raw['cols']):
        for k in range(raw['start'][j], raw['start'][j + 1]):
            result[raw['row'][k]][j] += Q(raw['value'][k])
    return result


def decimal(value):
    value = Q(value)
    denominator = value.denominator
    twos = fives = 0
    while denominator % 2 == 0:
        twos += 1
        denominator //= 2
    while denominator % 5 == 0:
        fives += 1
        denominator //= 5
    assert denominator == 1
    places = max(twos, fives)
    digits = str(abs(value.numerator) * 2 ** (places - twos) * 5 ** (places - fives)).zfill(places + 1)
    text = digits if places == 0 else digits[:-places] + '.' + digits[-places:]
    return ('-' if value < 0 else '') + text


def fixture_mps(record):
    """Own-generated source matrix with the recorded box; no global premise claim."""
    a = matrix(record['source'])
    text = ['NAME exported_fixture', 'ROWS', ' N objective']
    rhs = []
    ranges = []
    for i, box in enumerate(record['source_row_box']):
        lo, hi = box['lower'], box['upper']
        if lo is not None and hi is not None and Q(lo) == Q(hi):
            kind, side = 'E', lo
        elif hi is not None:
            kind, side = 'L', hi
            if lo is not None:
                ranges.append(f' rng r{i} {decimal(Q(hi) - Q(lo))}')
        else:
            assert lo is not None
            kind, side = 'G', lo
        text.append(f' {kind} r{i}')
        rhs.append(f' rhs r{i} {decimal(side)}')
    text.append('COLUMNS')
    for j in range(len(record['box'])):
        if record['integer'][j]:
            text.append(f" mark{j}a 'MARKER' 'INTORG'")
        text.append(f' x{j} objective 0')
        for i in range(len(a)):
            if a[i][j]:
                text.append(f' x{j} r{i} {decimal(a[i][j])}')
        if record['integer'][j]:
            text.append(f" mark{j}b 'MARKER' 'INTEND'")
    text.extend(['RHS', *rhs])
    if ranges:
        text.extend(['RANGES', *ranges])
    text.append('BOUNDS')
    for j, box in enumerate(record['box']):
        text.append(f' FR box x{j}')
        if box['lower'] is not None:
            text.append(f" LO box x{j} {decimal(box['lower'])}")
        if box['upper'] is not None:
            text.append(f" UP box x{j} {decimal(box['upper'])}")
    text.append('ENDATA')
    return '\n'.join(text) + '\n'


def prefix(record, row):
    text = ['MILPLAB_GMI_1', 'SCOPE ORIGINAL_MIP_IN_BOX', f"BASIC x{row['basic']}", f"BETA {row['beta']}"]
    for j, box in enumerate(record['box']):
        text.append(f"BOX x{j} {box['lower'] if box['lower'] is not None else '*'} {box['upper'] if box['upper'] is not None else '*'}")
    for k, term in enumerate(row['terms']):
        kind = ('VAR', 'ROW', 'BOUND')[term['kind']]
        name = ('r' if term['kind'] == 1 else 'x') + str(term['index'])
        text.append(f"TERM {k} {term['coefficient']} {'I' if term['integer'] else 'C'} {kind} {name} {'U' if term['upper'] else 'L'} {term['anchor']}")
    text += [f"EQ r{item['index']} {item['multiplier']}" for item in row['equalities']]
    text += [f"FIX x{item['index']} {item['multiplier']}" for item in row['fixed']]
    # PUB_RHS/PUB/END intentionally absent: this is NOT an accepted proof.
    return '\n'.join(text) + '\n'


def verify(record):
    a = matrix(record['source'])
    b = matrix(record['basis'])
    canonical = matrix(record['canonical'])
    n = record['dimension']
    z = list(map(Q, record['basic_values']))
    rhs = list(map(Q, record['rhs']))
    assert all(sum(b[i][j] * z[j] for j in range(n)) == rhs[i] for i in range(n))
    assert all(value >= 0 for value in z)
    for row in record['tableaux']:
        lam = list(map(Q, row['weights']))
        p = row['basis_position']
        assert all(sum(lam[i] * b[i][j] for i in range(n)) == (1 if p == j else 0) for j in range(n))
        assert sum(lam[i] * rhs[i] for i in range(n)) == z[p]
        coef = [Q(0)] * len(record['box'])
        coef[row['basic']] = Q(1)
        constant = -Q(row['beta'])
        beta = Q(row['basic_shift']) + Q(row['signed_basic_scale']) * z[p]
        for term in row['terms']:
            alpha = sum(lam[i] * canonical[i][term['canonical_column']] for i in range(n))
            assert Q(term['coefficient']) == Q(row['signed_basic_scale']) * alpha * Q(term['canonical_scale'])
            beta -= Q(row['signed_basic_scale']) * alpha * Q(term['canonical_offset'])
            sign = -1 if term['upper'] else 1
            factor = sign * Q(term['coefficient'])
            index = term['index']
            if term['kind'] == 1:
                assert Q(term['anchor']) == Q(record['source_row_box'][index]['upper' if term['upper'] else 'lower'])
                for j in range(len(coef)):
                    coef[j] += factor * a[index][j]
            else:
                coef[index] += factor
                assert Q(term['anchor']) == Q(record['box'][index]['upper' if term['upper'] else 'lower'])
            constant -= factor * Q(term['anchor'])
            if term['integer']:
                assert term['kind'] == 0 and record['integer'][index] and Q(term['anchor']).denominator == 1
        assert beta == Q(row['beta'])
        for item in row['equalities']:
            index, factor = item['index'], Q(item['multiplier'])
            box = record['source_row_box'][index]
            assert box['lower'] is not None and box['upper'] is not None and Q(box['lower']) == Q(box['upper'])
            for j in range(len(coef)):
                coef[j] -= factor * a[index][j]
            constant += factor * Q(box['lower'])
        for item in row['fixed']:
            index, factor = item['index'], Q(item['multiplier'])
            box = record['box'][index]
            assert box['lower'] is not None and box['upper'] is not None and Q(box['lower']) == Q(box['upper'])
            coef[index] -= factor
            constant += factor * Q(box['lower'])
        assert not any(coef) and constant == 0
    return len(record['tableaux'])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('batches', type=pathlib.Path)
    parser.add_argument('--out', required=True, type=pathlib.Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    records = [json.loads(line) for line in args.batches.read_text().splitlines()]
    rows = 0
    for index, record in enumerate(records):
        rows += verify(record)
        directory = args.out / f'batch-{index:04}'
        directory.mkdir()
        (directory / 'source-with-recorded-box.mps').write_text(fixture_mps(record))
        for j, row in enumerate(record['tableaux']):
            (directory / f'row-{j:02}.gmi-prefix').write_text(prefix(record, row))
    summary = {'passed': True, 'scope': 'Independent Fraction test oracle verifies exported algebra only. Prefix files are incomplete derivations with no published cut and are not validity certificates.',
               'exported_batches': len(records), 'exported_rows': rows,
               'input_sha256': hashlib.sha256(args.batches.read_bytes()).hexdigest(),
               'checker_milestone_required': True, 'cut_insertion': False,
               'independent_cut_validity_claim': False, 'global_completion_claim': False}
    (args.out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary))


if __name__ == '__main__':
    main()
