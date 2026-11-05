#!/usr/bin/env python3
"""Replay individual LP-box Farkas proofs; this does not certify a search tree."""
import argparse
import fractions
import gzip
import hashlib
import json
import math
import pathlib
import re
import subprocess
import sys

SCOPE = 'individual LP boxes; not a branch-and-bound tree proof'
RATIONAL = re.compile(r'[+-]?[0-9]+(?:/[0-9]+)?\Z')
SECTIONS = {'NAME', 'ROWS', 'COLUMNS', 'RHS', 'RANGES', 'BOUNDS',
            'ENDATA', 'OBJSENSE', 'OBJNAME'}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f'duplicate JSON key {key!r}')
        result[key] = value
    return result


def strict_json(text):
    def invalid_constant(value):
        raise ValueError(f'nonfinite JSON number {value}')
    def finite_float(value):
        parsed = float(value)
        if not math.isfinite(parsed):
            invalid_constant(value)
        return parsed
    return json.loads(text, object_pairs_hook=unique_object,
                      parse_constant=invalid_constant, parse_float=finite_float)


def rational(value):
    if not isinstance(value, str) or len(value) > 4096 or not RATIONAL.fullmatch(value):
        raise ValueError('expected an integer or integer/positive-integer rational string')
    try:
        return fractions.Fraction(value)
    except (ValueError, ZeroDivisionError) as error:
        raise ValueError('invalid rational string') from error


def finite_decimal(value):
    """Convert without rounding; MPS decimals cannot express other denominators."""
    value = rational(value)
    denominator = value.denominator
    twos = fives = 0
    while denominator % 2 == 0:
        denominator //= 2
        twos += 1
    while denominator % 5 == 0:
        denominator //= 5
        fives += 1
    if denominator != 1:
        raise ValueError('LP bound has no exact finite decimal (denominator contains other than 2/5)')
    places = max(twos, fives)
    scaled = abs(value.numerator) * 2**(places - twos) * 5**(places - fives)
    digits = str(scaled).zfill(places + 1)
    text = digits if not places else digits[:-places] + '.' + digits[-places:]
    return ('-' if value < 0 else '') + text


def model_without_bounds(path):
    """Retain every non-BOUNDS source line, including original decimal tokens."""
    raw = path.read_bytes()
    if raw.startswith(b'\x1f\x8b'):
        raw = gzip.decompress(raw)
    text = raw.decode('utf-8')
    if '\x00' in text:
        raise ValueError('NUL byte in source model')
    columns = {}
    retained = []
    section, last_column, ended = None, None, False
    for line in text.splitlines(keepends=True):
        fields = line.split()
        fields = fields[:next((i for i, field in enumerate(fields)
                               if field.startswith('$')), len(fields))]
        if not fields or line.startswith('*'):
            if section != 'BOUNDS':
                retained.append(line)
            continue
        header = not line[0].isspace() and fields[0] in SECTIONS
        if header:
            section = fields[0]
            if section == 'ENDATA':
                ended = True
                break
        elif section == 'COLUMNS':
            if len(fields) >= 3 and fields[1] == "'MARKER'":
                pass
            elif len(fields) in (3, 5):
                last_column = fields[0]
                columns.setdefault(last_column, None)
            elif len(fields) in (2, 4) and last_column is not None:
                columns.setdefault(last_column, None)
            else:
                raise ValueError('cannot identify source COLUMNS record')
        elif section == 'BOUNDS':
            if len(fields) not in (3, 4):
                raise ValueError('cannot identify source BOUNDS record')
            columns.setdefault(fields[2], None)
        if section != 'BOUNDS':
            retained.append(line)
    if not ended:
        raise ValueError('source model has no ENDATA')
    prefix = ''.join(retained)
    if prefix and not prefix.endswith('\n'):
        prefix += '\n'
    return prefix, list(columns)


def local_model(prefix, columns, bounds):
    if not isinstance(bounds, list):
        raise ValueError('bounds must be the complete list of source columns')
    box = {}
    for record in bounds:
        if not isinstance(record, dict) or set(record) != {'column', 'lower', 'upper'}:
            raise ValueError('bound records require exactly column, lower and upper')
        name = record['column']
        if not isinstance(name, str) or name in box:
            raise ValueError('invalid or duplicate bound column')
        box[name] = record
    missing, extra = set(columns) - set(box), set(box) - set(columns)
    if missing or extra:
        raise ValueError(f'trace box differs from source columns: missing={sorted(missing)}, extra={sorted(extra)}')
    lines = [prefix, 'BOUNDS\n']
    for name in columns:
        # Reset implicit nonnegativity and integer-marker upper bounds first.
        lines.append(f' FR replay {name}\n')
        for side, code in [('lower', 'LO'), ('upper', 'UP')]:
            value = box[name][side]
            if value is not None:
                lines.append(f' {code} replay {name} {finite_decimal(value)}\n')
    lines.append('ENDATA\n')
    return ''.join(lines)


def certificate_text(multipliers):
    if not isinstance(multipliers, dict):
        raise ValueError('multipliers must be a row-name map')
    lines = ['MILPLAB_FARKAS_1\n']
    for name, value in multipliers.items():
        if (not isinstance(name, str) or not name or any(ch.isspace() for ch in name)
                or name.startswith(('#', '$'))):
            raise ValueError('invalid certificate row name')
        lines.append(f'{name} {rational(value)}\n')
    return ''.join(lines)


def check(command, prefix, timeout):
    error = None
    try:
        result = subprocess.run(command, text=True, errors='replace', capture_output=True,
                                timeout=timeout)
        code, stdout, stderr = result.returncode, result.stdout, result.stderr
    except subprocess.TimeoutExpired as failure:
        code, stdout, stderr = None, failure.stdout or '', failure.stderr or ''
        stdout = stdout.decode(errors='replace') if isinstance(stdout, bytes) else stdout
        stderr = stderr.decode(errors='replace') if isinstance(stderr, bytes) else stderr
        error = 'checker timeout'
    pathlib.Path(str(prefix) + '.stdout.txt').write_text(stdout)
    pathlib.Path(str(prefix) + '.stderr.txt').write_text(stderr)
    try:
        parsed = strict_json(stdout)
        if not isinstance(parsed, dict):
            raise ValueError('checker JSON is not an object')
    except (ValueError, json.JSONDecodeError):
        parsed = None
        error = error or 'checker did not emit a valid JSON object'
    pathlib.Path(str(prefix) + '.json').write_text(json.dumps(
        {'exit_code': code, 'checker': parsed, 'error': error}, indent=2) + '\n')
    return code, parsed, error


def run(args, summary):
    code, inspection, error = check([str(args.checker), '--model', str(args.model), '--inspect'],
                                   args.out / 'original-inspect', args.timeout)
    if code != 0 or error or inspection.get('status') != 'supported':
        raise ValueError('original model is not supported by the independent checker')
    prefix, columns = model_without_bounds(args.model)
    summary['source_columns'] = len(columns)
    for line_number, line in enumerate(args.trace.read_text().splitlines(), 1):
        if not line.strip():
            continue
        try:
            event = strict_json(line)
            if not isinstance(event, dict):
                raise ValueError('trace record is not an object')
        except ValueError as error:
            summary['errors'].append(f'trace line {line_number}: {error}')
            continue
        if event.get('event') != 'lp_infeasibility_certificate':
            continue
        index = len(summary['proofs']) + 1
        directory = args.out / f'lp-{index:06d}'
        directory.mkdir()
        (directory / 'event.raw.json').write_text(line + '\n')
        (directory / 'event.json').write_text(json.dumps(event, indent=2) + '\n')
        record = {'trace_line': line_number, 'lp_solve': event.get('lp_solve'),
                  'directory': directory.name, 'independently_valid': False,
                  'claimed_margin': event.get('contradiction_margin'), 'checked_margin': None,
                  'margin_kind': event.get('margin_kind', 'exact'),
                  'claimed_margin_lower_bound': event.get('contradiction_margin_lower_bound'),
                  'margin_matches': False, 'margin_lower_bound_valid': False, 'error': None}
        summary['proofs'].append(record)
        try:
            model = local_model(prefix, columns, event.get('bounds'))
            certificate = certificate_text(event.get('multipliers'))
            kind = record['margin_kind']
            if kind not in ('exact', 'source_enclosure'):
                raise ValueError('unknown contradiction margin kind')
            if kind == 'source_enclosure':
                if 'contradiction_margin' in event:
                    raise ValueError('source enclosure cannot claim an exact margin')
                claimed = rational(event.get('contradiction_margin_lower_bound'))
            else:
                if 'contradiction_margin_lower_bound' in event:
                    raise ValueError('exact margin cannot carry an untyped lower bound')
                claimed = rational(event.get('contradiction_margin'))
            (directory / 'model.mps').write_text(model)
            (directory / 'certificate.farkas').write_text(certificate)
            record['model_sha256'] = digest(directory / 'model.mps')
            code, checked, error = check([str(args.checker), '--model', str(directory / 'model.mps'),
                '--farkas', str(directory / 'certificate.farkas')], directory / 'checker', args.timeout)
            record['checker_exit_code'], record['checker'] = code, checked
            if error:
                raise ValueError(error)
            record['checked_margin'] = checked.get('contradiction_margin')
            margin = (rational(record['checked_margin'])
                      if record['checked_margin'] is not None else None)
            record['margin_matches'] = kind == 'exact' and margin == claimed
            record['margin_lower_bound_valid'] = (kind == 'source_enclosure' and
                                                  margin is not None and 0 < claimed <= margin)
            valid_comparison = record['margin_matches'] or record['margin_lower_bound_valid']
            record['independently_valid'] = (code == 0 and checked.get('certificate_status') == 'valid'
                and checked.get('certifies_infeasibility') is True and margin is not None
                and margin > 0 and valid_comparison)
            if not record['independently_valid']:
                record['error'] = 'independent proof rejected or contradiction margin claim is invalid'
        except (ValueError, OSError) as error:
            record['error'] = str(error)
    if not summary['proofs'] and not args.allow_empty:
        summary['errors'].append('no LP infeasibility certificates found')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('model', 'trace', 'checker', 'out'):
        parser.add_argument('--' + name, required=True, type=pathlib.Path)
    parser.add_argument('--timeout', type=float, default=120)
    parser.add_argument('--allow-empty', action='store_true',
                        help='allow zero recorded proofs; never implies a tree proof')
    args = parser.parse_args()
    if not math.isfinite(args.timeout) or args.timeout <= 0:
        parser.error('--timeout must be finite and positive')
    args.model, args.trace, args.checker, args.out = (getattr(args, name).resolve()
                                                   for name in ('model', 'trace', 'checker', 'out'))
    try:
        args.out.mkdir(parents=True, exist_ok=False)
    except OSError as error:
        parser.error(f'--out must be a fresh directory: {error}')
    summary = {'schema': 1, 'scope': SCOPE, 'certifies_tree_infeasibility': False,
               'model': str(args.model), 'trace': str(args.trace), 'checker': str(args.checker),
               'all_lp_certificates_valid': False, 'proofs': [], 'errors': []}
    setup_failed = False
    try:
        summary.update(model_sha256=digest(args.model), trace_sha256=digest(args.trace),
                       checker_sha256=digest(args.checker), driver_sha256=digest(pathlib.Path(__file__)))
        run(args, summary)
    except (ValueError, OSError, EOFError) as error:
        summary['errors'].append(str(error))
        setup_failed = True
    summary['proof_count'] = len(summary['proofs'])
    summary['valid_proof_count'] = sum(p['independently_valid'] for p in summary['proofs'])
    summary['all_lp_certificates_valid'] = (not summary['errors'] and
        summary['valid_proof_count'] == summary['proof_count'])
    (args.out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary))
    return 2 if setup_failed else 0 if summary['all_lp_certificates_valid'] else 1


if __name__ == '__main__':
    sys.exit(main())
