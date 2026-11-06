#!/usr/bin/env python3
"""Replay independent LP-box and cut-dependent MIP-box proofs, never a tree."""
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

SCOPE = ('individual LP boxes and original-MIP-in-recorded-box cut proofs; '
         'not a branch-and-bound tree proof')
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


def box_records(columns, bounds):
    if not isinstance(bounds, list):
        raise ValueError('bounds must be the complete list of source columns')
    box = {}
    for record in bounds:
        if not isinstance(record, dict) or set(record) != {'column', 'lower', 'upper'}:
            raise ValueError('bound records require exactly column, lower and upper')
        name = record['column']
        if not isinstance(name, str) or name in box:
            raise ValueError('invalid or duplicate bound column')
        lo = None if record['lower'] is None else rational(record['lower'])
        hi = None if record['upper'] is None else rational(record['upper'])
        if lo is not None and hi is not None and lo > hi:
            raise ValueError('contradictory trace box')
        box[name] = (lo, hi)
    missing, extra = set(columns) - set(box), set(box) - set(columns)
    if missing or extra:
        raise ValueError(f'trace box differs from source columns: missing={sorted(missing)}, extra={sorted(extra)}')
    return box


def local_model(prefix, columns, bounds):
    box = box_records(columns, bounds)
    lines = [prefix, 'BOUNDS\n']
    for name in columns:
        # Reset implicit nonnegativity and integer-marker upper bounds first.
        lines.append(f' FR replay {name}\n')
        for value, code in zip(box[name], ('LO', 'UP')):
            if value is not None:
                lines.append(f' {code} replay {name} {finite_decimal(str(value))}\n')
    lines.append('ENDATA\n')
    return ''.join(lines)


def valid_name(name):
    return (isinstance(name, str) and bool(name) and len(name) <= 4096
            and not any(ch.isspace() or ch == '\x00' for ch in name)
            and not name.startswith(('#', '$')))


def source_row_names(prefix):
    result, section = set(), None
    for line in prefix.splitlines():
        fields = line.split()
        if not fields or line.startswith('*'):
            continue
        if not line[0].isspace() and fields[0] in SECTIONS:
            section = fields[0]
        elif section == 'ROWS':
            result.add(fields[1])
    return result


def proof_number(text):
    # Only called on fields of a certificate already accepted independently.
    # The checker also supports exact decimal/scientific strings, including D.
    if not isinstance(text, str) or len(text) > 4096:
        raise ValueError('invalid validated proof number')
    try:
        return fractions.Fraction(text.replace('D', 'E').replace('d', 'e'))
    except (ValueError, ZeroDivisionError) as error:
        raise ValueError('invalid validated proof number') from error


def checked_cut_fields(text, columns):
    box, coefficient, rhs = {}, {}, None
    for line in text.splitlines():
        fields = line.split()
        if not fields or fields[0].startswith('#'):
            continue
        if fields[0] == 'BOX':
            if len(fields) != 4 or fields[1] in box:
                raise ValueError('invalid/duplicate validated BOX record')
            box[fields[1]] = tuple(None if t == '*' else proof_number(t) for t in fields[2:])
        elif fields[0] == 'PUB':
            if len(fields) != 3 or fields[1] in coefficient:
                raise ValueError('invalid/duplicate validated PUB record')
            coefficient[fields[1]] = proof_number(fields[2])
        elif fields[0] == 'PUB_RHS':
            if len(fields) != 2 or rhs is not None:
                raise ValueError('invalid/duplicate validated PUB_RHS record')
            rhs = proof_number(fields[1])
    if set(box) != set(columns) or set(coefficient) != set(columns) or rhs is None:
        raise ValueError('validated proof BOX/PUB differs from complete source columns')
    for q in [rhs, *coefficient.values()]:
        try:
            value = float(q)
        except OverflowError as error:
            raise ValueError('nonfinite published cut value') from error
        if not math.isfinite(value) or fractions.Fraction(value) != q:
            raise ValueError('published cut value is not an exact finite binary64')
    return {'box': box, 'coefficient': coefficient, 'rhs': rhs}


def contained_box(local, premise):
    for name, (lo, hi) in local.items():
        outer_lo, outer_hi = premise[name]
        if outer_lo is not None and (lo is None or lo < outer_lo):
            return False
        if outer_hi is not None and (hi is None or hi > outer_hi):
            return False
    return True


def augmented_prefix(prefix, columns, cuts):
    """Insert cuts while preserving source decimals and integer-marker modes."""
    if not cuts:
        return prefix
    lines = prefix.splitlines(keepends=True)
    rows = source_row_names(prefix)
    for cut in cuts:
        if not valid_name(cut['cut_id']) or cut['cut_id'] in rows:
            raise ValueError('cut name collides with an original or added source row')
        rows.add(cut['cut_id'])
    last, headers, section, previous, rhs_vector = {}, {}, None, None, None
    for i, line in enumerate(lines):
        fields = line.split()
        fields = fields[:next((k for k, token in enumerate(fields)
                               if token.startswith('$')), len(fields))]
        if not fields or line.startswith('*'):
            continue
        if not line[0].isspace() and fields[0] in SECTIONS:
            section = fields[0]
            headers[section] = i
        elif section == 'COLUMNS':
            if len(fields) >= 3 and fields[1] == "'MARKER'":
                continue
            if len(fields) in (3, 5):
                previous = fields[0]
            elif len(fields) not in (2, 4) or previous is None:
                raise ValueError('cannot locate original column for exact cut insertion')
            last[previous] = i
        elif section == 'RHS' and rhs_vector is None:
            rhs_vector = fields[0]
    if 'COLUMNS' not in headers:
        raise ValueError('source has no COLUMNS section')
    column_end = min([i for i in headers.values() if i > headers['COLUMNS']] + [len(lines)])
    before, after = {}, {}
    before[headers['COLUMNS']] = [f" G {cut['cut_id']}\n" for cut in cuts]
    extras = []
    for name in columns:
        entries = [f" {name} {cut['cut_id']} {finite_decimal(str(cut['coefficient'][name]))}\n"
                   for cut in cuts if cut['coefficient'][name]]
        if name in last:
            after.setdefault(last[name], []).extend(entries)
        else:
            extras.extend(entries)  # Bound-only columns, outside integer markers.
    before.setdefault(column_end, []).extend(extras)
    rhs_entries = [f" {rhs_vector or 'replaycuts'} {cut['cut_id']} {finite_decimal(str(cut['rhs']))}\n"
                   for cut in cuts]
    if 'RHS' in headers:
        rhs_end = min([i for i in headers.values() if i > headers['RHS']] + [len(lines)])
        before.setdefault(rhs_end, []).extend(rhs_entries)
    else:
        before.setdefault(column_end, []).extend(['RHS\n', *rhs_entries])
    output = []
    for i in range(len(lines) + 1):
        output.extend(before.get(i, []))
        if i < len(lines):
            output.append(lines[i])
            output.extend(after.get(i, []))
    return ''.join(output)


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


def margin_claim(event, record):
    kind = record['margin_kind']
    if kind not in ('exact', 'source_enclosure'):
        raise ValueError('unknown contradiction margin kind')
    if kind == 'source_enclosure':
        if 'contradiction_margin' in event:
            raise ValueError('source enclosure cannot claim an exact margin')
        return rational(event.get('contradiction_margin_lower_bound'))
    if 'contradiction_margin_lower_bound' in event:
        raise ValueError('exact margin cannot carry an untyped lower bound')
    return rational(event.get('contradiction_margin'))


def run(args, summary):
    code, inspection, error = check([str(args.checker), '--model', str(args.model), '--inspect'],
                                   args.out / 'original-inspect', args.timeout)
    if code != 0 or error or inspection.get('status') != 'supported':
        raise ValueError('original model is not supported by the independent checker')
    prefix, columns = model_without_bounds(args.model)
    original_rows = source_row_names(prefix)
    summary['source_columns'] = len(columns)
    cuts, declared_ids, seen_events = {}, set(), set()
    types = {'lp_infeasibility_certificate', 'root_gmi_cut',
             'mip_in_box_infeasibility_certificate'}
    for line_number, line in enumerate(args.trace.read_text().splitlines(), 1):
        if not line.strip():
            continue
        try:
            event = strict_json(line)
            if not isinstance(event, dict):
                raise ValueError('trace record is not an object')
            event_type = event.get('event')
            if event_type is not None and not isinstance(event_type, str):
                raise ValueError('trace event type must be a string')
            if event_type not in types:
                if isinstance(event_type, str) and ('certificate' in event_type or
                                                   event_type.startswith('root_gmi')):
                    raise ValueError(f'unknown proof event {event_type!r}')
                continue
            identity = json.dumps(event, sort_keys=True, separators=(',', ':'))
            if identity in seen_events:
                raise ValueError('duplicate proof event')
            seen_events.add(identity)
        except ValueError as error:
            summary['errors'].append(f'trace line {line_number}: {error}')
            continue
        if event_type == 'root_gmi_cut':
            directory = args.out / f"cut-{len(summary['cuts']) + 1:06d}"
            directory.mkdir()
            (directory / 'event.raw.json').write_text(line + '\n')
            record = {'trace_line': line_number, 'cut_id': event.get('cut_id'),
                      'directory': directory.name, 'independently_valid': False,
                      'certificate_scope': 'original_mip_in_recorded_box', 'error': None}
            summary['cuts'].append(record)
            try:
                if set(event) != {'event', 'cut_id', 'scope', 'certificate'}:
                    raise ValueError('root cut event has unexpected or missing fields')
                cut_id = event.get('cut_id')
                if not valid_name(cut_id) or cut_id in original_rows or cut_id in declared_ids:
                    raise ValueError('invalid, colliding or duplicate cut ID')
                declared_ids.add(cut_id)
                if event.get('scope') != 'ORIGINAL_MIP_IN_BOX':
                    raise ValueError('root cut event has an unsupported scope')
                text = event.get('certificate')
                if not isinstance(text, str) or len(text.encode()) > (1 << 24):
                    raise ValueError('root cut certificate must be bounded text')
                path = directory / 'certificate.gmi'
                path.write_text(text)
                code, checked, error = check([str(args.checker), '--model', str(args.model),
                    '--gmi', str(path)], directory / 'checker', args.timeout)
                record['checker_exit_code'], record['checker'] = code, checked
                record['certificate_sha256'] = digest(path)
                if error:
                    raise ValueError(error)
                if (code != 0 or checked.get('status') != 'valid' or
                    checked.get('certifies_cut_validity') is not True or
                    checked.get('certificate_scope') != 'original_mip_in_recorded_box' or
                    checked.get('certifies_original_lp_infeasibility') is not False or
                    checked.get('certifies_original_mip_infeasibility') is not False or
                    checked.get('certifies_full_tree') is not False or
                    checked.get('model_sha256') != summary['model_sha256'] or
                    checked.get('certificate_sha256') != record['certificate_sha256']):
                    raise ValueError('independent root cut derivation rejected or scope/hash mismatch')
                cut = checked_cut_fields(text, columns)
                cut['cut_id'] = cut_id
                cuts[cut_id] = cut
                record['independently_valid'] = True
            except (ValueError, OSError) as error:
                record['error'] = str(error)
            continue
        index = len(summary['proofs']) + 1
        directory = args.out / f'lp-{index:06d}'
        directory.mkdir()
        (directory / 'event.raw.json').write_text(line + '\n')
        (directory / 'event.json').write_text(json.dumps(event, indent=2) + '\n')
        mip = event_type == 'mip_in_box_infeasibility_certificate'
        record = {'trace_line': line_number, 'lp_solve': event.get('lp_solve'),
                  'directory': directory.name, 'independently_valid': False,
                  'claimed_margin': event.get('contradiction_margin'), 'checked_margin': None,
                  'margin_kind': event.get('margin_kind', 'exact'),
                  'claimed_margin_lower_bound': event.get('contradiction_margin_lower_bound'),
                  'margin_matches': False, 'margin_lower_bound_valid': False, 'error': None,
                  'proof_type': event_type,
                  'certificate_scope': ('original_mip_in_recorded_box' if mip else 'lp_in_recorded_box'),
                  'certifies_original_lp_infeasibility': False,
                  'certifies_original_mip_infeasibility': False, 'certifies_full_tree': False}
        summary['proofs'].append(record)
        try:
            selected = []
            if mip:
                ids = event.get('cut_ids')
                if not isinstance(ids, list) or not ids or not all(valid_name(i) for i in ids):
                    raise ValueError('cut-dependent proof requires a nonempty cut ID list')
                if len(set(ids)) != len(ids):
                    raise ValueError('duplicate cut IDs in cut-dependent proof')
                if any(i not in cuts for i in ids):
                    raise ValueError('unknown or independently rejected cut ID')
                selected = [cuts[i] for i in ids]
                local_box = box_records(columns, event.get('bounds'))
                for cut in selected:
                    if not contained_box(local_box, cut['box']):
                        raise ValueError(f"proof box is outside cut {cut['cut_id']!r} validity box")
                record['cut_ids'] = ids
            elif event.get('cut_ids'):
                raise ValueError('original LP proof cannot carry cut IDs')
            model = local_model(augmented_prefix(prefix, columns, selected), columns,
                                event.get('bounds'))
            certificate = certificate_text(event.get('multipliers'))
            claimed = margin_claim(event, record)
            kind = record['margin_kind']
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
    if not summary['proofs'] and not summary['cuts'] and not args.allow_empty:
        summary['errors'].append('no LP infeasibility or root cut certificates found')


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
               'all_lp_certificates_valid': False, 'all_certificates_valid': False,
               'certifies_original_lp_infeasibility': False,
               'certifies_original_mip_infeasibility': False,
               'proofs': [], 'cuts': [], 'errors': []}
    setup_failed = False
    try:
        summary.update(model_sha256=digest(args.model), trace_sha256=digest(args.trace),
                       checker_sha256=digest(args.checker), driver_sha256=digest(pathlib.Path(__file__)))
        run(args, summary)
    except (ValueError, OSError, EOFError) as error:
        summary['errors'].append(str(error))
        setup_failed = True
    for name in ('model', 'trace', 'checker'):
        try:
            if digest(getattr(args, name)) != summary.get(name + '_sha256'):
                summary['errors'].append(f'{name} input changed during replay')
        except OSError as error:
            summary['errors'].append(f'{name} input unavailable after replay: {error}')
    summary['proof_count'] = len(summary['proofs'])
    summary['valid_proof_count'] = sum(p['independently_valid'] for p in summary['proofs'])
    summary['cut_count'] = len(summary['cuts'])
    summary['valid_cut_count'] = sum(p['independently_valid'] for p in summary['cuts'])
    summary['mip_in_box_proof_count'] = sum(p['proof_type'] ==
        'mip_in_box_infeasibility_certificate' for p in summary['proofs'])
    summary['all_certificates_valid'] = (not summary['errors'] and
        summary['valid_proof_count'] == summary['proof_count'] and
        summary['valid_cut_count'] == summary['cut_count'])
    # Keep the historical all-LP field conservative when cuts were involved.
    summary['all_lp_certificates_valid'] = (summary['all_certificates_valid'] and
        not summary['mip_in_box_proof_count'] and not summary['cut_count'])
    (args.out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary))
    return 2 if setup_failed else 0 if summary['all_certificates_valid'] else 1


if __name__ == '__main__':
    sys.exit(main())
