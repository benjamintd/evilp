#!/usr/bin/env python3
"""Proved point intervals, source getter controls and large native cuts."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import tempfile


def checked(argv):
    run = subprocess.run([str(x) for x in argv], text=True, capture_output=True, timeout=120)
    if run.returncode:
        raise AssertionError((argv, run.returncode, run.stdout, run.stderr))
    return run.stdout


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    for name in ['interval-driver', 'getters', 'producer', 'checker']:
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    here = Path(__file__).resolve().parent
    intervals = json.loads(checked([sys.executable, here / 'source_interval_oracle.py',
                                   args.interval_driver.resolve()]))
    assert intervals['source_number_cases'] == 5015 and not intervals['false_point_claims']
    with tempfile.TemporaryDirectory(prefix='milplab-point-enclosures-') as name:
        root = Path(name)
        pipeline = root / 'pipeline'
        large = json.loads(checked([sys.executable, here / 'source_enclosure_pipeline.py',
                                    args.producer.resolve(), args.checker.resolve(), pipeline]))
        assert large['passed'] and large['independently_valid_cut_proofs'] == 15
        proof = root / 'masked.gmi'
        checked([args.getters.resolve(), pipeline / 'integer.mps', proof])
        result = json.loads(checked([args.checker.resolve(), '--model',
                                     pipeline / 'integer.mps', '--gmi', proof]))
        assert result['certifies_cut_validity']
        assert not result['certifies_original_lp_infeasibility']
        assert not result['certifies_original_mip_infeasibility']
        assert not result['certifies_full_tree']
    print(json.dumps({'passed': True, 'interval_cases': 5015,
                      'large_source_proofs': 15, 'masked_source_proofs': 1}))
