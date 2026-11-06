"""Replay each actual exported tableau with independent Fraction algebra."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from native_tableau_oracle import verify

with tempfile.TemporaryDirectory(prefix="native-tableau-contract-") as directory:
    output = Path(directory) / "batches.jsonl"
    run = subprocess.run([sys.argv[1], str(output)], text=True, capture_output=True, timeout=120)
    if run.returncode:
        raise AssertionError(run.stdout + run.stderr)
    rows = 0
    for line in output.read_text().splitlines():
        rows += verify(json.loads(line))
    summary = json.loads(run.stdout)
    assert rows == summary["checked_exports"] == 34
    print(json.dumps({"passed": True, "independent_exact_rows": rows,
                      "allocation_rejections": summary["allocation_rejections"],
                      "stale_rejections": summary["stale_rejections"]}))
