#!/usr/bin/env python3
"""Optional native LP integration: exact-known fixtures, four real engines.

This is a correctness crosscheck, not a speed comparison. Models, identities,
native/lab run records and a summary remain in a fresh --out directory.
"""
import argparse
import collections
import csv
import datetime
import fractions
import hashlib
import importlib.util
import json
import pathlib
import platform
import shutil
import subprocess

ROOT = pathlib.Path(__file__).resolve().parents[1]
F = fractions.Fraction


def digest(path):
    h = hashlib.sha256()
    with pathlib.Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def absolute(path):
    path = pathlib.Path(path)
    return (path if path.is_absolute() else ROOT / path).resolve()


def helpers():
    # Import only the independently written model writer and rational oracle.
    spec = importlib.util.spec_from_file_location("milplab_lp_test_helpers",
                                                ROOT / "tests/test_lp_solver.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def cases():
    return [
        ("min_offset", {"cost": [2, 3], "bounds": [(0, 4), (0, 4)],
                        "rows": [([1, 1], F(3, 2), None)], "offset": F(5, 4)}, F(17, 4)),
        ("max_offset_fixed", {"cost": [3, 2, F(1, 10)],
                              "bounds": [(0, 1), (0, 1), (F(3, 10), F(3, 10))],
                              "rows": [([2, 1, 0], None, 2)],
                              "maximize": True, "offset": 5}, F(853, 100)),
        ("nonterminating", {"cost": [1, 2], "bounds": [(0, 3), (0, 3)],
                            "rows": [([3, 0], 1, None), ([0, 7], 2, None)]}, F(19, 21)),
        ("ranged_free", {"cost": [1, -1], "bounds": [(None, None), (0, 2)],
                         "rows": [([1, 1], 1, 3), ([1, -1], -1, None)]}, F(-1)),
        ("upper_only", {"cost": [-1], "bounds": [(None, -2)], "rows": []}, F(2)),
        ("fixed_redundant", {"cost": [1, 2, 3],
                             "bounds": [(0, 3), (0, 3), (F(-1, 2), F(-1, 2))],
                             "rows": [([1, 1, 0], 2, 2), ([2, 2, 0], 4, 4),
                                      ([0, 0, 0], 0, 0)]}, F(1, 2)),
        ("negative_rhs", {"cost": [3, F(-1, 2)], "bounds": [(-1, 2), (-4, 2)],
                          "rows": [([0, -1], 3, None), ([0, F(3, 2)], F(-9, 2), None),
                                   ([F(3, 2), F(3, 2)], None, F(-9, 4))],
                          "offset": F(7, 2)}, F(2)),
        ("max_ranged", {"cost": [1, 3], "bounds": [(0, 2), (0, 2)],
                        "rows": [([1, 2], F(1, 2), F(5, 2))],
                        "maximize": True, "offset": F(-1, 4)}, F(7, 2)),
    ]


def known_highs_row_name_warning(record, log):
    """Accept one output-only warning while retaining its nonzero process exit."""
    if not (record["solver"] == "highs" and record["exit_code"] == 1
            and record["signal"] == 0 and not record["watchdog_timeout"]
            and record["claimed_status"] == "process_error"
            and record.get("reported_solver_status") == "optimal"):
        return None
    warnings = [" ".join(line.split()) for line in log.splitlines() if "WARNING" in line.upper()]
    expected = ('WARNING: Row names are not present, or contain duplicates: using names '
                'with prefix "r", beginning with suffix 0')
    if warnings != [expected] or any("ERROR" in line for line in log.splitlines()):
        return None
    return expected


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", required=True, help="fresh artifact directory")
    ap.add_argument("--runner", default="build/milpbench")
    ap.add_argument("--checker", default="build/milpcheck")
    ap.add_argument("--solver", default="build/milpsolve")
    ap.add_argument("--native-solvers", default="config/solvers.tsv")
    args = ap.parse_args()
    out = absolute(args.out)
    if out.exists():
        ap.error("--out must not exist; historical artifacts are never replaced")
    runner, checker, solver = map(absolute, (args.runner, args.checker, args.solver))
    native_config = absolute(args.native_solvers)
    for path in (runner, checker, solver, native_config):
        if not path.is_file():
            ap.error(f"required input is missing: {path}")
    with native_config.open() as f:
        native_solvers = list(csv.DictReader(f, delimiter="\t"))
    if {x["solver"] for x in native_solvers} != {"highs", "scip", "cbc"} or len(native_solvers) != 3:
        ap.error("native configuration must contain exactly HiGHS, SCIP and Cbc")
    helper = helpers()
    fixtures = cases()
    for name, model, expected in fixtures:
        if helper.oracle(model) != expected:
            raise AssertionError(f"incorrect exact reference for {name}")
    out.mkdir(parents=True)
    # Keep native evaluator binaries fixed while developers may rebuild the
    # working tree. Upstream solvers remain identified by their locked hashes.
    evaluator = out / "evaluator"
    evaluator.mkdir()
    origins = {"runner": str(runner), "checker": str(checker), "lab": str(solver)}
    for source, name in ((runner, "milpbench"), (checker, "milpcheck"), (solver, "milpsolve")):
        shutil.copy2(source, evaluator / name)
    runner, checker, solver = (evaluator / name for name in ("milpbench", "milpcheck", "milpsolve"))
    model_dir = out / "models"
    model_dir.mkdir()
    suite = out / "suite.tsv"
    expected_models = {}
    with suite.open("x") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(["id", "path", "sha256", "reference_status", "reference_objective"])
        for name, model, expected in fixtures:
            path = model_dir / (name + ".mps")
            path.write_text(helper.mps(model))
            sha = digest(path)
            writer.writerow([name, path, sha, "=opt=", str(expected)])
            expected_models[name] = {"objective": str(expected), "maximize": bool(model.get("maximize")),
                                     "sha256": sha}
    lab_config = out / "lab-solvers.tsv"
    with lab_config.open("x") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(["solver", "binary", "sha256", "version"])
        writer.writerow(["lab", solver, digest(solver), "milestone-2A-dense-reference"])
    # Copy the exact native configuration into the result package. Its relative
    # executable paths continue to resolve from the documented project root.
    pinned_native = out / "native-solvers.tsv"
    pinned_native.write_bytes(native_config.read_bytes())
    commands = {}
    failures = []
    for group, config in (("native", pinned_native), ("lab", lab_config)):
        command = [str(runner), "--suite", str(suite), "--out", str(out / group),
                   "--solvers", str(config), "--checker", str(checker), "--limit", "5", "--seed", "1"]
        if group == "lab":
            command.append("--relax")
        proc = subprocess.run(command, cwd=ROOT, capture_output=True, text=True, timeout=300)
        commands[group] = {"argv": command, "exit_code": proc.returncode}
        (out / (group + "-driver.log")).write_text(proc.stdout + proc.stderr)
        if proc.returncode:
            failures.append(f"{group} harness exit {proc.returncode}")
    records = []
    for group in commands:
        runs = out / group / "runs.jsonl"
        if not runs.exists():
            failures.append(f"{group} has no run records")
            continue
        for line in runs.read_text().splitlines():
            r = json.loads(line)
            name, engine = r["instance"], r["solver"]
            check = r["checker"]
            prefix = f"{name}/{engine}"
            expected = F(expected_models[name]["objective"])
            log = (out / group / name / engine / "seed-1/solver.log").read_text()
            process_warning = known_highs_row_name_warning(r, log)
            process_clean = r["exit_code"] == 0 and not r["signal"] and not r["watchdog_timeout"]
            if not process_clean and not process_warning:
                failures.append(f"{prefix}: process exit {r['exit_code']}, signal {r['signal']}, watchdog {r['watchdog_timeout']}")
            if r["claimed_status"] != "optimal" and not process_warning:
                failures.append(f"{prefix}: solver status {r['claimed_status']}")
            feasible = check.get("checked_feasible") and check.get("status") == "valid"
            if not feasible:
                failures.append(f"{prefix}: checker rejects primal witness ({check.get('status')})")
            actual = F(check["objective"]) if check.get("objective") is not None else None
            error = abs(actual - expected) if actual is not None else None
            tolerance = F(1, 10**7) + F(1, 10**9) * abs(expected)
            objective_matches = error is not None and error <= tolerance
            if not objective_matches:
                failures.append(f"{prefix}: objective {actual}, exact reference {expected}")
            if engine == "lab":
                status = check.get("bound_status")
                bound = check.get("certified_bound")
                direction = "upper" if expected_models[name]["maximize"] else "lower"
                if status == "valid":
                    if bound is None or check.get("bound_direction") != direction:
                        failures.append(f"{prefix}: malformed certified bound")
                    elif (F(bound) < expected if direction == "upper" else F(bound) > expected):
                        failures.append(f"{prefix}: certified {direction} bound {bound} is unsafe")
                elif status != "unavailable" or bound is not None:
                    failures.append(f"{prefix}: unexpected bound status {status}")
            records.append({"instance": name, "solver": engine, "claimed_status": r["claimed_status"],
                            "reported_solver_status": r.get("reported_solver_status"),
                            "passed_checker": bool(feasible), "objective": str(actual) if actual is not None else None,
                            "objective_matches": objective_matches, "process_clean": process_clean,
                            "accepted_process_warning": process_warning,
                            "exit_code": r["exit_code"], "signal": r["signal"], "watchdog_timeout": r["watchdog_timeout"],
                            "exact_reference": str(expected), "objective_error": str(error) if error is not None else None,
                            "bound_status": check.get("bound_status"),
                            "certified_bound": check.get("certified_bound"),
                            "bound_direction": check.get("bound_direction"),
                            "model_sha256": r["model_sha256"], "solver_sha256": r["solver_sha256"],
                            "checker_sha256": r["checker_sha256"], "runner_sha256": r["runner_sha256"],
                            "record": f"{group}/{name}/{engine}/seed-1/run.json"})
    expected_pairs = {(name, engine) for name, _, _ in fixtures for engine in ("highs", "scip", "cbc", "lab")}
    actual_pairs = [(r["instance"], r["solver"]) for r in records]
    if len(actual_pairs) != len(set(actual_pairs)) or set(actual_pairs) != expected_pairs:
        failures.append("missing, duplicate or unexpected solver–model pairs")
    revision = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    dirty = bool(subprocess.run(["git", "status", "--porcelain"], cwd=ROOT, capture_output=True, text=True).stdout.strip())
    summary = {"schema": 1, "kind": "native LP correctness crosscheck; not a speed ranking",
               "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
               "instances": len(fixtures), "expected_runs": len(expected_pairs), "actual_runs": len(records),
               "passed": not failures, "failures": failures,
               "objective_tolerance": "absolute 1e-7 + relative 1e-9 against exact reference",
               "bound_policy": "certified lower/upper bounds must be on the exact safe side; unavailable permitted",
               "process_policy": "zero exit required except HiGHS exit1 with terminal optimal status and exactly the known output-only row-name warning; all raw exits/statuses retained",
               "process_warnings": [{"instance": r["instance"], "solver": r["solver"],
                                     "exit_code": r["exit_code"], "claimed_status": r["claimed_status"],
                                     "reported_solver_status": r["reported_solver_status"],
                                     "message": r["accepted_process_warning"]}
                                    for r in records if r["accepted_process_warning"]],
               "native_solvers": native_solvers, "reference_models": expected_models,
               "identities": {"git_revision": revision, "working_tree_dirty": dirty,
                              "binary_origins": origins,
                              "script_sha256": digest(__file__), "fixture_helper_sha256": digest(ROOT / "tests/test_lp_solver.py"),
                              "runner_sha256": digest(runner), "checker_sha256": digest(checker),
                              "lab_sha256": digest(solver), "native_config_sha256": digest(pinned_native),
                              "lab_config_sha256": digest(lab_config), "suite_sha256": digest(suite)},
               "platform": {"system": platform.system(), "machine": platform.machine(), "release": platform.release()},
               "counts": {"passed_checker": sum(r["passed_checker"] for r in records),
                          "matched_exact_reference": sum(r["objective_matches"] for r in records),
                          "claimed_optimal": sum(r["claimed_status"] == "optimal" for r in records),
                          "clean_process_runs": sum(r["process_clean"] for r in records),
                          "accepted_process_warnings": sum(bool(r["accepted_process_warning"]) for r in records),
                          "lab_bound_statuses": dict(collections.Counter(r["bound_status"] for r in records if r["solver"] == "lab"))},
               "commands": commands, "records": records}
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: summary[k] for k in ("passed", "instances", "actual_runs", "counts", "failures")}, indent=2))
    return 0 if summary["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
