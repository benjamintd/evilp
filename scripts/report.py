#!/usr/bin/env python3
"""Summarize C runner artifacts. Claimed optimality is never certified here."""
import argparse,collections,fractions,html,json,pathlib,statistics
F=fractions.Fraction
def report(root,runs_file='runs.jsonl',prefix='summary'):
    root=pathlib.Path(root);rows=[json.loads(x) for x in (root/runs_file).read_text().splitlines() if x.strip()]
    groups=collections.defaultdict(list);details=[];contradictions=[]
    for r in rows:
        c=r['checker'];valid=c.get('checked_feasible',False) and c.get('status')=='valid'
        objective=F(c['objective']) if valid else None
        ref=F(r['reference_objective']) if r['reference_objective']!='-' else None
        tol=F('0.0001')+F('0.000000001')*abs(ref or 0)
        match=valid and ref is not None and r['reference_status']=='=opt=' and abs(objective-ref)<=tol
        better=valid and ref is not None and r['reference_status']=='=opt=' and not match and (objective<ref if c.get('objective_sense')=='min' else objective>ref)
        contradiction=(r['claimed_status']=='infeasible' and r['reference_status'] in ('=opt=','=best=')) or (valid and r['reference_status']=='=inf=') or (r['claimed_status']=='optimal' and valid and ref is not None and r['reference_status']=='=opt=' and not match) or better
        d={**r,'validated_objective_decimal':float(objective) if objective is not None else None,'reached_reference':bool(match),'reference_contradiction':contradiction}
        groups[r['solver']].append(d);details.append(d)
        if contradiction:contradictions.append({'instance':r['instance'],'solver':r['solver'],'claimed_status':r['claimed_status']})
    summary={}
    for name,runs in sorted(groups.items()):
        summary[name]={'attempts':len(runs),'claimed_statuses':dict(collections.Counter(x['claimed_status'] for x in runs)),
          'checked_feasible':sum(x['checker'].get('checked_feasible',False) and x['checker'].get('status')=='valid' for x in runs),
          'invalid_candidates':sum(x['checker'].get('status')=='invalid' for x in runs),
          'checker_errors':sum(x['checker'].get('status') in ('error','unsupported') for x in runs),
          'certified_lp_bounds':sum(x['checker'].get('bound_status')=='valid' for x in runs),
          'certified_optima':sum(x['checker'].get('certifies_optimality',False) for x in runs),
          'reached_reference':sum(x['reached_reference'] for x in runs),
          'reference_contradictions':sum(x['reference_contradiction'] for x in runs),
          'total_solver_wall_s':sum(x['wall_s'] for x in runs),'total_checker_wall_s':sum(x['checker_wall_s'] for x in runs),
          'median_solver_wall_s':statistics.median(x['wall_s'] for x in runs)}
    data={'schema':1,'kind':'short-budget pipeline baseline; not a solver ranking','reference_match_policy':'absolute 1e-4 + relative 1e-9 against pinned =opt= reference; not an optimality proof',
          'summary':summary,'contradictions':contradictions,'attempts':len(rows),'instances':len({r['instance'] for r in rows}),'environment':json.loads((root/'environment.json').read_text())}
    data['run_records']=runs_file
    data['relaxation']=bool(data['environment'].get('relaxation',False))
    if data['relaxation']:data['kind']='LP relaxation correctness baseline; not a solver ranking'
    if (root/'host.json').exists():data['host']=json.loads((root/'host.json').read_text())
    if any('revalidation' in r for r in rows):data['revalidation']={'checker_sha256':sorted({r['checker_sha256'] for r in rows}),'driver_sha256':sorted({r['revalidation']['driver_sha256'] for r in rows if 'revalidation' in r}),'solver_times_unchanged':True}
    (root/(prefix+'.json')).write_text(json.dumps(data,indent=2)+'\n')
    E=lambda x:html.escape(str(x))
    parts=['<!doctype html><meta charset="utf-8"><title>MILP baseline report</title><style>body{font:16px system-ui;max-width:1200px;margin:40px auto;padding:0 24px;color:#102a43}table{border-collapse:collapse;width:100%;font-size:13px;margin:24px 0}th,td{text-align:left;padding:9px;border-bottom:1px solid #ccd}th{background:#102a43;color:white}tr:nth-child(even){background:#eff5f7}code{overflow-wrap:anywhere}.bad{color:#b22}a{color:#087f8c}</style><h1>MILP baseline report</h1>']
    env=data['environment'];parts.append(f'<p>{data["instances"]} instances, {len(rows)} runs. Requested one thread; seed {env["seed"]}; internal solve limit {env["time_limit_s"]} seconds. Process wall time includes startup, input parsing and output; checking is separate.</p>')
    if data['relaxation']:parts.append('<p><strong>Continuous LP relaxations: integrality is explicitly ignored and suite reference objectives must be LP values.</strong></p>')
    parts.append('<p><strong>This short-budget run validates the laboratory. It does not establish a performance ranking.</strong> A checked feasible solution satisfies the model under the recorded tolerances. Reference agreement is a comparison with versioned metadata, not a proof. Only an exact matching primal/dual certificate can independently certify optimality.</p>')
    parts.append('<p>The checker requires absolute row/bound tolerance 1e-4, integer tolerance 1e-5, and claimed-objective consistency 1e-4. A rejection can involve feasibility or objective consistency. With large coefficients, a point accepted under solver scaling/tolerances can fail this absolute policy. Such a failure alone does not establish a solver bug.</p>')
    if 'revalidation' in data:parts.append('<p>Independent validation was rerun from retained solution files. Original solver timings and raw run records are unchanged. Each row links to the newer validation record.</p>')
    parts.append('<table><tr><th>Solver</th><th>Attempts</th><th>Passed checker</th><th>Reached reference</th><th>Certified LP bounds</th><th>Certified optima</th><th>Invalid candidates</th><th>Checker errors</th><th>Claimed statuses</th><th>Solver wall seconds</th></tr>')
    for name,s in summary.items():parts.append('<tr>'+''.join(f'<td>{E(x)}</td>' for x in [name,s['attempts'],s['checked_feasible'],s['reached_reference'],s['certified_lp_bounds'],s['certified_optima'],s['invalid_candidates'],s['checker_errors'],s['claimed_statuses'],round(s['total_solver_wall_s'],2)])+'</tr>')
    parts.append('</table><p>Reference matching: absolute 1e-4 + relative 1e-9, for pinned =opt= records only. No comparisons use a solver unvalidated dual bound as truth. All failures and missing candidates remain in the table.</p>')
    parts.append('<h2>Individual runs</h2><p>Run status includes process outcomes and explicit terminal solution-header refinements. Solver report preserves the parsed log status and any separate solution status. Raw claims, exits and signals remain in each record.</p><table><tr><th>Instance</th><th>Solver</th><th>Run status</th><th>Solver report</th><th>Checker</th><th>Checked objective</th><th>Certified LP bound</th><th>Reference reached</th><th>Wall seconds</th><th>Evidence</th></tr>')
    for r in details:
        relative_run=f'{r["instance"]}/{r["solver"]}/seed-{r["seed"]}'
        bound=r['checker'].get('certified_bound') if r['checker'].get('bound_status')=='valid' else None
        reported=r.get('reported_solver_status')
        if r.get('reported_solution_status'):reported=str(reported)+' / solution: '+r['reported_solution_status']
        values=[r['instance'],r['solver'],r['claimed_status'],reported,r['checker']['status'],r['validated_objective_decimal'],bound,r['reached_reference'],round(r['wall_s'],4)]
        validation=f' / <a href="{E(relative_run)}/{E(r["revalidation"].get("record","checker-revalidated.json"))}">validation</a>' if 'revalidation' in r else ''
        if r.get('dual_file'):validation+=f' / <a href="{E(relative_run)}/dual.sol">dual</a>'
        parts.append('<tr>'+''.join(f'<td>{E(x)}</td>' for x in values)+f'<td><a href="{E(relative_run)}/run.json">record</a> / <a href="{E(relative_run)}/solver.log">log</a>{validation}</td></tr>')
    parts.append('</table><h2>Original timed evaluator</h2><pre><code>'+E(json.dumps(env,indent=2))+'</code></pre>')
    if 'revalidation' in data:parts.append('<h2>Final independent validation</h2><pre><code>'+E(json.dumps(data['revalidation'],indent=2))+'</code></pre>')
    if 'host' in data:parts.append('<h2>Host and measurement conditions</h2><pre><code>'+E(json.dumps(data['host'],indent=2))+'</code></pre>')
    (root/(prefix+'.html')).write_text('\n'.join(parts)+'\n')
    print(json.dumps(data['summary'],indent=2))
    return data
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('directory');ap.add_argument('--runs',default='runs.jsonl');ap.add_argument('--prefix',default='summary');args=ap.parse_args();report(args.directory,args.runs,args.prefix)
