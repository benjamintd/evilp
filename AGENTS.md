# Working on MILP Lab

- Keep the numerical engine, independent checker, and process runner in C17.
  Python is for setup, test oracles and reporting, not optimization inside a run.
- The checker must not link an optimization solver or reuse solver feasibility
  APIs. Preserve rational arithmetic and validate original, unpresolved models.
- Dataset and solver lock files identify experiments. Version intentional updates
  and create a fresh output directory; never rewrite historical results.
- A reported bound/status, a checked feasible witness, reference-objective
  agreement, and an independent proof are distinct. Preserve that distinction.
- Build and run the contract tests after relevant C changes. Numerical changes
  require exact/adversarial cases; run sanitizers for memory/control-flow changes.
- Preserve minimized regressions and record unsupported formats explicitly.
  Do not weaken tolerances or skip failing cases to make a benchmark look better.
- Do not specialize algorithms using benchmark filenames, hashes or stored answers.
- Keep fetched datasets, binaries and build trees outside the tracked source set.
  Imported artifacts retain their own licenses; do not relabel them as MIT.
