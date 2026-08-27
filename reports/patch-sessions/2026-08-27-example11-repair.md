# Example 11 gradient-search repair

- Status: `completed_verified`
- Date: 2026-08-27
- Mode: discovery
- Scope: `workshopScriptsFilIn/example11Optimization2DGradientKoFillIn.py`
- Patch count: 2 of 10; 8 remaining
- Defect: `searchNextPoint` accepted a GP, acquisition type, and best-observed value, but its optimizer callback read module-global variables instead. Independent calls could therefore optimize the wrong model or acquisition configuration.
- Repair: moved the minimization objective into `searchNextPoint`, where it closes over the arguments supplied to that call. Removed the stale fill-in marker for the completed function.
- Path defect: `getExperiment` compiled and ran KO using paths relative to the caller's current working directory, and `postProcessKo` did the same for `ko.dat` and `tv.txt`. Launching Example 11 from the repository root therefore failed to find `KOv13.f`.
- Path repair: anchored KO subprocesses and experiment-data paths to the directory containing `getExperiment.py`; added backward-compatible optional input/output paths to `postProcessKo`.
- Contract baseline: slides 109-110 of `AFRLRegionalNet-Workshop-2026.pdf` define the parameter bounds and require a search using the GP, acquisition function, and gradient search. No repository architecture or interface declaration governs this standalone workshop function.
- Contract result: no new drift; the implementation preserves the slide bounds and multistart L-BFGS-B acquisition search.
- Targeted validation: Python compilation passed; an isolated regression smoke test called `searchNextPoint` with two different dummy GP models and UCB/EI settings, returned distinct expected optima, and kept both results within bounds.
- End-to-end validation: the headless Example 11 run was launched from the repository root, completed three KO compile/run/postprocess evaluations, exited 0, and reported optimum `[1.75, 0.01]` in 68.319 seconds.
- Intermediary audit: the patch changes no simulation, plotting, data format, or public feature behavior. The optimizer callback now uses only the values belonging to its invocation.
- Generated side effects: the full run refreshed `ko.in`, `ke.dat`, and `tv.txt`, as expected for the existing KO workflow.
- Skipped broader lint gate: the file contains pre-existing workshop-style lint findings unrelated to this one behavioral repair; targeted `git diff --check` passed.
- Tracking records: no prior issue ledger, TODO, or known-issues record named this defect; this session log is the durable closure record.
