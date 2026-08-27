# Example 10 candidate-search repair

- Status: `completed_verified`
- Date: 2026-08-27
- Scope: `workshopScriptsFilIn/example10Optimization2DFillIn.py`
- Defect: `searchNextPoint` had no executable function body, so the script failed with an `IndentationError` after candidate-search code was pasted at module scope.
- Repair: moved the intended Sobol candidate generation, scaling, acquisition scoring, and maximum-score selection into `searchNextPoint`; removed the accidental test-fixture import and discarded acquisition recalculation from the working-tree edit.
- Contract check: no repository architecture or interface specification governs this standalone workshop exercise; the function's existing comments and callers require one two-element point, which the repaired function returns.
- Validation: `uv run python -m py_compile workshopScriptsFilIn/example10Optimization2DFillIn.py` passed. Headless full-script runs passed for UCB, PI, and EI. PI emitted a non-fatal scikit-learn optimizer convergence warning.
- Regression check: candidate acquisition shapes are `(128,)`, the returned point has shape `(2,)`, and `git diff --check` passed.
- Broader gate: `./scripts/check.sh` remains blocked by eight pre-existing Ruff findings outside the repaired workshop script.
- Tracking records: no prior defect ledger, TODO, or known-issues entry named this defect.
