# SCENT Gut Multiome Runtime Report

## Snapshot
- Prepared cell types: 29
- Batches per cell type: 100
- Peak-gene pairs per batch: 99,702
- Total peak-gene pairs per cell type: 9,970,188
- Total batch jobs for all cell types: 2,900
- Largest cell types: Th1-17 cells (20,113), CD8 T cells (16,680), Epithelial (12,464)

## Measured Runs
- Completed benchmark: `Th1-17 cells`, `batch_001`, downsampled to 5,000 cells, finished in 189.61 minutes (3.2 h).
- Current full-cell run: `Th1-17 cells`, `batch_001`, started 2026-05-04T06:16:17+00:00, successful `END` marker present: False.
- Current full-cell run lower bound: at least 14995.6 minutes (10.4 days) from start to the last monitor sample at 2026-05-14T16:11:55+00:00.
- Current full-cell run memory: last sampled RSS 40.6 GB; peak sampled RSS 42.9 GB across the process tree.

## What This Implies
- If every batch behaved like the completed 5k-cell benchmark, the whole project would still take about 1.0 years serially.
- If each cell type were downsampled and capped at 5,000 cells, and runtime scaled roughly with retained cells, the whole project drops to about 151.6 days serially.
- If the stalled all-cell `Th1-17` batch is representative, then the full project lower bound is already about 82.7 years serially.
- One large-cell-type lower bound alone is brutal: `Th1-17 cells` would need at least 2.9 years for 100 batches with the current full-cell behavior.

## Parallel Wall-Time Scenarios
| Concurrent jobs | Flat 5k benchmark | Scaled, cap at 5k | Current full-run lower bound |
| --- | ---: | ---: | ---: |
| 1 | 381.9 | 151.6 | 30199.5 |
| 4 | 95.5 | 37.9 | 7549.9 |
| 8 | 47.7 | 19.0 | 3774.9 |
| 16 | 23.9 | 9.5 | 1887.5 |

## Output Files
- `celltype_runtime_scaled_cap.png`: estimated serial days per cell type under a 5,000-cell cap.
- `project_runtime_scenarios.png`: whole-project wall time across several concurrency levels and runtime assumptions.
- `full_run_memory.png`: memory growth for the incomplete all-cell `Th1-17` batch.

## Bottom Line
- The pair count alone is not the whole story; the all-cell run looks much worse than the 5,000-cell benchmark and appears to have gone pathological.
- If you want this project to finish in a sane amount of time, the numbers argue for either cell downsampling, a much smaller candidate pair set, or both.
