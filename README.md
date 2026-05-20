# Gut Multiome SCENT Runner

This directory contains a command-line SCENT pipeline for:

- RNA: `/data/pinello/PROJECTS/Gut_multiome/rna_subset_matched.h5ad`
- ATAC: `/data/pinello/PROJECTS/Gut_multiome/atac_subset_single_v2.h5ad`
- Genes: `/data/pinello/PROJECTS/Gut_multiome/selected_genes.txt`
- GTF: `/data/pinello/PROJECTS/Gut_multiome/gencode.v32.annotation.gtf.gz`

The wrapper prepares gene-cis-ATAC-bin pairs within 512 kb of each TSS, writes a selected-gene RNA matrix, and creates per-batch ATAC matrices on demand so R does not need to load the full 6M-bin ATAC object.

## Commands

Prepare shared inputs only:

```bash
bash wenkai/SCENT_gut_multiome/run_SCENT_gut_multiome.sh --prepare-only
```

Run a small smoke test using the first selected gene split into 2 batches:

```bash
bash wenkai/SCENT_gut_multiome/run_SCENT_gut_multiome.sh --smoke-test
```

Run all cell types and all batches:

```bash
bash wenkai/SCENT_gut_multiome/run_SCENT_gut_multiome.sh --celltypes all --batches all --cores 6
```

Run one cell type and one batch:

```bash
bash wenkai/SCENT_gut_multiome/run_SCENT_gut_multiome.sh --celltype "Th1-17 cells" --batch 1 --cores 6
```

Run one cell type and one batch with optional cell downsampling:

```bash
bash wenkai/SCENT_gut_multiome/run_SCENT_gut_multiome.sh \
  --celltype "Th1-17 cells" \
  --batch 1 \
  --cores 6 \
  --downsample-cells 5000 \
  --downsample-seed 1
```

You can also label a run explicitly so the log filenames are distinct:

```bash
bash wenkai/SCENT_gut_multiome/run_SCENT_gut_multiome.sh \
  --celltype "Th1-17 cells" \
  --batch 1 \
  --cores 6 \
  --downsample-cells 5000 \
  --downsample-seed 1 \
  --run-label ds5000_try1
```

Use a shorter resource sampling interval while profiling:

```bash
bash wenkai/SCENT_gut_multiome/run_SCENT_gut_multiome.sh \
  --celltype "Th1-17 cells" \
  --batch 1 \
  --cores 6 \
  --monitor-interval 30
```

If `Rscript` is not on `PATH`, the wrapper defaults to:

```bash
/data/pinello/SHARED_SOFTWARE/anaconda_latest/envs/SCENT/bin/Rscript
```

You can override binaries and input paths:

```bash
R_BIN=/path/to/Rscript PY_BIN=/path/to/python bash wenkai/SCENT_gut_multiome/run_SCENT_gut_multiome.sh --smoke-test
```

## Outputs

- `prepared/`: metadata, selected-gene RNA matrix, gene TSS table, and peak-info batches.
- `prepared/batches/batch_###/`: batch-specific ATAC matrices.
- `results/<celltype_slug>/batch_###.csv`: SCENT output per cell type and batch.
- `results/combined/scent_results_<celltype_slug>.csv`: concatenated outputs.
- `logs/`: per-run stdout/stderr logs.

Failed SCENT batches also write `results/<celltype_slug>/batch_###.failed.txt`.

## Runtime and Resource Logs

Every heavy command writes three files into `logs/`:

- `*.log`: normal stdout/stderr plus R-side timestamps for loading, object construction, and SCENT runtime.
  The R log now also records whether all cells were used or the target cell type was downsampled.
- `*.time.txt`: `/usr/bin/time -v` output, including elapsed time and maximum resident set size.
- `*.monitor.tsv`: sampled process-tree resource usage while the command is running.

If `--run-label` is omitted and `--downsample-cells` is used, the wrapper automatically appends a label like `ds5000_seed1` to the log filenames.

The monitor table columns are:

```text
timestamp root_pid pid ppid etime pcpu pmem rss_kb vsz_kb comm
```

Each sample also includes a `TOTAL` row that sums RSS and VSZ across the observed process tree. Set `--monitor-interval 0` to disable live sampling, or use `--monitor-interval 30` for more frequent profiling.
