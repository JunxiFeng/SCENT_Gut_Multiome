#!/usr/bin/env python3
"""Prepare Gut multiome inputs for batch-wise SCENT runs.

The Gut ATAC h5ad has millions of 500 bp bins.  This helper builds the
gene-cis-bin batches once, writes a selected-gene RNA matrix once, and writes
batch-specific ATAC matrices on demand so the R runner only loads what it needs.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite


DEFAULT_DATA_DIR = Path("/data/pinello/PROJECTS/Gut_multiome")
DEFAULT_RNA = DEFAULT_DATA_DIR / "rna_subset_matched.h5ad"
DEFAULT_ATAC = DEFAULT_DATA_DIR / "atac_subset_single_v2.h5ad"
DEFAULT_GENES = DEFAULT_DATA_DIR / "selected_genes.txt"
DEFAULT_GTF = DEFAULT_DATA_DIR / "gencode.v32.annotation.gtf.gz"
WINDOW_BP = 512_000


def log(msg: str) -> None:
    print(f"[prepare_gut_multiome_scent] {msg}", flush=True)


def open_text(path: Path, mode: str = "rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode, newline="")


def write_lines_gz(path: Path, values: list[str]) -> None:
    with gzip.open(path, "wt") as handle:
        for value in values:
            handle.write(f"{value}\n")


def read_lines_gz(path: Path) -> list[str]:
    with gzip.open(path, "rt") as handle:
        return [line.rstrip("\n") for line in handle]


def slugify(value: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_")
    return slug or "celltype"


def decode_array(values) -> list[str]:
    return [x.decode() if isinstance(x, bytes) else str(x) for x in values]


def read_h5ad_index(path: Path, axis: str) -> list[str]:
    with h5py.File(path, "r") as handle:
        group = handle[axis]
        index_key = group.attrs.get("_index")
        if isinstance(index_key, bytes):
            index_key = index_key.decode()
        if index_key is None:
            index_key = "_index"
        return decode_array(group[index_key][:])


def read_obs_column(path: Path, column: str):
    with h5py.File(path, "r") as handle:
        obj = handle["obs"][column]
        if isinstance(obj, h5py.Group):
            categories = decode_array(obj["categories"][:])
            codes = obj["codes"][:]
            return [categories[int(code)] if int(code) >= 0 else None for code in codes]
        values = obj[:]
        if values.dtype.kind in {"S", "O"}:
            return decode_array(values)
        return values.tolist()


def make_unique_cell_ids(cells: list[str]) -> list[str]:
    return [f"{cell}__row{i + 1:06d}" for i, cell in enumerate(cells)]


def load_selected_genes(path: Path) -> list[str]:
    genes = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    duplicates = [gene for gene, count in Counter(genes).items() if count > 1]
    if duplicates:
        raise ValueError(f"selected gene file contains duplicate genes, e.g. {duplicates[:5]}")
    return genes


GTF_ATTR_RE = re.compile(r'(\S+) "([^"]+)"')


def parse_gtf_tss(path: Path, selected: set[str]) -> dict[str, dict[str, object]]:
    tss = {}
    duplicates = 0
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attrs = dict(GTF_ATTR_RE.findall(fields[8]))
            gene_name = attrs.get("gene_name")
            if gene_name not in selected:
                continue
            if gene_name in tss:
                duplicates += 1
                continue
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            tss_pos = start if strand == "+" else end
            tss[gene_name] = {
                "gene": gene_name,
                "gene_id": attrs.get("gene_id", ""),
                "chrom": fields[0],
                "start": start,
                "end": end,
                "strand": strand,
                "tss": tss_pos,
            }
    log(f"Parsed GTF TSS records for {len(tss)} selected genes; skipped {duplicates} duplicate gene records.")
    return tss


PEAK_RE = re.compile(r"^([^:]+):([0-9]+)-([0-9]+)$")


def normalize_peak_name(peak: str) -> str:
    return re.sub(r"[:-]", "_", peak)


def parse_atac_peaks(atac_h5ad: Path):
    peaks = read_h5ad_index(atac_h5ad, "var")
    by_chrom = defaultdict(lambda: {"starts": [], "ends": [], "indices": []})
    normalized = []
    for idx, peak in enumerate(peaks):
        match = PEAK_RE.match(peak)
        if not match:
            raise ValueError(f"ATAC feature does not match chr:start-end format: {peak}")
        chrom, start, end = match.group(1), int(match.group(2)), int(match.group(3))
        by_chrom[chrom]["starts"].append(start)
        by_chrom[chrom]["ends"].append(end)
        by_chrom[chrom]["indices"].append(idx)
        normalized.append(normalize_peak_name(peak))

    for chrom, vals in by_chrom.items():
        order = np.argsort(np.asarray(vals["starts"], dtype=np.int64), kind="mergesort")
        vals["starts"] = np.asarray(vals["starts"], dtype=np.int64)[order]
        vals["ends"] = np.asarray(vals["ends"], dtype=np.int64)[order]
        vals["indices"] = np.asarray(vals["indices"], dtype=np.int64)[order]

    return peaks, normalized, by_chrom


def iter_cis_pairs(genes: list[str], tss: dict[str, dict[str, object]], peaks: list[str], normalized: list[str], by_chrom):
    for gene in genes:
        rec = tss[gene]
        chrom = rec["chrom"]
        vals = by_chrom.get(chrom)
        if vals is None:
            continue
        left = int(rec["tss"]) - WINDOW_BP
        right = int(rec["tss"]) + WINDOW_BP
        starts = vals["starts"]
        lo = int(np.searchsorted(starts, left, side="left"))
        hi = int(np.searchsorted(starts, right, side="right"))
        for pos in range(lo, hi):
            peak_idx = int(vals["indices"][pos])
            if int(vals["ends"][pos]) >= left:
                yield gene, normalized[peak_idx], peak_idx, peaks[peak_idx]


def write_metadata(prepared_dir: Path, rna_h5ad: Path, atac_h5ad: Path) -> list[str]:
    rna_cells = read_h5ad_index(rna_h5ad, "obs")
    atac_cells = read_h5ad_index(atac_h5ad, "obs")
    if rna_cells != atac_cells:
        raise ValueError("RNA and ATAC obs indices are not identical and in the same order.")

    cell_ids = make_unique_cell_ids(rna_cells)
    level2 = read_obs_column(rna_h5ad, "level_2_annotation")
    n_counts = np.asarray(read_obs_column(rna_h5ad, "n_counts"), dtype=float)
    n_fragment = np.asarray(read_obs_column(atac_h5ad, "n_fragment"), dtype=float)
    tsse = np.asarray(read_obs_column(atac_h5ad, "tsse"), dtype=float)
    frac_mito = np.asarray(read_obs_column(atac_h5ad, "frac_mito"), dtype=float)
    frac_dup = np.asarray(read_obs_column(atac_h5ad, "frac_dup"), dtype=float)

    metadata_path = prepared_dir / "metadata.tsv.gz"
    with gzip.open(metadata_path, "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["cell_id", "original_cell", "celltype", "n_counts", "log1p_n_fragment", "tsse", "frac_mito", "frac_dup"])
        for row in zip(cell_ids, rna_cells, level2, n_counts, np.log1p(n_fragment), tsse, frac_mito, frac_dup):
            writer.writerow(row)

    counts = Counter(level2)
    celltype_path = prepared_dir / "celltypes.tsv"
    with open(celltype_path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["celltype", "slug", "n_cells"])
        for celltype in sorted(counts):
            writer.writerow([celltype, slugify(celltype), counts[celltype]])

    write_lines_gz(prepared_dir / "cells.tsv.gz", cell_ids)
    log(f"Wrote metadata for {len(cell_ids)} matched cells to {metadata_path}.")
    return cell_ids


def write_rna_matrix(prepared_dir: Path, rna_h5ad: Path, genes: list[str]) -> None:
    genes_path = prepared_dir / "genes.tsv.gz"
    matrix_path = prepared_dir / "rna_selected_genes.mtx.gz"
    if matrix_path.exists() and genes_path.exists():
        log(f"RNA matrix already exists at {matrix_path}; reusing it.")
        return

    log(f"Writing RNA matrix for {len(genes)} selected genes.")
    rna = ad.read_h5ad(rna_h5ad, backed="r")
    try:
        var_names = list(map(str, rna.var_names))
        gene_to_idx = {gene: idx for idx, gene in enumerate(var_names)}
        indices = [gene_to_idx[gene] for gene in genes]
        x = rna[:, indices].X
        if not sparse.issparse(x):
            x = sparse.csr_matrix(x)
        x = x.T.tocsr().astype(np.float64)
        with gzip.open(matrix_path, "wb") as handle:
            mmwrite(handle, x)
        write_lines_gz(genes_path, genes)
    finally:
        rna.file.close()
    log(f"Wrote RNA matrix to {matrix_path}.")


def write_gene_tss(prepared_dir: Path, genes: list[str], tss: dict[str, dict[str, object]]) -> None:
    path = prepared_dir / "gene_tss.tsv.gz"
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["gene", "gene_id", "chrom", "start", "end", "strand", "tss"], delimiter="\t")
        writer.writeheader()
        for gene in genes:
            writer.writerow(tss[gene])


def write_batches(prepared_dir: Path, genes: list[str], tss, peaks, normalized, by_chrom, nbatches: int) -> dict[str, object]:
    batches_dir = prepared_dir / "peak_info_batches"
    batches_dir.mkdir(parents=True, exist_ok=True)
    batch_paths = [batches_dir / f"peak_info_batch_{i:03d}.tsv.gz" for i in range(1, nbatches + 1)]
    if all(path.exists() for path in batch_paths):
        log(f"All {nbatches} peak-info batches already exist; reusing them.")
        total_pairs = sum(1 for path in batch_paths for _ in gzip.open(path, "rt")) - nbatches
        return {"total_pairs": total_pairs, "nbatches": nbatches}

    log("Counting cis gene-peak pairs.")
    total_pairs = sum(1 for _ in iter_cis_pairs(genes, tss, peaks, normalized, by_chrom))
    if total_pairs == 0:
        raise ValueError("No cis gene-peak pairs were generated.")

    log(f"Writing {total_pairs:,} cis pairs into {nbatches} batches.")
    handles = [gzip.open(path, "wt", newline="") for path in batch_paths]
    writers = []
    try:
        for handle in handles:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["Gene", "Peak", "atac_index", "PeakOriginal"])
            writers.append(writer)

        counts = [0] * nbatches
        unique_peak_indices = set()
        for pair_idx, (gene, peak_norm, peak_idx, peak_original) in enumerate(iter_cis_pairs(genes, tss, peaks, normalized, by_chrom)):
            batch_idx = min((pair_idx * nbatches) // total_pairs, nbatches - 1)
            writers[batch_idx].writerow([gene, peak_norm, peak_idx, peak_original])
            counts[batch_idx] += 1
            unique_peak_indices.add(peak_idx)
    finally:
        for handle in handles:
            handle.close()

    manifest = {
        "nbatches": nbatches,
        "total_pairs": total_pairs,
        "unique_atac_peaks": len(unique_peak_indices),
        "batch_pair_counts": counts,
    }
    with open(prepared_dir / "batch_manifest.json", "w") as handle:
        json.dump(manifest, handle, indent=2)
    return manifest


def prepare(args) -> None:
    prepared_dir = Path(args.prepared_dir)
    prepared_dir.mkdir(parents=True, exist_ok=True)

    selected_all = load_selected_genes(Path(args.genes))
    log(f"Loaded {len(selected_all)} selected genes from {args.genes}.")

    rna_genes = read_h5ad_index(Path(args.rna), "var")
    missing_rna = sorted(set(selected_all) - set(rna_genes))
    if missing_rna:
        raise ValueError(f"{len(missing_rna)} selected genes are missing from RNA; first examples: {missing_rna[:10]}")

    tss_all = parse_gtf_tss(Path(args.gtf), set(selected_all))
    missing_gtf = sorted(set(selected_all) - set(tss_all))
    if missing_gtf:
        raise ValueError(f"{len(missing_gtf)} selected genes are missing from GTF; first examples: {missing_gtf[:10]}")

    genes = selected_all[: args.gene_limit] if args.gene_limit else selected_all
    if args.gene_limit:
        log(f"Smoke/test mode: using first {len(genes)} genes after validating the full selected-gene list.")

    write_metadata(prepared_dir, Path(args.rna), Path(args.atac))
    write_gene_tss(prepared_dir, genes, tss_all)
    write_rna_matrix(prepared_dir, Path(args.rna), genes)

    log("Parsing ATAC bin coordinates.")
    peaks, normalized, by_chrom = parse_atac_peaks(Path(args.atac))
    manifest = write_batches(prepared_dir, genes, tss_all, peaks, normalized, by_chrom, args.nbatches)

    config = {
        "rna_h5ad": str(Path(args.rna)),
        "atac_h5ad": str(Path(args.atac)),
        "genes_file": str(Path(args.genes)),
        "gtf": str(Path(args.gtf)),
        "window_bp": WINDOW_BP,
        "nbatches": args.nbatches,
        "gene_count": len(genes),
        "validated_selected_gene_count": len(selected_all),
        **manifest,
    }
    with open(prepared_dir / "config.json", "w") as handle:
        json.dump(config, handle, indent=2)
    log(f"Preparation complete in {prepared_dir}.")


def read_config(prepared_dir: Path) -> dict[str, object]:
    with open(prepared_dir / "config.json") as handle:
        return json.load(handle)


def batch_matrices(args) -> None:
    prepared_dir = Path(args.prepared_dir)
    config = read_config(prepared_dir)
    batch = int(args.batch)
    batch_dir = prepared_dir / "batches" / f"batch_{batch:03d}"
    batch_dir.mkdir(parents=True, exist_ok=True)

    peak_info_path = prepared_dir / "peak_info_batches" / f"peak_info_batch_{batch:03d}.tsv.gz"
    if not peak_info_path.exists():
        raise FileNotFoundError(f"Missing peak-info batch: {peak_info_path}")

    atac_matrix_path = batch_dir / "atac.mtx.gz"
    peaks_path = batch_dir / "atac_peaks.tsv.gz"
    if atac_matrix_path.exists() and peaks_path.exists() and not args.force:
        log(f"ATAC matrix for batch {batch} already exists; reusing it.")
        return

    peak_df = pd.read_csv(peak_info_path, sep="\t", usecols=["Peak", "atac_index"])
    peak_df = peak_df.drop_duplicates("Peak", keep="first")
    peak_df = peak_df.sort_values("atac_index", kind="mergesort")
    peak_indices = peak_df["atac_index"].to_numpy(dtype=int)
    peak_names = peak_df["Peak"].astype(str).tolist()
    if len(peak_indices) == 0:
        raise ValueError(f"Batch {batch} has no ATAC peaks.")

    log(f"Writing ATAC matrix for batch {batch}: {len(peak_indices)} unique bins.")
    atac = ad.read_h5ad(config["atac_h5ad"], backed="r")
    try:
        x = atac[:, peak_indices].X
        if not sparse.issparse(x):
            x = sparse.csr_matrix(x)
        x = x.T.tocsr().astype(np.float64)
        with gzip.open(atac_matrix_path, "wb") as handle:
            mmwrite(handle, x)
        write_lines_gz(peaks_path, peak_names)
    finally:
        atac.file.close()
    log(f"Wrote ATAC matrix to {atac_matrix_path}.")


def list_celltypes(args) -> None:
    path = Path(args.prepared_dir) / "celltypes.tsv"
    df = pd.read_csv(path, sep="\t")
    if args.celltype:
        wanted = set(args.celltype)
        missing = wanted - set(df["celltype"])
        if missing:
            raise ValueError(f"Requested cell types are not present: {sorted(missing)}")
        df = df[df["celltype"].isin(args.celltype)]
    for celltype in df["celltype"].tolist():
        print(celltype)


def combine(args) -> None:
    results_dir = Path(args.results_dir)
    combined_dir = results_dir / "combined"
    combined_dir.mkdir(parents=True, exist_ok=True)

    celltypes = args.celltype
    if not celltypes:
        df = pd.read_csv(Path(args.prepared_dir) / "celltypes.tsv", sep="\t")
        celltypes = df["celltype"].tolist()

    for celltype in celltypes:
        slug = slugify(celltype)
        ct_dir = results_dir / slug
        files = sorted(ct_dir.glob("batch_*.csv"))
        frames = []
        for path in files:
            if path.stat().st_size == 0:
                continue
            frame = pd.read_csv(path)
            if frame.empty:
                continue
            frame.insert(0, "celltype", celltype)
            frames.append(frame)
        if not frames:
            log(f"No non-empty result batches found for {celltype}; skipping combine.")
            continue
        out = combined_dir / f"scent_results_{slug}.csv"
        pd.concat(frames, ignore_index=True).to_csv(out, index=False)
        log(f"Wrote combined results for {celltype} to {out}.")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    p_prepare = sub.add_parser("prepare", help="Validate inputs and write shared SCENT preparation files.")
    p_prepare.add_argument("--rna", default=str(DEFAULT_RNA))
    p_prepare.add_argument("--atac", default=str(DEFAULT_ATAC))
    p_prepare.add_argument("--genes", default=str(DEFAULT_GENES))
    p_prepare.add_argument("--gtf", default=str(DEFAULT_GTF))
    p_prepare.add_argument("--prepared-dir", required=True)
    p_prepare.add_argument("--nbatches", type=int, default=100)
    p_prepare.add_argument("--gene-limit", type=int, default=0, help="Use first N genes after validation; intended for smoke tests.")
    p_prepare.set_defaults(func=prepare)

    p_batch = sub.add_parser("batch-matrices", help="Write batch-specific ATAC matrix.")
    p_batch.add_argument("--prepared-dir", required=True)
    p_batch.add_argument("--batch", type=int, required=True)
    p_batch.add_argument("--force", action="store_true")
    p_batch.set_defaults(func=batch_matrices)

    p_list = sub.add_parser("list-celltypes", help="Print prepared cell types, one per line.")
    p_list.add_argument("--prepared-dir", required=True)
    p_list.add_argument("--celltype", action="append")
    p_list.set_defaults(func=list_celltypes)

    p_combine = sub.add_parser("combine", help="Combine per-batch SCENT results.")
    p_combine.add_argument("--prepared-dir", required=True)
    p_combine.add_argument("--results-dir", required=True)
    p_combine.add_argument("--celltype", action="append")
    p_combine.set_defaults(func=combine)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        args.func(args)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
