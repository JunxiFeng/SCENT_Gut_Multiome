#!/usr/bin/env python3

from __future__ import annotations

import csv
import gzip
import json
import math
import re
from collections import Counter
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt


ROOT = Path("/data/pinello/PROJECTS/2023_09_JF_SIMBAvariant/wenkai/SCENT_gut_multiome")
PREPARED = ROOT / "prepared"
LOGS = ROOT / "logs"
OUTDIR = ROOT / "runtime_report"

CONFIG_PATH = PREPARED / "config.json"
METADATA_PATH = PREPARED / "metadata.tsv.gz"
BENCHMARK_LOG = LOGS / "Th1_17_cells_batch_001_ds5000_seed1.log"
BENCHMARK_TIME = LOGS / "Th1_17_cells_batch_001_ds5000_seed1.time.txt"
FULL_LOG = LOGS / "Th1_17_cells_batch_001.log"
FULL_MONITOR = LOGS / "Th1_17_cells_batch_001.monitor.tsv"


def read_json(path: Path) -> dict:
    return json.loads(path.read_text())


def read_metadata_counts(path: Path) -> list[tuple[str, int]]:
    counts: Counter[str] = Counter()
    with gzip.open(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            counts[row["celltype"]] += 1
    return sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))


def parse_time_elapsed_minutes(path: Path) -> float | None:
    if not path.exists():
        return None
    text = path.read_text()
    match = re.search(r"Elapsed \(wall clock\) time .*:\s+([0-9:]+)", text)
    if not match:
        return None
    parts = [int(x) for x in match.group(1).split(":")]
    if len(parts) == 2:
        minutes, seconds = parts
        return minutes + seconds / 60
    if len(parts) == 3:
        hours, minutes, seconds = parts
        return hours * 60 + minutes + seconds / 60
    return None


def parse_log_runtime_minutes(path: Path) -> float | None:
    if not path.exists():
        return None
    text = path.read_text()
    match = re.search(r"Total R script runtime minutes=([0-9.]+)", text)
    if match:
        return float(match.group(1))
    return None


def parse_log_start(path: Path) -> datetime | None:
    if not path.exists():
        return None
    first_line = path.read_text().splitlines()[0]
    match = re.match(r"\[([^\]]+)\]", first_line)
    return datetime.strptime(match.group(1), "%Y-%m-%dT%H:%M:%S%z") if match else None


def parse_monitor_totals(path: Path) -> list[tuple[datetime, int]]:
    totals: list[tuple[datetime, int]] = []
    with path.open() as handle:
        next(handle)
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10 or parts[2] != "TOTAL":
                continue
            ts = datetime.strptime(parts[0], "%Y-%m-%dT%H:%M:%S%z")
            totals.append((ts, int(parts[7])))
    return totals


def fmt_days(days: float) -> str:
    if days >= 365:
        return f"{days / 365:.1f} years"
    return f"{days:.1f} days"


def fmt_hours(hours: float) -> str:
    if hours >= 24:
        return f"{hours / 24:.1f} days"
    return f"{hours:.1f} h"


def ensure_outdir() -> None:
    OUTDIR.mkdir(exist_ok=True)


def plot_celltype_runtime(celltype_counts: list[tuple[str, int]], benchmark_batch_min: float, nbatches: int) -> Path:
    labels = [name for name, _ in reversed(celltype_counts)]
    counts = [count for _, count in reversed(celltype_counts)]
    days = [min(count, 5000) / 5000 * benchmark_batch_min * nbatches / (60 * 24) for count in counts]

    fig, ax = plt.subplots(figsize=(10, 9))
    bars = ax.barh(labels, days, color="#3b82f6")
    ax.set_xlabel("Estimated serial runtime per cell type (days)")
    ax.set_title("Projected runtime if each cell type is capped at 5,000 cells")
    ax.grid(axis="x", alpha=0.25)

    for bar, count in zip(bars, counts):
        ax.text(bar.get_width() + 0.05, bar.get_y() + bar.get_height() / 2, f"{count:,} cells", va="center", fontsize=8)

    fig.tight_layout()
    out = OUTDIR / "celltype_runtime_scaled_cap.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out


def plot_parallel_scenarios(parallelism: list[int], scenario_days: dict[str, list[float]]) -> Path:
    fig, ax = plt.subplots(figsize=(10, 6))
    width = 0.25
    x = list(range(len(parallelism)))
    colors = ["#2563eb", "#16a34a", "#dc2626"]

    for idx, (label, values) in enumerate(scenario_days.items()):
        xpos = [v + (idx - 1) * width for v in x]
        ax.bar(xpos, values, width=width, label=label, color=colors[idx])

    ax.set_xticks(x)
    ax.set_xticklabels([str(p) for p in parallelism])
    ax.set_xlabel("Concurrent SCENT jobs")
    ax.set_ylabel("Projected project wall time (days, log scale)")
    ax.set_yscale("log")
    ax.set_title("Whole-project runtime under different assumptions")
    ax.grid(axis="y", which="both", alpha=0.25)
    ax.legend(frameon=False)

    fig.tight_layout()
    out = OUTDIR / "project_runtime_scenarios.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out


def plot_full_run_memory(totals: list[tuple[datetime, int]]) -> Path:
    start = totals[0][0]
    x_days = [((ts - start).total_seconds() / 86400) for ts, _ in totals]
    y_gb = [rss_kb / (1024 * 1024) for _, rss_kb in totals]

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(x_days, y_gb, color="#b91c1c", linewidth=1.5)
    ax.set_xlabel("Elapsed days since run start")
    ax.set_ylabel("Process-tree RSS (GB)")
    ax.set_title("Current full-cell Th1-17 batch: memory growth over time")
    ax.grid(alpha=0.25)

    fig.tight_layout()
    out = OUTDIR / "full_run_memory.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out


def main() -> None:
    ensure_outdir()

    config = read_json(CONFIG_PATH)
    celltype_counts = read_metadata_counts(METADATA_PATH)
    n_celltypes = len(celltype_counts)
    nbatches = int(config["nbatches"])
    total_pairs = int(config["total_pairs"])
    batch_pairs = int(config["batch_pair_counts"][0])

    benchmark_batch_min = parse_log_runtime_minutes(BENCHMARK_LOG) or parse_time_elapsed_minutes(BENCHMARK_TIME)
    if benchmark_batch_min is None:
        raise RuntimeError("Could not parse the completed downsampled benchmark runtime.")

    full_start = parse_log_start(FULL_LOG)
    full_totals = parse_monitor_totals(FULL_MONITOR)
    if full_start is None or not full_totals:
        raise RuntimeError("Could not parse the current full-run monitor log.")

    full_last_ts = full_totals[-1][0]
    full_lower_bound_min = (full_last_ts - full_start).total_seconds() / 60
    full_max_rss_kb = max(rss_kb for _, rss_kb in full_totals)
    full_last_rss_kb = full_totals[-1][1]
    full_has_end = "END scent_Th1_17_cells_batch_001 status=0" in FULL_LOG.read_text()

    fixed5k_batch_days = benchmark_batch_min / (60 * 24)
    fixed5k_project_days = fixed5k_batch_days * nbatches * n_celltypes

    scaled_cap_per_celltype_days = [
        min(count, 5000) / 5000 * benchmark_batch_min * nbatches / (60 * 24)
        for _, count in celltype_counts
    ]
    scaled_cap_project_days = sum(scaled_cap_per_celltype_days)

    full_lb_batch_days = full_lower_bound_min / (60 * 24)
    full_lb_project_days = full_lb_batch_days * nbatches * n_celltypes

    parallelism = [1, 4, 8, 16]
    scenario_days = {
        "Flat 5k benchmark": [fixed5k_project_days / p for p in parallelism],
        "Scaled, cap at 5k": [scaled_cap_project_days / p for p in parallelism],
        "Current full-run lower bound": [full_lb_project_days / p for p in parallelism],
    }

    plot1 = plot_celltype_runtime(celltype_counts, benchmark_batch_min, nbatches)
    plot2 = plot_parallel_scenarios(parallelism, scenario_days)
    plot3 = plot_full_run_memory(full_totals)

    top3 = celltype_counts[:3]
    top3_text = ", ".join(f"{name} ({count:,})" for name, count in top3)

    report_lines = [
        "# SCENT Gut Multiome Runtime Report",
        "",
        "## Snapshot",
        f"- Prepared cell types: {n_celltypes}",
        f"- Batches per cell type: {nbatches}",
        f"- Peak-gene pairs per batch: {batch_pairs:,}",
        f"- Total peak-gene pairs per cell type: {total_pairs:,}",
        f"- Total batch jobs for all cell types: {nbatches * n_celltypes:,}",
        f"- Largest cell types: {top3_text}",
        "",
        "## Measured Runs",
        f"- Completed benchmark: `Th1-17 cells`, `batch_001`, downsampled to 5,000 cells, finished in {benchmark_batch_min:.2f} minutes ({fmt_hours(benchmark_batch_min / 60)}).",
        f"- Current full-cell run: `Th1-17 cells`, `batch_001`, started {full_start.isoformat()}, successful `END` marker present: {full_has_end}.",
        f"- Current full-cell run lower bound: at least {full_lower_bound_min:.1f} minutes ({fmt_days(full_lb_batch_days)}) from start to the last monitor sample at {full_last_ts.isoformat()}.",
        f"- Current full-cell run memory: last sampled RSS {full_last_rss_kb / (1024 * 1024):.1f} GB; peak sampled RSS {full_max_rss_kb / (1024 * 1024):.1f} GB across the process tree.",
        "",
        "## What This Implies",
        f"- If every batch behaved like the completed 5k-cell benchmark, the whole project would still take about {fmt_days(fixed5k_project_days)} serially.",
        f"- If each cell type were downsampled and capped at 5,000 cells, and runtime scaled roughly with retained cells, the whole project drops to about {fmt_days(scaled_cap_project_days)} serially.",
        f"- If the stalled all-cell `Th1-17` batch is representative, then the full project lower bound is already about {fmt_days(full_lb_project_days)} serially.",
        f"- One large-cell-type lower bound alone is brutal: `Th1-17 cells` would need at least {fmt_days(full_lb_batch_days * nbatches)} for 100 batches with the current full-cell behavior.",
        "",
        "## Parallel Wall-Time Scenarios",
        "| Concurrent jobs | Flat 5k benchmark | Scaled, cap at 5k | Current full-run lower bound |",
        "| --- | ---: | ---: | ---: |",
    ]
    for idx, p in enumerate(parallelism):
        report_lines.append(
            f"| {p} | {scenario_days['Flat 5k benchmark'][idx]:.1f} | "
            f"{scenario_days['Scaled, cap at 5k'][idx]:.1f} | "
            f"{scenario_days['Current full-run lower bound'][idx]:.1f} |"
        )

    report_lines.extend(
        [
            "",
            "## Output Files",
            f"- `{plot1.name}`: estimated serial days per cell type under a 5,000-cell cap.",
            f"- `{plot2.name}`: whole-project wall time across several concurrency levels and runtime assumptions.",
            f"- `{plot3.name}`: memory growth for the incomplete all-cell `Th1-17` batch.",
            "",
            "## Bottom Line",
            "- The pair count alone is not the whole story; the all-cell run looks much worse than the 5,000-cell benchmark and appears to have gone pathological.",
            "- If you want this project to finish in a sane amount of time, the numbers argue for either cell downsampling, a much smaller candidate pair set, or both.",
        ]
    )

    (OUTDIR / "runtime_report.md").write_text("\n".join(report_lines) + "\n")


if __name__ == "__main__":
    main()
