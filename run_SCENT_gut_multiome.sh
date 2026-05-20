#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY_SCRIPT="${SCRIPT_DIR}/prepare_gut_multiome_scent.py"
R_SCRIPT="${SCRIPT_DIR}/run_SCENT_gut_multiome.R"

DATA_DIR="${DATA_DIR:-/data/pinello/PROJECTS/Gut_multiome}"
RNA_H5AD="${RNA_H5AD:-${DATA_DIR}/rna_subset_matched.h5ad}"
ATAC_H5AD="${ATAC_H5AD:-${DATA_DIR}/atac_subset_single_v2.h5ad}"
GENES_FILE="${GENES_FILE:-${DATA_DIR}/selected_genes.txt}"
GTF_FILE="${GTF_FILE:-${DATA_DIR}/gencode.v32.annotation.gtf.gz}"

PY_BIN="${PY_BIN:-python}"
R_BIN="${R_BIN:-}"
if [[ -z "${R_BIN}" ]]; then
  if command -v Rscript >/dev/null 2>&1; then
    R_BIN="$(command -v Rscript)"
  elif [[ -x /data/pinello/SHARED_SOFTWARE/anaconda_latest/envs/SCENT/bin/Rscript ]]; then
    R_BIN="/data/pinello/SHARED_SOFTWARE/anaconda_latest/envs/SCENT/bin/Rscript"
  else
    echo "ERROR: Rscript is not on PATH and the shared SCENT Rscript was not found." >&2
    echo "Set R_BIN=/path/to/Rscript and rerun." >&2
    exit 1
  fi
fi

PREPARED_DIR="${PREPARED_DIR:-${SCRIPT_DIR}/prepared}"
RESULTS_DIR="${RESULTS_DIR:-${SCRIPT_DIR}/results}"
LOGS_DIR="${LOGS_DIR:-${SCRIPT_DIR}/logs}"
MONITOR_INTERVAL="${MONITOR_INTERVAL:-60}"
NBATCHES=100
CORES=6
REGR="poisson"
BIN="TRUE"
DOWNSAMPLE_CELLS=0
DOWNSAMPLE_SEED=1
RUN_LABEL=""
PREPARE_ONLY=0
FORCE_PREPARE=0
FORCE_BATCH=0
FORCE_SCENT=0
SMOKE_TEST=0
GENE_LIMIT=0
COMBINE_ONLY=0
CELLTYPES_MODE="all"
BATCHES_MODE="all"
declare -a REQUESTED_CELLTYPES=()
declare -a REQUESTED_BATCHES=()

usage() {
  cat <<'EOF'
Usage:
  bash run_SCENT_gut_multiome.sh [options]

Common examples:
  bash run_SCENT_gut_multiome.sh --prepare-only
  bash run_SCENT_gut_multiome.sh --celltypes all --batches all --cores 6
  bash run_SCENT_gut_multiome.sh --celltype "Th1-17 cells" --batch 1 --cores 6
  bash run_SCENT_gut_multiome.sh --smoke-test

Options:
  --prepare-only             Validate inputs and write shared prepared files, then stop.
  --smoke-test               Use first selected gene, 2 batches, and one default cell type.
  --celltypes all            Run all prepared level_2_annotation cell types. Default.
  --celltype NAME            Run one cell type. Can be repeated.
  --batches all              Run all batches. Default.
  --batch N                  Run one batch. Can be repeated.
  --nbatches N               Number of peak-info batches for preparation. Default: 100.
  --cores N                  Cores passed to SCENT_algorithm. Default: 6.
  --regr NAME                SCENT regression model. Default: poisson.
  --bin TRUE|FALSE           Binarize ATAC counts for SCENT. Default: TRUE.
  --downsample-cells N       Randomly keep at most N cells from the requested cell type. Default: off.
  --downsample-seed N        Seed used for downsampling. Default: 1.
  --run-label LABEL          Append LABEL to log filenames for this run.
  --combine-only             Combine existing per-batch results only.
  --force-prepare            Recreate the prepared directory.
  --force-batch              Recreate batch-specific ATAC matrices.
  --force-scent              Rerun SCENT even if a batch output exists.
  --prepared-dir DIR         Prepared output directory.
  --results-dir DIR          SCENT result directory.
  --logs-dir DIR             Log directory.
  --monitor-interval N       Seconds between resource monitor samples. Default: 60.
  --help                     Show this help.

Environment overrides:
  PY_BIN, R_BIN, DATA_DIR, RNA_H5AD, ATAC_H5AD, GENES_FILE, GTF_FILE, MONITOR_INTERVAL
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --prepare-only) PREPARE_ONLY=1; shift ;;
    --smoke-test) SMOKE_TEST=1; shift ;;
    --celltypes) CELLTYPES_MODE="$2"; shift 2 ;;
    --celltype) REQUESTED_CELLTYPES+=("$2"); CELLTYPES_MODE="selected"; shift 2 ;;
    --batches) BATCHES_MODE="$2"; shift 2 ;;
    --batch) REQUESTED_BATCHES+=("$2"); BATCHES_MODE="selected"; shift 2 ;;
    --nbatches) NBATCHES="$2"; shift 2 ;;
    --cores) CORES="$2"; shift 2 ;;
    --regr) REGR="$2"; shift 2 ;;
    --bin) BIN="$2"; shift 2 ;;
    --downsample-cells) DOWNSAMPLE_CELLS="$2"; shift 2 ;;
    --downsample-seed) DOWNSAMPLE_SEED="$2"; shift 2 ;;
    --run-label) RUN_LABEL="$2"; shift 2 ;;
    --combine-only) COMBINE_ONLY=1; shift ;;
    --force-prepare) FORCE_PREPARE=1; shift ;;
    --force-batch) FORCE_BATCH=1; shift ;;
    --force-scent) FORCE_SCENT=1; shift ;;
    --prepared-dir) PREPARED_DIR="$2"; shift 2 ;;
    --results-dir) RESULTS_DIR="$2"; shift 2 ;;
    --logs-dir) LOGS_DIR="$2"; shift 2 ;;
    --monitor-interval) MONITOR_INTERVAL="$2"; shift 2 ;;
    --help|-h) usage; exit 0 ;;
    *) echo "ERROR: unknown option: $1" >&2; usage >&2; exit 1 ;;
  esac
done

if [[ "${SMOKE_TEST}" -eq 1 ]]; then
  PREPARED_DIR="${PREPARED_DIR%/}_smoke"
  RESULTS_DIR="${RESULTS_DIR%/}_smoke"
  LOGS_DIR="${LOGS_DIR%/}_smoke"
  NBATCHES=2
  GENE_LIMIT=1
  if [[ "${#REQUESTED_CELLTYPES[@]}" -eq 0 ]]; then
    REQUESTED_CELLTYPES=("Th1-17 cells")
    CELLTYPES_MODE="selected"
  fi
  if [[ "${#REQUESTED_BATCHES[@]}" -eq 0 ]]; then
    BATCHES_MODE="all"
  fi
fi

mkdir -p "${PREPARED_DIR}" "${RESULTS_DIR}" "${LOGS_DIR}"

slugify() {
  python - "$1" <<'PY'
import re
import sys
slug = re.sub(r"[^A-Za-z0-9]+", "_", sys.argv[1]).strip("_") or "run"
print(slug)
PY
}

if [[ -z "${RUN_LABEL}" && "${DOWNSAMPLE_CELLS}" -gt 0 ]]; then
  RUN_LABEL="ds${DOWNSAMPLE_CELLS}_seed${DOWNSAMPLE_SEED}"
fi

RUN_LABEL_SLUG=""
if [[ -n "${RUN_LABEL}" ]]; then
  RUN_LABEL_SLUG="$(slugify "${RUN_LABEL}")"
fi

timestamp() {
  date '+%Y-%m-%dT%H:%M:%S%z'
}

sample_process_tree() {
  local root_pid="$1"
  local monitor_file="$2"
  local pids="${root_pid}"
  local children=""

  while true; do
    children="$(ps -eo pid=,ppid= | awk -v parents="${pids}" '
      BEGIN {
        split(parents, parent_arr, " ")
        for (i in parent_arr) {
          if (parent_arr[i] != "") seen[parent_arr[i]] = 1
        }
      }
      seen[$2] && !seen[$1] { print $1 }
    ')"
    if [[ -z "${children}" ]]; then
      break
    fi
    pids="${pids} ${children}"
  done

  local pid_csv
  pid_csv="$(printf '%s\n' ${pids} | paste -sd, -)"
  local ts
  ts="$(timestamp)"
  ps -o pid=,ppid=,etime=,%cpu=,%mem=,rss=,vsz=,comm= -p "${pid_csv}" 2>/dev/null | \
    awk -v ts="${ts}" -v root="${root_pid}" '
      {
        rss += $6
        vsz += $7
        print ts "\t" root "\t" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8
      }
      END {
        if (NR > 0) {
          print ts "\t" root "\tTOTAL\tNA\tNA\tNA\tNA\t" rss "\t" vsz "\tprocess_tree"
        }
      }
    ' >> "${monitor_file}"
}

run_monitored() {
  local label="$1"
  local stdout_log="$2"
  local monitor_file="$3"
  local time_file="$4"
  shift 4

  mkdir -p "$(dirname "${stdout_log}")" "$(dirname "${monitor_file}")" "$(dirname "${time_file}")"
  {
    echo "[$(timestamp)] START ${label}"
    printf 'Command:'
    printf ' %q' "$@"
    echo
  } >> "${stdout_log}"
  echo -e "timestamp\troot_pid\tpid\tppid\tetime\tpcpu\tpmem\trss_kb\tvsz_kb\tcomm" > "${monitor_file}"

  local had_errexit=0
  if [[ $- == *e* ]]; then
    had_errexit=1
  fi
  set +e

  if [[ -x /usr/bin/time ]]; then
    /usr/bin/time -v -o "${time_file}" "$@" >> "${stdout_log}" 2>&1 &
  else
    echo "[$(timestamp)] /usr/bin/time not found; writing live monitor only." >> "${stdout_log}"
    "$@" >> "${stdout_log}" 2>&1 &
  fi
  local cmd_pid=$!

  if [[ "${MONITOR_INTERVAL}" != "0" ]]; then
    (
      while kill -0 "${cmd_pid}" 2>/dev/null; do
        sample_process_tree "${cmd_pid}" "${monitor_file}"
        sleep "${MONITOR_INTERVAL}"
      done
    ) &
    local monitor_pid=$!
  else
    local monitor_pid=""
  fi

  wait "${cmd_pid}"
  local status=$?
  if [[ -n "${monitor_pid}" ]]; then
    kill "${monitor_pid}" 2>/dev/null
    wait "${monitor_pid}" 2>/dev/null
  fi

  echo "[$(timestamp)] END ${label} status=${status}" >> "${stdout_log}"
  if [[ "${had_errexit}" -eq 1 ]]; then
    set -e
  fi
  return "${status}"
}

prepare_inputs() {
  if [[ "${FORCE_PREPARE}" -eq 1 ]]; then
    rm -rf "${PREPARED_DIR}"
    mkdir -p "${PREPARED_DIR}"
  fi

  if [[ -s "${PREPARED_DIR}/config.json" && "${FORCE_PREPARE}" -eq 0 ]]; then
    echo "[run_SCENT_gut_multiome] Reusing prepared files in ${PREPARED_DIR}"
    return
  fi

  echo "[run_SCENT_gut_multiome] Preparing Gut multiome SCENT inputs."
  local prep_log="${LOGS_DIR}/prepare.log"
  echo "[run_SCENT_gut_multiome] Prepare log: ${prep_log}"
  run_monitored "prepare_inputs" \
    "${prep_log}" \
    "${LOGS_DIR}/prepare.monitor.tsv" \
    "${LOGS_DIR}/prepare.time.txt" \
    "${PY_BIN}" "${PY_SCRIPT}" prepare \
    --rna "${RNA_H5AD}" \
    --atac "${ATAC_H5AD}" \
    --genes "${GENES_FILE}" \
    --gtf "${GTF_FILE}" \
    --prepared-dir "${PREPARED_DIR}" \
    --nbatches "${NBATCHES}" \
    --gene-limit "${GENE_LIMIT}"
}

config_value() {
  "${PY_BIN}" - "$PREPARED_DIR/config.json" "$1" <<'PY'
import json
import sys
with open(sys.argv[1]) as handle:
    data = json.load(handle)
print(data[sys.argv[2]])
PY
}

expand_batches() {
  local nbatches="$1"
  if [[ "${BATCHES_MODE}" == "all" ]]; then
    seq 1 "${nbatches}"
  elif [[ "${BATCHES_MODE}" == "selected" ]]; then
    printf '%s\n' "${REQUESTED_BATCHES[@]}"
  else
    echo "ERROR: --batches currently supports only 'all'." >&2
    exit 1
  fi
}

prepare_inputs

if [[ "${PREPARE_ONLY}" -eq 1 ]]; then
  echo "[run_SCENT_gut_multiome] Preparation complete. Stopping because --prepare-only was requested."
  exit 0
fi

if [[ "${COMBINE_ONLY}" -eq 0 ]]; then
  if [[ "${CELLTYPES_MODE}" == "all" ]]; then
    mapfile -t CELLTYPES < <("${PY_BIN}" "${PY_SCRIPT}" list-celltypes --prepared-dir "${PREPARED_DIR}")
  elif [[ "${CELLTYPES_MODE}" == "selected" ]]; then
    list_args=()
    for celltype in "${REQUESTED_CELLTYPES[@]}"; do
      list_args+=(--celltype "${celltype}")
    done
    mapfile -t CELLTYPES < <("${PY_BIN}" "${PY_SCRIPT}" list-celltypes --prepared-dir "${PREPARED_DIR}" "${list_args[@]}")
  else
    echo "ERROR: --celltypes currently supports only 'all'." >&2
    exit 1
  fi

  NBATCHES_PREPARED="$(config_value nbatches)"
  mapfile -t BATCHES < <(expand_batches "${NBATCHES_PREPARED}")

  echo "[run_SCENT_gut_multiome] Running ${#CELLTYPES[@]} cell type(s) across ${#BATCHES[@]} batch(es)."
  for batch in "${BATCHES[@]}"; do
    batch_args=()
    if [[ "${FORCE_BATCH}" -eq 1 ]]; then
      batch_args+=(--force)
    fi
    batch_label="$(printf 'batch_%03d' "${batch}")"
    batch_log_base="${LOGS_DIR}/${batch_label}_matrix"
    if [[ -n "${RUN_LABEL_SLUG}" ]]; then
      batch_log_base="${batch_log_base}_${RUN_LABEL_SLUG}"
    fi
    batch_log="${batch_log_base}.log"
    echo "[run_SCENT_gut_multiome] Preparing ATAC matrix for ${batch_label}; log=${batch_log}"
    run_monitored "${batch_label}_matrix" \
      "${batch_log}" \
      "${batch_log_base}.monitor.tsv" \
      "${batch_log_base}.time.txt" \
      "${PY_BIN}" "${PY_SCRIPT}" batch-matrices \
      --prepared-dir "${PREPARED_DIR}" \
      --batch "${batch}" \
      "${batch_args[@]}"

    for celltype in "${CELLTYPES[@]}"; do
      slug="$(slugify "${celltype}")"
      log_base="${LOGS_DIR}/${slug}_batch_$(printf '%03d' "${batch}")"
      if [[ -n "${RUN_LABEL_SLUG}" ]]; then
        log_base="${log_base}_${RUN_LABEL_SLUG}"
      fi
      log_file="${log_base}.log"
      scent_args=()
      if [[ "${FORCE_SCENT}" -eq 1 ]]; then
        scent_args+=(--force)
      fi
      echo "[run_SCENT_gut_multiome] SCENT celltype='${celltype}' batch=${batch}; log=${log_file}"
      set +e
      run_monitored "scent_${slug}_batch_$(printf '%03d' "${batch}")" \
        "${log_file}" \
        "${log_base}.monitor.tsv" \
        "${log_base}.time.txt" \
        "${R_BIN}" "${R_SCRIPT}" \
        --prepared-dir "${PREPARED_DIR}" \
        --results-dir "${RESULTS_DIR}" \
        --logs-dir "${LOGS_DIR}" \
        --batch "${batch}" \
        --celltype "${celltype}" \
        --cores "${CORES}" \
        --regr "${REGR}" \
        --bin "${BIN}" \
        --downsample-cells "${DOWNSAMPLE_CELLS}" \
        --downsample-seed "${DOWNSAMPLE_SEED}" \
        "${scent_args[@]}" >"${log_file}" 2>&1
      status=$?
      set -e
      if [[ "${status}" -ne 0 ]]; then
        echo "[run_SCENT_gut_multiome] WARNING: SCENT failed for celltype='${celltype}' batch=${batch}. See ${log_file}" >&2
      fi
    done
  done
fi

if [[ "${CELLTYPES_MODE}" == "all" ]]; then
  run_monitored "combine_all" \
    "${LOGS_DIR}/combine.log" \
    "${LOGS_DIR}/combine.monitor.tsv" \
    "${LOGS_DIR}/combine.time.txt" \
    "${PY_BIN}" "${PY_SCRIPT}" combine --prepared-dir "${PREPARED_DIR}" --results-dir "${RESULTS_DIR}"
else
  combine_args=()
  for celltype in "${REQUESTED_CELLTYPES[@]}"; do
    combine_args+=(--celltype "${celltype}")
  done
  run_monitored "combine_selected" \
    "${LOGS_DIR}/combine.log" \
    "${LOGS_DIR}/combine.monitor.tsv" \
    "${LOGS_DIR}/combine.time.txt" \
    "${PY_BIN}" "${PY_SCRIPT}" combine --prepared-dir "${PREPARED_DIR}" --results-dir "${RESULTS_DIR}" "${combine_args[@]}"
fi

echo "[run_SCENT_gut_multiome] Done."
