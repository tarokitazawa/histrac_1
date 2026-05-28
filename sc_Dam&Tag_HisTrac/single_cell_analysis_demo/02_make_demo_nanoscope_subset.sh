#!/usr/bin/env bash
set -Eeuo pipefail

trap 'echo "[ERROR] line ${LINENO}: ${BASH_COMMAND}" >&2' ERR

ROOT="${1:?Usage: bash 02_make_demo_nanoscope_subset.fixed.sh /path/to/original/Day7R_MERGED /path/to/Day7R_MERGED_demo 5000}"
OUT="${2:-.}"
NPEAKS="${3:-5000}"

SAMPLES=("rep1" "rep2")
MODS=("Day7R.Dam_bb" "Day7R.K27ac_aa")

SEL="${OUT}/selected_cells.tsv"
TMP="${OUT}/_tmp_subset"

mkdir -p "${TMP}"

for cmd in bedtools bgzip tabix awk sort zcat cut wc; do
  command -v "${cmd}" >/dev/null 2>&1 || {
    echo "ERROR: ${cmd} is required but not found in PATH." >&2
    exit 1
  }
done

if [[ ! -f "${SEL}" ]]; then
  echo "ERROR: selected cell file not found: ${SEL}" >&2
  echo "Run 01_choose_demo_cells.R first." >&2
  exit 1
fi

for MOD in "${MODS[@]}"; do
  PREFIX="${MOD%_*}"  # Day7R.Dam_bb -> Day7R.Dam; Day7R.K27ac_aa -> Day7R.K27ac

  echo "============================================================"
  echo "Processing modality: ${MOD}"
  echo "Peak file prefix:     ${PREFIX}"
  echo "============================================================"

  CAND="${TMP}/${MOD}.candidate_merged_peaks.bed"
  FRAGPOOL="${TMP}/${MOD}.selected_cells.fragments.tsv"
  COUNTS="${TMP}/${MOD}.candidate_peak_fragment_counts.tsv"
  SUPPORTED="${TMP}/${MOD}.supported_peak_fragment_counts.tsv"
  SORTED="${TMP}/${MOD}.supported_peak_fragment_counts.sorted.tsv"
  SELPEAK="${OUT}/${MOD}.selected_${NPEAKS}_peaks.bed"

  # ------------------------------------------------------------
  # 1. Build candidate peak universe from both replicates
  # ------------------------------------------------------------
  {
    for SMP in "${SAMPLES[@]}"; do
      BP="${ROOT}/${SMP}/${MOD}/peaks/macs_broad/${PREFIX}_peaks.broadPeak"

      if [[ ! -f "${BP}" ]]; then
        echo "ERROR: missing broadPeak file: ${BP}" >&2
        exit 1
      fi

      awk 'BEGIN{OFS="\t"} !/^#/ && NF>=3 && $1 ~ /^chr/ {print $1,$2,$3}' "${BP}"
    done
  } | sort -k1,1 -k2,2n -k3,3n \
    | bedtools merge -i - \
    > "${CAND}"

  echo "Candidate merged peaks: $(wc -l < "${CAND}")"

  # ------------------------------------------------------------
  # 2. Pool fragments from selected cells for this modality
  # ------------------------------------------------------------
  {
    for SMP in "${SAMPLES[@]}"; do
      INFRAG="${ROOT}/${SMP}/${MOD}/cellranger/outs/fragments.tsv.gz"

      if [[ ! -f "${INFRAG}" ]]; then
        echo "ERROR: missing fragments file: ${INFRAG}" >&2
        exit 1
      fi

      zcat "${INFRAG}" \
        | awk -v smp="${SMP}" 'BEGIN{FS=OFS="\t"}
            FNR==NR {
              if (FNR > 1 && $1 == smp) keep[$2] = 1
              next
            }
            !/^#/ && ($4 in keep) {print $0}
          ' "${SEL}" -
    done
  } | sort -k1,1 -k2,2n -k3,3n \
    > "${FRAGPOOL}"

  echo "Fragments from selected cells: $(wc -l < "${FRAGPOOL}")"

  if [[ ! -s "${FRAGPOOL}" ]]; then
    echo "ERROR: no fragments were retained for ${MOD}." >&2
    echo "Check whether selected_cells.tsv barcodes match ${ROOT}/${SMP}/${MOD}/cellranger/outs/fragments.tsv.gz" >&2
    exit 1
  fi

  # ------------------------------------------------------------
  # 3. Select NPEAKS supported peaks without sort|head pipefail issue
  # ------------------------------------------------------------
  bedtools intersect -a "${CAND}" -b "${FRAGPOOL}" -c > "${COUNTS}"

  awk 'BEGIN{OFS="\t"} $4 > 0 {print $1,$2,$3,$4}' "${COUNTS}" > "${SUPPORTED}"

  NSUPPORTED=$(wc -l < "${SUPPORTED}")
  echo "Supported peaks with >=1 demo fragment: ${NSUPPORTED}"

  if [[ "${NSUPPORTED}" -eq 0 ]]; then
    echo "ERROR: no supported peaks found for ${MOD}." >&2
    exit 1
  fi

  sort -k4,4nr "${SUPPORTED}" > "${SORTED}"

  awk -v n="${NPEAKS}" 'BEGIN{OFS="\t"} NR <= n {print $1,$2,$3}' "${SORTED}" \
    | sort -k1,1 -k2,2n -k3,3n \
    > "${SELPEAK}"

  NSEL=$(wc -l < "${SELPEAK}")
  echo "Selected peaks for ${MOD}: ${NSEL}"

  if [[ "${NSEL}" -lt "${NPEAKS}" ]]; then
    echo "WARNING: selected fewer than ${NPEAKS} peaks for ${MOD}." >&2
  fi

  # ------------------------------------------------------------
  # 4. Write selected peak files into each replicate
  # ------------------------------------------------------------
  for SMP in "${SAMPLES[@]}"; do
    mkdir -p "${OUT}/${SMP}/${MOD}/peaks/macs_broad"
    mkdir -p "${OUT}/${SMP}/${MOD}/cellranger/outs"
    mkdir -p "${OUT}/${SMP}/${MOD}/barcode_metrics"

    # broadPeak-like format.
    # Your R code reads only columns 1:3, but this keeps the file valid enough.
    awk 'BEGIN{OFS="\t"} {
      print $1,$2,$3,"demo_peak_"NR,1000,".",0,0,0
    }' "${SELPEAK}" \
      > "${OUT}/${SMP}/${MOD}/peaks/macs_broad/${PREFIX}_peaks.broadPeak"

    cp "${SELPEAK}" "${OUT}/${SMP}/${MOD}/cellranger/outs/peaks.bed"
  done

  # ------------------------------------------------------------
  # 5. Filter fragments to selected cells and selected peaks
  # ------------------------------------------------------------
  for SMP in "${SAMPLES[@]}"; do
    echo "Filtering fragments: ${SMP} / ${MOD}"

    INFRAG="${ROOT}/${SMP}/${MOD}/cellranger/outs/fragments.tsv.gz"
    TMPFRAG="${TMP}/${MOD}.${SMP}.selected_cells.fragments.tsv"
    OUTFRAG="${OUT}/${SMP}/${MOD}/cellranger/outs/fragments.tsv.gz"

    zcat "${INFRAG}" \
      | awk -v smp="${SMP}" 'BEGIN{FS=OFS="\t"}
          FNR==NR {
            if (FNR > 1 && $1 == smp) keep[$2] = 1
            next
          }
          !/^#/ && ($4 in keep) {print $0}
        ' "${SEL}" - \
      > "${TMPFRAG}"

    if [[ ! -s "${TMPFRAG}" ]]; then
      echo "ERROR: no selected-cell fragments for ${SMP} / ${MOD}" >&2
      exit 1
    fi

    bedtools intersect -a "${TMPFRAG}" -b "${SELPEAK}" -u \
      | sort -k1,1 -k2,2n -k3,3n \
      | bgzip -c \
      > "${OUTFRAG}"

    if [[ ! -s "${OUTFRAG}" ]]; then
      echo "ERROR: output fragment file is empty: ${OUTFRAG}" >&2
      exit 1
    fi

    tabix -f -p bed "${OUTFRAG}"

    # Keep both files for nanoscope-like folder structure.
    awk -v smp="${SMP}" 'BEGIN{FS=OFS="\t"} FNR > 1 && $1 == smp {print $2}' "${SEL}" \
      > "${OUT}/${SMP}/${MOD}/barcode_metrics/all_barcodes.txt"

    zcat "${OUTFRAG}" | cut -f4 | sort -u \
      > "${OUT}/${SMP}/${MOD}/barcode_metrics/peaks_barcodes.txt"

    echo "  Output fragments: ${OUTFRAG}"
    echo "  Retained fragments: $(zcat "${OUTFRAG}" | wc -l)"
    echo "  Cells with retained fragments: $(wc -l < "${OUT}/${SMP}/${MOD}/barcode_metrics/peaks_barcodes.txt")"
  done
done

echo "============================================================"
echo "Done."
echo "============================================================"