#!/bin/bash
# set -euo pipefail  # Optional: enable strict error checking

script_name="./repetitive_element_motif_search.sh"

usage() {
  echo "-h Help documentation for $script_name"
  echo "-l  -- Full path to file containing list of repetitive elements."
  echo "-g  -- Genome (human or mouse)."
  echo "-o  -- Full path to output directory (optional, otherwise inferred from list file)."
  echo "Example: $script_name -l /optional/path/to/myrepeatlist.list -g human"
  exit 1
}

main() {
  OPTIND=1
  while getopts :l:g:o:h opt; do
    case $opt in
      l) listfile="$OPTARG" ;;
      g) genome="$OPTARG" ;;
      o) outdir="$OPTARG" ;;
      h|*) usage ;;
    esac
  done
  shift $((OPTIND - 1))

  if [[ -z "${listfile:-}" || -z "${genome:-}" ]]; then usage; fi
  [[ -f "$listfile" ]] || { echo "List file not found: $listfile"; exit 1; }

  outdir="${outdir:-$(dirname "$listfile")}"
  mkdir -p "$outdir/temp" "$outdir/motif_search" "$outdir/heatmaps"

  cputhread=$((SLURM_CPUS_ON_NODE - 3))
  echo "Using $cputhread CPUs"

  case "$genome" in
    mouse)
      fasta="/project/GCRB/Banaszynski_lab/shared/repetitive_element_files/dfam/mm10/mm10_dfam_families.fa"
      beddir="/project/GCRB/Banaszynski_lab/shared/jaspar_motif_files/beds/mm10_motif_bedfiles/bw"
      build="mm10"
      repdir="/project/GCRB/Banaszynski_lab/shared/repetitive_element_files/dfam/mm10/individual_family_beds"
      ;;
    human)
      fasta="/project/GCRB/Banaszynski_lab/shared/repetitive_element_files/dfam/hg38/hg38_dfam_families.fa"
      beddir="/project/GCRB/Banaszynski_lab/shared/jaspar_motif_files/beds/hg38_motif_bedfiles/bw"
      build="hg38"
      repdir="/project/GCRB/Banaszynski_lab/shared/repetitive_element_files/dfam/hg38/individual_family_beds"
      ;;
    *) usage ;;
  esac

  module load meme/5.3.3
  module load deeptools/2.5.0.1

  motifs="/project/GCRB/Banaszynski_lab/shared/jaspar_motif_files/jaspar_motifs.meme"

  while IFS= read -r repeat; do
    [[ -z "$repeat" ]] && continue
    echo "Processing: $repeat"

    bedfile="$repdir/${repeat}.bed"
    [[ -f "$bedfile" ]] || { echo "Missing BED file: $bedfile"; continue; }

    grep -w "$repeat" -A 1 "$fasta" > "$outdir/temp/${repeat}.fa"

    fimo --text --max-strand --verbosity 1 "$motifs" "$outdir/temp/${repeat}.fa" \
      > "$outdir/temp/${repeat}.hitlist"

    tail -n +2 "$outdir/temp/${repeat}.hitlist" | sort -k7,7nr -k2 | \
      awk -v r="$repeat" 'BEGIN{OFS="\t"}{print $2,$4,$5,$9,$8,$6}' \
      > "$outdir/motif_search/${repeat}_motifs.tab"

    cut -f 1 "$outdir/motif_search/${repeat}_motifs.tab" | awk '!seen[$0]++' | head -n 10 \
      > "$outdir/temp/${repeat}.top10"

    readarray -t motifs_array < "$outdir/temp/${repeat}.top10"
    for i in {1..10}; do
      eval "motif$i='${motifs_array[$((i-1))]}'"
    done

    # Subset top 2000 longest elements (or all if <2000)
    topbed="$outdir/temp/${repeat}_top2000.bed"
    awk '{print $0 "\t" ($3 - $2)}' "$bedfile" | sort -k5,5nr | cut -f1-4 | head -n 2000 > "$topbed"

    numbed=$(wc -l < "$topbed")

    computeMatrix scale-regions -S \
      "$beddir/${build}_${motif1}.bw" \
      "$beddir/${build}_${motif2}.bw" \
      "$beddir/${build}_${motif3}.bw" \
      "$beddir/${build}_${motif4}.bw" \
      "$beddir/${build}_${motif5}.bw" \
      "$beddir/${build}_${motif6}.bw" \
      "$beddir/${build}_${motif7}.bw" \
      "$beddir/${build}_${motif8}.bw" \
      "$beddir/${build}_${motif9}.bw" \
      "$beddir/${build}_${motif10}.bw" \
      -R "$topbed" \
      --averageTypeBins max -b 250 -a 250 -p "$cputhread" --binSize 10 \
      --missingDataAsZero --sortRegions descend \
      -o "$outdir/heatmaps/${repeat}_motifs_heatmap.matrix.gz"

    plotHeatmap -m "$outdir/heatmaps/${repeat}_motifs_heatmap.matrix.gz" \
      --samplesLabel "${motifs_array[@]}" \
      -T "" --sortRegions no \
      --startLabel "" --endLabel "" --xAxisLabel "" \
      --regionsLabel "${repeat} n=${numbed}" \
      --colorList "white,black" --boxAroundHeatmaps yes \
      --zMin 0 --zMax 1 \
      -o "$outdir/heatmaps/${repeat}_motifs_heatmap.pdf"

    rm -f "$outdir/heatmaps/${repeat}_motifs_heatmap.matrix.gz"
  done < "$listfile"

  rm -rf "$outdir/temp"

  if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    mv "serialJob.${SLURM_JOB_ID}.out" "serialJob.${SLURM_JOB_ID}.repeat_motif_search.out"
    mv "serialJob.${SLURM_JOB_ID}.time" "serialJob.${SLURM_JOB_ID}.repeat_motif_search.err"
  fi
}

main "$@"
