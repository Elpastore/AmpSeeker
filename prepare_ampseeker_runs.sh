#!/bin/bash
#SBATCH -A asadio
#SBATCH -J AmpSeeker-prep
#SBATCH --output=logs/prepare.%J.out
#SBATCH --error=logs/prepare.%J.error
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#
# Script: prepare_ampseeker_runs.sh
# Purpose: Go from one or more raw MiSeq run folders to a ready-to-run
#          AmpSeeker pipeline, in one command:
#            1. convert each run's SampleSheet.csv (convert_samplesheet.py;
#               original kept as backup_SampleSheet.csv)
#            2. demultiplex with bcl2fastq (or link on-instrument FASTQs) and
#               build one combined metadata TSV (build_multirun_metadata.py)
#            3. download the reference genome if the config points to a
#               FASTA/GFF that doesn't exist yet
#            4. fix conda channel order in workflow/envs/*.yaml
#               (conda-forge before bioconda; needed with
#               `channel_priority: strict`, otherwise htslib/libdeflate
#               can't be solved)
#            5. build all Snakemake conda envs (--conda-create-envs-only)
#            6. Snakemake dry run (-n)
#            7. optionally submit run_AmpSeeker_pipeline.sh with sbatch
#          Every step is idempotent: re-running skips work already done.
#
# Usage (from anywhere; paths may be relative to where you run it):
#   bash AmpSeeker/prepare_ampseeker_runs.sh [options] RUN_DIR [RUN_DIR ...]
#   sbatch AmpSeeker/prepare_ampseeker_runs.sh [options] RUN_DIR ...   (on a compute node)
#
# Options:
#   -c, --config FILE      AmpSeeker config (default: config/moussa_agam.yaml)
#   -o, --metadata FILE    metadata TSV to write (default: the config's 'metadata:' value)
#       --fastq-mode MODE  rebuild   = bcl2fastq into <run>/fastq_files (default)
#                          instrument = symlink the MiSeq Alignment_*/Fastq files
#       --force-fastq      re-run bcl2fastq even if fastq_files/ is complete
#       --threads N        bcl2fastq threads (default: 8 or $SLURM_CPUS_PER_TASK)
#       --on-duplicate X   warn|error|suffix for sample_id seen in >1 run (default: error)
#       --skip-envs        don't build conda envs (step 5)
#       --submit           sbatch run_AmpSeeker_pipeline.sh if the dry run passes
#       --from-step N      start at step N (1-7), e.g. 5 to only rebuild envs + dry run
#   -h, --help
#
# Example (what was run for MACARD-1 + MACARD-2 on 2026-11-20):
#   cd /mnt/nas2/Moussa_Amplicon
#   bash AmpSeeker/prepare_ampseeker_runs.sh \
#       MACARD-1/261001_M08382_0057_000000000-MHMNW \
#       MACARD-2/261109_M08382_0059_000000000-MHMR4
#
# Note: run_AmpSeeker_pipeline.sh hard-codes its --configfile; --submit
# refuses to submit if that differs from --config.

set -euo pipefail

AMPSEEKER_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Under sbatch, BASH_SOURCE points to a spool copy; fall back to the submit dir.
if [[ ! -f "$AMPSEEKER_DIR/build_multirun_metadata.py" && -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    for d in "$SLURM_SUBMIT_DIR" "$SLURM_SUBMIT_DIR/AmpSeeker"; do
        [[ -f "$d/build_multirun_metadata.py" ]] && AMPSEEKER_DIR="$(cd "$d" && pwd)"
    done
fi
PREP_ENV="AmpSeeker-python"     # pandas + bcl2fastq
SMK_ENV="ampseeker"             # snakemake

CONFIG="config/moussa_agam.yaml"
METADATA=""
FASTQ_MODE="rebuild"
FORCE_FASTQ=""
THREADS="${SLURM_CPUS_PER_TASK:-8}"
ON_DUP="error"
SKIP_ENVS=0
SUBMIT=0
FROM_STEP=1
RUNS=()

usage() { sed -n '/^# Usage/,/^#   -h/p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }
log()  { printf '\n\033[1m[%s] %s\033[0m\n' "$(date +%T)" "$*"; }
die()  { printf '\nERROR: %s\n' "$*" >&2; exit 1; }

while [[ $# -gt 0 ]]; do
    case "$1" in
        -c|--config)     CONFIG="$2"; shift 2 ;;
        -o|--metadata)   METADATA="$2"; shift 2 ;;
        --fastq-mode)    FASTQ_MODE="$2"; shift 2 ;;
        --force-fastq)   FORCE_FASTQ="--force-fastq"; shift ;;
        --threads)       THREADS="$2"; shift 2 ;;
        --on-duplicate)  ON_DUP="$2"; shift 2 ;;
        --skip-envs)     SKIP_ENVS=1; shift ;;
        --submit)        SUBMIT=1; shift ;;
        --from-step)     FROM_STEP="$2"; shift 2 ;;
        -h|--help)       usage 0 ;;
        -*)              die "Unknown option: $1 (see --help)" ;;
        *)               RUNS+=("$(cd "$1" 2>/dev/null && pwd || echo "$1")"); shift ;;
    esac
done

[[ "$FASTQ_MODE" == rebuild || "$FASTQ_MODE" == instrument ]] || die "--fastq-mode must be rebuild or instrument"
if (( FROM_STEP <= 2 )); then
    [[ ${#RUNS[@]} -gt 0 ]] || die "Give at least one run directory (see --help)"
    for r in "${RUNS[@]}"; do
        [[ -d "$r" ]] || die "Run directory not found: $r"
        [[ -f "$r/SampleSheet.csv" || -f "$r/backup_SampleSheet.csv" ]] || die "No SampleSheet.csv in $r"
    done
fi

cd "$AMPSEEKER_DIR"
mkdir -p logs
[[ -f "$CONFIG" ]] || die "Config not found: $AMPSEEKER_DIR/$CONFIG"

# shellcheck disable=SC1091
source "$HOME/miniforge3/etc/profile.d/conda.sh"

# Read a top-level key from the YAML config (simple 'key: value' lines).
cfg() { grep -E "^$1:" "$CONFIG" | head -1 | sed -E "s/^$1:[[:space:]]*//; s/[[:space:]]+#.*$//; s/^[\"']//; s/[\"']$//"; }

[[ -n "$METADATA" ]] || METADATA="$(cfg metadata)"
[[ -n "$METADATA" ]] || die "No --metadata given and no 'metadata:' in $CONFIG"
[[ "$(cfg from-bcl)" == "False" ]] || echo "WARNING: $CONFIG has from-bcl != False; the FASTQ metadata won't be used."

log "AmpSeeker dir: $AMPSEEKER_DIR | config: $CONFIG | metadata: $METADATA | from step $FROM_STEP"

# ---------------------------------------------------------------- step 1
if (( FROM_STEP <= 1 )); then
    log "Step 1/7: convert SampleSheets"
    conda activate "$PREP_ENV"
    for r in "${RUNS[@]}"; do
        echo "-- $r"
        python convert_samplesheet.py -i "$r/SampleSheet.csv" -m resources/exampleSampleSheet.csv
    done
    conda deactivate
fi

# ---------------------------------------------------------------- step 2
if (( FROM_STEP <= 2 )); then
    log "Step 2/7: FASTQs + combined metadata ($FASTQ_MODE mode)"
    conda activate "$PREP_ENV"
    fq_args=()
    if [[ "$FASTQ_MODE" == rebuild ]]; then
        fq_args=(--rebuild-fastq-dir fastq_files)
    else
        fq_args=(--reads-dir resources/reads)
    fi
    PYTHONUNBUFFERED=1 python build_multirun_metadata.py "${RUNS[@]}" \
        -m config/example-metadata.tsv -o "$METADATA" \
        "${fq_args[@]}" $FORCE_FASTQ --threads "$THREADS" --on-duplicate "$ON_DUP"
    removed="${METADATA%.tsv}.removed_samples.tsv"
    n_removed=$(( $(wc -l < "$removed") - 1 ))
    (( n_removed == 0 )) || echo "WARNING: $n_removed sample(s) had no FASTQs — see $removed"
    conda deactivate
fi

# ---------------------------------------------------------------- step 3
if (( FROM_STEP <= 3 )); then
    log "Step 3/7: reference genome"
    fasta="$(cfg reference-fasta)"; gff="$(cfg reference-gff3)"
    if [[ -s "$fasta" && -s "$gff" ]]; then
        echo "Reference present: $fasta, $gff"
    else
        tmp="$(mktemp -d "$AMPSEEKER_DIR/resources/reference/.download.XXXX")"
        ( cd "$tmp" && bash "$AMPSEEKER_DIR/resources/reference/download-vectorbase-reference.sh" )
        mkdir -p "$(dirname "$fasta")" "$(dirname "$gff")"
        [[ -s "$fasta" ]] || mv "$tmp/resources/reference/AgamP4.fa" "$fasta"
        [[ -s "$gff"   ]] || mv "$tmp/resources/reference/AgamP4.gff" "$gff"
        rm -rf "$tmp"
        echo "Downloaded: $fasta ($(grep -c '>' "$fasta") sequences), $gff"
    fi
    # Target BED contigs must exist in the FASTA.
    targets="$(cfg targets)"
    if [[ -f "$targets" ]]; then
        missing=$(comm -23 <(cut -f1 "$targets" | sort -u) \
                           <(grep '>' "$fasta" | sed -E 's/^>([^ ]+).*/\1/' | sort -u) || true)
        [[ -z "$missing" ]] || die "Contigs in $targets not in $fasta: $missing"
        echo "All contigs in $targets found in the reference."
    fi
fi

# ---------------------------------------------------------------- step 4
if (( FROM_STEP <= 4 )); then
    log "Step 4/7: conda channel order in workflow/envs"
    python3 - <<'EOF'
import re, glob, shutil, os
os.makedirs("workflow/envs/orig_backup", exist_ok=True)
for p in sorted(glob.glob("workflow/envs/*.yaml")):
    s = open(p).read()
    m = re.search(r"channels:\n((?:[ \t]+-[ \t]*\S+\n)+)", s)
    if not m:
        continue
    chans = re.findall(r"-\s*(\S+)", m.group(1))
    # Only touch files that actually have the problem: conda-forge listed
    # after bioconda, or anaconda/defaults mixed in. Files without
    # conda-forge (e.g. bioconda-only snpeff, dranew bcl2fastq) build fine
    # and are left alone.
    bad_order = ("conda-forge" in chans and "bioconda" in chans
                 and chans.index("bioconda") < chans.index("conda-forge"))
    bad_extra = "conda-forge" in chans and any(c in chans for c in ("anaconda", "defaults"))
    if not (bad_order or bad_extra):
        print(f"ok       {p}: {chans}")
        continue
    # conda-forge first, then bioconda, then any custom channels; drop
    # anaconda/defaults (they mix incompatible builds under strict priority).
    new = ["conda-forge"] + (["bioconda"] if "bioconda" in chans else []) + \
          [c for c in chans if c not in ("conda-forge", "bioconda", "anaconda", "defaults")]
    if new == chans:
        print(f"ok       {p}: {chans}")
        continue
    bak = os.path.join("workflow/envs/orig_backup", os.path.basename(p))
    if not os.path.exists(bak):
        shutil.copy2(p, bak)
    s = s[:m.start(1)] + "".join(f"  - {c}\n" for c in new) + s[m.end(1):]
    open(p, "w").write(s)
    print(f"FIXED    {p}: {chans} -> {new}")
EOF
fi

# ---------------------------------------------------------------- step 5
conda activate "$SMK_ENV"
SMK=(snakemake --cores "$THREADS" --use-conda --configfile "$CONFIG")
if (( FROM_STEP <= 5 )) && (( ! SKIP_ENVS )); then
    log "Step 5/7: build conda envs (first time can take 30+ min)"
    "${SMK[@]}" --conda-create-envs-only
fi

# ---------------------------------------------------------------- step 6
log "Step 6/7: Snakemake dry run"
dry_log="logs/dryrun.$(date +%Y%m%d_%H%M%S).log"
if "${SMK[@]}" -n --resources mem_mb=90000 --rerun-incomplete --keep-going > "$dry_log" 2>&1; then
    awk '/^Job stats/{p=1} p{print} p&&/^total/{exit}' "$dry_log"
    echo "Dry run OK (full log: $dry_log)"
else
    tail -30 "$dry_log"
    die "Dry run failed — see $dry_log"
fi

# ---------------------------------------------------------------- step 7
if (( SUBMIT )); then
    log "Step 7/7: submit pipeline"
    run_cfg="$(grep -oE -- '--configfile[[:space:]]+[^[:space:]\\]+' run_AmpSeeker_pipeline.sh | tail -1 | awk '{print $2}')"
    [[ "$run_cfg" == "$CONFIG" ]] || die "run_AmpSeeker_pipeline.sh uses '$run_cfg' but you prepared '$CONFIG' — edit it first."
    command -v sbatch >/dev/null || die "sbatch not found"
    sbatch run_AmpSeeker_pipeline.sh
else
    log "Done. Submit with:  cd $AMPSEEKER_DIR && sbatch run_AmpSeeker_pipeline.sh"
fi
