#!/usr/bin/env python3
# Script: build_multirun_metadata.py
# Purpose: For each of several Illumina run folders, reuse on-instrument FASTQs
#          if present (else run bcl2fastq), clean/rename the FASTQ files, and
#          build one combined metadata.tsv across all runs for AmpSeeker's
#          FASTQ-input mode (AmpSeeker can't take multiple BCL folders at once,
#          only a single metadata sheet with fq1/fq2 paths).
# Input:   One or more run directories, each containing a SampleSheet.csv
#          (or backup_SampleSheet.csv) and either raw BCL data or an
#          on-instrument Alignment_*/*/Fastq/ folder already populated.
#          A model TSV (config/example-metadata.tsv) defining the target
#          metadata columns.
# Output:  resources/reads/<run_name>/<sample>_1.fastq.gz + _2.fastq.gz
#          (symlinks when reusing existing FASTQs, real files when bcl2fastq
#          is run) and one combined metadata TSV.
# Dependencies: pandas, bcl2fastq (only needed for runs without existing FASTQs)
# python3 build_multirun_metadata.py   resources/251211_M08382_0055_000000000-M783V/ --model config/example-metadata.tsv  --reads-dir ressources/cleaned_reads  --force-fastq -o config/metadata.tsv

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

# Illumina raw-sheet column -> intermediate name (before mapping to the model)
COLUMN_RENAME = {
    "Sample_ID": "sample_id_raw",
    "Sample_Name": "sample_name",
    "Sample_Plate": "plate",
    "Sample_Well": "well",
}

FASTQ_PATTERN = re.compile(
    r"^(?P<sample>.+)_S\d+_L\d{3}_(?P<read>R[12]|I[12])_001\.fastq\.gz$"
)


# ---------------------------------------------------------------------------
# SampleSheet reading
# ---------------------------------------------------------------------------

def find_raw_samplesheet(run_dir: Path) -> Path:
    """Prefer backup_SampleSheet.csv (the untouched original) over
    SampleSheet.csv, since an earlier conversion step may have overwritten
    the latter with AmpSeeker-formatted column names bcl2fastq won't
    recognise.
    """
   
    plain = run_dir / "SampleSheet.csv"
    if plain.exists():
        return plain
    backup = run_dir / "backup_SampleSheet.csv"
    if backup.exists():
            return backup
    raise FileNotFoundError(
        f"No SampleSheet.csv or backup_SampleSheet.csv found in {run_dir}"
    )


def read_samplesheet(path: Path):
    """Read an Illumina SampleSheet and return (header_lines, DataFrame)."""
    with open(path, encoding="utf-8-sig") as f:
        lines = f.readlines()

    data_start = None
    for i, line in enumerate(lines):
        if line.strip().lower().startswith("[data"):
            data_start = i + 1
            break
    if data_start is None:
        raise ValueError(f"Could not find a '[Data]' section in {path}")

    header = lines[:data_start]
    df = pd.read_csv(path, skiprows=data_start, encoding="utf-8-sig")
    return header, df


def read_model_columns(model_path: Path) -> list:
    """Read the target column list from the tab-separated example metadata."""
    df = pd.read_csv(model_path, sep="\t", nrows=0)
    return df.columns.tolist()


def find_col(columns, target_lower: str):
    for c in columns:
        if c.lower() == target_lower:
            return c
    return None


# ---------------------------------------------------------------------------
# FASTQ discovery / generation / staging
# ---------------------------------------------------------------------------

def find_existing_fastq_dir(run_dir: Path) -> Path | None:
    """Look for an on-instrument Alignment_*/*/Fastq/ folder already
    populated with FASTQs (e.g. from a MiSeq GenerateFASTQ run). Returns the
    most recent match, or None if the run needs bcl2fastq run against it.
    """
    candidates = sorted(run_dir.glob("Alignment_*/*/Fastq"))
    for d in reversed(candidates):  # most recent timestamp folder last
        if any(d.glob("*.fastq.gz")):
            return d
    # Fallback: a looser search in case the folder layout differs slightly.
    hits = list(run_dir.rglob("*_R1_*.fastq.gz"))
    if hits:
        return hits[0].parent
    return None


def run_bcl2fastq(run_dir: Path, samplesheet: Path, out_dir: Path,
                   threads: int, no_lane_splitting: bool, dry_run: bool):
    cmd = [
        "bcl2fastq",
        "--runfolder-dir", str(run_dir),
        "--sample-sheet", str(samplesheet),
        "--output-dir", str(out_dir),
        "-r", str(threads), "-p", str(threads), "-w", str(threads),
    ]
    if no_lane_splitting:
        cmd.append("--no-lane-splitting")

    print(f"  $ {' '.join(cmd)}")
    if dry_run:
        return

    if shutil.which("bcl2fastq") is None:
        raise RuntimeError(
            "bcl2fastq not found on PATH. On HPC it's usually provided via "
            "a module rather than conda — try `module avail bcl2fastq` "
            "and `module load <the module name>`."
        )
    out_dir.mkdir(parents=True, exist_ok=True)
    subprocess.run(cmd, check=True)


def cleaned_name(fastq_filename: str) -> str | None:
    """Map a raw bcl2fastq-style filename to the AmpSeeker-friendly form,
    or None if it's an index read (I1/I2) that should be dropped.
    e.g. 'SAMPLE_S1_L001_R1_001.fastq.gz' -> 'SAMPLE_1.fastq.gz'
    """
    m = FASTQ_PATTERN.match(fastq_filename)
    if not m:
        return None  # doesn't match the expected bcl2fastq naming; skip it
    read = m.group("read")
    if read.startswith("I"):
        return None
    read_num = read[1]  # 'R1' -> '1', 'R2' -> '2'
    return f"{m.group('sample')}_{read_num}.fastq.gz"


def stage_fastqs(source_dir: Path, dest_dir: Path, symlink: bool,
                  force: bool, dry_run: bool) -> dict:
    """Populate dest_dir with cleanly-named R1/R2 fastqs (dropping I1/I2),
    either by renaming in place (source_dir == dest_dir, e.g. right after
    bcl2fastq wrote there) or by symlinking from source_dir into dest_dir
    (when reusing an existing, read-only run folder — never renames files
    inside the original run directory).

    Returns {sample: [fq1_path, fq2_path]}.
    """
    if dest_dir.exists() and any(dest_dir.glob("*_1.fastq.gz")) and not force:
        print(f"  Reusing already-staged FASTQs in {dest_dir} (use --force-fastq to redo)")
        return collect_fastq_dict(dest_dir)

    if dry_run:
        print(f"  [dry-run] would stage FASTQs from {source_dir} -> {dest_dir}")
        return {}

    dest_dir.mkdir(parents=True, exist_ok=True)
    raw_files = sorted(source_dir.glob("*.fastq.gz"))
    if not raw_files:
        raise FileNotFoundError(f"No .fastq.gz files found in {source_dir}")

    same_dir = source_dir.resolve() == dest_dir.resolve()
    for raw in raw_files:
        new_name = cleaned_name(raw.name)
        if new_name is None:
            if same_dir:
                raw.unlink()  # drop I1/I2 (or unrecognised) files we created ourselves
            continue  # never delete files in a source dir we don't own
        target = dest_dir / new_name
        if same_dir:
            if raw.name != new_name:
                raw.rename(target)
        else:
            if target.exists() or target.is_symlink():
                target.unlink()
            target.symlink_to(raw.resolve())

    return collect_fastq_dict(dest_dir)


def collect_fastq_dict(fastq_dir: Path) -> dict:
    fastq_dict = {}
    for f in sorted(fastq_dir.glob("*_[12].fastq.gz")):
        sample = f.name.rsplit("_", 1)[0]
        fastq_dict.setdefault(sample, []).append(str(f))
    for sample, files in fastq_dict.items():
        files.sort()  # ensures _1 before _2
    return fastq_dict


# ---------------------------------------------------------------------------
# Per-run metadata
# ---------------------------------------------------------------------------

def build_run_metadata(run_dir: Path, fastq_dict: dict, ref_cols: list) -> pd.DataFrame:
    samplesheet_path = find_raw_samplesheet(run_dir)
    _, df = read_samplesheet(samplesheet_path)
    df = df.rename(columns=COLUMN_RENAME)

    # sample_id comes from Sample_Name (matches the fastq-dict key produced
    # by cleaned_name(), which uses bcl2fastq's Sample_Name-derived prefix).
    name_col = find_col(df.columns, "sample_id")
    if name_col is None:
        raise ValueError(
            f"{samplesheet_path} has no Sample_Name column to use as sample_id."
        )
    df["sample_id"] = df[name_col]

    for col in ["taxon", "location", "country", "latitude", "longitude"]:
        if col not in df.columns:
            df[col] = ""

    if "well" in df.columns:
        extracted = df["well"].astype(str).str.extract(r"([A-Za-z]+)(\d+)")
        df["well_letter"], df["well_number"] = extracted[0], extracted[1]
    else:
        df["well_letter"], df["well_number"] = "", ""

    df["fq1"] = ""
    df["fq2"] = ""
    missing = []
    for idx, row in df.iterrows():
        pair = fastq_dict.get(row["sample_id"])
        if pair and len(pair) == 2:
            df.at[idx, "fq1"], df.at[idx, "fq2"] = pair[0], pair[1]
        else:
            missing.append(row["sample_id"])

    if missing:
        print(f"  Warning: {len(missing)} sample(s) in {samplesheet_path.name} "
              f"had no matching FASTQ pair: {missing}")

    for col in ref_cols:
        if col.lower() == "sample_id":
            continue
        if col not in df.columns:
            df[col] = ""
    print(f"before final df: \n {df.head()}")
    out_cols = ["sample_id" if c.lower() == "sample_id" else c for c in ref_cols]
    out_cols += ["fq1", "fq2", "plate_name"]

    final_df = df[out_cols].copy()

    # Identify rows with missing or empty fq1/fq2
    remove_mask = (
        final_df["fq1"].isna() |
        final_df["fq2"].isna() |
        final_df["fq1"].astype(str).str.strip().eq("") |
        final_df["fq2"].astype(str).str.strip().eq("")
    )

    # Save removed samples before filtering
    removed_df = final_df.loc[remove_mask, ["sample_id", "fq1", "fq2", "plate_name"]].copy()

    # Print removed sample IDs
    print(f"Removing {len(removed_df)} samples:")
    for sample_id in removed_df["sample_id"]:
        print(sample_id)

    # Save removal log
    removed_df.to_csv("removed_samples.log", sep="\t", index=False)

    # Keep valid rows
    final_df = final_df.loc[~remove_mask].copy()

    print(f"Remaining samples: {len(final_df)}")
    print(f"Removal log saved to: removed_samples.log")

    return final_df


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------

def process_run(run_dir: Path, reads_dir: Path, ref_cols: list, threads: int,
                 no_lane_splitting: bool, force_fastq: bool, dry_run: bool,
                 rebuild_dirname: str | None) -> pd.DataFrame:
    run_name = run_dir.name
    print(f"\n=== {run_name} ===")

    if rebuild_dirname:
        # Always rebuild from BCL, regardless of what already exists, and
        # write into its own folder inside the run directory so any
        # previous FASTQs (on-instrument or a prior staged run) are left
        # untouched for comparison.
        dest_dir = run_dir / rebuild_dirname
        print(f"  --rebuild-fastq-dir set — forcing bcl2fastq into {dest_dir} "
              f"(ignoring any existing FASTQs)")
        samplesheet_path = find_raw_samplesheet(run_dir)
        run_bcl2fastq(run_dir, samplesheet_path, dest_dir, threads,
                      no_lane_splitting, dry_run)
        fastq_dict = stage_fastqs(dest_dir, dest_dir, symlink=False,
                                   force=True, dry_run=dry_run)
    else:
        dest_dir = reads_dir / run_name
        existing = None if force_fastq else find_existing_fastq_dir(run_dir)
        if existing:
            print(f"  Found existing on-instrument FASTQs: {existing}")
            fastq_dict = stage_fastqs(existing, dest_dir, symlink=True,
                                       force=force_fastq, dry_run=dry_run)
        else:
            print(f"  No existing FASTQs found — running bcl2fastq")
            samplesheet_path = find_raw_samplesheet(run_dir)
            run_bcl2fastq(run_dir, samplesheet_path, dest_dir, threads,
                          no_lane_splitting, dry_run)
            fastq_dict = stage_fastqs(dest_dir, dest_dir, symlink=False,
                                       force=force_fastq, dry_run=dry_run)
    if dry_run:
        return pd.DataFrame(columns=ref_cols)

    return build_run_metadata(run_dir, fastq_dict, ref_cols)


def main():
    parser = argparse.ArgumentParser(
        description="Convert multiple Illumina run folders to FASTQ (reusing "
                     "on-instrument FASTQs where available) and build one "
                     "combined metadata TSV for AmpSeeker."
    )
    parser.add_argument("runs", nargs="+", type=Path,
                         help="One or more run directories.")
    parser.add_argument("-m", "--model", type=Path,
                         default=Path("config/example-metadata.tsv"),
                         help="Model TSV defining the target metadata columns "
                              "(default: config/example-metadata.tsv).")
    parser.add_argument("--reads-dir", type=Path, default=Path("resources/reads"),
                         help="Base directory to stage cleaned FASTQs into; "
                              "each run gets its own subfolder (default: resources/reads).")
    parser.add_argument("-o", "--output", type=Path, default=Path("config/metadata.tsv"),
                         help="Path to write the combined metadata TSV "
                              "(default: config/metadata.tsv).")
    parser.add_argument("--threads", type=int, default=4,
                         help="Threads passed to bcl2fastq -r/-p/-w (default: 4).")
    parser.add_argument("--no-lane-splitting", action="store_true",
                         help="Pass --no-lane-splitting to bcl2fastq.")
    parser.add_argument("--force-fastq", action="store_true",
                         help="Re-run bcl2fastq / re-stage FASTQs even if "
                              "already present, ignoring any existing "
                              "on-instrument or previously-staged output. "
                              "Writes into the same resources/reads/<run>/ "
                              "location, overwriting what's there.")
    parser.add_argument("--rebuild-fastq-dir", metavar="DIRNAME", default=None,
                         help="Always rebuild FASTQs from BCL with bcl2fastq, "
                              "even if on-instrument FASTQs already exist, "
                              "and write them into <run_dir>/DIRNAME (e.g. "
                              "'Fastq_new') instead of resources/reads/. "
                              "Existing FASTQs (on-instrument or previously "
                              "staged) are left completely untouched. "
                              "Applies to every run passed in this invocation.")
    parser.add_argument("--on-duplicate", choices=["warn", "error", "suffix"],
                         default="warn",
                         help="What to do if the same sample_id appears in "
                              "more than one run: warn (default, keep as-is), "
                              "error (abort), or suffix (append the run name "
                              "to disambiguate).")
    parser.add_argument("-n", "--dry-run", action="store_true",
                         help="Print what would be done without running "
                              "bcl2fastq, staging files, or writing output.")
    args = parser.parse_args()

    for run_dir in args.runs:
        if not run_dir.is_dir():
            sys.exit(f"Run directory not found: {run_dir}")
    if not args.model.exists():
        sys.exit(f"Model file not found: {args.model}")

    ref_cols = read_model_columns(args.model)
    if not any(c.lower() == "sample_id" for c in ref_cols):
        sys.exit(f"Model file {args.model} has no sample_id column.")

    all_dfs = []
    for run_dir in args.runs:
        try:
            df = process_run(run_dir, args.reads_dir, ref_cols, args.threads,
                              args.no_lane_splitting, args.force_fastq, args.dry_run,
                              args.rebuild_fastq_dir)
        except (FileNotFoundError, ValueError, RuntimeError,
                subprocess.CalledProcessError) as e:
            sys.exit(f"\nFailed on run {run_dir.name}: {e}")
        df["_source_run"] = run_dir.name
        all_dfs.append(df)

    if args.dry_run:
        print("\n[dry-run] Stopping before writing combined metadata.")
        return

    combined = pd.concat(all_dfs, ignore_index=True)

    dupes = combined["sample_id"][combined["sample_id"].duplicated(keep=False)]
    if not dupes.empty:
        dupe_report = (
            combined.loc[dupes.index, ["sample_id", "_source_run"]]
            .groupby("sample_id")["_source_run"].apply(list)
        )
        msg = f"{len(dupe_report)} sample_id(s) appear in more than one run:\n{dupe_report}"
        if args.on_duplicate == "error":
            sys.exit(f"Conversion failed: {msg}")
        elif args.on_duplicate == "suffix":
            combined.loc[dupes.index, "sample_id"] = (
                combined.loc[dupes.index, "sample_id"] + "_" + combined.loc[dupes.index, "_source_run"]
            )
            print(f"Note: disambiguated duplicate sample_id(s) with run suffix:\n{dupe_report}")
        else:
            print(f"Warning: {msg}")

    combined = combined.drop(columns=["_source_run"])

    args.output.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(args.output, sep="\t", index=False, lineterminator="\n")
    print(f"final metadata:\n {combined.head()}")

    print(f"\nWrote combined metadata for {len(combined)} sample(s) across "
          f"{len(args.runs)} run(s) to: {args.output}")


if __name__ == "__main__":
    main()
