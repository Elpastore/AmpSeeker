#!/usr/bin/env python3
# Script: build_multirun_metadata.py
# Purpose: For each of several Illumina run folders, get per-sample FASTQs
#          (reuse on-instrument FASTQs, or run bcl2fastq), give them clean
#          AmpSeeker names (<sample_id>_1/_2.fastq.gz), and build one
#          combined metadata.tsv across all runs for AmpSeeker's
#          FASTQ-input mode (from-bcl: False).
# Input:   One or more run directories, each containing SampleSheet.csv
#          (raw Illumina or already converted by convert_samplesheet.py;
#          backup_SampleSheet.csv is used if SampleSheet.csv is missing).
#          A model TSV (config/example-metadata.tsv) defining the first,
#          mandatory metadata columns.
# Output:  FASTQs staged under --reads-dir/<run_name>/ (symlinks to
#          on-instrument FASTQs, or real files from bcl2fastq), or under
#          <run_dir>/<DIRNAME> with --rebuild-fastq-dir; one combined
#          metadata TSV with absolute fq1/fq2 paths; and a
#          <output>.removed_samples.tsv log of samples without usable FASTQs.
# Dependencies: pandas, bcl2fastq (only for runs that need demultiplexing)
#
# Example:
#   conda activate AmpSeeker-python
#   cd /mnt/nas2/Moussa_Amplicon
#   python AmpSeeker/build_multirun_metadata.py MACARD-1/<run> MACARD-2/<run> \
#       -m AmpSeeker/config/example-metadata.tsv \
#       --rebuild-fastq-dir fastq_files -o AmpSeeker/config/metadata.tsv

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from convert_samplesheet import (  # noqa: E402
    DROP_EXTRA,
    METADATA_COLUMNS,
    normalise_columns,
    read_samplesheet,
    snake,
)

FASTQ_PATTERN = re.compile(
    r"^(?P<sample>.+)_S(?P<snum>\d+)(?:_L\d{3})?_(?P<read>R[12]|I[12])_001\.fastq\.gz$"
)
STAGED_PATTERN = re.compile(r"^(?P<sample>.+)_(?P<read>[12])\.fastq\.gz$")


# ---------------------------------------------------------------------------
# SampleSheet reading
# ---------------------------------------------------------------------------

def find_samplesheet(run_dir: Path) -> Path:
    """SampleSheet.csv may be the raw Illumina sheet or the AmpSeeker
    conversion — both are normalised the same way, so either works.
    backup_SampleSheet.csv is only a fallback."""
    for name in ("SampleSheet.csv", "backup_SampleSheet.csv"):
        p = run_dir / name
        if p.exists():
            return p
    raise FileNotFoundError(
        f"No SampleSheet.csv or backup_SampleSheet.csv found in {run_dir}"
    )


def load_sheet(run_dir: Path):
    """Return (header_lines, normalised DataFrame, sheet path).

    Columns are lowercased/renamed exactly as convert_samplesheet.py does
    (Country -> country, Sample_Plate -> plate_name, ...). sample_name is
    set to sample_id (the original kept as sample_name_original) because
    that's what bcl2fastq names FASTQ files after.
    """
    path = find_samplesheet(run_dir)
    header, df = read_samplesheet(path)
    df = normalise_columns(df)
    if "sample_id" not in df.columns:
        raise ValueError(f"{path} has no Sample_ID / sample_id column.")
    if "sample_name" in df.columns and "sample_name_original" not in df.columns:
        if df["sample_name"].ne(df["sample_id"]).any():
            df["sample_name_original"] = df["sample_name"]
    df["sample_name"] = df["sample_id"]
    dup = df["sample_id"].duplicated(keep=False)
    if dup.any():
        raise ValueError(f"Duplicate sample_id(s) in {path}: "
                         f"{sorted(df.loc[dup, 'sample_id'].unique())}")
    return header, df, path


def read_model_columns(model_path: Path) -> list:
    """Read the target column list from the tab-separated example metadata."""
    cols = pd.read_csv(model_path, sep="\t", nrows=0).columns.tolist()
    return ["sample_id" if snake(c) == "sample_id" else c for c in cols]


# ---------------------------------------------------------------------------
# FASTQ discovery / generation / staging
# ---------------------------------------------------------------------------

def find_existing_fastq_dir(run_dir: Path) -> Path | None:
    """Most recent on-instrument Alignment_*/<timestamp>/Fastq/ folder
    containing non-empty FASTQs, or None."""
    candidates = sorted(run_dir.glob("Alignment_*/*/Fastq"))
    for d in reversed(candidates):
        if any(f.stat().st_size > 0 for f in d.glob("*_R1_*.fastq.gz")):
            return d
    return None


def write_bcl2fastq_sheet(header: list, sheet: pd.DataFrame, path: Path):
    """Write a minimal sheet for bcl2fastq: the original [Header]/[Reads]/
    [Settings] sections + [Data] with only Sample_ID/Sample_Name/index/index2.

    - Sample_Name = sample_id, so output files are <sample_id>_S<n>_... .
    - No Sample_Project column, so output is flat (no project subfolder).
    - Index-read FASTQs are not requested (they'd be deleted anyway).
    """
    out_header = []
    for line in header:
        stripped = line.rstrip("\r\n").rstrip(",")
        if stripped.lower().startswith("createfastqforindexreads"):
            stripped = "CreateFastqForIndexReads,0"
        if stripped.lower().startswith("[data"):
            stripped = "[Data]"
        out_header.append(stripped + "\n")

    cols = {"Sample_ID": sheet["sample_id"], "Sample_Name": sheet["sample_id"]}
    for c in ("index", "index2"):
        if c in sheet.columns:
            cols[c] = sheet[c]
    data = pd.DataFrame(cols)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as f:
        f.writelines(out_header)
        data.to_csv(f, index=False, lineterminator="\n")


def run_bcl2fastq(run_dir: Path, samplesheet: Path, out_dir: Path,
                  threads: int, no_lane_splitting: bool, dry_run: bool,
                  log_name: str = "bcl2fastq.log"):
    io_threads = max(1, min(4, threads))
    cmd = [
        "bcl2fastq",
        "--runfolder-dir", str(run_dir),
        "--sample-sheet", str(samplesheet),
        "--output-dir", str(out_dir),
        "-r", str(io_threads), "-p", str(threads), "-w", str(io_threads),
    ]
    if no_lane_splitting:
        cmd.append("--no-lane-splitting")

    print(f"  $ {' '.join(cmd)}")
    if dry_run:
        return

    if shutil.which("bcl2fastq") is None:
        raise RuntimeError(
            "bcl2fastq not found on PATH. Activate the env that provides it: "
            "`conda activate AmpSeeker-python`."
        )
    out_dir.mkdir(parents=True, exist_ok=True)
    log = out_dir / log_name
    print(f"  bcl2fastq log: {log}", flush=True)
    with open(log, "w") as fh:
        res = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT)
    if res.returncode != 0:
        tail = log.read_text().splitlines()[-8:]
        hint = ""
        if "double free" in "\n".join(tail):
            hint = ("\nThis bcl2fastq crash usually means it ran out of open file "
                    "descriptors (`ulimit -n`); lower --samples-per-pass.")
        raise RuntimeError(f"bcl2fastq failed (exit {res.returncode}). Last log lines:\n"
                           + "\n".join(tail) + hint)


def name_key(s: str) -> str:
    """Instrument software rewrites '_' as '-' in FASTQ names; compare
    names ignoring case and '_' vs '-'."""
    return str(s).lower().replace("_", "-")


def resolve_sample(prefix: str, snum: int, sheet: pd.DataFrame) -> str | None:
    """Map a raw FASTQ (name prefix + S number) to a sample_id.

    1. Exact sample_id match (bcl2fastq output from our generated sheet).
    2. S number = row order in the sheet (Illumina convention), accepted
       only if the name also matches sample_id or the original
       Sample_Name (modulo '_'/'-' and case).
    """
    ids = sheet["sample_id"].tolist()
    if prefix in ids:
        return prefix
    if 1 <= snum <= len(sheet):
        row = sheet.iloc[snum - 1]
        names = {name_key(row["sample_id"])}
        if "sample_name_original" in sheet.columns:
            names.add(name_key(row["sample_name_original"]))
        if name_key(prefix) in names:
            return row["sample_id"]
    return None


def staged_complete(dest_dir: Path, sheet: pd.DataFrame) -> bool:
    """True only if dest_dir has a non-empty staged pair for EVERY sample in
    the sheet (a partial folder, e.g. from a crashed bcl2fastq, is redone)."""
    if not dest_dir.is_dir():
        return False
    d = collect_fastq_dict(dest_dir)
    missing = [s for s in sheet["sample_id"] if s not in d]
    if missing and d:
        print(f"  {dest_dir} has staged FASTQs but {len(missing)} sample(s) are "
              f"missing/empty, e.g. {missing[:5]} — redoing this run.")
    return not missing


def clear_previous_output(dest_dir: Path):
    """Delete FASTQs a previous (possibly crashed) bcl2fastq/staging left in
    dest_dir. Stale raw files carry different _S<n> numbers than a new run,
    so without this they get renamed over the fresh FASTQs. Only files
    matching bcl2fastq/staged naming are touched."""
    if not dest_dir.is_dir():
        return
    n = 0
    for f in list(dest_dir.rglob("*.fastq.gz")):
        if (FASTQ_PATTERN.match(f.name) or STAGED_PATTERN.match(f.name)
                or f.name.startswith("Undetermined_")):
            f.unlink()
            n += 1
    if n:
        print(f"  Removed {n} FASTQ file(s) left in {dest_dir} by a previous run")


def stage_fastqs(source_dir: Path, dest_dir: Path, sheet: pd.DataFrame,
                 dry_run: bool) -> dict:
    """Put cleanly named <sample_id>_1/_2.fastq.gz files in dest_dir.

    - source_dir == dest_dir (bcl2fastq just wrote there): rename in
      place and delete index-read / unrecognised files we produced.
    - otherwise (on-instrument folder): symlink, never touching the
      original run folder.
    Searches source_dir recursively (handles Sample_Project subfolders)
    but ignores Undetermined reads.
    """
    if dry_run:
        print(f"  [dry-run] would stage FASTQs from {source_dir} -> {dest_dir}")
        return {}

    dest_dir.mkdir(parents=True, exist_ok=True)
    same_dir = source_dir.resolve() == dest_dir.resolve()
    raw_files = sorted(p for p in source_dir.rglob("*.fastq.gz")
                       if not p.name.startswith("Undetermined")
                       and FASTQ_PATTERN.match(p.name))
    if not raw_files and not collect_fastq_dict(dest_dir):
        raise FileNotFoundError(f"No raw .fastq.gz files found in {source_dir}")

    unmatched = set()
    plan = {}
    for raw in raw_files:
        m = FASTQ_PATTERN.match(raw.name)
        read = m.group("read")
        if read.startswith("I"):
            continue
        sample = resolve_sample(m.group("sample"), int(m.group("snum")), sheet)
        if sample is None:
            unmatched.add(raw.name)
            continue
        target = dest_dir / f"{sample}_{read[1]}.fastq.gz"
        if target in plan:
            raise ValueError(
                f"Two FASTQs map to {target.name}: {plan[target].name} and {raw.name}. "
                f"Stale files from an earlier run are probably in {source_dir}; "
                "re-run with --force-fastq (clears old bcl2fastq output) or remove them.")
        plan[target] = raw

    n_staged = 0
    for raw in raw_files:
        if FASTQ_PATTERN.match(raw.name).group("read").startswith("I") and same_dir:
            raw.unlink()
    for target, raw in plan.items():
        if target.exists() or target.is_symlink():
            target.unlink()
        if same_dir:
            raw.rename(target)
        else:
            target.symlink_to(raw.resolve())
        n_staged += 1

    if same_dir:
        # Remove now-empty project subfolders bcl2fastq may have created.
        for d in sorted((p for p in dest_dir.rglob("*") if p.is_dir()), reverse=True):
            if d.name not in ("Stats", "Reports") and not any(d.iterdir()):
                d.rmdir()

    print(f"  Staged {n_staged} FASTQ file(s) into {dest_dir}")
    if unmatched:
        print(f"  Warning: {len(unmatched)} FASTQ file(s) matched no sample in the "
              f"sheet and were skipped, e.g. {sorted(unmatched)[:5]}")
    return collect_fastq_dict(dest_dir)


def collect_fastq_dict(fastq_dir: Path) -> dict:
    """{sample_id: [fq1, fq2]} for non-empty staged files (absolute paths)."""
    pairs = {}
    for f in fastq_dir.glob("*_[12].fastq.gz"):
        m = STAGED_PATTERN.match(f.name)
        if not m or not f.exists() or f.stat().st_size == 0:
            continue
        pairs.setdefault(m.group("sample"), {})[m.group("read")] = str(f.absolute())
    return {s: [r["1"], r["2"]] for s, r in pairs.items() if "1" in r and "2" in r}


# ---------------------------------------------------------------------------
# Per-run metadata
# ---------------------------------------------------------------------------

def build_run_metadata(sheet: pd.DataFrame, fastq_dict: dict, ref_cols: list,
                       run_name: str, keep_extra: bool):
    """Return (kept_df, removed_df) for one run."""
    df = sheet.copy()

    for col in METADATA_COLUMNS:
        if col not in df.columns:
            df[col] = ""

    if "plate" not in df.columns or df["plate"].astype(str).str.strip().eq("").all():
        df["plate"] = df["plate_name"] if "plate_name" in df.columns else ""

    if "well" in df.columns:
        ext = df["well"].astype(str).str.extract(r"^\s*([A-Za-z]+)0*(\d+)\s*$")
        df["well_letter"] = ext[0].fillna("").str.upper()
        df["well_number"] = ext[1].fillna("")
    for col in ("well_letter", "well_number"):
        if col not in df.columns:
            df[col] = ""

    for col in ref_cols:
        if col not in df.columns:
            df[col] = ""

    df["run_id"] = run_name
    df["fq1"] = df["sample_id"].map(lambda s: fastq_dict.get(s, ["", ""])[0])
    df["fq2"] = df["sample_id"].map(lambda s: fastq_dict.get(s, ["", ""])[1])

    out_cols = list(ref_cols) + ["fq1", "fq2", "run_id"]
    if keep_extra:
        skip = set(out_cols) | DROP_EXTRA | {"sample_name"}
        out_cols += [c for c in df.columns if c not in skip
                     and df[c].astype(str).str.strip().ne("").any()]
    df = df[out_cols]

    remove_mask = df["fq1"].eq("") | df["fq2"].eq("")
    removed = df.loc[remove_mask, ["sample_id", "run_id", "fq1", "fq2"]].copy()
    removed["reason"] = "no non-empty FASTQ pair found"
    kept = df.loc[~remove_mask].copy()
    print(f"  {len(kept)} sample(s) with FASTQs, {len(removed)} removed"
          + (f": {removed['sample_id'].tolist()}" if len(removed) else ""))
    return kept, removed


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------

def open_file_budget() -> int:
    """Raise the soft open-file limit to the hard limit and return it."""
    try:
        import resource
        soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
        if hard == resource.RLIM_INFINITY or hard > soft:
            target = hard if hard != resource.RLIM_INFINITY else max(soft, 65536)
            resource.setrlimit(resource.RLIMIT_NOFILE, (target, hard))
            soft = target
        return soft
    except (ImportError, ValueError, OSError):
        return 1024


def auto_samples_per_pass(n_samples: int) -> int:
    """bcl2fastq keeps one open file per sample per read (2 per sample
    here) plus a few hundred for BCL/filter input, and when it runs out of
    descriptors it crashes with 'free(): double free detected' instead of
    a clear error (observed with a 500-file limit: crash after ~175
    output files). Split the samples into passes that fit."""
    limit = open_file_budget()
    if limit >= 2 * n_samples + 1024:
        return n_samples
    return max(8, (limit - 372) // 2)


def demultiplex(run_dir, header, sheet, dest_dir, threads, no_lane_splitting,
                dry_run, samples_per_pass: int = 0):
    if not dry_run:
        clear_previous_output(dest_dir)
    n = len(sheet)
    per_pass = samples_per_pass or auto_samples_per_pass(n)
    chunks = [sheet.iloc[i:i + per_pass] for i in range(0, n, per_pass)]
    if len(chunks) > 1:
        print(f"  Open-file limit is {open_file_budget()}: demultiplexing in "
              f"{len(chunks)} passes of <= {per_pass} samples "
              f"(override with --samples-per-pass)")
    for k, chunk in enumerate(chunks, 1):
        suffix = f".part{k}" if len(chunks) > 1 else ""
        bcl_sheet = dest_dir / f"SampleSheet.bcl2fastq{suffix}.csv"
        if not dry_run:
            write_bcl2fastq_sheet(header, chunk, bcl_sheet)
        if len(chunks) > 1:
            print(f"  pass {k}/{len(chunks)}: {len(chunk)} samples")
        run_bcl2fastq(run_dir, bcl_sheet, dest_dir, threads, no_lane_splitting,
                      dry_run, log_name=f"bcl2fastq{suffix}.log")
        if not dry_run and len(chunks) > 1:
            # Each pass rewrites Stats/Reports and Undetermined; keep each
            # pass's reports apart so nothing is silently overwritten.
            for sub in ("Stats", "Reports"):
                src = dest_dir / sub
                if src.is_dir():
                    dst = dest_dir / f"{sub}{suffix}"
                    if dst.exists():
                        shutil.rmtree(dst)
                    src.rename(dst)
    return stage_fastqs(dest_dir, dest_dir, sheet, dry_run)


def process_run(run_dir: Path, reads_dir: Path, ref_cols: list, threads: int,
                no_lane_splitting: bool, force_fastq: bool, dry_run: bool,
                rebuild_dirname: str | None, keep_extra: bool,
                samples_per_pass: int = 0):
    run_dir = run_dir.resolve()
    run_name = run_dir.name
    print(f"\n=== {run_name} ===", flush=True)
    header, sheet, sheet_path = load_sheet(run_dir)
    print(f"  Sample sheet: {sheet_path} ({len(sheet)} samples)")

    if rebuild_dirname:
        dest_dir = run_dir / rebuild_dirname.strip("/")
        if staged_complete(dest_dir, sheet) and not force_fastq:
            print(f"  Reusing FASTQs already rebuilt in {dest_dir} "
                  f"(add --force-fastq to re-run bcl2fastq)")
            fastq_dict = collect_fastq_dict(dest_dir)
        else:
            print(f"  Rebuilding FASTQs from BCL with bcl2fastq into {dest_dir}")
            fastq_dict = demultiplex(run_dir, header, sheet, dest_dir, threads,
                                     no_lane_splitting, dry_run, samples_per_pass)
    else:
        dest_dir = reads_dir.resolve() / run_name
        existing = None if force_fastq else find_existing_fastq_dir(run_dir)
        if staged_complete(dest_dir, sheet) and not force_fastq:
            print(f"  Reusing already-staged FASTQs in {dest_dir} (use --force-fastq to redo)")
            fastq_dict = collect_fastq_dict(dest_dir)
        elif existing:
            print(f"  Linking on-instrument FASTQs from {existing}")
            fastq_dict = stage_fastqs(existing, dest_dir, sheet, dry_run)
        else:
            print("  No usable on-instrument FASTQs (or --force-fastq) — running bcl2fastq")
            fastq_dict = demultiplex(run_dir, header, sheet, dest_dir, threads,
                                     no_lane_splitting, dry_run, samples_per_pass)

    if dry_run:
        return pd.DataFrame(columns=ref_cols), pd.DataFrame()
    return build_run_metadata(sheet, fastq_dict, ref_cols, run_name, keep_extra)


def main():
    parser = argparse.ArgumentParser(
        description="Get FASTQs for multiple Illumina run folders and build one "
                    "combined metadata TSV for AmpSeeker (from-bcl: False)."
    )
    parser.add_argument("runs", nargs="+", type=Path, help="One or more run directories.")
    parser.add_argument("-m", "--model", type=Path,
                        default=Path(__file__).resolve().parent / "config/example-metadata.tsv",
                        help="Model TSV defining the mandatory metadata columns "
                             "(default: config/example-metadata.tsv next to this script).")
    parser.add_argument("--reads-dir", type=Path, default=Path("resources/reads"),
                        help="Base directory for staged FASTQs when --rebuild-fastq-dir "
                             "is not used; each run gets a subfolder (default: resources/reads).")
    parser.add_argument("-o", "--output", type=Path,
                        default=Path(__file__).resolve().parent / "config/metadata.tsv",
                        help="Combined metadata TSV (default: config/metadata.tsv next to this script).")
    parser.add_argument("--threads", type=int, default=4,
                        help="bcl2fastq processing threads (-p); -r/-w use min(4, threads). Default 4.")
    parser.add_argument("--no-lane-splitting", action="store_true",
                        help="Pass --no-lane-splitting to bcl2fastq.")
    parser.add_argument("--samples-per-pass", type=int, default=0,
                        help="Max samples per bcl2fastq pass (0 = auto from the "
                             "open-file limit, `ulimit -n`).")
    parser.add_argument("--force-fastq", action="store_true",
                        help="Always (re)run bcl2fastq, even if FASTQs were already "
                             "staged/rebuilt or on-instrument FASTQs exist.")
    parser.add_argument("--rebuild-fastq-dir", metavar="DIRNAME", default=None,
                        help="Demultiplex with bcl2fastq into <run_dir>/DIRNAME instead "
                             "of --reads-dir (ignores on-instrument FASTQs). Re-uses a "
                             "previous rebuild there unless --force-fastq is given.")
    parser.add_argument("--on-duplicate", choices=["warn", "error", "suffix"], default="warn",
                        help="Same sample_id in more than one run: warn (default), "
                             "error, or suffix (append _<run_name>).")
    parser.add_argument("--no-keep-extra", dest="keep_extra", action="store_false",
                        help="Only write the model columns + fq1/fq2/run_id "
                             "(by default extra sheet columns such as project, "
                             "description, sample_name_original, index, index2, "
                             "plate_name, well are kept).")
    parser.add_argument("-n", "--dry-run", action="store_true",
                        help="Print what would be done without running bcl2fastq, "
                             "staging files, or writing output.")
    args = parser.parse_args()

    for run_dir in args.runs:
        if not run_dir.is_dir():
            sys.exit(f"Run directory not found: {run_dir}")
    if not args.model.exists():
        sys.exit(f"Model file not found: {args.model}")

    ref_cols = read_model_columns(args.model)
    if "sample_id" not in ref_cols:
        sys.exit(f"Model file {args.model} has no sample_id column.")

    kept_all, removed_all = [], []
    for run_dir in args.runs:
        try:
            kept, removed = process_run(
                run_dir, args.reads_dir, ref_cols, args.threads,
                args.no_lane_splitting, args.force_fastq, args.dry_run,
                args.rebuild_fastq_dir, args.keep_extra, args.samples_per_pass)
        except (FileNotFoundError, ValueError, KeyError, RuntimeError) as e:
            sys.exit(f"\nFailed on run {run_dir.name}: {e}")
        kept_all.append(kept)
        removed_all.append(removed)

    if args.dry_run:
        print("\n[dry-run] Stopping before writing combined metadata.")
        return

    combined = pd.concat(kept_all, ignore_index=True).fillna("")
    if combined.empty:
        sys.exit("No samples with FASTQs in any run — nothing written.")

    dup_mask = combined["sample_id"].duplicated(keep=False)
    if dup_mask.any():
        report = combined.loc[dup_mask].groupby("sample_id")["run_id"].apply(list)
        msg = f"{len(report)} sample_id(s) appear in more than one run:\n{report}"
        if args.on_duplicate == "error":
            sys.exit(f"Failed: {msg}")
        elif args.on_duplicate == "suffix":
            combined.loc[dup_mask, "sample_id"] = (
                combined.loc[dup_mask, "sample_id"] + "_" + combined.loc[dup_mask, "run_id"])
            print(f"Note: disambiguated duplicate sample_id(s) with run suffix:\n{report}")
        else:
            print(f"Warning: {msg}\nAmpSeeker needs unique sample_id values — "
                  "re-run with --on-duplicate suffix or error.")

    # Put model columns first, then fq1/fq2/run_id, then extras in first-seen order.
    lead = list(ref_cols) + ["fq1", "fq2", "run_id"]
    combined = combined[lead + [c for c in combined.columns if c not in lead]]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(args.output, sep="\t", index=False, lineterminator="\n")

    removed = pd.concat(removed_all, ignore_index=True)
    removed_log = args.output.with_name(args.output.stem + ".removed_samples.tsv")
    removed.to_csv(removed_log, sep="\t", index=False, lineterminator="\n")

    print(f"\nWrote combined metadata for {len(combined)} sample(s) across "
          f"{len(args.runs)} run(s) to: {args.output}")
    print(f"Removed {len(removed)} sample(s) without FASTQs; log: {removed_log}")
    print(f"Columns: {combined.columns.tolist()}")


if __name__ == "__main__":
    main()
