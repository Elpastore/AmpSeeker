#!/usr/bin/env python3
# Script: convert_samplesheet.py
# Purpose: Convert an Illumina SampleSheet.csv into the flat column format
#          AmpSeeker expects for its sample metadata input.
# Input:   Illumina SampleSheet.csv (with [Header]/[Data] sections) +
#          a model/example samplesheet defining the target column set
# Output:  Reformatted CSV in AmpSeeker's expected format; original is backed up
# Dependencies: pandas

import argparse
import shutil
import sys
from pathlib import Path

import pandas as pd

# Illumina column name -> AmpSeeker column name
COLUMN_RENAME = {
    "Sample_ID": "sample_id",
    "Sample_Name": "sample_name",
    "Sample_Plate": "plate_name",
    "Sample_Well": "well",
}

# Columns AmpSeeker expects but that Illumina sheets don't provide;
# added blank so the user can fill them in by hand.
METADATA_COLUMNS = ["taxon", "location", "country", "latitude", "longitude"]


def read_samplesheet(path: Path):
    """Read an Illumina SampleSheet and return (header_lines, DataFrame).

    Uses utf-8-sig so a leading UTF-8 BOM (common in Illumina-exported
    sheets) is stripped rather than silently embedded in the first header
    line, where it can break exact-match parsing downstream.
    """
    with open(path, encoding="utf-8-sig") as f:
        lines = f.readlines()

    data_start = None
    for i, line in enumerate(lines):
        if line.strip().lower().startswith("[data"):
            data_start = i + 1
            break

    if data_start is None:
        raise ValueError(
            f"Could not find a '[Data]' section in {path}. "
            "Is this a standard Illumina SampleSheet?"
        )

    header = lines[:data_start]
    df = pd.read_csv(path, skiprows=data_start, encoding="utf-8-sig")
    return header, df


def pad_line_to_width(line: str, width: int) -> str:
    """Pad a raw header line's comma-separated fields with trailing empty
    fields until it has `width` fields total.

    Mirrors what Excel/Numbers does automatically when it opens a ragged
    CSV and saves it back out: it treats the sheet as a rectangular grid
    sized to its widest row and pads every shorter row with empty cells.
    Only pads (never truncates), so a header line that's already wider
    than `width` is left untouched.
    """
    stripped = line.rstrip("\n").rstrip("\r")
    fields = stripped.split(",")
    if len(fields) < width:
        fields = fields + [""] * (width - len(fields))
    return ",".join(fields) + "\n"


def write_samplesheet(path: Path, header, df: pd.DataFrame, pad_width: int = None):
    """Write the SampleSheet back out with original header + new data rows.

    encoding='utf-8' (no BOM) keeps the output byte-for-byte plain ASCII/UTF-8.
    newline='' hands full control of line endings to the csv writer inside
    df.to_csv(); without it, Python's own text-mode newline translation and
    the csv writer's newline insertion can disagree, producing a file with
    mixed \\r\\n / \\n line endings that some CSV parsers reject outright
    (this is a well-known pandas gotcha, not specific to this script).
    lineterminator='\\n' then makes every line in the file consistently
    LF-only, matching a standard Linux/HPC text file.

    If pad_width is given, every preserved header line (including blank
    lines and the '[Data]' marker itself) is padded with trailing commas
    out to that many fields, so the whole file is a rectangular grid —
    matching what AmpSeeker's parser expects and what Excel produces when
    it re-saves a ragged CSV.
    """
    if pad_width:
        header = [pad_line_to_width(line, pad_width) for line in header]
    with open(path, "w", encoding="utf-8", newline="") as f:
        f.writelines(header)
        df.to_csv(f, index=False, lineterminator="\n")


def parse_rename_arg(rename_str: str) -> dict:
    """Parse '--rename Old=New,Old2=New2' into a dict, raising on bad format."""
    mapping = {}
    for pair in rename_str.split(","):
        pair = pair.strip()
        if not pair:
            continue
        if "=" not in pair:
            raise ValueError(
                f"Malformed --rename entry '{pair}' — expected 'OldName=NewName'."
            )
        old, new = pair.split("=", 1)
        mapping[old.strip()] = new.strip()
    return mapping


def dedupe_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Collapse columns that share a name after renaming, if their values agree.

    Raises if two same-named columns disagree, since silently picking one
    would risk writing the wrong sample metadata.
    """
    dupe_names = df.columns[df.columns.duplicated()].unique().tolist()
    if not dupe_names:
        return df

    for name in dupe_names:
        cols = df.loc[:, df.columns == name]
        first = cols.iloc[:, 0]
        for j in range(1, cols.shape[1]):
            other = cols.iloc[:, j]
            if not first.astype(str).fillna("").equals(other.astype(str).fillna("")):
                mismatches = first.astype(str) != other.astype(str)
                raise ValueError(
                    f"Column '{name}' appears {cols.shape[1]} times after renaming "
                    f"with conflicting values in {mismatches.sum()} row(s). "
                    "Fix the source sheet or adjust --rename so these map to distinct names."
                )
    # All duplicate groups agreed — keep just the first occurrence of each.
    df = df.loc[:, ~df.columns.duplicated(keep="first")]
    print(f"Note: collapsed duplicate column(s) with matching values: {dupe_names}")
    return df


def match_case_insensitive(df_columns, ref_cols):
    """Map each ref_col to the actual df column name, ignoring case.

    Returns (rename_map, missing) where rename_map maps the df's current
    column name -> the model sheet's exact spelling, and missing lists any
    ref_cols with no match at all (case-insensitive) in df_columns.
    """
    lower_lookup = {c.lower(): c for c in df_columns}
    rename_map = {}
    missing = []
    for ref_col in ref_cols:
        match = lower_lookup.get(ref_col.lower())
        if match is None:
            missing.append(ref_col)
        else:
            rename_map[match] = ref_col
    return rename_map, missing


def find_col(df_columns, target_lower: str):
    """Return the actual column name matching target_lower case-insensitively, or None."""
    for c in df_columns:
        if c.lower() == target_lower:
            return c
    return None


def apply_sample_name_from_id(df: pd.DataFrame) -> pd.DataFrame:
    """Overwrite the sample_name column's values with the sample_id column's values.

    Used when Sample_Name in the raw sheet is blank/unreliable and the
    meaningful identifier actually lives in Sample_ID.
    """
    id_col = find_col(df.columns, "sample_id")
    if id_col is None:
        raise ValueError(
            "--sample-name-from-id was requested but no 'sample_id' column "
            f"was found. Columns present: {df.columns.tolist()}"
        )
    name_col = find_col(df.columns, "sample_name")
    if name_col is None:
        # No sample_name column yet — create one from sample_id.
        df["sample_name"] = df[id_col]
    else:
        df[name_col] = df[id_col]
    return df


def force_lowercase_sample_id(df: pd.DataFrame, ref_cols: list) -> tuple:
    """Force whatever the model calls the ID column (e.g. 'sample_ID') to be
    output as literal lowercase 'sample_id', regardless of the model sheet's
    own casing. Returns (df, updated_ref_cols).
    """
    id_col = find_col(df.columns, "sample_id")
    if id_col and id_col != "sample_id":
        df = df.rename(columns={id_col: "sample_id"})
    ref_cols = ["sample_id" if c.lower() == "sample_id" else c for c in ref_cols]
    return df, ref_cols


def convert(
    input_path: Path,
    model_path: Path,
    output_path: Path,
    backup: bool,
    extra_rename: dict,
    sample_name_from_id: bool,
    pad_header: bool,
):
    header, df = read_samplesheet(input_path)

    # Rename Illumina columns to AmpSeeker's expected names.
    # CLI-supplied --rename entries take priority over the defaults.
    rename_map = {**COLUMN_RENAME, **extra_rename}
    df = df.rename(columns=rename_map)
    df = dedupe_columns(df)

    # Add any missing metadata columns, blank.
    for col in METADATA_COLUMNS:
        if col not in df.columns:
            df[col] = ""

    # Load the model/example sheet to get the exact target column set.
    _, sheet_model = read_samplesheet(model_path)
    ref_cols = sheet_model.columns.tolist()

    # Match case-insensitively (e.g. model's 'sample_ID' vs our 'sample_id'),
    # then rename df's columns to the model's exact spelling.
    col_rename, missing = match_case_insensitive(df.columns, ref_cols)
    if missing:
        raise KeyError(
            f"Input sheet is missing columns required by the model sheet: {missing}. "
            f"Input columns after renaming: {df.columns.tolist()}"
        )
    df = df.rename(columns=col_rename)

    if sample_name_from_id:
        df = apply_sample_name_from_id(df)

    # Always output the ID column as lowercase 'sample_id', regardless of
    # how the model sheet happens to spell it (e.g. 'sample_ID').
    df, ref_cols = force_lowercase_sample_id(df, ref_cols)

    df = df[ref_cols]  # keep only + reorder to match the model exactly

    # Defensive check: don't silently write a sheet with missing sample IDs.
    if df["sample_id"].isna().any() or (df["sample_id"].astype(str).str.strip() == "").any():
        bad_rows = df[df["sample_id"].isna() | (df["sample_id"].astype(str).str.strip() == "")]
        raise ValueError(
            f"{len(bad_rows)} row(s) have a missing 'sample_id' after conversion:\n{bad_rows}"
        )

    if backup and output_path == input_path:
        backup_path = input_path.with_name(f"backup_{input_path.name}")
        shutil.copy2(input_path, backup_path)
        print(f"Backed up original to: {backup_path}")

    write_samplesheet(
        output_path, header, df,
        pad_width=len(df.columns) if pad_header else None,
    )
    print(f"Wrote converted samplesheet to: {output_path}")
    print(f"Rows: {len(df)}  Columns: {df.columns.tolist()}")


def main():
    parser = argparse.ArgumentParser(
        description="Convert an Illumina SampleSheet.csv to AmpSeeker's expected format."
    )
    parser.add_argument(
        "-i", "--input", required=True, type=Path,
        help="Path to the raw Illumina SampleSheet.csv",
    )
    parser.add_argument(
        "-m", "--model", required=True, type=Path,
        help="Path to the model/example samplesheet defining the target columns "
             "(e.g. resources/exampleSampleSheet.csv)",
    )
    parser.add_argument(
        "-o", "--output", type=Path, default=None,
        help="Path to write the converted samplesheet. "
             "Defaults to overwriting --input (with a backup made first).",
    )
    parser.add_argument(
        "--no-backup", action="store_true",
        help="Skip creating a backup_<name> copy when writing in place.",
    )
    parser.add_argument(
        "--rename", type=str, default="",
        help="Extra/override column renames as 'Old=New,Old2=New2' "
             "(applied on top of, and taking priority over, the built-in "
             "Illumina -> AmpSeeker defaults). Useful for sheets with "
             "nonstandard column names.",
    )
    parser.add_argument(
        "--sample-name-from-id", dest="sample_name_from_id",
        action="store_true", default=True,
        help="Set sample_name values equal to sample_id values (default: on). "
             "Use when Sample_Name in the raw sheet is blank/unreliable and "
             "the meaningful identifier lives in Sample_ID.",
    )
    parser.add_argument(
        "--no-sample-name-from-id", dest="sample_name_from_id",
        action="store_false",
        help="Keep the original Sample_Name values instead of overwriting "
             "them with sample_id.",
    )
    parser.add_argument(
        "--pad-header", dest="pad_header",
        action="store_true", default=True,
        help="Pad every preserved header line with trailing commas so the "
             "whole file is a rectangular grid matching the data section's "
             "column count (default: on) — this is what AmpSeeker expects "
             "and what Excel produces when it re-saves a ragged CSV.",
    )
    parser.add_argument(
        "--no-pad-header", dest="pad_header",
        action="store_false",
        help="Leave header lines at their original (ragged) width.",
    )
    args = parser.parse_args()

    if not args.input.exists():
        sys.exit(f"Input file not found: {args.input}")
    if not args.model.exists():
        sys.exit(f"Model file not found: {args.model}")

    try:
        extra_rename = parse_rename_arg(args.rename)
    except ValueError as e:
        sys.exit(str(e))

    output_path = args.output or args.input

    try:
        convert(
            args.input, args.model, output_path,
            backup=not args.no_backup,
            extra_rename=extra_rename,
            sample_name_from_id=args.sample_name_from_id,
            pad_header=args.pad_header,
        )
    except (ValueError, KeyError) as e:
        sys.exit(f"Conversion failed: {e}")


if __name__ == "__main__":
    main()
