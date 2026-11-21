#!/usr/bin/env python3
# Script: convert_samplesheet.py
# Purpose: Convert an Illumina SampleSheet.csv into the flat column format
#          AmpSeeker expects for its sample metadata input.
# Input:   Illumina SampleSheet.csv (with [Header]/[Data] sections) +
#          a model/example samplesheet defining the target column set
# Output:  Reformatted CSV in AmpSeeker's expected format; original is backed up
#          to backup_<name> (never overwritten once it exists).
# Dependencies: pandas
#
# Example:
#   python convert_samplesheet.py -i ../MACARD-1/<run>/SampleSheet.csv \
#       -m resources/exampleSampleSheet.csv
#
# The normalisation logic (normalise_columns / to_ampseeker_frame) is also
# imported by build_multirun_metadata.py so both scripts read sheets the
# same way.

import argparse
import re
import shutil
import sys
from pathlib import Path

import pandas as pd

# Illumina column name (any case) -> AmpSeeker column name.
COLUMN_RENAME = {
    "sample_id": "sample_id",
    "sample_name": "sample_name",
    "sample_plate": "plate_name",
    "sample_well": "well",
    # bcl2fastq matches column names case-insensitively, so a lowercase
    # 'sample_project' would still make it write FASTQs into a
    # per-project subfolder. 'project' is ignored by bcl2fastq.
    "sample_project": "project",
}

# Columns AmpSeeker expects but that Illumina sheets may not provide;
# added blank so the user can fill them in by hand.
METADATA_COLUMNS = ["taxon", "location", "country", "latitude", "longitude"]

# Extra (non-model) columns that are dropped even with --keep-extra because
# they only duplicate another column.
DROP_EXTRA = {"i7_index_id", "i5_index_id"}


def read_samplesheet(path: Path):
    """Read an Illumina SampleSheet and return (header_lines, DataFrame).

    Uses utf-8-sig so a leading UTF-8 BOM (common in Illumina-exported
    sheets) is stripped. Every value is read as a string so IDs such as
    plate '03' or well 'A01' are not silently turned into numbers.
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
    df = pd.read_csv(path, skiprows=data_start, encoding="utf-8-sig",
                     dtype=str, keep_default_na=False)
    # Drop fully-empty trailing rows/columns (Excel padding).
    df = df.loc[:, [c for c in df.columns if not str(c).startswith("Unnamed:")
                    or df[c].str.strip().ne("").any()]]
    df = df[df.apply(lambda r: r.str.strip().ne("").any(), axis=1)]
    df = df.apply(lambda s: s.str.strip())
    return header, df.reset_index(drop=True)


def snake(name: str) -> str:
    """'Sample_ID' -> 'sample_id', 'Collection Date' -> 'collection_date'."""
    return re.sub(r"[^0-9a-z]+", "_", str(name).strip().lower()).strip("_")


def pad_line_to_width(line: str, width: int) -> str:
    """Pad a raw header line with trailing empty fields up to `width`
    fields (never truncates). Mirrors what Excel does when it re-saves a
    ragged CSV."""
    stripped = line.rstrip("\n").rstrip("\r")
    fields = stripped.split(",")
    if len(fields) < width:
        fields = fields + [""] * (width - len(fields))
    return ",".join(fields) + "\n"


def write_samplesheet(path: Path, header, df: pd.DataFrame, pad_width: int = None):
    """Write the SampleSheet back out with original header + new data rows,
    LF line endings, no BOM. With pad_width, header lines are padded to a
    rectangular grid (what AmpSeeker's [Data] regex expects)."""
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
    """Collapse columns that share a name after renaming.

    Values are merged cell by cell: a blank cell takes the other column's
    value. Raises only if two non-blank values actually disagree, since
    silently picking one would risk writing the wrong sample metadata.
    """
    dupe_names = df.columns[df.columns.duplicated()].unique().tolist()
    if not dupe_names:
        return df

    merged = {}
    for name in dupe_names:
        cols = df.loc[:, df.columns == name].astype(str)
        result = cols.iloc[:, 0].copy()
        for j in range(1, cols.shape[1]):
            other = cols.iloc[:, j]
            conflict = result.ne("") & other.ne("") & result.ne(other)
            if conflict.any():
                raise ValueError(
                    f"Column '{name}' appears {cols.shape[1]} times after renaming "
                    f"with conflicting values in {conflict.sum()} row(s). "
                    "Fix the source sheet or adjust --rename so these map to distinct names."
                )
            result = result.where(result.ne(""), other)
        merged[name] = result

    df = df.loc[:, ~df.columns.duplicated(keep="first")].copy()
    for name, values in merged.items():
        df[name] = values
    print(f"Note: merged duplicate column(s): {dupe_names}")
    return df


def find_col(df_columns, target_lower: str):
    """Return the actual column name matching target_lower case-insensitively, or None."""
    for c in df_columns:
        if c.lower() == target_lower:
            return c
    return None


def normalise_columns(df: pd.DataFrame, extra_rename: dict = None) -> pd.DataFrame:
    """Rename every column to AmpSeeker's lowercase snake_case naming.

    - Illumina names (Sample_ID, Sample_Plate, ...) are mapped through
      COLUMN_RENAME, case-insensitively.
    - Every other column is lowercased (Country -> country,
      Latitude -> latitude, ...). Previously these kept their capital
      letter, a blank lowercase twin was added, and the blank one won —
      silently wiping country/location/taxon/latitude/longitude.
    - --rename entries (matched case-insensitively) win over both.
    - Already-converted sheets (sample_ID / plate_name / well) pass
      through unchanged, so the function is idempotent.
    """
    extra = {snake(k): v for k, v in (extra_rename or {}).items()}
    new_cols = []
    for c in df.columns:
        key = snake(c)
        if key in extra:
            new_cols.append(extra[key])
        else:
            new_cols.append(COLUMN_RENAME.get(key, key))
    df = df.copy()
    df.columns = new_cols
    return dedupe_columns(df)


def to_ampseeker_frame(
    df: pd.DataFrame,
    ref_cols: list,
    extra_rename: dict = None,
    sample_name_from_id: bool = True,
    keep_extra: bool = True,
) -> pd.DataFrame:
    """Turn a raw (or already converted) [Data] frame into AmpSeeker's
    column layout: the model's columns first, in the model's order, then
    (with keep_extra) any additional metadata columns from the sheet.
    """
    df = normalise_columns(df, extra_rename)
    ref_cols = ["sample_id" if snake(c) == "sample_id" else snake(c) for c in ref_cols]

    for col in METADATA_COLUMNS:
        if col not in df.columns:
            df[col] = ""

    missing = [c for c in ref_cols if c not in df.columns]
    if missing:
        raise KeyError(
            f"Input sheet is missing columns required by the model sheet: {missing}. "
            f"Input columns after renaming: {df.columns.tolist()}"
        )

    if sample_name_from_id:
        # Keep the original short name; bcl2fastq names FASTQs after
        # sample_name, so sample_name is set to sample_id to make FASTQ
        # file names match the sample IDs.
        if "sample_name" in df.columns and "sample_name_original" not in df.columns:
            if df["sample_name"].ne(df["sample_id"]).any():
                df["sample_name_original"] = df["sample_name"]
        df["sample_name"] = df["sample_id"]

    bad = df["sample_id"].isna() | df["sample_id"].astype(str).str.strip().eq("")
    if bad.any():
        raise ValueError(
            f"{bad.sum()} row(s) have a missing 'sample_id' after conversion:\n{df[bad]}"
        )
    dup = df["sample_id"].duplicated(keep=False)
    if dup.any():
        raise ValueError(
            f"Duplicate sample_id(s) within one sheet: {sorted(df.loc[dup, 'sample_id'].unique())}"
        )

    for col in ("latitude", "longitude"):
        vals = df[col].astype(str).str.strip()
        bad_num = vals.ne("") & pd.to_numeric(vals, errors="coerce").isna()
        if bad_num.any():
            print(f"Warning: {bad_num.sum()} non-numeric value(s) in '{col}': "
                  f"{sorted(vals[bad_num].unique())[:10]}")

    extras = []
    if keep_extra:
        extras = [c for c in df.columns
                  if c not in ref_cols and c not in DROP_EXTRA and df[c].astype(str).str.strip().ne("").any()]
    out = df[ref_cols + extras].copy()

    # AmpSeeker's own SampleSheet parser splits [Data] lines on bare commas
    # and doesn't understand CSV quoting, so commas inside values would
    # shift columns. Replace them.
    for c in out.columns:
        s = out[c].astype(str)
        if s.str.contains(",").any():
            print(f"Note: replaced ',' with ';' inside values of column '{c}'.")
            out[c] = s.str.replace(",", ";", regex=False)
    return out


def read_model_columns(model_path: Path) -> list:
    """Model may be an Illumina-style sheet (.csv with [Data]) or a flat TSV."""
    if model_path.suffix.lower() == ".tsv":
        return pd.read_csv(model_path, sep="\t", nrows=0).columns.tolist()
    _, model_df = read_samplesheet(model_path)
    return model_df.columns.tolist()


def convert(
    input_path: Path,
    model_path: Path,
    output_path: Path,
    backup: bool,
    extra_rename: dict,
    sample_name_from_id: bool,
    pad_header: bool,
    keep_extra: bool = True,
):
    backup_path = input_path.with_name(f"backup_{input_path.name}")
    in_place = output_path.resolve() == input_path.resolve()

    source = input_path
    if in_place and backup_path.exists():
        # Already converted once: re-convert from the untouched original
        # instead of the already-converted copy, and never overwrite the
        # backup (doing so would destroy the only original).
        source = backup_path
        print(f"Found existing {backup_path.name}; converting from it (backup left untouched).")

    header, df = read_samplesheet(source)
    ref_cols = read_model_columns(model_path)
    df = to_ampseeker_frame(df, ref_cols, extra_rename, sample_name_from_id, keep_extra)

    if backup and in_place and not backup_path.exists():
        shutil.copy2(input_path, backup_path)
        print(f"Backed up original to: {backup_path}")

    write_samplesheet(
        output_path, header, df,
        pad_width=len(df.columns) if pad_header else None,
    )
    print(f"Wrote converted samplesheet to: {output_path}")
    print(f"Rows: {len(df)}  Columns: {df.columns.tolist()}")
    for col in METADATA_COLUMNS:
        n_blank = df[col].astype(str).str.strip().eq("").sum()
        if n_blank:
            print(f"  {col}: {n_blank}/{len(df)} blank")


def main():
    parser = argparse.ArgumentParser(
        description="Convert an Illumina SampleSheet.csv to AmpSeeker's expected format."
    )
    parser.add_argument("-i", "--input", required=True, type=Path,
                        help="Path to the raw Illumina SampleSheet.csv")
    parser.add_argument("-m", "--model", required=True, type=Path,
                        help="Model/example samplesheet defining the target columns "
                             "(e.g. resources/exampleSampleSheet.csv)")
    parser.add_argument("-o", "--output", type=Path, default=None,
                        help="Where to write the converted samplesheet. "
                             "Defaults to overwriting --input (with a backup made first).")
    parser.add_argument("--no-backup", action="store_true",
                        help="Skip creating a backup_<name> copy when writing in place.")
    parser.add_argument("--rename", type=str, default="",
                        help="Extra/override column renames as 'Old=New,Old2=New2' "
                             "(case-insensitive; takes priority over the defaults).")
    parser.add_argument("--sample-name-from-id", dest="sample_name_from_id",
                        action="store_true", default=True,
                        help="Set sample_name = sample_id (default: on), so bcl2fastq "
                             "names FASTQs after sample_id. The original value is kept "
                             "in 'sample_name_original'.")
    parser.add_argument("--no-sample-name-from-id", dest="sample_name_from_id",
                        action="store_false",
                        help="Keep the original Sample_Name values.")
    parser.add_argument("--keep-extra", dest="keep_extra", action="store_true", default=True,
                        help="Append sheet columns not in the model (e.g. description, "
                             "project, sample_name_original) after the model columns (default: on).")
    parser.add_argument("--no-keep-extra", dest="keep_extra", action="store_false",
                        help="Output exactly the model's columns.")
    parser.add_argument("--pad-header", dest="pad_header", action="store_true", default=True,
                        help="Pad header lines with trailing commas to a rectangular grid (default: on).")
    parser.add_argument("--no-pad-header", dest="pad_header", action="store_false",
                        help="Leave header lines at their original (ragged) width.")
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
            keep_extra=args.keep_extra,
        )
    except (ValueError, KeyError) as e:
        sys.exit(f"Conversion failed: {e}")


if __name__ == "__main__":
    main()
