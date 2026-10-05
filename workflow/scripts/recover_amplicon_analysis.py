#!/usr/bin/env python
"""Standalone recovery path for the AmpSeeker amplicon-scale analysis stage.

Why this exists
----------------
`mpileup_call_amplicons` (workflow/rules/map-call-illumina.smk) used to run
bcftools mpileup with no region restriction, so it called genotypes at every
covered base genome-wide instead of just the panel. For this run's 1631
samples that produced a 25GB, ~2.5M-record VCF
(results/vcfs/amplicons/<dataset>.annot.vcf). Both the sample_quality_control
and snp_dataframe notebooks load that file in full via
`allel.read_vcf(..., fields="*")` and get OOM-killed by the SLURM cgroup
(96GB allocation, ~123GB observed peak) -- see logs/AmpSeeker.118240.error.

The rule itself is now fixed to pass -R going forward, but that requires
re-running bcftools mpileup plus a two-tier merge across 1631 samples, which
is expensive. This script instead recovers usable results from the VCF that
has *already* been generated, without needing that re-run:

  1. Stream the existing oversized VCF once (never loading it fully into
     memory) and keep only records at the panel's target positions
     (+/- --pad bases), writing a small filtered VCF.
  2. Run the same PCA-outlier-detection logic as
     workflow/notebooks/sample-quality-control.ipynb, and the same
     vcf_to_excel export as workflow/notebooks/snp-dataframe.ipynb, against
     that filtered VCF instead of the original 25GB one.

Each stage has its own try/except and its own logging, so a failure in one
stage doesn't take down the others, and you get a real Python traceback
in results/recovery/recovery.log instead of papermill's silent
"Kernel died".

Usage (run inside the AmpSeeker-python conda env, which has
pandas/numpy/scikit-allel):

    conda activate AmpSeeker-python
    python workflow/scripts/recover_amplicon_analysis.py

All arguments default to this run's config/eniyou_agam.yaml values; override
as needed. --pad defaults to 0 (i.e. it matches the panel's exact SNP
positions in config/ag-vampir.bed -- the same scope as the "targets" VCF)
because no amplicon-insert BED exists in this repo. Pass --pad <bases> if
you know the real amplicon insert length and want the wider "whole-amplicon"
scope this stage is meant to have.
"""
import argparse
import logging
import os
import sys
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--wkdir", default=".", help="AmpSeeker working directory")
    p.add_argument("--dataset", default="agvampir-eniyou-agam")
    p.add_argument("--platform", default="illumina")
    p.add_argument("--cohort-col", default="location")
    p.add_argument("--bed", default="config/ag-vampir.bed", help="Panel target BED")
    p.add_argument("--metadata", default="results/config/metadata.tsv")
    p.add_argument(
        "--amplicons-vcf",
        default=None,
        help="Oversized amplicons VCF (default: results/vcfs/amplicons/<dataset>.annot.vcf)",
    )
    p.add_argument(
        "--targets-vcf",
        default=None,
        help="Targets VCF (default: results/vcfs/targets/<dataset>.annot.vcf)",
    )
    p.add_argument(
        "--pad",
        type=int,
        default=0,
        help="Bases to pad each target position by before filtering the amplicons VCF",
    )
    p.add_argument(
        "--force-filter",
        action="store_true",
        help="Redo the streaming filter pass even if a filtered VCF already exists in --outdir",
    )
    p.add_argument(
        "--outdir",
        default="results/recovery",
        help="Output directory, kept separate from the pipeline's own results/",
    )
    p.add_argument("--missing-threshold", type=float, default=0.2)
    p.add_argument("--zscore-threshold", type=float, default=4)
    return p.parse_args()


def resolve(wkdir, path):
    return path if os.path.isabs(path) else os.path.join(wkdir, path)


def setup_logging(outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    log_path = outdir / "recovery.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler(sys.stdout)],
    )
    return log_path


def load_target_positions(bed_path, pad):
    """Return {contig: sorted [(start, end), ...]} 1-based inclusive intervals."""
    import pandas as pd

    bed = pd.read_csv(bed_path, sep="\t", header=None)
    by_contig = {}
    for _, row in bed.iterrows():
        contig = str(row[0])
        pos = int(row[2])  # bed 'end' == 1-based VCF POS for these single-base rows
        by_contig.setdefault(contig, []).append((pos - pad, pos + pad))
    for contig in by_contig:
        by_contig[contig].sort()
    return by_contig


def position_in_intervals(pos, intervals_for_contig):
    for start, end in intervals_for_contig:
        if pos < start:
            break
        if start <= pos <= end:
            return True
    return False


def filter_vcf_by_position(src_path, dst_path, intervals):
    """Stream src_path -> dst_path, keeping header lines and only data lines
    whose CHROM/POS fall inside `intervals`. Never loads the file into memory."""
    kept = 0
    seen = 0
    with open(src_path) as src, open(dst_path, "w") as dst:
        for line in src:
            if line.startswith("#"):
                dst.write(line)
                continue
            seen += 1
            chrom, pos = line.split("\t", 2)[:2]
            contig_intervals = intervals.get(chrom)
            if contig_intervals and position_in_intervals(int(pos), contig_intervals):
                dst.write(line)
                kept += 1
            if seen % 500_000 == 0:
                logging.info("  ...scanned %d records, kept %d so far", seen, kept)
    logging.info("Filtered %s -> %s: kept %d/%d records", src_path, dst_path, kept, seen)
    return kept


def stage_filter_amplicons_vcf(args, outdir):
    logging.info("STAGE 1: filtering oversized amplicons VCF down to panel positions")
    amplicons_vcf = args.amplicons_vcf or resolve(
        args.wkdir, f"results/vcfs/amplicons/{args.dataset}.annot.vcf"
    )
    if not os.path.exists(amplicons_vcf):
        raise FileNotFoundError(f"amplicons VCF not found: {amplicons_vcf}")

    filtered_path = outdir / f"{args.dataset}.amplicons.filtered.vcf"
    if filtered_path.exists() and not args.force_filter:
        kept = sum(1 for line in open(filtered_path) if not line.startswith("#"))
        logging.info(
            "  reusing existing %s (%d records) -- pass --force-filter to redo the "
            "streaming pass over the 25GB source VCF",
            filtered_path, kept,
        )
        return str(filtered_path)

    bed_path = resolve(args.wkdir, args.bed)
    intervals = load_target_positions(bed_path, args.pad)
    kept = filter_vcf_by_position(amplicons_vcf, filtered_path, intervals)
    if kept == 0:
        raise RuntimeError(
            "No records survived filtering -- check that --bed contig names "
            "match the VCF (e.g. '2L' vs 'chr2L') and that --pad is not negative."
        )
    return str(filtered_path)


def stage_snp_dataframes(args, outdir, filtered_amplicons_vcf):
    logging.info("STAGE 2: SNP dataframe / Excel export (mirrors snp-dataframe.ipynb)")
    import ampseeker as amp

    targets_vcf = args.targets_vcf or resolve(
        args.wkdir, f"results/vcfs/targets/{args.dataset}.annot.vcf"
    )

    try:
        df = amp.vcf_to_excel(
            vcf_path=targets_vcf,
            excel_path=str(outdir / f"{args.dataset}-targets-snps.xlsx"),
            convert_genotypes=True,
            split_multiallelic=True,
        )
        logging.info("  targets SNP dataframe: %d rows", len(df))
    except Exception:
        logging.exception("  targets SNP dataframe failed")

    if not filtered_amplicons_vcf:
        logging.warning("  skipping amplicons SNP dataframe: no filtered VCF available")
        return

    try:
        df = amp.vcf_to_excel(
            vcf_path=filtered_amplicons_vcf,
            excel_path=str(outdir / f"{args.dataset}-amplicons-snps.xlsx"),
            convert_genotypes=True,
            split_multiallelic=True,
        )
        logging.info("  amplicons SNP dataframe: %d rows", len(df))
    except Exception:
        logging.exception("  amplicons SNP dataframe failed")


def find_pca_outliers(pca_df, zscore_threshold=4):
    """Ported from workflow/notebooks/sample-quality-control.ipynb."""
    import numpy as np
    import pandas as pd
    from scipy import stats

    pca_df = pca_df.filter(like="PC")
    if pca_df.shape[1] == 0:
        # Too few informative sites for this cohort to produce any PCA
        # components at all -- nothing to score as an outlier.
        return pd.DataFrame({"max_zscore": [], "is_outlier": [], "outlier_components": []})
    zscores = pd.DataFrame(
        np.abs(stats.zscore(pca_df)), columns=pca_df.columns, index=pca_df.index
    )
    max_zscores = zscores.max(axis=1)
    outlier_components = zscores.apply(lambda x: x > zscore_threshold)
    outlier_component_lists = outlier_components.apply(lambda x: list(x.index[x]), axis=1)
    results = pd.DataFrame(
        {
            "max_zscore": max_zscores,
            "is_outlier": max_zscores > zscore_threshold,
            "outlier_components": outlier_component_lists,
        }
    )
    return results.sort_values("max_zscore", ascending=False)


def stage_pca_outliers(args, outdir, filtered_amplicons_vcf):
    logging.info("STAGE 3: PCA outlier detection (mirrors sample-quality-control.ipynb)")
    import pandas as pd
    import ampseeker as amp

    metadata_path = resolve(args.wkdir, args.metadata)
    metadata = pd.read_csv(metadata_path, sep="\t")
    geno, pos, contigs, metadata, ref, alt, ann = amp.load_variants(
        filtered_amplicons_vcf, metadata, platform=args.platform, filter_indel=True
    )

    pca_exclude_samples = []
    for coh in metadata[args.cohort_col].unique():
        # amp.pca() runs `metadata.eval(query)` in its own stack frame, so a
        # pandas "@coh" local-variable reference (resolved against the
        # caller's frame) doesn't reach it -- interpolate the value directly
        # into the query string instead, as the source notebook does.
        cohort_query = f"{args.cohort_col} == '{coh}'"
        cohort_meta = metadata.query(cohort_query)
        if cohort_meta.shape[0] < 5:
            logging.info(
                "  skipping cohort %r (only %d samples, need >=5 for PCA)",
                coh, cohort_meta.shape[0],
            )
            continue
        try:
            pca_df, model = amp.pca(
                geno, metadata, query=cohort_query, n_components=3,
                missing_threshold=args.missing_threshold,
            )
            df_outliers = find_pca_outliers(
                pca_df.set_index("sample_id"), zscore_threshold=args.zscore_threshold
            )
            n_outliers = int(df_outliers["is_outlier"].sum())
            logging.info("  %s: %d PCA outliers in %d samples", coh, n_outliers, df_outliers.shape[0])
            pca_exclude_samples.extend(df_outliers[df_outliers["is_outlier"]].index.tolist())
        except Exception:
            logging.exception("  PCA failed for cohort %r", coh)

    out_path = outdir / f"{args.dataset}.pca_exclude_samples.txt"
    out_path.write_text("\n".join(pca_exclude_samples) + ("\n" if pca_exclude_samples else ""))
    logging.info("  wrote %d PCA-outlier sample IDs to %s", len(pca_exclude_samples), out_path)
    return pca_exclude_samples


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    log_path = setup_logging(outdir)
    logging.info("Log file: %s", log_path)
    logging.info("Args: %s", vars(args))

    sys.path.append(resolve(args.wkdir, "workflow/lib"))

    stage_results = {}

    try:
        filtered_amplicons_vcf = stage_filter_amplicons_vcf(args, outdir)
        stage_results["filter"] = "ok"
    except Exception:
        logging.exception("STAGE 1 failed")
        filtered_amplicons_vcf = None
        stage_results["filter"] = "failed"

    try:
        stage_snp_dataframes(args, outdir, filtered_amplicons_vcf)
        stage_results["snp_dataframe"] = "ok"
    except Exception:
        logging.exception("STAGE 2 failed")
        stage_results["snp_dataframe"] = "failed"

    if filtered_amplicons_vcf:
        try:
            stage_pca_outliers(args, outdir, filtered_amplicons_vcf)
            stage_results["pca_outliers"] = "ok"
        except Exception:
            logging.exception("STAGE 3 failed")
            stage_results["pca_outliers"] = "failed"
    else:
        stage_results["pca_outliers"] = "skipped (no filtered VCF)"

    logging.info("Summary: %s", stage_results)
    if any(v == "failed" for v in stage_results.values()):
        sys.exit(1)


if __name__ == "__main__":
    main()
