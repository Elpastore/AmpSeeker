"""
Batched (N-way) VCF merge — an alternative to the {call_type}.merged.vcf rules
in utilities.smk, meant to handle sample counts the hardcoded 2-way split in
common.smk (large_sample_size / samples1 / samples2) cannot. bcftools merge
has an internal ceiling of roughly 1000 input files, and the OS's open-file
descriptor limit (ulimit -n, commonly 1024) imposes a similar cap
independently of bcftools itself. A fixed 50/50 split only postpones the
problem: at ~2000 samples each half is already at/over that ceiling, and at
~10,000 samples each half (5000 files) blows well past it.

NOT YET WIRED IN. This file's rules (bcftools_merge_batch, tabix_merge_batch,
bcftools_merge_final) produce results/vcfs/{call_type}/{dataset}.merged.vcf,
the same path as utilities.smk's bcftools_merge / bcftools_merge1-3, so do
not `include:` both at once — Snakemake will raise an ambiguous/duplicate
output error.

To test this instead of the existing merge logic:
    1. In the Snakefile, replace `include: "rules/utilities.smk"`'s merge
       rules — comment out (or temporarily remove) bcftools_merge,
       bcftools_merge1, bcftools_merge2, bgzip2, tabix2 and bcftools_merge3
       in utilities.smk (or the whole `if large_sample_size:` block).
    2. Add `include: "rules/scalable-merge.smk"` to the Snakefile, after
       `include: "rules/common.smk"` — this file relies on the `samples`
       and `dataset` globals that common.smk/Snakefile already define, and
       does not use large_sample_size/samples1/samples2 at all.
    3. common.smk's ampseeker_outputs() only requires the
       ".complete.{dataset}.merge_vcfs" touch file when large_sample_size is
       True (see common.smk around line 139). This file's
       bcftools_merge_final always produces that touch file, so once wired
       in for real, ampseeker_outputs() should require it unconditionally
       (or the touch() should be dropped) rather than gating on
       large_sample_size — left as-is here since that's a separate file.

How it works
------------
Samples are split into batches of at most MERGE_BATCH_SIZE (config key
"merge-batch-size", default 500 — comfortably under both bcftools's ~1000
file ceiling and a default ulimit -n of 1024, leaving headroom for logs,
tabix indices, etc. held open at the same time). Each batch is merged
independently (bcftools_merge_batch) into a compressed, tabixed intermediate
VCF under results/vcfs/{call_type}/.batches/, and then those (far fewer)
batch VCFs are merged into the final {dataset}.merged.vcf
(bcftools_merge_final).

This two-level fan-in comfortably covers datasets well beyond the original
1000-sample limit, e.g.:
  - 2,000 samples / batch size 500  ->   4 batches; final merge opens 4 files.
  - 10,000 samples / batch size 500 ->  20 batches; final merge opens 20 files.
  - 100,000 samples / batch size 500 -> 200 batches; final merge opens 200 files.
and scales up to roughly MERGE_BATCH_SIZE**2 samples (~250,000 at the
default) before the *final* merge step itself would need a third batching
level (not implemented here, since it's not needed at the sizes above —
flag it if you're testing something past ~250k samples).

A batch containing only a single sample has nothing to merge, so both merge
rules special-case n==1 (copy / reformat instead of calling bcftools merge).
"""

MERGE_BATCH_SIZE = int(config.get("merge-batch-size", 500))


wildcard_constraints:
    batch="[0-9]+",


def _make_batches(sample_list, batch_size):
    sample_list = list(sample_list)
    return [
        sample_list[i : i + batch_size] for i in range(0, len(sample_list), batch_size)
    ]


sample_batches = _make_batches(samples, MERGE_BATCH_SIZE)
n_batches = len(sample_batches)
batch_ids = [str(i) for i in range(n_batches)]
batch_samples = {str(i): batch for i, batch in enumerate(sample_batches)}


def get_batch_vcfs(wildcards):
    return expand(
        "results/vcfs/{call_type}/{sample}.calls.vcf.gz",
        call_type=wildcards.call_type,
        sample=batch_samples[wildcards.batch],
    )


def get_batch_tbis(wildcards):
    return expand(
        "results/vcfs/{call_type}/{sample}.calls.vcf.gz.tbi",
        call_type=wildcards.call_type,
        sample=batch_samples[wildcards.batch],
    )


rule bcftools_merge_batch:
    input:
        vcfs=get_batch_vcfs,
        idx=get_batch_tbis,
    output:
        vcf="results/vcfs/{call_type}/.batches/{dataset}.batch{batch}.vcf.gz",
    log:
        "logs/bcftools/{call_type}/merge_batch{batch}_{dataset}.log",
    conda:
        "../envs/AmpSeeker-cli.yaml"
    threads: 12
    shell:
        """
        n=$(echo {input.vcfs} | wc -w)
        if [ "$n" -eq 1 ]; then
            cp {input.vcfs} {output.vcf} 2> {log}
        else
            bcftools merge --threads {threads} -O z -o {output.vcf} {input.vcfs} --force-samples 2> {log}
        fi
        """


rule tabix_merge_batch:
    input:
        "results/vcfs/{call_type}/.batches/{dataset}.batch{batch}.vcf.gz",
    output:
        "results/vcfs/{call_type}/.batches/{dataset}.batch{batch}.vcf.gz.tbi",
    log:
        "logs/tabix/{call_type}/merge_batch{batch}_{dataset}.log",
    conda:
        "../envs/AmpSeeker-cli.yaml"
    shell:
        """
        tabix {input} 2> {log}
        """


rule bcftools_merge_final:
    input:
        vcfs=expand(
            "results/vcfs/{{call_type}}/.batches/{{dataset}}.batch{batch}.vcf.gz",
            batch=batch_ids,
        ),
        idx=expand(
            "results/vcfs/{{call_type}}/.batches/{{dataset}}.batch{batch}.vcf.gz.tbi",
            batch=batch_ids,
        ),
    output:
        vcf="results/vcfs/{call_type}/{dataset}.merged.vcf",
        holder=touch("results/vcfs/{call_type}/.complete.{dataset}.merge_vcfs"),
    log:
        "logs/bcftools/{call_type}/merge_final_{dataset}.log",
    conda:
        "../envs/AmpSeeker-cli.yaml"
    threads: 12
    shell:
        """
        n=$(echo {input.vcfs} | wc -w)
        if [ "$n" -eq 1 ]; then
            bcftools view -O v -o {output.vcf} {input.vcfs} 2> {log}
        else
            bcftools merge --threads {threads} -o {output.vcf} -O v {input.vcfs} --force-samples 2> {log}
        fi
        """
