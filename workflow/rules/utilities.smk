"""
A snakemake rule file which includes rules for bgzipping, tabix indexing, and merging bcf/vcfs.
Kept in a separate rule file due to maintain visibility of the main, analysis rules in analysis.smk
"""

rule reference_index:
    input:
        ref=config["reference-fasta"],
    output:
        idx=config["reference-fasta"] + ".fai",
    conda:
        "../envs/AmpSeeker-cli.yaml"
    log:
        "logs/reference_index.log",
    shell:
        """
        samtools faidx {input.ref} 2> {log}
        """

rule bam_index:
    input:
        "results/alignments/{sample}.bam",
    output:
        "results/alignments/{sample}.bam.bai",
    conda:
        "../envs/AmpSeeker-cli.yaml"
    log:
        "logs/index_bams/{sample}_index.log",
    shell:
        "samtools index {input} {output} 2> {log}"


rule bgzip:
    input:
        vcf="results/vcfs/{call_type}/{sample}.calls.vcf",
    output:
        vcfgz="results/vcfs/{call_type}/{sample}.calls.vcf.gz",
    log:
        log="logs/bgzip/{call_type}/{sample}.log",
    wrapper:
        "v1.17.4-17-g62b55d45/bio/bgzip"


rule tabix:
    input:
        calls="results/vcfs/{call_type}/{sample}.calls.vcf.gz",
    output:
        calls_tbi="results/vcfs/{call_type}/{sample}.calls.vcf.gz.tbi",
    conda:
        "../envs/AmpSeeker-cli.yaml"
    log:
        "logs/tabix/{call_type}/{sample}.log",
    shell:
        """
        tabix {input.calls} 2> {log}
        """


# bcftools merge cannot take >1000 input files at once (see common.smk, where
# `large_sample_size` and the samples1/samples2 split are computed), so datasets
# above that threshold are merged in two rounds (merge1 + merge2, each bgzipped/
# tabixed via bgzip2/tabix2, then combined in merge3) instead of the single-round
# bcftools_merge used for smaller datasets. Both branches produce the same output
# path (results/vcfs/{call_type}/{dataset}.merged.vcf), so they are defined inside
# this if/else rather than side by side, to avoid an AmbiguousRuleException.
#
# Fixed 2026-09-15 (previously broken whenever large_sample_size was True):
#   - bcftools_merge1/2 used a single-braced {call_type} inside expand(), which
#     expand() resolves via str.format() before Snakemake ever sees the rule. With
#     no call_type= argument supplied to expand(), this raised
#     `WildcardError: No values given for wildcard 'call_type'`. Fixed by doubling
#     the braces ({{call_type}}) so it survives expand() as a real wildcard, the
#     same way bcftools_merge already did it.
#   - tabix2's output used the wildcard constraint {n, [/d]}, a character class
#     matching a literal '/' or 'd' (probably meant to be the regex `\d`). It could
#     never match the actual values ("1", "2"), so bcftools_merge3's .tbi inputs
#     were unreachable. Fixed to {n,[0-9]}, matching bgzip2's constraint.
#   - --force-samples was added to merge1/2/3, to match bcftools_merge, so that
#     duplicate sample names fail (or are handled) the same way on both paths.
if not large_sample_size:

    rule bcftools_merge:
        input:
            vcfs=expand(
                "results/vcfs/{{call_type}}/{sample}.calls.vcf.gz", sample=samples
            ),
            idx=expand(
                "results/vcfs/{{call_type}}/{sample}.calls.vcf.gz.tbi", sample=samples
            ),
        output:
            vcf="results/vcfs/{call_type}/{dataset}.merged.vcf",
        log:
            "logs/bcftools/{call_type}/merge_{dataset}.log",
        conda:
            "../envs/AmpSeeker-cli.yaml"
        threads: 12
        shell:
            """
            bcftools merge --threads {threads} -o {output.vcf} -O v {input.vcfs} --force-samples 2> {log}
            """

else:

    rule bcftools_merge1:
        input:
            vcfs=expand(
                "results/vcfs/{{call_type}}/{sample}.calls.vcf.gz", sample=samples1
            ),
            idx=expand(
                "results/vcfs/{{call_type}}/{sample}.calls.vcf.gz.tbi", sample=samples1
            ),
        output:
            vcf="results/vcfs/{call_type}/{dataset}.1.vcf",
        log:
            "logs/bcftools/{call_type}/merge1_{dataset}.log",
        conda:
            "../envs/AmpSeeker-cli.yaml"
        threads: 12
        shell:
            """
            bcftools merge --threads {threads} -o {output.vcf} -O v {input.vcfs} --force-samples 2> {log}
            """

    rule bcftools_merge2:
        input:
            vcfs=expand(
                "results/vcfs/{{call_type}}/{sample}.calls.vcf.gz", sample=samples2
            ),
            idx=expand(
                "results/vcfs/{{call_type}}/{sample}.calls.vcf.gz.tbi", sample=samples2
            ),
        output:
            vcf="results/vcfs/{call_type}/{dataset}.2.vcf",
        log:
            "logs/bcftools/{call_type}/merge2_{dataset}.log",
        conda:
            "../envs/AmpSeeker-cli.yaml"
        threads: 12
        shell:
            """
            bcftools merge --threads {threads} -o {output.vcf} -O v {input.vcfs} --force-samples 2> {log}
            """

    rule bgzip2:
        input:
            vcf="results/vcfs/{call_type}/{dataset}.{n}.vcf",
        output:
            vcfgz="results/vcfs/{call_type}/{dataset}.{n,[0-9]}.vcf.gz",
        log:
            "logs/bgzip/{call_type}/{dataset}.{n}.log",
        wrapper:
            "v1.17.4-17-g62b55d45/bio/bgzip"

    rule tabix2:
        input:
            vcfgz="results/vcfs/{call_type}/{dataset}.{n}.vcf.gz",
        output:
            tbi="results/vcfs/{call_type}/{dataset}.{n,[0-9]}.vcf.gz.tbi",
        log:
            "logs/tabix/{call_type}/{dataset}.{n}.log",
        conda:
            "../envs/AmpSeeker-cli.yaml"
        shell:
            """
            tabix {input.vcfgz} 2> {log}
            """

    rule bcftools_merge3:
        input:
            vcf=expand("results/vcfs/{{call_type}}/{{dataset}}.{n}.vcf.gz", n=[1, 2]),
            tbi=expand(
                "results/vcfs/{{call_type}}/{{dataset}}.{n}.vcf.gz.tbi", n=[1, 2]
            ),
        output:
            vcf="results/vcfs/{call_type}/{dataset}.merged.vcf",
            holder=touch("results/vcfs/{call_type}/.complete.{dataset}.merge_vcfs"),
        log:
            "logs/bcftools/{call_type}/merge3_{dataset}.log",
        conda:
            "../envs/AmpSeeker-cli.yaml"
        shell:
            """
            bcftools merge -o {output.vcf} -Ov {input.vcf} --force-samples 2> {log}
            """
