# AmpSeeker pipeline failures (Sept 2026) — root causes, fixes, and how to verify them

This documents everything that broke across SLURM jobs `118236` → `118264`
of `sbatch run_AmpSeeker_pipeline.sh` (dataset `agvampir-eniyou-agam`,
1631 samples), what was actually wrong in each case, and what was changed to
fix it. Job `118264` (2026-09-18) completed **15/15 steps, 100%**, so the
fixes below are confirmed working end-to-end, not just theorized.

Read this top to bottom if you're picking this up cold — each section
depends on understanding why the previous "fix" wasn't actually enough.

## Summary table

| # | Symptom | Root cause | Fix |
|---|---|---|---|
| 1 | `sample_quality_control`, `snp_dataframe`, `coverage`, `multi_qc` all died together | 4 memory-heavy notebook rules had no (or too-low) `mem_mb` declared, so Snakemake ran them concurrently on a 96GB SLURM allocation | Declared/raised `mem_mb` on all four rules |
| 2 | `coverage` notebook: `pandas.errors.EmptyDataError` | 6 samples (incl. a negative control) had zero aligned reads → empty `per-base.bed.gz` | `coverage.ipynb` now skips samples with empty per-base files |
| 3 | `sample_quality_control` / `snp_dataframe` still OOM-killed (~123GB against a 96GB allocation) even after fix #1 | `mpileup_call_amplicons` had no `-R` region restriction → pileup-called every covered base genome-wide → 25GB, 2.5M-record VCF for 1631 samples, loaded whole into memory by `allel.read_vcf(fields="*")` | Added `-R` to `mpileup_call_amplicons`, same as its `mpileup_call_targets` sibling |
| 4 | PCA crashed for small/low-diversity cohorts (`IndexError`) | `amp.pca()` had an indentation bug: the axis-flip `if` had fallen outside its `for` loop, and it also assumed exactly `n_components` PCs came back | Fixed indentation; now uses however many components `allel.pca` actually returns |
| 5 | PCA outlier step crashed for cohorts with **zero** usable PCA components (e.g. negative controls) | `find_pca_outliers` didn't handle a 0-column input | Added a guard: return "no outliers" instead of crashing |
| 6 | `multi_qc` `TypeError: unsupported operand type(s) for |: 'NoneType' and 'type'` | `multi_qc` is a `wrapper:`-based rule; Snakemake runs the wrapper script using **the executing rule's own conda env's Python**, and that env (`AmpSeeker-qc.yaml`) pinned `python=3.9`, which cannot evaluate the `X | Y` union-type syntax Snakemake 9.x's own source uses internally (needs Python ≥3.10) | Bumped `python=3.9` → `python=3.12` in `workflow/envs/AmpSeeker-qc.yaml` |
| 7 | (Side investigation, not the real cause of #6) | The **orchestrating** `ampseeker` conda env (the one that runs `snakemake` itself) had drifted to Python 3.14, unusually new | Downgraded to Python 3.12 anyway for stability — harmless, but note this did **not** fix `multi_qc`; #6 was the real fix |

## Detailed fixes, in the order they were found

### 1. Memory: declare real `mem_mb` on the heavy notebook/QC rules

`run_AmpSeeker_pipeline.sh` submits with `--cores 4` and `--resources mem_mb=90000`
against a 96GB SLURM allocation (`#SBATCH --mem=96G`). Snakemake only avoids
running memory-heavy jobs concurrently if their `mem_mb` is declared
accurately — undeclared defaults to 0, meaning "free to run any time."

Two of the four notebook rules had `mem_mb=32000` (too low for 1631 samples),
and two had no `mem_mb` at all. Fixed in:

- `workflow/rules/qc.smk:111` — `multi_qc`: `mem_mb=8000`
- `workflow/rules/qc-notebooks.smk:135` — `coverage`: `mem_mb=16000`
- `workflow/rules/qc-notebooks.smk:168` — `sample_quality_control`: `mem_mb=32000` → `48000`
- `workflow/rules/analysis.smk:116` — `snp_dataframe`: `mem_mb=32000` → `48000`

**This alone was not sufficient** — see fix #3. `seff <jobid>` after this fix
still showed ~123GB peak usage against the 96GB allocation, because a single
cell in each notebook (loading the amplicons VCF) used more memory by itself
than any amount of `mem_mb` tuning could schedule around.

### 2. `coverage.ipynb`: handle samples with zero aligned reads

`workflow/notebooks/coverage.ipynb`, the cell building `cov_list`, called
`pd.read_csv` on every sample's `results/coverage/{sample}.per-base.bed.gz`.
Six samples (one negative control, five low-yield samples) had genuinely
empty per-base files (0 rows — mosdepth had nothing to report), which raises
`pandas.errors.EmptyDataError`, not something a `try/except` was there to
catch.

**Fix:** the cell now does:
```python
try:
    cov_df = pd.read_csv(f"{wkdir}/results/coverage/{sample_id}.per-base.bed.gz", sep="\t", header=None)
except pd.errors.EmptyDataError:
    continue  # sample had no aligned reads
```
Their `.regions.bed.gz` files are unaffected (mosdepth always emits one row
per target region even at zero depth), so these samples still appear
everywhere else in the notebook.

### 3. The real OOM cause: unrestricted `mpileup_call_amplicons`

After fix #1, `sample_quality_control` and `snp_dataframe` **still** died
with `nbclient.exceptions.DeadKernelError: Kernel died`, at the exact same
cell every run, regardless of what else was running concurrently. `seff`
showed ~123GB peak against the 96GB SLURM allocation — a single process, not
contention between jobs.

Traced to: `results/vcfs/amplicons/agvampir-eniyou-agam.annot.vcf` was
**25GB** (2,513,464 records, 1631 samples), versus the sibling
`results/vcfs/targets/...annot.vcf` at **3.5MB** for the same samples.
`workflow/lib/shared.py`'s `load_vcf()` calls
`allel.read_vcf(vcf_path, fields="*")`, which loads the entire file into
memory — no chunking, no field filtering. For 25GB of VCF text that becomes
far more than 25GB of numpy arrays.

Root cause of the 25GB file: `workflow/rules/map-call-illumina.smk`,
`mpileup_call_amplicons`, ran:
```
bcftools mpileup -Ov -I -f {params.ref} -a AD --max-depth {params.depth} {input.bam}
```
with **no `-R` region restriction** — unlike its sibling
`mpileup_call_targets` (line 106), which has `-R {params.regions}`. Without
it, `bcftools mpileup` calls genotypes at *every* covered base genome-wide,
not just the panel's positions. Confirmed by inspecting the VCF body:
sequential single-base invariant records (`2L 13807`, `13808`, `13809`, ...)
far outside the ~93-position `ag-vampir` panel.

**Fix** (`workflow/rules/map-call-illumina.smk:140`, `mpileup_call_amplicons`):
```
regions=config["targets"],   # added to params
...
bcftools mpileup -Ov -I -f {params.ref} -R {params.regions} -a AD --max-depth {params.depth} {input.bam}
```

**Caveat:** `config["targets"]` (`config/ag-vampir.bed`) has the panel's 93
*exact* SNP positions, not a padded "whole amplicon insert" region. No
broader amplicon-insert BED exists anywhere in this repo. This means
`amplicons` and `targets` VCFs are now the same genomic scope, until/unless
someone supplies a proper amplicon-insert BED with real padding (the
`snp-dataframe.ipynb` notebook markdown calls the amplicons output
"whole-amplicon SNP data" — that's the scope this rule is *meant* to have).
If you get real primer/insert coordinates, swap them in as `regions=`.

**This required regenerating data**, not just a code fix: it means
re-running `mpileup_call_amplicons` (1631 samples) + the two-tier
`bcftools_merge` + `snp_eff` annotation. Snakemake's default
`--rerun-triggers` (`code input mtime params software-env` — i.e. all of
them) correctly detected the rule's code/params changed and rescheduled all
1631 jobs on the next `sbatch`, with no manual file deletion needed. This
was confirmed via `snakemake -n` (dry run) before resubmitting, and the file
size after job `118264` dropped from 25GB to 2.3MB.

### 4 & 5. Two latent bugs in `amp.pca()` / `find_pca_outliers`

These were only discovered because the standalone recovery script (below)
kept going past individual cohort failures instead of dying on the first
one — the real notebook would have hit both of these too, once it got past
the OOM crash.

**Bug 4** — `workflow/lib/shared.py`, function `pca()`: the axis-flip logic
had an indentation bug (the `if` was outside the `for` loop it belonged to),
and it assumed `allel.pca()` always returns exactly `n_components` columns —
it doesn't, for cohorts with too few informative sites. Fixed at
`workflow/lib/shared.py:239-251`:
```python
coords, model = allel.pca(gn_var, n_components=n_components)
n_components_out = coords.shape[1]   # use what we actually got

for i in range(n_components_out):
    c = coords[:, i]
    if np.abs(c.min()) > np.abs(c.max()):    # now correctly inside the loop
        coords[:, i] = c * -1

pca_df.columns = [f"PC{pc+1}" for pc in range(n_components_out)]
```

**Bug 5** — `find_pca_outliers()`, defined both in
`workflow/notebooks/sample-quality-control.ipynb` (cell with `import allel`)
and ported into `workflow/scripts/recover_amplicon_analysis.py`: for
cohorts where PCA produced **zero** components (e.g. `NEG` — negative
controls, essentially no real genetic signal), `pca_df.filter(like='PC')`
returns a 0-column frame, and `.apply(..., axis=1)` on that raises
`ValueError: Length of values (0) does not match length of index (N)`. Fixed
in both places by returning an empty "no outliers" result when there are no
PC columns to score.

### 6. `multi_qc`: the real fix (Python version in its own conda env)

`multi_qc` uses a snakemake **wrapper** (`v2.2.1/bio/multiqc`), not a plain
`shell:` block. Wrapper rules work like this: Snakemake generates a small
`wrapper.py`, and runs it with **the target rule's own conda env's Python**
(after activating that env) — but injects the orchestrating Snakemake
installation onto that interpreter's import path, because the wrapper needs
`from snakemake.shell import shell` and `snakemake.params`/`.input`/`.output`
to work. The target env does **not** have its own `snakemake` package
installed (confirmed: `AmpSeeker-qc.yaml`'s env has no `snakemake` in
`site-packages`).

`workflow/envs/AmpSeeker-qc.yaml` pinned `python=3.9`. Snakemake 9.x's own
source (`snakemake/io/__init__.py`) uses `algorithm: None | str | Callable`
— PEP 604 union-type syntax, which requires the *executing* interpreter to
be Python ≥3.10, regardless of which env installed the `.py` file. Under
3.9, evaluating that annotation raises exactly the reported error:
```
TypeError: unsupported operand type(s) for |: 'NoneType' and 'type'
```
Reproduced directly: `.snakemake/conda/<multiqc-env-hash>/bin/python3 -c "x: None | str = None"`
raises the identical error on that env's own Python 3.9.

**Fix** (`workflow/envs/AmpSeeker-qc.yaml:10`):
```yaml
dependencies:
  - python=3.12   # was python=3.9
  - fastqc
  ...
```
`--use-conda`'s default rerun triggers include `software-env`, so Snakemake
rebuilds this env automatically once the yaml changes — no manual
`--conda-cleanup-envs` or similar needed.

**A dead end worth recording:** the *first* attempt at this fix downgraded
the orchestrating `ampseeker` env (the one `snakemake` itself runs from)
from Python 3.14 → 3.12, based on seeing the same `TypeError` with file
paths pointing at `ampseeker`'s site-packages in an earlier failed job. That
was a reasonable environment hygiene fix (3.14 is very new and plenty of
packages don't fully support it yet) but it **did not fix `multi_qc`** — a
follow-up `sbatch` still failed identically. The actual executing
interpreter for this specific wrapper is the target rule's own env
(`AmpSeeker-qc.yaml`), not the orchestrator. If you ever see this error
again for a *different* wrapper-based rule, check that rule's own
`envs/*.yaml` Python pin first.

## New tool: standalone recovery script

`workflow/scripts/recover_amplicon_analysis.py` — runs outside
Snakemake/papermill entirely, for recovering usable results without waiting
for a full expensive re-run (e.g. if you hit the 25GB-VCF OOM again before
the rule fix has been applied / re-run). It:

1. Streams a `results/vcfs/amplicons/*.annot.vcf` once (never loads it fully
   into memory) and keeps only records at the panel's BED positions
   (`--pad` to widen the window).
2. Reproduces `snp-dataframe.ipynb`'s Excel export.
3. Reproduces `sample-quality-control.ipynb`'s PCA-outlier detection.

Each stage has its own try/except and logs to `<outdir>/recovery.log`, so
one cohort or one stage failing doesn't take the rest down with it — this is
precisely how bugs #4 and #5 above were found, instead of the whole thing
dying on the first bad cohort.

```bash
conda activate AmpSeeker-python
python workflow/scripts/recover_amplicon_analysis.py   # see --help for overrides
```

## Verification performed

- `snakemake -n` (dry run) after the `-R` fix: confirmed 1631×
  `mpileup_call_amplicons` + `bcftools_merge1/2/3` + `snp_eff` scheduled,
  reason `code has changed since last execution`.
- Direct reproduction of the `multi_qc` `TypeError` using the exact failing
  env's own Python, before and after the `AmpSeeker-qc.yaml` fix.
- Manual `snakemake -R multi_qc ...` run after the fix: completed
  (`Finished jobid: 0 (Rule: multi_qc)`), cascaded through to a full
  jupyterbook rebuild (`build succeeded, 23 warnings`).
- **Real `sbatch` resubmission, job `118264` (2026-09-18 10:00:11):
  15 of 15 steps, 100% done.** `results/vcfs/amplicons/agvampir-eniyou-agam.annot.vcf`
  is now 2.3MB (was 25GB). `results/qc/multiqc/multiqc_report.html` (21MB)
  and `results/ampseeker-results/_build/html/index.html` both exist and are
  current.

## Fixed: cosmetic Sphinx error in `intro.md`

Job `118264`'s log showed one Sphinx build error during the `jupyterbook`
step (non-fatal — jupyter-book still produced `index.html`):
```
docs/ampseeker-results/intro.md:7: ERROR: Document or section may not begin with a transition.
```
`docs/ampseeker-results/intro.md` had a `---` horizontal rule directly after
the "### User guide" heading with no content between them, which docutils
treats as an invalid section start. **Fixed** by removing the stray `---`.
Verified with a direct rebuild (`jupyter-book build --all
docs/ampseeker-results --path-output results/ampseeker-results` from the
`AmpSeeker-jupyterbook` env): the error is gone, warning count dropped from
23 to 22 (the remainder are pre-existing "heading level" style warnings in
the generated notebooks, unrelated to this fix and not addressed here).

## Notes for whoever touches this pipeline next

- **`results/config/metadata.tsv` is rewritten unconditionally on every
  `snakemake` invocation**, including `--unlock` and `--dry-run` — see
  `workflow/Snakefile:57-58`, which calls `load_metadata(config['metadata'],
  write=True)` at parse time, not inside a rule. Editing `config/metadata.tsv`
  or `config/eniyou_agam.yaml` (or any `--configfile`) takes effect on the
  very next invocation automatically; no cache to clear. Side effect: since
  the file's mtime bumps every run regardless of content, some downstream
  rules may be re-flagged under the default `mtime` trigger even when
  nothing meaningful changed.
- Snakemake's default `--rerun-triggers` here is `code input mtime params
  software-env` (i.e. all of them) — editing a rule's `shell:`/`params:`/
  `conda:` env spec is enough to get it to regenerate stale outputs on the
  next run. You do not need to manually delete old output files after a
  rule-code fix.
- If you get a real amplicon-insert BED (broader than
  `config/ag-vampir.bed`'s exact SNP positions), pass it via `regions=` in
  `mpileup_call_amplicons` to restore the intended "whole-amplicon" scope
  described in `snp-dataframe.ipynb`.
