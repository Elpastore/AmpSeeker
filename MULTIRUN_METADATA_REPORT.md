# MACARD-1 + MACARD-2: SampleSheet conversion, FASTQ rebuild and metadata — report

Date: 2026-11-20 · Environment: `conda activate AmpSeeker-python` (Python 3.12, pandas 2.2.3, bcl2fastq v2.19.0)

## Result

| Run | Samples in sheet | FASTQ pairs | Removed | bcl2fastq passes |
|---|---|---|---|---|
| MACARD-1 `261001_M08382_0057_000000000-MHMNW` | 384 | 384 | 0 | 6/6, 0 errors, 0 warnings |
| MACARD-2 `261109_M08382_0059_000000000-MHMR4` | 384 | 384 | 0 | 6/6, 0 errors, 0 warnings |
| **Total** | **768** | **768** | **0** | |

- Combined metadata: [`config/metadata.tsv`](config/metadata.tsv) (768 rows, unique `sample_id`, absolute `fq1`/`fq2` paths).
- Removal log: [`config/metadata.removed_samples.tsv`](config/metadata.removed_samples.tsv) (empty).
- FASTQs: `<run>/fastq_files/<sample_id>_1.fastq.gz` / `_2.fastq.gz` (MACARD-1 7.1 GB, MACARD-2 5.2 GB).
- Read counts checked against the on-instrument FASTQs (`Alignment_1/*/Fastq`): identical (e.g. MP_102 28 899 = 28 899, MP_141 4 392 = 4 392, Blank2_Plate8 17 = 17). `gzip -t` passes.
- Metadata filled for every field except 16 control samples (positives, negatives, blanks; `project = GSM_Phase97`), whose location/country/latitude/longitude are blank in the original sheets as well.

## Commands that were run

```bash
conda activate AmpSeeker-python
cd /mnt/nas2/Moussa_Amplicon

# 1. Convert each SampleSheet in place (original kept as backup_SampleSheet.csv)
python AmpSeeker/convert_samplesheet.py -i MACARD-1/261001_M08382_0057_000000000-MHMNW/SampleSheet.csv -m AmpSeeker/resources/exampleSampleSheet.csv
python AmpSeeker/convert_samplesheet.py -i MACARD-2/261109_M08382_0059_000000000-MHMR4/SampleSheet.csv -m AmpSeeker/resources/exampleSampleSheet.csv

# 2. Rebuild FASTQs with bcl2fastq and build the combined metadata
python AmpSeeker/build_multirun_metadata.py \
  MACARD-1/261001_M08382_0057_000000000-MHMNW \
  MACARD-2/261109_M08382_0059_000000000-MHMR4 \
  -m AmpSeeker/config/example-metadata.tsv \
  --rebuild-fastq-dir fastq_files --threads 8 \
  -o AmpSeeker/config/metadata.tsv
```

Running step 2 again takes seconds: complete `fastq_files/` folders are reused. Add `--force-fastq` to demultiplex again.

---

## Problems found and fixed

### `convert_samplesheet.py`

| # | Problem | Effect | Fix |
|---|---|---|---|
| C1 | Only `Sample_ID/Sample_Name/Sample_Plate/Sample_Well` were renamed. `Country`, `Location`, `Taxon`, `Latitude`, `Longitude` kept their capital letter, so the script added **blank** lowercase `country`, `location`, ... columns, and the case-insensitive matcher then picked the blank ones. | All geographic and taxon metadata was **silently wiped** (`,,,,,` at the end of every row). | Every column is now lowercased/snake_cased (`Country → country`). Illumina names are mapped case-insensitively. If duplicate columns appear, they are merged cell by cell, and the script raises an error only on a real conflict. |
| C2 | Columns not in the model sheet (`Description`, `Sample_Project`, the original `Sample_Name`) were dropped. | Useful metadata was lost, e.g. `project` tells controls (`GSM_Phase97`) apart from samples (`GSM_Phase2`). | They are now kept after the model columns (`--no-keep-extra` turns this off). New columns: `description`, `project`, `sample_name_original`. `Sample_Project` is renamed to `project` because bcl2fastq would otherwise write FASTQs into a per-project subfolder. |
| C3 | A second run overwrote `backup_SampleSheet.csv` with the already-converted sheet. | The only original copy would be destroyed. | If a backup exists, the script converts **from** the backup and never overwrites it. Running it twice gives the same result. |
| C4 | All values were read with pandas type inference. | IDs such as `03` or `A01` could be changed. | All values are read as text and stripped of whitespace. |
| C5 | No sanity checks. | Problems went unnoticed. | The script now warns on non-numeric latitude/longitude, errors on duplicate `sample_id`, replaces commas inside values with `;` (AmpSeeker's SampleSheet parser splits on bare commas), and prints how many values are blank per metadata column. |

The converted sheet was checked with AmpSeeker's own `metadata_from_sample_sheet()` ([workflow/rules/metadata.smk](workflow/rules/metadata.smk)): it parses all 384 rows with country, location, taxon, latitude, longitude and plate filled.

### `build_multirun_metadata.py`

| # | Problem | Effect | Fix |
|---|---|---|---|
| B1 | **bcl2fastq crashed** with `free(): double free detected in tcache 2` (exit −6). Cause: this machine allows only **500 open files** per process (`ulimit -n`, hard limit 500), and bcl2fastq keeps one output file open per sample per read (768 for 384 samples) plus its BCL inputs. | This crash is why the earlier `MACARD-1/.../fastq_files` contained only **0-byte** FASTQs. | The script raises the soft limit to the hard limit. If that is still too low, it demultiplexes in **passes** (automatically 64 samples per pass here; `--samples-per-pass N` overrides). Each pass has its own `bcl2fastq.partN.log`, `Stats.partN/` and `Reports.partN/`. bcl2fastq errors now show the end of the log and a hint about the file limit. |
| B2 | Leftover FASTQs from a crashed run (same sample name, different `_S<n>`) were renamed **over** the fresh FASTQs during staging. | The first rebuild lost 20 MACARD-1 samples (MP_141, MP_205, ...). | Before demultiplexing, old bcl2fastq/staged FASTQs in the destination folder are deleted. Two raw files that map to one sample now cause an error instead of an overwrite. A folder is reused only if **every** sample has a non-empty pair. |
| B3 | Metadata required a `sample_id` column and a `plate_name` column, so it only worked on converted sheets. The docstring claimed the backup sheet was preferred while the code did the opposite. | The script crashed with a raw Illumina sheet, and the docstring was misleading. | The script imports `convert_samplesheet`'s normalisation, so raw and converted sheets work the same way. |
| B4 | The bcl2fastq sample sheet was the run's own sheet. | With a raw sheet, `Sample_Project` gives subfolder output that staging could not find, FASTQ names come from the short `Sample_Name` (MP_102, not MP_102_AgVampIR), and index-read FASTQs are written only to be deleted. | A minimal sheet `SampleSheet.bcl2fastq*.csv` is generated with `Sample_ID = Sample_Name = sample_id` plus `index`/`index2`, no project column, and `CreateFastqForIndexReads,0`. Staging also searches subfolders. |
| B5 | Model column `plate` was always blank (the sheet has `plate_name`). Wells like `A01` gave `well_number = 01`. | Plate and well metadata were wrong for AmpSeeker's plate and well plots. | `plate ← plate_name`, and `well_number` loses leading zeros (`A01 → A`, `1`), matching AmpSeeker's example. |
| B6 | `country`, `latitude`, `longitude`, `taxon`, `location` came out blank (same capital-letter problem as C1). | Maps and cohorts were empty. | Fixed by the shared normalisation. |
| B7 | `fq1`/`fq2` were relative to the directory where the script was run (`MACARD-1/...`). | Snakemake runs from `AmpSeeker/`, so it would not find any FASTQ. | All paths are absolute. |
| B8 | `-o` defaulted to `./config/metadata.tsv` relative to the current directory. | Running from `/mnt/nas2/Moussa_Amplicon` wrote a stray `config/` folder outside AmpSeeker. | `-o` and `-m` now default to the `AmpSeeker/config/` folder next to the script. |
| B9 | `removed_samples.log` went to the current directory and was overwritten for each run. | Only the last run's removals were kept. | There is now one combined `<output>.removed_samples.tsv` with `run_id` and `reason` columns. |
| B10 | `--force-fastq` was ignored together with `--rebuild-fastq-dir`. | The flag did nothing. | `--rebuild-fastq-dir` now reuses a complete earlier rebuild, and `--force-fastq` forces a new one. |
| B11 | On-instrument FASTQs could never be matched: MiSeq Local Run Manager renames `MP_102` to `MP-102`. | Without `--rebuild-fastq-dir`, every sample was "missing". | Matching uses the `_S<n>` number (= row order in the sheet) and is accepted only when the name also matches (ignoring `_`/`-` and case). Tested: 768/768 linked as symlinks with no bcl2fastq. This is a fast alternative: leave out `--rebuild-fastq-dir` and use `--reads-dir`. |
| B12 | Extra metadata was dropped. | Useful fields were lost. | New columns after the model columns: `fq1, fq2, run_id, index, index2, well, plate_name, description, project, sample_name_original` (`--no-keep-extra` turns this off). |
| B13 | Python output was buffered when redirected to a log. | Progress did not appear in SLURM logs. | Key progress lines are flushed. |

## Output columns (`config/metadata.tsv`)

`sample_id, taxon, location, country, plate, well_letter, well_number, latitude, longitude` (model order), then `fq1, fq2, run_id, index, index2, well, plate_name, description, project, sample_name_original`.

## Before running AmpSeeker

1. **Config:** use `from-bcl: False`, `fastq: auto: False`, `metadata: config/metadata.tsv`.
2. **Cohort columns:** `config.yaml` (location, taxon, country) and `config_ms.yaml` work with this metadata. **`eniyou_agam.yaml` lists `population`, which does not exist in these sheets.** AmpSeeker's `load_metadata` will fail its assertion. Either remove `population` from `cohort-columns` or add a `population` column.
3. **Controls:** the 16 positive, negative and blank controls are included (`project = GSM_Phase97`, blank location). Drop them from `metadata.tsv` if you don't want them in population analyses. The QC step removes low-read samples anyway.
4. **Paths:** `run_AmpSeeker_pipeline.sh` `cd`s into `/mnt/beegfs/projects/.../AmpSeeker`. The FASTQ paths are absolute under `/mnt/nas2/Moussa_Amplicon/...`, so that node must be able to see `/mnt/nas2`.
5. `demul.sh` still calls the script with old relative paths (`../251211_...`, `ressources/cleaned_reads`). Update it if you reuse it.

## Files changed or created

- Modified: [`convert_samplesheet.py`](convert_samplesheet.py), [`build_multirun_metadata.py`](build_multirun_metadata.py).
- Per run: `SampleSheet.csv` (converted), `backup_SampleSheet.csv` (untouched original), `fastq_files/` (FASTQs, `SampleSheet.bcl2fastq.partN.csv`, `bcl2fastq.partN.log`, `Stats.partN/`, `Reports.partN/`, `Undetermined_S0_*` from the last pass only).
- Created: `config/metadata.tsv`, `config/metadata.removed_samples.tsv`, this report.

---

## Update (2026-11-20, afternoon): pipeline run, conda fix, one-command script

### SLURM job 124625 failed: conda environment could not be created
`run_AmpSeeker_pipeline.sh` stopped before any analysis with
`CreateCondaEnvironmentException ... AmpSeeker-cli.yaml ... LibMambaUnsatisfiableError: htslib requires libdeflate ...`.

**Cause:** `~/.condarc` has `channel_priority: strict`, and `workflow/envs/AmpSeeker-{cli,python,qc,nanopore}.yaml` listed `bioconda` before `conda-forge` (plus `anaconda`/`defaults`). Under strict priority, conda takes each package only from the highest channel that has it, so it picked bioconda's old `libdeflate`, which no `htslib` accepts.

**Fix:** those 4 files now list `conda-forge`, then `bioconda`, the order bioconda requires. Originals are in `workflow/envs/orig_backup/`. `AmpSeeker-snpeff.yaml` (bioconda only) and `AmpSeeker-bcl2fastq.yaml` (dranew, bioconda) are unchanged.

### Reference genome
`moussa_agam.yaml` expects `resources/reference/resources/reference/AgamP4.{fa,gff}`, and these files were missing in this copy. They were downloaded with `resources/reference/download-vectorbase-reference.sh`, run from `resources/reference/`. The FASTA is 286 MB with 153 sequences, and all contigs in `ag-vampir.bed` are present.

### Dry run
`snakemake -n --configfile config/moussa_agam.yaml` plans **8,476 jobs** for all 768 samples, with no errors.

### Re-doing all of this later: `prepare_ampseeker_runs.sh`
One script runs every step above in order. Each step skips work that is already done, so re-running it is safe and fast:

1. convert each SampleSheet (`convert_samplesheet.py`)
2. bcl2fastq in batches + combined metadata (`build_multirun_metadata.py`)
3. download the reference if the config's FASTA/GFF is missing, and check that the BED contigs are in the FASTA
4. fix conda channel order in `workflow/envs/*.yaml` where needed
5. build all conda envs (`snakemake --conda-create-envs-only`)
6. Snakemake dry run (log in `logs/dryrun.<time>.log`)
7. with `--submit`: `sbatch run_AmpSeeker_pipeline.sh` (refuses if that script uses a different `--configfile`)

```bash
cd /mnt/nas2/Moussa_Amplicon
# new runs (bash, or sbatch to run it on a compute node)
bash AmpSeeker/prepare_ampseeker_runs.sh MACARD-3/<run_folder> MACARD-4/<run_folder>
# then prepare + submit in one go
bash AmpSeeker/prepare_ampseeker_runs.sh --submit MACARD-1/<run> MACARD-2/<run>
# only re-check envs + dry run
bash AmpSeeker/prepare_ampseeker_runs.sh --from-step 4
```

Options: `-c/--config` (default `config/moussa_agam.yaml`), `-o/--metadata` (default: the config's `metadata:`), `--fastq-mode rebuild|instrument`, `--force-fastq`, `--threads`, `--on-duplicate warn|error|suffix` (default `error`), `--skip-envs`, `--from-step N`, `--submit`. Run `--help` for details.
