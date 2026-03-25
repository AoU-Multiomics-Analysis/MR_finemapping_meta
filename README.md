# MR Finemapping Meta

A [WDL](https://openwdl.org/) workflow for **Mendelian Randomization (MR)** using fine-mapped QTL credible sets as genetic instruments. The workflow integrates GWAS summary statistics with SuSiE fine-mapping results to estimate the causal effect of molecular traits (e.g., gene expression, protein levels) on complex traits or diseases.

---

## Repository Contents

| File | Description |
|------|-------------|
| `MR.wdl` | WDL workflow definition |
| `mendelian_randomization.R` | Core MR analysis script |
| `Dockerfile` | Docker image for the runtime environment |
| `.dockstore.yml` | Dockstore workflow registration |
| `.github/workflows/docker-image.yml` | CI/CD pipeline to build and push the Docker image |

---

## Workflow Overview

The workflow (`MR.wdl`) runs a single WDL task (`MR`) that calls `mendelian_randomization.R` inside a Docker container. The analysis proceeds as follows:

1. **Load GWAS summary statistics** (munged format) and **SuSiE fine-mapping results** (Parquet format).
2. **Filter variants** to those present in both datasets and adjust posterior inclusion probabilities (PIPs) to account for missing variants.
3. **Format and harmonise** exposure (QTL) and outcome (GWAS) data using [TwoSampleMR](https://mrcieu.github.io/TwoSampleMR/).
4. **Compute per-credible-set composite MR estimates** by weighting SNP-level ratio estimates by their PIPs.
5. **Meta-analyse** composite estimates across credible sets (inverse-variance weighting) to produce a gene-level causal effect estimate.
6. **Run single-SNP MR** for genes with only one credible set containing one variant.
7. **Write results** to a TSV file.

---

## WDL Workflow (`MR.wdl`)

### Inputs

| Parameter | Type | Description |
|-----------|------|-------------|
| `MungedSumstats` | `File` | Munged GWAS summary statistics (see [Data Preparation](#data-preparation) below) |
| `SusieFinemapping` | `File` | SuSiE fine-mapping results in Parquet format (see [Data Preparation](#data-preparation) below) |
| `OutputPrefix` | `String` | Prefix for output file names |
| `Memory` | `Int` | Memory to allocate in GB |
| `QTLSampleSize` | `Int` | Sample size of the QTL study (used to compute standard errors from posterior SDs) |
| `QTLGroup` | `String` | Label for the QTL group (e.g., tissue or cell type name) |

### Outputs

| File | Description |
|------|-------------|
| `{OutputPrefix}_MR.tsv` | MR results table (one row per molecular trait) |

### Runtime

- **Docker image:** `ghcr.io/aou-multiomics-analysis/mr_finemapping_meta:main`
- **Memory:** `{Memory}GB`
- **Disk:** `500 GB SSD`
- **Boot disk:** `25 GB`
- **CPUs:** `1`
- **Zone:** `us-central1-c`

---

## R Script (`mendelian_randomization.R`)

The script is invoked by the WDL task and accepts the same parameters via command-line arguments:

```bash
Rscript /mendelian_randomization.R \
  --MungedSumstats <path> \
  --SusieFinemapping <path> \
  --OutputPrefix <prefix> \
  --QTLSampleSize <n> \
  --QTLGroup <group_label>
```

### Key Functions

| Function | Description |
|----------|-------------|
| `load_gwas_data()` | Reads and standardises GWAS summary statistics; handles missing `FRQ`, `SE`, `OR`, and `BETA` columns |
| `load_finemapping_data()` | Reads a SuSiE Parquet file and renames the `position` column to `pos` |
| `filter_variants_adjust_pips()` | Filters to variants shared with the GWAS data and rescales PIPs within each credible set |
| `format_fm_data_MR()` | Converts fine-mapping output to TwoSampleMR exposure format; computes `se.exposure = posterior_sd / sqrt(sample_size)` |
| `create_MR_input()` | Builds per-SNP ratio estimates and per-credible-set composite causal estimates weighted by PIPs |
| `run_MR()` | Inverse-variance meta-analysis of composite estimates across credible sets; reports effect size, SE, p-value, Q statistic, Q p-value, and I² |
| `LOO_analysis()` / `run_LOO_analysis()` | Leave-one-credible-set-out sensitivity analysis (currently commented out) |

### Output Columns (`{OutputPrefix}_MR.tsv`)

| Column | Description |
|--------|-------------|
| `molecular_trait_id` | Molecular trait identifier (e.g., gene/protein ID) |
| `num_CS` | Number of credible sets used as instruments |
| `num_IV` | Total number of SNPs across all credible sets |
| `cpip` | Sum of adjusted PIPs across all credible sets |
| `meta_eff` | Meta-analysed causal effect estimate |
| `se_meta_eff` | Standard error of the meta-analysed effect |
| `meta_pval` | P-value for the causal effect |
| `Q` | Cochran's Q heterogeneity statistic |
| `Q_pval` | P-value for the Q statistic |
| `I2` | I² heterogeneity measure |
| `group` | QTL group label (from `--QTLGroup`) |
| `trait` | Output prefix (from `--OutputPrefix`) |
| `analysis_type` | `"meta"` (multi-CS/SNP) or `"singlesnp"` (single-CS/SNP) |

---

## Data Preparation

### GWAS Summary Statistics (`--MungedSumstats`)

The GWAS file should be a tab-separated or space-separated text file (read by `data.table::fread`) with the following columns:

| Column | Required | Description |
|--------|----------|-------------|
| `CHR` | Yes | Chromosome (integer, no `chr` prefix) |
| `BP` | Yes | Base-pair position |
| `SNP` | Yes | Variant identifier (e.g., rsID) |
| `A1` | Yes | Effect allele |
| `A2` | Yes | Other allele |
| `P` | Yes | P-value |
| `outcome` | Yes | Trait/outcome label |
| `BETA` | Recommended | Effect size (beta coefficient); used if present |
| `OR` | Alternative | Odds ratio; used if `BETA` is absent |
| `SE` | Recommended | Standard error; computed from `OR`/`BETA` and `P` if absent |
| `FRQ` | Optional | Effect allele frequency; filled with `NA` if absent |

Variants are matched to the QTL fine-mapping data using a composite key of the form `CHR_BP_A1_A2`.

> **Tip:** GWAS summary statistics can be munged using tools such as [ldsc](https://github.com/bulik/ldsc) (`munge_sumstats.py`) or similar pipelines, provided the output is reformatted to the column schema above.

### SuSiE Fine-Mapping Results (`--SusieFinemapping`)

The fine-mapping file must be in **Apache Parquet** format (e.g., produced by [eQTL Catalogue](https://www.ebi.ac.uk/eqtl/) or a custom SuSiE run) with the following columns:

| Column | Description |
|--------|-------------|
| `variant` | Variant identifier in the format `chrCHR_BP_REF_ALT` (e.g., `chr1_925952_G_A`); the `chr` prefix is stripped automatically |
| `position` | Base-pair position (renamed to `pos` internally) |
| `molecular_trait_id` | Molecular trait identifier (e.g., Ensembl gene ID) |
| `cs_id` | Credible set identifier |
| `pip` | Posterior inclusion probability |
| `posterior_mean` | Posterior mean effect size from SuSiE |
| `posterior_sd` | Posterior standard deviation from SuSiE |
| `ref` | Reference allele |
| `alt` | Alternative allele |

> **Tip:** Fine-mapping results can be generated using [SuSiE](https://stephenslab.github.io/susieR/) or downloaded from resources such as the [eQTL Catalogue](https://www.ebi.ac.uk/eqtl/). Only variants with `cs_id` assigned (i.e., belonging to a credible set) are used as instruments.

---

## Docker Image

The runtime environment is built from `mambaorg/micromamba:1.5.3` and includes the following R packages:

- `tidyverse`
- `data.table`
- `optparse`
- `TwoSampleMR`
- `enrichR`
- `gprofiler2`
- `R.utils`
- `arrow`

The Docker image is automatically built and pushed to the GitHub Container Registry (`ghcr.io/aou-multiomics-analysis/mr_finemapping_meta:main`) on every push to `main` via the GitHub Actions workflow in `.github/workflows/docker-image.yml`.

---

## Dockstore

This workflow is registered on [Dockstore](https://dockstore.org/) as a WDL workflow. The registration is configured in `.dockstore.yml` with the primary descriptor path set to `/MR.wdl`.
