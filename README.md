# 📦 SequoiaR

**sequoiar** is an R package designed to streamline common tasks in genomic data analysis. It provides:

- Encapsulated utility functions (e.g., mlogp calculation, harmonisation, SNP-to-gene mapping)
- A template system for quickly generating various analysis scripts (R and Python)
- Optional automatic setup of external tools (e.g., PLINK and LD reference panel)
- Minimal namespace pollution via qualified function usage (e.g., `dplyr::filter()`)

---

## 🚀 Installation

Install directly from GitHub:

```r
# install.packages("devtools")  # if needed
devtools::install_github("SequoiaGenetics/sequoiar")
```

---

## 🧰 Functionality

The package provides a suite of modular tools for working with GWAS and other summary statistics:

| Function            | Purpose                                  |
| ------------------- | ---------------------------------------- |
| `cal_mlogp()`       | Compute `-log10(p)` using high precision |
| `clump_local()`     | Locally LD-clump SNPs using PLINK        |
| `create_template()` | Copy a template script to your workspace |
| `harmonise()`       | Harmonise alleles across datasets        |
| `map_snp_to_gene()` | Map SNPs to nearest or overlapping gene  |
| `complete_stats()`  | Calculate missing GWAS stats (beta/se/p) |
| `link_LD_local()`   | Add LD information to a data.frame       |

All internal functions use safe, encapsulated calls to avoid namespace conflicts.

---

## 🧪 Templates Included (TODO)

You can generate pre-written scripts using `create_template()`. Templates available include:

| Template Name             | Type   | Description                               |
| ------------------------- | ------ | ----------------------------------------- |
| `manhattan_plot`          | R      | Code for generating a Manhattan plot      |
| `forest_plot`             | R      | Forest plot template                      |
| `meta_analysis`           | R      | Fixed/random effects meta-analysis        |
| `liftover`                | Python | GRCh37/38 position liftOver               |
| `SNP_lookup`              | Python | Find and filter SNPs from reference files |
| `target_analysis`         | Python | PheWAS/PWAS target analysis               |
| `mendelian_randomisation` | R      | Basic two-sample MR pipeline              |

To use:

```r
create_template("manhattan_plot", type = "R", out_dir = "my_analysis/")
```

---

## ⚙️ PLINK and LD Reference Setup

The function `clump_local()` uses PLINK and EUR reference data for LD-based clumping. If not found, these files will be:

- Automatically downloaded to: `~/.yourpkg/plink/`
- Made executable on install
- Used transparently by internal functions

You can override this behavior by setting a custom path:

```r
options(YourPackageName.plink_path = "/custom/path/plink")
```

---

## 📁 Project Structure

```
YourPackageName/
├── R/                   # R function definitions
├── inst/
│   ├── templates/       # R and Python templates
│   └── scripts/         # setup or shell scripts
├── man/                 # Auto-generated documentation
├── DESCRIPTION          # Package metadata
├── NAMESPACE            # Function exports/imports
```

---

## 👥 Authors

- **Lin** – Creator and maintainer (`aut`, `cre`)

---

## 📜 License

MIT © 2025 SequoiaGenetics

---

## ✨ Coming Soon (hopefully/or not)

- Optional Docker image
- Advanced QC workflows
- Full VCF support for SNP/SV merging
