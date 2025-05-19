# 📦 SequoiaR

**SequoiaR** is an R package designed to streamline common genomic data tasks with minimal setup and maximum flexibility. It provides:

- Utility functions for SNP annotation, harmonisation, clumping, and statistical completion
- Support for local LD calculations using PLINK
- A lightweight template system for generating R/Python analysis scripts
- Minimal namespace pollution — fully qualified, encapsulated functions

---


## 🚀 Installation

Install directly from GitHub:

```r
# install.packages("devtools")  # if needed
devtools::install_github("SequoiaGenetics/sequoiar")
```

---

## 🔄 Upgrade

To upgrade **SequoiaR** to the latest version from GitHub, simply reinstall it:

```r
devtools::install_github("SequoiaGenetics/sequoiar", force = TRUE)


---

After installation, if you're on Windows, you must manually provide a local copy of PLINK and a reference panel:

```r
sg._copy_plink_win("C:/path/to/your/plink_dir")
```

## 🧰 Functionality

The following high-level functions are exported and ready to use:

| Function                  | Purpose                                                 |
| ------------------------- | ------------------------------------------------------- |
| `sg.clump()`              | LD-based SNP clumping using mlogp col and local PLINK + EUR panel     |
| `sg.harmonise()`          | Harmonise alleles across multiple datasets              |
| `sg.link_LD_local()`      | Annotate SNPs with LD (r²) relative to a lead SNP       |
| `sg.map_snps_to_genes()`  | Map SNPs to overlapping or nearest gene(s)              |
| `sg.complete_stats()`     | Compute missing GWAS stats (beta, se, or p)             |
| `sg._copy_plink_win()`    | Helper to copy PLINK and reference files to package dir |
| `sg.create_gene_region()` | Get gene coordinates from RACER reference               |

All functions are fully encapsulated — use them via sg. prefix to avoid namespace conflicts.

---

## 🧪 Templates Included (Coming Soon)

create_template() will allow you （in the future) to scaffold common analyses from predefined templates:

| Template Name             | Type   | Description                               |
| ------------------------- | ------ | ----------------------------------------- |
| `manhattan_plot`          | R      | Manhattan plot generation                 |
| `forest_plot`             | R      | Forest plot for MR/meta-analysis          |
| `meta_analysis`           | R      | Fixed/random-effects GWAS meta-analysis   |
| `liftover`                | Python | GRCh37 ↔ GRCh38 SNP coordinate conversion |
| `SNP_lookup`              | Python | Reference-based SNP filtering             |
| `target_analysis`         | Python | Target prioritisation (PheWAS/PWAS)       |
| `mendelian_randomisation` | R      | 2-sample MR pipeline                      |


To use:

```r
create_template("manhattan_plot", type = "R", out_dir = "my_analysis/")
```

---

## ⚙️ PLINK and LD Reference Setup

The sg.clump() and sg.link_LD_local() functions use a local PLINK binary and reference panel.

To configure:

- Download plink.exe (Windows) and EUR reference panel files (BED/BIM/FAM)
- Run:

```r
sg._copy_plink_win("C:/path/to/your/plink_dir")
```

This copies the files into:

```
<package install location>/sequoiar/plink/
```

---

## 📁 Project Structure

```
sequoiar/
├── R/                   # Function scripts
├── man/                 # Auto-generated help files
├── inst/
│   └── templates/       # R + Python templates (optional)
├── NAMESPACE            # Auto-generated exports
├── DESCRIPTION          # Package metadata

```

---

## 👥 Authors

- **Lin** – Creator and maintainer (`aut`, `cre`)

---

## 📜 License

MIT © 2025 SequoiaGenetics

---

## ✨ Roadmap
- Docker image with PLINK pre-installed
- VCF integration for SNP/SV merging
- Template engine for automated workflow scaffolding
