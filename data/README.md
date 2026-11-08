# Data

**Source:** the `ALL` dataset (Chiaretti et al., 2004), distributed as the Bioconductor R package `ALL` (version 1.52.0, chip `hgu95av2`).

Chiaretti, S., Li, X., Gentleman, R., et al. (2004). Gene expression profile of adult T-cell acute lymphocytic leukemia identifies distinct subsets of patients with different response to therapy and survival. *Blood*, 103(7), 2771-2778.

**Used here:** the 94 B-cell patients whose subtype is BCR/ABL, ALL1/AF4, E2A/PBX1 or NEG, with 12,625 gene measurements each. Class counts: NEG 42, BCR/ABL 37, ALL1/AF4 10, E2A/PBX1 5. The subset is selected as in Chapter 7 of Torgo, *Data Mining with R: Learning with Case Studies* (2nd ed.): `BT` in its first five levels (B, B1 to B4) and `mol.biol` in its first four levels.

## Files

| File | Content |
| --- | --- |
| `expression.csv` | 94 rows (patients) by 12,625 gene columns, plus a `patient` id column. Values are the expression levels distributed with the package, unchanged. About 20 MB. |
| `genes_kept.csv` | Probe ids kept by the label-free first filter (4,025). |
| `labels.csv` | One row per patient, same order as `expression.csv`: `patient`, `class` (the four subtypes), `BT` (B-cell stage). |

Patient ids are strings with leading zeros (for example `01005`); read them as text, for example `pd.read_csv(..., index_col='patient', dtype={'patient': str})`.

## Export steps

1. Install R (4.5.3 was used) and the Bioconductor package:

   ```r
   install.packages("BiocManager")
   BiocManager::install("ALL")
   ```

2. From the project root, run:

   ```
   Rscript src/export_all.R
   ```

   The script loads `ALL`, applies the subsetting above, transposes the matrix so that patients are rows, and writes `data/expression.csv` and `data/labels.csv`. It prints the package versions used (ALL 1.52.0, Biobase 2.70.0, R 4.5.3).

3. Install the packages for the first filter and run it:

   ```r
   BiocManager::install(c("genefilter", "hgu95av2.db"))
   ```

   ```
   Rscript src/export_filter.R
   ```

   This applies the book's label-free filter (`nsFilter`: control probes, probes without a gene id, duplicate probes, and probes with an IQR below one fifth of the global IQR) and writes the kept probe ids to `data/genes_kept.csv`. With genefilter 1.92.0 and hgu95av2.db 3.13.0 it keeps 4,025 probes; the book reports 3,943. The difference comes from the annotation version (fewer probes lack a gene id in the current database: 942 against 1,151 in the book).

4. Check the export by running `notebooks/01_data_check.ipynb`. It confirms the shape, the class counts, the value range and the absence of missing values.

The data files are small enough for GitHub (under 100 MB). If they are ever removed from the repository, step 2 regenerates them.
