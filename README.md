# Leukaemia Mining

*Can gene activity reveal which type of leukaemia a patient has?*

**[Read the project write-up](https://haqiqak.github.io/leukaemia-mining/site/)**: a one-page explanation of the problem, the method and the results, with charts.

A Python reproduction and extension of the microarray case study in Torgo's *Data Mining with R: Learning with Case Studies* (2nd ed., Chapter 7). The project rebuilds the original classification workflow, then measures how much accuracy is inflated when genes are chosen with information from the test patients.

## The problem

Acute lymphoblastic leukaemia comes in subtypes, each linked to a different genetic change. A lab chip measures the activity of 12,625 genes per patient. From 94 B-cell patients, the task is to name one of four classes: BCR/ABL, ALL1/AF4, E2A/PBX1 or NEG (none of these). There are far more genes than patients, and one class has only five patients.

## What was done

1. **Data.** The `ALL` data set was exported from R and checked against the book: 94 patients, 12,625 genes, class counts 42, 37, 10 and 5. A label-free filter (control, duplicate and low-variance probes) keeps 4,025 genes.
2. **Baseline.** The book's workflow, rebuilt with scikit-learn. In every round of resampling, on the training patients only, an ANOVA filter keeps the genes whose mean differs between classes. Three strategies follow: **s1** uses all of them, **s2** keeps the 30 with the highest random-forest importance, and **s3** groups them into 30 clusters and lets 50 models, each trained on one random gene per cluster, vote. Random forest, SVM and k-NN are then fitted: 34 workflows, 100 rounds each.
3. **Extension.** The same workflows with gene selection done once on all patients, before splitting. Every test patient then helps choose the genes, and the gap to the baseline shows how much accuracy that leak adds.
4. **Checks.** Shuffled labels (one shuffle for all workflows, then a 200-shuffle permutation test for one workflow) show that the honest pipeline finds no signal where there is none. A cell-stage check tests whether ALL1/AF4 is recognised only through its cell stage.
5. **Number of genes.** The book fixes s2 at 30 genes. A sensitivity study scores 5 to 200 genes.

Besides accuracy, results report balanced accuracy (the mean of the per-class recalls), per-class recall and a confusion matrix, because the smallest class is easy to miss.

## Results at a glance

Means over the resampling rounds. Always answering NEG scores 0.447.

| Result | Value |
| --- | --- |
| Best accuracy | 0.852 (`s2_knn3_z`) |
| Best balanced accuracy | 0.885 (`s1_knn3_iqr`), with E2A/PBX1 recall 0.95 |
| Best of the cluster ensemble (s3) | 0.831 (`s3_svm_C100_g0.001`) |
| Accuracy gain when genes are selected outside the loop | +0.039 on average over 34 workflows (s1 +0.014, s3 +0.041, s2 +0.063) |
| Best workflow outside the loop | 0.928 (`s2_knn5_z`, 0.839 inside) |
| Shuffled labels, best workflow, selection inside / outside | 0.448 / 0.578 |
| Permutation test, `s1_knn3_iqr`, 200 shuffles | balanced accuracy 0.877, against at most 0.387 with shuffled labels (p = 0.005) |
| ALL1/AF4 against NEG within one cell stage (18 patients) | accuracy 0.957 |
| k-NN with 100 genes instead of 30 | balanced accuracy 0.899 (0.847 with 30) |

In short: the ranking of methods follows the book (k-NN and the better SVM settings lead, no random forest in the top ten, s2 ahead of s3 and s1), the best model depends on the metric, and choosing genes before splitting the data inflates results and exaggerates the differences between selection strategies.

## Repository layout

```
docs/        problem statement, chapter walkthrough, project page (docs/site)
data/        exported data set (see data/README.md)
notebooks/   01 data check, 02 baseline, 03 leakage, 04 label shuffle, 05 cell stage, 06 gene count
src/         pipeline code, R export scripts, page data builder
results/     tables and figures written by the notebooks
```

Background reading: [docs/02_problem_statement.md](docs/02_problem_statement.md) explains the biology and the task in plain language, and [docs/01_chapter7_walkthrough.md](docs/01_chapter7_walkthrough.md) walks through the book's chapter.

## Run it

Tested with Python 3.13 and R 4.5.3. Random draws use fixed seeds, and a clean rerun of every notebook reproduces the saved results. The one exception is the 18-patient check in notebook 05, which can differ in the third decimal between machines.

```
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

The exported data is already in `data/`. To regenerate it, install the Bioconductor packages `ALL`, `genefilter` and `hgu95av2.db`, then run `Rscript src/export_all.R` and `Rscript src/export_filter.R` (details in [data/README.md](data/README.md)).

Run the notebooks in order. Run times are from a two-core machine.

| Notebook | Content | Run time |
| --- | --- | --- |
| `01_data_check` | Shape, class counts and value range against the book | seconds |
| `02_baseline` | 34 workflows, 100 rounds, selection inside every round | about 17 min |
| `03_leakage_extension` | Selection fitted outside the evaluation loop | about 15 min |
| `04_label_shuffle` | Shuffled-label control and permutation test | about 20 min |
| `05_b1_confound` | Cell-stage check for ALL1/AF4 | about 8 min |
| `06_gene_count` | Number of genes kept by s2 | about 3 min |

## Project page

The write-up is the file `docs/site/index.html`. It opens in any browser, with no server and no internet connection. Every number on it comes from `docs/site/assets/data.js`; after rerunning the notebooks, refresh it with `python src/build_site_data.py`.

It is published with GitHub Pages from the `main` branch and the `/docs` folder (repository Settings, then Pages).

## Differences from the book

- **Resampling.** Repeated stratified hold-out (30% per round) instead of the bootstrap, so accuracies are not directly comparable; the ranking pattern is what is compared. Stratification puts exactly 2 of the 5 E2A/PBX1 patients in every test set, while a bootstrap test set can contain none. It also avoids the duplicated training patients of the bootstrap, which count more than once among the k nearest neighbours.
- **Gene counts.** The first filter keeps 4,025 genes (book: 3,943) and ANOVA keeps 751 (book: 746), because the current annotation package maps more probes.
- **Gene ranking.** The 30-gene ranking uses impurity importance, where the book uses permutation importance.
- **Grid.** 34 workflows instead of 54. Random forests are not run in s3 (50 forests per round is very slow, and a forest is already an ensemble), random forests use 500 trees only, and the SVM grid has four of the book's six settings.
- **k-NN scaling.** The book compares standardised against unstandardised k-NN. Here both variants are standardised, one with mean and standard deviation (`_z`), one with median and IQR (`_iqr`, the book's choice).
- **Low-variance filter.** The book describes it as removing genes whose IQR is below a fifth of the global IQR. The code (kept as in the book) actually removes the least variable 54% of the remaining genes; see [data/README.md](data/README.md).

## Limits

- 94 patients, one seed, one resampling scheme, and a five-patient class whose recall is coarse.
- Paired significance tests are optimistic because rounds share training patients.
- The best workflow is picked after comparing 34 means on the same splits, so its score slightly overstates performance on new patients. A nested evaluation would correct this and was not run.
- This is a learning exercise and not a medical tool.

## References

- Torgo, L. (2017). *Data Mining with R: Learning with Case Studies* (2nd ed.). CRC Press. Chapter 7, pp. 353-381.
- Chiaretti, S., et al. (2004). Gene expression profile of adult T-cell acute lymphocytic leukemia identifies distinct subsets of patients with different response to therapy and survival. *Blood*, 103(7), 2771-2778. Data via the Bioconductor `ALL` package.
- Ambroise, C. and McLachlan, G. J. (2002). Selection bias in gene extraction on the basis of microarray gene-expression data. *PNAS*, 99(10), 6562-6566.
