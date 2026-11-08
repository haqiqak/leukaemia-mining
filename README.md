# Leukaemia Mining

Can gene activity reveal which type of leukaemia a patient has?

A Python reproduction and extension of the microarray case study from Torgo's *Data Mining with R: Learning with Case Studies* (2nd ed., Chapter 7). The project reproduces the original classification workflow and adds one evaluated extension.

## The problem in one paragraph

Acute lymphoblastic leukaemia (ALL) comes in subtypes, each linked to a different genetic abnormality. A lab chip measures the activity of about 12,000 genes per patient. Using 94 B-cell patients, the goal is to predict which of four classes (BCR/ABL, ALL1/AF4, E2A/PBX1, NEG) a patient belongs to. The difficulty is many genes and very few patients.

## Documentation

| File | Contents |
| --- | --- |
| [docs/02_problem_statement.md](docs/02_problem_statement.md) | Plain-language explanation of the biology, the task, and why it is hard |
| [docs/01_chapter7_walkthrough.md](docs/01_chapter7_walkthrough.md) | Section-by-section walkthrough of the book's chapter (original wording, not a copy) |

## Repository layout

```
docs/        problem statement and chapter walkthrough
data/        the dataset (see data/README.md)
notebooks/   analysis notebooks 01 to 05
src/         reusable code
results/     tables and figures produced by the notebooks
slides/      presentation material
```

## Reproducing the results

Tested with Python 3.13.5 and R 4.5.3 on Windows. All random draws use fixed seeds; a clean rerun of every notebook reproduced the saved result files exactly.

1. Create an environment and install the Python packages (versions are pinned):

   ```
   python -m venv .venv
   .venv\Scripts\activate        # Linux and macOS: source .venv/bin/activate
   pip install -r requirements.txt
   ```

2. Export the data from R (one time; the exported files are also included in `data/`). Install the Bioconductor packages `ALL`, `genefilter` and `hgu95av2.db`, then run, from the project root:

   ```
   Rscript src/export_all.R
   Rscript src/export_filter.R
   ```

   Details and package versions are in [data/README.md](data/README.md).

3. Run the notebooks in `notebooks/` in numeric order (from the `notebooks/` folder, or with `jupyter nbconvert --to notebook --execute --inplace <notebook>`):

   | Notebook | Content | Run time |
   | --- | --- | --- |
   | `01_data_check.ipynb` | Shape, class counts and value range against the book | seconds |
   | `02_baseline.ipynb` | Baseline: 24 workflows, 100 rounds, selection inside every round | about 13 min |
   | `03_leakage_extension.ipynb` | Extension: selection fitted outside the evaluation loop | about 13 min |
   | `04_label_shuffle.ipynb` | Control: shuffled labels, selection inside and outside | about 12 min |
   | `05_b1_confound.ipynb` | Check of the B-cell stage confound for ALL1/AF4 | about 7 min |

   Notebook 03 reads the scores written by notebook 02, and notebook 04 reuses its saved scores if they exist. Reusable code is in `src/pipeline.py`. Figures and score tables are written to `results/`.

## Method in brief

Per round of repeated stratified resampling (100 rounds, 30% of patients held out, seed 0), using the training patients only: ANOVA filter (p < 0.01) on the 4,025 probes left by the label-free first filter, then either all ANOVA genes (s1) or the 30 genes with the highest random-forest importance (s2). Random forest (2 settings), SVM (4 settings, standardised) and k-NN (k = 3, 5, 7; z-score or median/IQR scaling) are fitted and scored on the held-out patients: 24 workflows in total. Accuracy, balanced accuracy, per-class recall and the confusion matrix are reported.

## Results summary

All values are means over the rounds, computed in the notebooks. The always-NEG accuracy is 0.447.

| Result | Value | Notebook |
| --- | --- | --- |
| Best accuracy (selection inside the loop) | 0.852, `s2_knn3_z` (balanced accuracy 0.847) | 02 |
| Best balanced accuracy | 0.885, `s1_knn3_iqr` (accuracy 0.834, E2A/PBX1 recall 0.95) | 02 |
| E2A/PBX1 recall of the best-accuracy workflow | 0.695 | 02 |
| Best workflow with selection outside the loop | 0.928, `s2_knn5_z` (0.839 inside) | 03 |
| Mean accuracy gain from selecting outside the loop, 24 workflows | +0.039 (+0.063 for s2, +0.014 for s1) | 03 |
| Shuffled labels, selection inside: mean accuracy, best workflow | 0.380, 0.448 | 04 |
| Shuffled labels, selection outside: mean accuracy, best workflow | 0.502, 0.578 | 04 |
| ALL1/AF4 against NEG within stage B1 only (18 patients) | accuracy 0.960 (`s2_knn3_z`), larger class 0.556 | 05 |
| Stage prediction from genes: balanced accuracy, five stages / B1 against rest | 0.420 (chance 0.20) / 0.774 (chance 0.50) | 05 |

Main findings: the workflow ranking follows the book (k-NN and SVM ahead of random forest, 30-gene selection mostly ahead of all genes); accuracy and balanced accuracy rank workflows differently because of the five-patient E2A/PBX1 class; selecting genes outside the loop inflates accuracy, and the shuffled-label control shows that selection inside the loop does not.

## Differences from the book

- Resampling is repeated stratified hold-out (30%), not the bootstrap, so accuracies are not directly comparable. The ranking pattern is what is compared.
- The first filter keeps 4,025 probes (book: 3,943) and ANOVA keeps 751 (book: 746), because the current annotation package maps more probes.
- The random-forest ranking uses impurity importance; the book uses permutation importance.
- 24 workflows instead of 54: the clustering ensemble (s3) is not reproduced, the SVM grid has four of six settings and the random-forest grid has two of six.
- Paired p-values are optimistic because the rounds share training patients.

## Status

Complete: data export and check, baseline, extension, label-shuffle control and stage-confound check.

## Reference

Torgo, L. (2017). *Data Mining with R: Learning with Case Studies* (2nd ed.). CRC Press. Chapter 7, pp. 353-381.
Original data: Chiaretti et al. (2004), via the Bioconductor `ALL` package.
