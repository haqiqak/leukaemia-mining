# Leukaemia Mining

*Can gene activity reveal which type of leukaemia a patient has?*

**[Read the project write-up](docs/site/index.html)**: a one-page explanation of the problem, the method and the results, with charts. Open the file in any browser, or publish it with GitHub Pages (see [Project page](#project-page)).

A Python reproduction and extension of the microarray case study in Torgo's *Data Mining with R: Learning with Case Studies* (2nd ed., Chapter 7). The project rebuilds the original classification workflow, then adds an experiment on how easily such results become too good to be true.

## The problem

Acute lymphoblastic leukaemia comes in subtypes, each linked to a different genetic change. A lab chip measures the activity of 12,625 genes per patient. From 94 B-cell patients, the task is to name one of four classes: BCR/ABL, ALL1/AF4, E2A/PBX1 or NEG (none of these). The difficulty is far more genes than patients, and one class with only five patients.

## What was done

1. **Data.** The `ALL` data set was exported from R and checked against the book: 94 patients, 12,625 genes, class counts 42, 37, 10 and 5.
2. **Baseline.** The book's workflow was rebuilt in Python: an ANOVA gene filter, then either all remaining genes or the 30 best by random-forest importance, then random forest, SVM and k-NN. 24 workflows were scored over 100 resampling rounds, with gene selection redone inside every round so the test patients never influence it.
3. **Extension.** The same workflows were run with gene selection done once on all patients, to measure how much accuracy that leak adds.
4. **Checks.** A shuffled-label control shows the honest pipeline invents no signal. A cell-stage check tests whether ALL1/AF4 is only being recognised through its cell stage.

Besides accuracy, results report balanced accuracy, per-class recall and a confusion matrix, because the smallest class is easy to miss.

## Results at a glance

Means over the resampling rounds. Always answering NEG scores 0.447.

| Result | Value |
| --- | --- |
| Best accuracy | 0.852 (`s2_knn3_z`) |
| Best balanced accuracy | 0.885 (`s1_knn3_iqr`), with E2A/PBX1 recall 0.95 |
| Best workflow when genes are selected outside the loop | 0.928 (0.839 inside) |
| Average accuracy gain from that leak, 24 workflows | +0.039 (+0.063 for the 30-gene strategy) |
| Shuffled labels, selection inside / outside: best workflow | 0.448 / 0.578 |
| ALL1/AF4 against NEG within one cell stage (18 patients) | accuracy 0.960 |

In short: the ranking of methods follows the book (k-NN and the better SVM settings lead, no random forest in the top ten), the best model depends on the metric, and choosing genes before splitting the data inflates results.

## Repository layout

```
docs/        problem statement, chapter walkthrough, project page (docs/site)
data/        exported data set (see data/README.md)
notebooks/   01 data check, 02 baseline, 03 leakage, 04 label shuffle, 05 cell stage
src/         pipeline code, R export scripts, page data builder
results/     tables and figures written by the notebooks
```

Background reading: [docs/02_problem_statement.md](docs/02_problem_statement.md) explains the biology and the task in plain language, and [docs/01_chapter7_walkthrough.md](docs/01_chapter7_walkthrough.md) walks through the book's chapter.

## Run it

Tested with Python 3.13.5 and R 4.5.3. Random draws use fixed seeds, and a clean rerun of every notebook reproduced the saved results exactly.

```
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

The exported data is already in `data/`. To regenerate it, install the Bioconductor packages `ALL`, `genefilter` and `hgu95av2.db`, then run `Rscript src/export_all.R` and `Rscript src/export_filter.R` (details in [data/README.md](data/README.md)).

Run the notebooks in order:

| Notebook | Content | Run time |
| --- | --- | --- |
| `01_data_check` | Shape, class counts and value range against the book | seconds |
| `02_baseline` | 24 workflows, 100 rounds, selection inside every round | about 13 min |
| `03_leakage_extension` | Selection fitted outside the evaluation loop | about 13 min |
| `04_label_shuffle` | Shuffled-label control | about 12 min |
| `05_b1_confound` | Cell-stage check for ALL1/AF4 | about 7 min |

## Project page

`docs/site/index.html` needs no server and no internet connection. Every number on it comes from `docs/site/assets/data.js`; after rerunning the notebooks, refresh it with `python src/build_site_data.py`.

To publish it, open the repository's Settings, then Pages, choose "Deploy from a branch", select the `main` branch and the `/docs` folder, and save. The page is then served at `https://<account>.github.io/<repository>/` (the root redirects to `/site/`).

## Differences from the book

- Resampling is repeated stratified hold-out (30% per round) instead of the bootstrap, so accuracies are not directly comparable. The ranking pattern is what is compared.
- The first filter keeps 4,025 genes (book: 3,943) and ANOVA keeps 751 (book: 746), because the current annotation package maps more probes.
- The 30-gene ranking uses impurity importance, where the book uses permutation importance.
- 24 workflows instead of 54: the clustering ensemble is not reproduced, and the SVM and random-forest grids are smaller.

## Limits

94 patients, one seed, one resampling scheme, and a five-patient class whose recall is coarse. Paired significance tests are optimistic because rounds share training patients. This is a learning exercise and not a medical tool.

## Reference

- Torgo, L. (2017). *Data Mining with R: Learning with Case Studies* (2nd ed.). CRC Press. Chapter 7, pp. 353-381.
- Chiaretti, S., et al. (2004). Gene expression profile of adult T-cell acute lymphocytic leukemia identifies distinct subsets of patients with different response to therapy and survival. *Blood*, 103(7), 2771-2778. Data via the Bioconductor `ALL` package.
