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
notebooks/   analysis notebooks (baseline, extension, comparison)
src/         reusable code
results/     tables and figures produced by the notebooks
slides/      presentation material
```

## Reproducing the results

1. `pip install -r requirements.txt`
2. Obtain the data as described in `data/README.md`.
3. Run the notebooks in `notebooks/` in numeric order.

Detailed steps will be added as the work is completed.

## Status

Problem understood and documented. Methodology, baseline and extension: not started.

## Reference

Torgo, L. (2017). *Data Mining with R: Learning with Case Studies* (2nd ed.). CRC Press. Chapter 7, pp. 353-381.
Original data: Chiaretti et al. (2004), via the Bioconductor `ALL` package.
