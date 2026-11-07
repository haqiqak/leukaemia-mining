# The Problem in Plain Words: Telling Leukaemia Subtypes Apart from Gene Activity

*A problem statement for the Leukaemia Mining project, based on a case study from Torgo's* Data Mining with R *(Chapter 7). It explains the biology in everyday language, states exactly what is being tested, what a data-mining solution involves, and the planned approach. Specific methods are deliberately left open.*

---

## 1. The whole thing in one paragraph

Some patients have a blood cancer called acute lymphoblastic leukaemia (ALL). ALL is not one single disease. It comes in subtypes, each caused by a different change in the patient's genetic material, and the subtype matters for how a patient is treated. A lab can measure how active about 12,000 genes are in a patient's sample. The task is to teach a computer to look at those 12,000 activity numbers and name the patient's subtype. Only 94 patients are available to learn from, which is the central difficulty.

---

## 2. The biology, decoded

No medical training is needed for this project. About ten words cover it.

| Term | Plain meaning |
| --- | --- |
| **Cell** | The tiny building block of the body. Blood contains several kinds. |
| **Leukaemia** | Cancer of the blood-making system, where abnormal blood cells multiply uncontrollably in the bone marrow and blood. |
| **Acute** | Fast-developing, not slow. |
| **Lymphoblastic** | The runaway cells are immature versions ("blasts") of a white blood cell type called lymphocytes, which fight infection. |
| **B-cell and T-cell** | Two kinds of lymphocyte. B-cell ALL starts in B-cells, T-cell ALL in T-cells. The original study has both; this project uses only B-cell ALL. |
| **Gene** | A stretch of DNA that holds the instructions for making one protein. Think of it as one recipe in a huge cookbook. |
| **Gene expression (activity level)** | How much a gene is "switched on" right now in a cell. A recipe being cooked a lot has high expression; one barely used has low expression. |
| **Microarray** | A lab chip that measures the expression of thousands of genes in one go, from one sample. The output is one number per gene. |
| **Mutation / genetic abnormality** | A change in the DNA. The book calls the ones used here *cytogenetic abnormalities*, meaning changes visible at the level of chromosomes. |
| **Chromosome** | A long packaged strand of DNA. Humans have 23 pairs. |
| **Translocation** | Two chromosomes swap pieces. The swap can stick two genes together into a new "fused" gene. |
| **Subtype** | A named category of the same disease. Here each subtype is named after the genes involved in the swap. |

### The four classes to predict

| Class | Meaning in plain terms |
| --- | --- |
| **BCR/ABL** | Two genes, BCR and ABL, have been fused by a translocation (the book's data labels this as t(9;22)). |
| **ALL1/AF4** | Two genes, ALL1 and AF4, fused (labelled t(4;11) in the book's data). |
| **E2A/PBX1** | Two genes, E2A and PBX1, fused. (The chromosome swap behind it, t(1;19), is general background and is not stated in the book.) |
| **NEG** | "Negative". Tests found none of these abnormalities. |

Why does the subtype matter? Different subtypes behave differently and can call for different treatment. That is general background, not something the book discusses in detail; it is the reason the question is worth asking.

A fused gene changes which genes in the cell are busy. So the idea behind the whole case is: **each subtype should leave a recognisable pattern in the gene activity numbers.** The book states this as its key hypothesis.

---

## 3. What the data looks like

Picture a spreadsheet.

- **Each row is one patient's sample** (94 rows).
- **Each column is one gene** (12,625 columns before any cleaning).
- **Each cell is a number**: how active that gene is in that patient.
- **One extra column holds the answer**: the patient's subtype (BCR/ABL, ALL1/AF4, E2A/PBX1 or NEG), as determined by lab tests.

One small oddity: in biology files this table is usually stored sideways, with genes as rows and patients as columns. The book's data arrives that way, so it is flipped before modelling.

### How many patients per class

| Class | Patients |
| --- | --- |
| NEG | 42 |
| BCR/ABL | 37 |
| ALL1/AF4 | 10 |
| E2A/PBX1 | 5 |
| **Total** | **94** |

These counts come from the book's cross-table of the B-cell patients. The classes are very uneven. E2A/PBX1 has only five patients.

---

## 4. What exactly is being tested

**The question:** given only a patient's 12,625 gene activity numbers, can the correct class out of four be predicted, for patients the model has never seen?

- **Input:** gene activity numbers (all of them, or a chosen subset).
- **Output:** one of four labels.
- **What counts as success:** the share of new patients labelled correctly (accuracy). The book uses plain accuracy because it has no information about which mistakes are worse.
- **What "never seen" means:** a model that only remembers its training patients is useless. Testing must use patients that were held out of training.

A useful reference point: always answering "NEG" would be right for 42 of 94 patients, about 45%. Any real model must clearly beat that.

---

## 5. Why this is hard

**5.1 Far more columns than rows.** There are 12,625 genes and 94 patients, about 130 columns per patient. With that many columns, some genes look predictive purely by luck. A model can latch onto such lucky genes and look excellent on the patients it studied, then fail on new ones. This is **overfitting**, and it is the heart of the problem.

**5.2 Most genes are irrelevant.** Many genes barely change between patients, and some are lab control probes. They add noise and no information.

**5.3 Very few patients.** 94 patients cannot be split into a large training set and a large test set without leaving both too small. The way accuracy is measured needs care.

**5.4 Uneven classes.** Five patients in one class is thin evidence. Overall accuracy can look fine while that class is almost always wrong.

**5.5 Information leakage (the quiet trap).** If the answers are used to pick the best genes using *all* patients, and the model is then tested on those same patients, the test patients have already influenced the choice. The reported accuracy is then too good to be true. The book is careful to avoid this.

**5.6 Different scales.** Some genes have larger numbers than others. Methods that compare patients using distances let big-number genes dominate unless scales are handled.

---

## 6. What the book did, in brief

The full detail is in [01_chapter7_walkthrough.md](01_chapter7_walkthrough.md). In short:

1. Removed genes that clearly cannot help (control probes, duplicates, genes that hardly vary): 12,625 down to 3,943.
2. Kept only genes whose average activity differs between the classes: down to 746.
3. Optionally cut further to a small set, in two ways: top 30 by a random forest ranking, or many models each using one gene from each of 30 gene clusters, voting together.
4. Trained random forests, SVMs and k-nearest neighbours.
5. Measured accuracy with 100 rounds of bootstrap resampling, doing the label-using gene selection inside each round.
6. Best result: about 81% accuracy, with several settings statistically indistinguishable.

---

## 7. What a data-mining solution involves

Whoever takes this on, in any language, goes through these stages and justifies each choice.

| Stage | What it involves |
| --- | --- |
| **Understand the problem** | State the biology question and what a good answer looks like. |
| **Understand the data** | Know the shape, the value range, the class balance, and what each column means. |
| **Prepare the data** | Orient the table correctly, restrict to the right patients, remove obviously useless columns, decide how to handle scale. |
| **Reduce the columns** | Get from thousands of genes to a number a model can handle, without leaking answers into the test. |
| **Model** | Choose learning methods and justify them. |
| **Evaluate honestly** | Measure performance only on patients not used for training, with a method suited to a tiny dataset, against a known baseline. |
| **Interpret** | Explain what the numbers mean, where the model fails, and why. |

---

## 8. Planned approach

The project works in two layers.

**Layer 1: Baseline.** Reproduce the book's workflow well enough to get comparable results. This confirms the problem, the data and the evaluation are understood.

**Layer 2: One extension.** Change one thing on purpose, measure the effect against the baseline on the same data splits, and explain what changed and why. The extension does not need to beat the baseline; it needs to be evaluated and interpreted properly.

### Order of work

1. Finish the description of the problem (this document and the chapter walkthrough).
2. Export the data from R to plain files usable in Python.
3. Explore the data and check it matches the book (94 patients, 12,625 genes, the class counts above).
4. Decide the baseline.
5. Decide the extension.
6. Run, compare, interpret.
7. Summarise the results.

### Decisions deliberately left open

Nothing below is chosen yet. Choices will be recorded here as they are made.

- Which parts of the book's workflow to reproduce, and how many of its 54 setups.
- How to get the data out of Bioconductor and into Python.
- How to scale the genes.
- How to select genes, and how many to keep.
- Which models to train.
- How to evaluate (bootstrap, cross-validation, leave-one-out), and which measures besides accuracy.
- What the extension will be.
- How to handle the five-patient class.

---

## 9. Things to keep honest about

- This is a learning exercise. A model built on 94 patients is not a medical tool and must not be presented as one.
- Results will not match the book's exactly: different language, random draws and library defaults. The comparison should look at which setups do well, not exact decimals.
- The biology background beyond the book (the role of subtypes in treatment, the chromosome numbers not given in the book) is general knowledge and should be checked against a reliable source before being cited as fact.

---

## 10. Data mining words used above

| Word | Meaning |
| --- | --- |
| **Supervised classification** | Learning to assign labels from examples where the right label is known. |
| **Features / predictors** | The input columns (here, genes). |
| **Target** | The thing to predict (here, the subtype). |
| **Overfitting** | Fitting noise and quirks of the training data so well that performance on new data suffers. |
| **Feature selection** | Choosing which columns to keep. |
| **Filter** | Feature selection that looks at columns on their own, without a model. |
| **Wrapper** | Feature selection that repeatedly tries a model on candidate column sets. |
| **Data leakage** | Information from the test data sneaking into training or selection, making results look better than they are. |
| **Bootstrap** | Resampling with replacement to estimate performance when data is small. |
| **Accuracy** | Share of predictions that are correct. |
| **Baseline** | A simple reference result to compare against. |
