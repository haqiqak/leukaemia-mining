# Chapter 7 Walkthrough: Classifying Microarray Samples

*Source: Luís Torgo, **Data Mining with R: Learning with Case Studies**, 2nd edition (CRC Press, 2017), Chapter 7, pages 353 to 381.*

This is a detailed walkthrough of the chapter, written in original wording, following the chapter's own section order and keeping all of its figures, decisions and results. It is not a copy of the book text. It is meant to be read alongside the book; the book itself is not included in this repository.

---

## 0. What the chapter is for

The chapter is the book's fourth and last case study. It deals with a bioinformatics problem: given a measurement of how active thousands of genes are in a leukaemia patient's sample, predict which genetic abnormality that patient has.

Its main teaching focus is **feature selection**: how to shrink a huge number of columns (genes) down to a useful few. The chapter also introduces three other things:

- the **k-nearest neighbours** classifier
- **bootstrap** estimates of accuracy
- some new variants of **ensemble** models (several models voting)

The tools come from **Bioconductor**, a large collection of R packages for biological data. Python equivalents replace these in this project, but the logic is the same.

---

## 7.1 Problem description and objectives

### 7.1.1 Background on microarray experiments

The author first translates biology jargon into data mining terms:

| Biology term | Data mining term |
| --- | --- |
| Sample (a patient's measured tissue) | Observation, case, row |
| Gene | Variable, predictor, column |
| Expression level (how active a gene is) | The value in a cell |
| Disease type | Target variable (the class to predict) |

Key points the chapter makes:

- **Microarray experiments** measure the expression level of a very large set of genes for each sample.
- If each sample has a known disease state, a function can be learned from gene expression to disease state. That is a **supervised classification** task.
- The **key hypothesis**: different disease types are associated with different gene expression patterns, and measuring those patterns should allow accurate prediction of a patient's disease type.
- The raw output of the lab technology (short oligonucleotide arrays) is an image. After several pre-processing steps (image analysis, normalization and so on) it becomes a table of expression levels. Those steps are done by other Bioconductor packages and are **not covered** in the chapter. The chapter starts from the finished table.
- Two typical features of this kind of data:
  1. There are **many more predictors (genes) than observations (patients)**.
  2. The table is stored **transposed** compared with the usual layout: rows are genes, columns are samples.
- Other information about patients (sex, age and so on) can exist as extra covariates.

### 7.1.2 The ALL dataset

- It comes from a study on **acute lymphoblastic leukaemia** (Chiaretti et al., 2004).
- It has microarray samples from **128 individuals**.
- Two tumour types: **T-cell ALL (33 samples)** and **B-cell ALL (95 samples)**.
- The chapter focuses on **B-cell ALL only**.
- Within B-cell ALL there are different **cytogenetic abnormalities** (mutation types): ALL1/AF4, BCR/ABL, E2A/PBX1, p15/p16, and patients with no abnormality (called NEG).
- **p15/p16 is discarded** because there is only one sample of it.
- The final target variable is nominal with **4 classes**. This is a multi-class classification task.

---

## 7.2 The available data

### Getting the data

- The ALL dataset is part of Bioconductor, so it is **not** shipped in the book's R package. It is installed through the basic Bioconductor packages, followed by the data package `ALL`.
- Loading it creates an `ExpressionSet` object, Bioconductor's container for a microarray experiment.

### What is inside the ALL object

| Part | Content |
| --- | --- |
| Assay data | **12,625 features (genes)** by **128 samples**, one expression matrix |
| Phenotype data | Information about samples: 21 covariates such as patient ID, sex, age, B-cell or T-cell type (`BT`), molecular biology result (`mol.biol`), remission, transplant, follow-up |
| Annotation | `hgu95av2`, the chip type, linking genes to biomedical databases |

### Two covariates that matter

- **`BT`**: the type of leukaemia. Values B, B1, B2, B3, B4 are B-cell; T, T1 to T4 are T-cell.
- **`mol.biol`**: the molecular result, which is the **target**. Counts over all 128 patients: ALL1/AF4 = 10, BCR/ABL = 37, E2A/PBX1 = 5, NEG = 74, NUP-98 = 1, p15/p16 = 1.

### Choosing the subset

- Keep patients whose `BT` is one of the first five levels (the B-cell stages) **and** whose `mol.biol` is one of the first four levels.
- This gives **94 samples**.
- Class counts in this subset (derived from the book's cross-table): **NEG = 42, BCR/ABL = 37, ALL1/AF4 = 10, E2A/PBX1 = 5**.
- After subsetting, the factor levels of `BT` and `mol.biol` are refreshed so unused levels are dropped.
- The author suggests saving this subset to a file so the pre-processing need not be repeated.

### 7.2.1 Exploring the dataset

- The expression matrix is **12,625 rows by 94 columns**.
- The main challenge is stated plainly: **far too many variables for the number of cases**. Most modelling techniques struggle with this.
- Summary of all expression values: minimum about 1.99, median about 5.47, mean about 5.62, maximum about 14.04. **Most values lie between 4 and 7.**
- A histogram with the quartiles, median and a robust centre estimate called the **shorth** (mean of the middle 50% of values) shows values packed around the centre with a few large ones.
- The distribution of expression values is **almost the same for each mutation group** and similar to the global one. So no single glance at overall values separates the classes. The signal must be in specific genes.

---

## 7.3 Gene (feature) selection

### The two families of feature selection

| Family | How it works | Cost |
| --- | --- | --- |
| **Filters** | Use statistical properties of the features alone, in one pass | Cheap |
| **Wrappers** | Include the model in the search, trying candidate gene sets repeatedly until some stop rule | Expensive |

The chapter uses **filters** only.

### 7.3.1 Simple filters based on distribution properties

Goal: remove genes that cannot possibly help.

- Some genes are not expressed at all or barely vary. They cannot distinguish samples.
- Some probes on the chip are **controls**. On this chip type (Affymetrix U95Av2) their names start with `AFFX` and they can be removed.
- To see variability, each gene is summarised by its **median** and **inter-quartile range (IQR)**, then plotted. The plot shows a large share of genes with an IQR near zero.
- **Caveat the author states:** judging genes one at a time risks discarding a gene that is useless alone but useful in combination with others. Methods that account for dependencies between features exist (for example RELIEF), but gene-by-gene filtering is the common practice for data this size.
- **Rule used:** remove genes whose variability (IQR) is below **one fifth of the global IQR** of all expression values.
- A function from the `genefilter` package applies this together with the control-probe removal.

**Result: 12,625 genes reduced to 3,943.** The filter log breaks it down:

| Removed because | Genes |
| --- | --- |
| Duplicates (several probes for the same gene) | 2,858 |
| Low variability | 4,654 |
| No gene identifier (Entrez ID) | 1,151 |
| Control probes (`AFFX`) | 19 |

The author notes 3,943 is still far too many for 94 observations.

### 7.3.2 ANOVA filters

Idea: a gene whose expression has **the same average in all four classes** cannot help tell the classes apart.

- For each gene, an **ANOVA test** compares the mean expression across the four mutation groups. This uses the labels, unlike step 7.3.1.
- Genes whose means are confidently equal across classes are discarded. The significance level used is **p = 0.01**.
- **Result: 746 genes remain.**
- The remaining genes sit on very different scales (different medians and IQRs). This matters for distance-based methods, since genes with larger values would dominate. The author therefore **leaves standardization to the modelling stage**, applying it only for models that need it.

### 7.3.3 Filtering using random forests

The 746 genes are manageable, but the author asks whether a more "standard" dimensionality does better. So a further cut is tried.

- A **random forest** is trained and asked for **variable importance**, here the average drop in classification accuracy when a gene is removed.
- Genes are ranked by that score and the **top 30** are kept.
- The median expression of these 30 genes, split by mutation, shows clear differences between classes for many of them, which supports their usefulness. A grid plot of median and IQR per gene per class shows the same.

### 7.3.4 Filtering using feature-clustering ensembles

A third option, aimed at building ensembles:

- Ensembles work best when the members differ. One way to create difference is to give each member **different predictors**.
- This works better if the gene pool is **redundant**, with groups of genes behaving alike.
- So the genes (after the ANOVA filter) are **clustered by similarity of their expression across the 94 samples**, using hierarchical clustering of variables, and the tree is cut into **30 clusters**. Cluster sizes range from 11 to 52 genes.
- A helper function builds a gene set by **picking one random gene from each of the 30 clusters**. Each call gives a different set of 30 genes.
- Later, many models are each trained on a different such set and **vote**.

---

## 7.4 Predicting cytogenetic abnormalities

### 7.4.1 The prediction task

- A **classification problem**: forecast a nominal target (four values) from predictors.
- Target: the cytogenetic abnormality type of a B-cell ALL sample. Values: ALL1/AF4, BCR/ABL, E2A/PBX1, NEG.
- Predictors: sets of selected genes, depending on which selection strategy is used.
- Observations: **94 B-cell ALL cases**.

### 7.4.2 The evaluation metric

- Alternatives exist (area under the ROC curve, precision and recall, Brier score), but the author uses **plain accuracy**.
- Reasons: **no information about misclassification costs**, so all errors are treated as equally serious, and multi-class generalizations of ROC analysis are less settled.
- Accuracy is the share of test cases where the predicted class equals the true class.

### 7.4.3 The experimental procedure

- With only 94 cases, the suitable options are **bootstrap** or **leave-one-out cross-validation (LOOCV)**.
- **LOOCV:** train N models, each leaving out one case and testing on it.
- **Bootstrap:** draw a random sample **with replacement** of the same size as the data. The drawn cases train the model; the cases not drawn form the test set. Repeat many times (for example 100 or 200) and average the scores.
- Two bootstrap versions exist (the "e0" and the ".632" bootstrap). A small demonstration uses the iris dataset and 100 repetitions.

### 7.4.4 The modelling techniques

Three models are compared:

| Model | Why it was chosen |
| --- | --- |
| **Random forest** | Handles many features well, because each tree uses random subsets of features |
| **Support vector machine (SVM)** | Can capture non-linear relations between gene expression and abnormality |
| **k-nearest neighbours (k-NN)** | New here. Assumes patients with the same mutation have similar gene "signatures" |

**How k-NN works:**

- It is a **lazy learner**: it builds no model, it just stores the training data.
- To predict a new case, it finds the **k most similar training cases** and takes a vote.
- Similarity is a distance between cases, usually Euclidean distance over the predictors.
- It is **sensitive to the distance used**, to irrelevant variables, and to **variable scale**, so variables should be standardized.
- The number of neighbours k matters. Common values are 1, 3, 5, 7, 11. Too large a k pulls in cases that are far away.
- An odd k is preferred for voting.
- The book wraps the standard function in a formula-style helper that can **standardize** the data first, and can use **user-supplied centre and spread statistics**. In the experiments it uses the **median and IQR** rather than mean and standard deviation.

### 7.4.5 Comparing the models

**A key methodological point.** The ANOVA, random-forest and clustering filters all use the labels. If genes are chosen using **all** the data and only then evaluated, the test cases have influenced which genes were selected, so the accuracy estimate is **biased (too optimistic)**. Therefore these filtering steps are **put inside the workflow being evaluated**, so they run on each training set only. The simple first filter (7.3.1) uses no labels and is done once beforehand.

**Three feature-selection strategies:**

| Code | Meaning |
| --- | --- |
| **s1** | All genes that pass the ANOVA filter |
| **s2** | Random-forest selection of the top 30 genes after ANOVA |
| **s3** | Ensemble of 50 models, each on 30 genes, one random gene from each of 30 clusters, after ANOVA |

**The workflow function**, in order:

1. ANOVA filter (p = 0.01) using the training set only.
2. Apply the chosen strategy (s1, s2 or s3).
3. Train the chosen learner and predict the test set. For s3, predictions come from a **majority vote** among the 50 members.
4. Return true classes, predicted classes and optionally the model.

**Variants compared (parameter grid):**

| Learner | Settings varied | Variants |
| --- | --- | --- |
| Random forest | trees = 500, 750, 1000; features tried per split = 5, 15 | 6 |
| SVM | cost = 1, 100; gamma = 0.01, 0.001, 0.0001 | 6 |
| k-NN | k = 3, 5, 7; standardize = yes or no | 6 |

Each of these 18 is combined with the 3 selection strategies, giving **54 workflows**. Each gets **100 bootstrap repetitions**. The author warns this takes a long time and supplies saved results.

**Results:**

- Best workflow: **svm.v8**, accuracy about **0.813**. It uses an SVM with cost = 100 and gamma = 0.01, with strategy **s2** (30 genes picked by random forest).
- Close behind: svm.v12 (about 0.811), then six k-NN variants at about **0.808**, svm.v10 (about 0.806) and svm.v6 (about 0.805).
- **No random-forest variant appears in the top 10**, which the author calls surprising.
- Of the top 10, **all but one use strategy s2**. The exception svm.v6 uses s1 (all ANOVA-selected genes).
- A paired **Wilcoxon signed-rank test** against the best workflow finds **no statistically significant differences among the top 10** (all p-values are large).
- A look at one bootstrap round for the best workflow, as a confusion matrix:

| True class | Predicted correctly | Where the errors went |
| --- | --- | --- |
| ALL1/AF4 | 3 of 3 | none |
| BCR/ABL | 12 of 15 | 3 predicted as NEG |
| E2A/PBX1 | 0 of 1 | 1 predicted as NEG |
| NEG | 14 of 16 | 1 as BCR/ABL, 1 as E2A/PBX1 |

Most errors are mutated patients labelled NEG (false negatives). The author notes this is not particularly interesting for this application, though the reverse also occurs.

---

## 7.5 Summary (the chapter's own list of topics)

- Feature selection for problems with very many predictors
- Classification methods
- Random forests
- k-nearest neighbours
- SVMs
- Ensembles using different subsets of predictors
- Bootstrap experiments

R skills covered: handling microarray data, ANOVA comparison across groups, clustering variables, ensembles with different predictors, k-NN models, bootstrap accuracy estimates.

---

## Notes and small inconsistencies spotted while reading

These are points to be aware of when reproducing the work. Check them against the book yourself.

1. The parameter list for the random-forest variants is spelled `preditor.pars` in the book code (apparently a typo for `predictor.pars`), so the intended setting probably does not take effect.
2. The text under Figure 7.4 refers to gene `X4045_at`, while the figure's genes include `X40454_at`. Likely a typo.
3. The book's accuracy figures come from R's random number generator and package versions of 2016. Python numbers will differ by a few points. Compare which setups rank well, not exact decimals.
4. The full 54-workflow grid is expensive. The book itself recommends not running it casually.
5. The chapter fixes the number of genes at 30 for s2 and s3 and never studies how accuracy changes with that number.
