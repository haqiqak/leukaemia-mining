# Export the B-cell ALL subset (94 patients) to CSV, patients in rows and genes in columns.
# Run from the project root: Rscript src/export_all.R
library(ALL)
data(ALL)

# Same subsetting as the book: B-cell stages and the first four molecular-biology levels
tgt.cases <- which(ALL$BT %in% levels(ALL$BT)[1:5] &
                   ALL$mol.biol %in% levels(ALL$mol.biol)[1:4])
ALLb <- ALL[, tgt.cases]
ALLb$mol.biol <- factor(ALLb$mol.biol)

expr <- data.frame(patient = sampleNames(ALLb), t(exprs(ALLb)), check.names = FALSE)
write.csv(expr, "data/expression.csv", row.names = FALSE)
labels <- data.frame(patient = sampleNames(ALLb), class = ALLb$mol.biol, BT = ALLb$BT)
write.csv(labels, "data/labels.csv", row.names = FALSE)

cat("Package versions: ALL", as.character(packageVersion("ALL")),
    "| Biobase", as.character(packageVersion("Biobase")),
    "|", R.version.string, "\n")
