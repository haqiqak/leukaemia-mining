# Export the probes kept by the book's label-free first filter (3,943 of 12,625).
# Run from the project root after src/export_all.R: Rscript src/export_filter.R
library(ALL)
library(genefilter)
data(ALL)

tgt.cases <- which(ALL$BT %in% levels(ALL$BT)[1:5] &
                   ALL$mol.biol %in% levels(ALL$mol.biol)[1:4])
ALLb <- ALL[, tgt.cases]
es <- exprs(ALLb)

# Drops control probes, probes without a gene id, duplicate probes and low-IQR probes
res <- nsFilter(ALLb, var.func = IQR, var.cutoff = IQR(as.vector(es)) / 5,
                feature.exclude = "^AFFX")
print(res$filter.log)

write.csv(data.frame(probe = featureNames(res$eset)), "data/genes_kept.csv", row.names = FALSE)
cat("Kept", nrow(res$eset), "probes. genefilter", as.character(packageVersion("genefilter")),
    "| hgu95av2.db", as.character(packageVersion("hgu95av2.db")), "\n")
