library('affy')

pData = read.table('phenodata_40.tab', row.names=1, header=TRUE, sep="\t")

metadata = data.frame(labelDescription=c("p53 Status", "DLDA Classifier", "DLDA Error", "Elston histologic grade", "ER status", "PgR status", "Age at diagnosis", "Tumor Size (mm)", "Lymph Node Status", "DSS Time", "DSS Event"))

phenoData = new("AnnotatedDataFrame", data=pData, varMetadata=metadata)

file_content = read.table('newAA')
A_files = as.character(file_content$V1)
A_batch = ReadAffy(filenames = A_files, phenoData=phenoData)

A_batch

A_batch_RMA = rma(A_batch)


library(limma)
library(genefilter)
library(hgu133a.db)

probenames = probeNames(A_batch)

IQRs = esApply(A_batch_RMA, 1, IQR)
#plot.ecdf(IQRs, pch=".")
#points(IQRs[probenames],rep(0,length(probenames)),pch=1, col="red")
#abline(v=quantile(IQRs, 0.6), lwd=2, col="orange")

# Non-specific filtering
filtered_A_batch_RMA <- nsFilter(A_batch_RMA, require.entrez=FALSE, remove.dupEntrez=FALSE, var.func=IQR, var.cutoff=0.6, var.filter=TRUE, filterByQuantile=TRUE, feature.exclude="^AFFX")
feset = filtered_A_batch_RMA$eset
# Functional annotation based on the GO annotation
library(org.Hs.eg.db)
library(GOstats)

design = model.matrix(~ 0 + feset$p53.seq.mut.status..p53..mutant..p53..wt., sep="")

colnames(design)= c("WT", "M")
cont.matrix = makeContrasts(WTvsM = (M - WT), levels = design)
cont.matrix

fit <- lmFit(feset, design)

fit2  <- contrasts.fit(fit, cont.matrix)
fit2  <- eBayes(fit2)

library(annotate)
library(hgu133a.db)
fit2$genes$Symbol <- getSYMBOL(fit2$genes$ID, "hgu133a")

topTable(fit2, coef=1, adjust.method="BH", sort.by="P")

results <- decideTests(fit2, p.value=0.1)
results
summary(results)

affyUniverse = featureNames(feset)
entrezGeneUniverse <- unlist(mget(affyUniverse, hgu133aENTREZID))
DE_WTvsM = names(results[results[, "WTvsM"]!=0,])
selectedEntrezGeneIDs = unlist(mget(DE_WTvsM, hgu133aENTREZID))
params <- new("GOHyperGParams", geneIds=selectedEntrezGeneIDs,
                 universeGeneIds=entrezGeneUniverse, annotation="hgu133a.db",
                 ontology="CC", pvalueCutoff=0.01, conditional=FALSE,
                 testDirection="over")

paramsCond <- params
conditional(paramsCond) <- TRUE
hgOverA <- hyperGTest(paramsCond)
summary(hgOverA)
