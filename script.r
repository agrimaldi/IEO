library('affy')

pData = read.table('phenodata_20.tab', row.names=1, header=TRUE, sep="\t")

metadata = data.frame(labelDescription=c("p53 Status", "DLDA Classifier", "DLDA Error", "Elston histologic grade", "ER status", "PgR status", "Age at diagnosis", "Tumor Size (mm)", "Lymph Node Status", "DSS Time", "DSS Event"))

phenoData = new("AnnotatedDataFrame", data=pData, varMetadata=metadata)

file_content = read.table('gsm_A')
A_files = as.character(file_content$V1)
A_batch = ReadAffy(filenames = A_files, phenoData=phenoData)

A_batch

#library(affyPLM)
#png('tex/A_MAplot.png', width=1920, height=1200)
#par(mfrow = c(4, 5))
#MAplot(A_batch)
#dev.off()

#Pset = fitPLM(A_batch)
#pdf('tex/A_NUSEplot.pdf')
#NUSE(Pset)
#dev.off()
#pdf('tex/A_RLEplot.pdf')
#RLE(Pset)
#dev.off()

#apms <- pm(A_batch)
#amms <- mm(A_batch)

#pdf('tex/A_SCATTERplot.pdf')
#smoothScatter(log2(amms[, 1]), log2(apms[, 1]), xlab=expression(log[2] * "(MM) values"), ylab=expression(log[2] * "(PM) values"), asp=1)
#dev.off()

#pdf('tex/A_INTENSITY_DISTRIBplot.pdf')
#boxplot(A_batch, names=1:length(sampleNames(A_batch)), xlab="Chip A", ylab=expression(log[2] * "(intensity)"), main="Hybridization intensities distributions across 20 subjects on the A chip")
#dev.off()

A_batch_RMA = rma(A_batch)

#png('tex/A_batch_RMA.png', width=1920, height=1200)
#par(mfrow = c(4, 5))
#MAplot(A_batch_RMA)
#dev.off()

library(limma)
library(genefilter)
library(hgu133a.db)

probenames = probeNames(A_batch)

IQRs <- esApply(A_batch_RMA, 1, IQR)
plot.ecdf(IQRs, pch=".")
points(IQRs[probenames],rep(0,length(probenames)),pch=1, col="red")
abline(v=quantile(IQRs, 0.2), lwd=2, col="orange")

filtered_A_batch_RMA <- nsFilter(A_batch_RMA, require.entrez=FALSE, remove.dupEntrez=FALSE, var.func=IQR, var.cutoff=0.3, var.filter=TRUE, filterByQuantile=TRUE, feature.exclue="^AFFX")

feset = filtered_A_batch_RMA$eset
combinedLevels = factor(paste(feset$p53.seq.mut.status..p53..mutant..p53..wt., feset$p53.DLDA.classifier.result..0.wt.like..1.mt.like., sep=""))
design <- model.matrix(~ 0 + combinedLevels)
colnames(design) = c("p53m0", "p53m1", "p53p0", "p53p1")

cont.matrix <- makeContrasts(
MT.WTlike=(p53p0 - p53m0),
MT.MTlike=(p53p1 - p53m1),
WT.WTlike=(p53m1 - p53m0),
WT.MTlike=(p53p1 - p53m0),
levels=design)

fit <- lmFit(feset, design)

fit2  <- contrasts.fit(fit, cont.matrix)
fit2  <- eBayes(fit2)

library(annotate)
library(hgu133a.db)
fit2$genes$Symbol <- getSYMBOL(fit2$genes$ID, "hgu133a")

MT_WTlike = topTable(fit2, coef=1, adjust.method="BH", sort.by="P")
MT_MTlike = topTable(fit2, coef=2, adjust.method="BH", sort.by="P")
WT_WTlike = topTable(fit2, coef=3, adjust.method="BH", sort.by="P")
WT_MTlike = topTable(fit2, coef=4, adjust.method="BH", sort.by="P")

results <- decideTests(fit2)
results
summary(results)