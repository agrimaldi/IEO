library('affy')

pData = read.table('phenodata_40.tab', row.names=1, header=TRUE, sep="\t")

metadata = data.frame(labelDescription=c("p53 Status", "DLDA Classifier", "DLDA Error", "Elston histologic grade", "ER status", "PgR status", "Age at diagnosis", "Tumor Size (mm)", "Lymph Node Status", "DSS Time", "DSS Event"))

phenoData = new("AnnotatedDataFrame", data=pData, varMetadata=metadata)

file_content = read.table('newAA')
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

IQRs = esApply(A_batch_RMA, 1, IQR)
#plot.ecdf(IQRs, pch=".")
#points(IQRs[probenames],rep(0,length(probenames)),pch=1, col="red")
#abline(v=quantile(IQRs, 0.6), lwd=2, col="orange")

# Non-specific filtering
filtered_A_batch_RMA <- nsFilter(A_batch_RMA, require.entrez=FALSE, remove.dupEntrez=FALSE, var.func=IQR, var.cutoff=0.6, var.filter=TRUE, filterByQuantile=TRUE, feature.exclude="^AFFX")

# Factorial design
feset = filtered_A_batch_RMA$eset
combinedLevels = factor(paste(feset$p53.seq.mut.status..p53..mutant..p53..wt., feset$p53.DLDA.classifier.result..0.wt.like..1.mt.like., sep=""))
design <- model.matrix(~ 0 + combinedLevels)
colnames(design) = c("p53m0", "p53m1", "p53p0", "p53p1")
# 
# # contrasts : look at differencially expressed genes due to the visible phenotype "p53 mutant-like" in p53-wt tumors, and DE genes due to p53 status in both visible phenotypes mutant-like and wt-like
 cont.matrix <- makeContrasts(
 #visible_phenotype_in_p53WT=(p53m1 - p53m0),
 #visible_phenotype_in_p53M=(p53p1 - p53p0),
p53_status_WTlike=(p53p0 - p53m0),
p53_status_Mlike=(p53p1 - p53m1),
 levels=design)
 
 fit <- lmFit(feset, design)
 
 fit2  <- contrasts.fit(fit, cont.matrix)
 fit2  <- eBayes(fit2)
 
 library(annotate)
 library(hgu133a.db)
 fit2$genes$Symbol <- getSYMBOL(fit2$genes$ID, "hgu133a")
 
 # ranking of probesets
 topTable(fit2, coef=1, adjust.method="BH", sort.by="P")
 topTable(fit2, coef=2, adjust.method="BH", sort.by="P")
 topTable(fit2, coef=3, adjust.method="BH", sort.by="P")
 
 results <- decideTests(fit2, p.value=0.1)
 results
 summary(results)
 
 
 # Functional annotation based on the GO annotation
 library(org.Hs.eg.db)
 library(GOstats)
 
 affyUniverse = featureNames(feset)
 entrezGeneUniverse <- unlist(mget(affyUniverse, hgu133aENTREZID))
 
 aa = names(results[results[, "visible_phenotype_in_p53WT"]!=0, ])
 bb = names(results[results[, "visible_phenotype_in_p53M"]!=0, ])
 # p53_status_Mlike_DE = rownames(results[results[, "p53_status_Mlike"]!=0, ])
 p53_status_DE = intersect(aa, bb)



 selectedEntrezGeneIDs = unlist(mget(p53_status_DE, hgu133aENTREZID))
 params <- new("GOHyperGParams", geneIds=selectedEntrezGeneIDs,
                 universeGeneIds=entrezGeneUniverse, annotation="hgu133a.db",
                 ontology="BP", pvalueCutoff=0.05, conditional=FALSE,
                 testDirection="over")
 paramsCond <- params
 conditional(paramsCond) <- TRUE
 hgOverA <- hyperGTest(paramsCond)
# 
# selectedEntrezGeneIDs = unlist(mget(Mlike_p53WT_DE, hgu133aENTREZID))
# params <- new("GOHyperGParams", geneIds=selectedEntrezGeneIDs,
#                 universeGeneIds=entrezGeneUniverse, annotation="hgu133a.db",
#                 ontology="BP", pvalueCutoff=0.05, conditional=FALSE,
#                 testDirection="over")
# paramsCond <- params
# conditional(paramsCond) <- TRUE
# hgOverB <- hyperGTest(paramsCond)
# 
# summary(hgOverA)
# summary(hgOverB)

#############################################################
#Look at DE genes due only to the p53 status of the tumors
 
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
