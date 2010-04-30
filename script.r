library('affy')

pData = read.table('phenoData_20.tab', row.names=1, header=TRUE, sep="\t")

metadata = data.frame(labelDescription=c("p53 Status", "DLDA Classifier", "DLDA Error", "Elston histologic grade", "ER status", "PgR status", "Age at diagnosis", "Tumor Size (mm)", "Lymph Node Status", "DSS Time", "DSS Event"))

phenoData = new("AnnotatedDataFrame", data=pData, varMetadata=metadata)

file_content = read.table('gsm_A')
A_files = as.character(file_content$V1)
A_batch = ReadAffy(filenames = A_files, phenoData=phenoData)

# file_content = read.table('gsm_B')
# B_files = as.character(file_content$V1)
# B_batch = ReadAffy(filenames = B_files, phenoData=phenoData)

A_batch

# B_batch

library(affyPLM)
png('tex/A_MAplot.png', width=1920, height=1200)
par(mfrow = c(4, 5))
MAplot(A_batch)
dev.off()

# png('tex/B_MAplot.png', width=1920, height=1200)
# par(mfrow = c(4, 5))
# MAplot(B_batch)
# dev.off()

pdf('tex/A_NUSEplot.pdf')
NUSE(Pset)
dev.off()
pdf('tex/A_RLEplot.pdf')
RLE(Pset)
dev.off()

# Pset <- fitPLM(B_batch)
# pdf('tex/B_NUSEplot.pdf')
# NUSE(Pset)
# dev.off()
# pdf('tex/B_RLEplot.pdf')
# RLE(Pset)
# dev.off()

apms <- pm(A_batch)
amms <- mm(A_batch)

# bpms <- pm(B_batch)
# bmms <- mm(B_batch)

pdf('tex/A_SCATTERplot.pdf')
smoothScatter(log2(amms[, 1]), log2(apms[, 1]), xlab=expression(log[2] * "(MM) values"), ylab=expression(log[2] * "(PM) values"), asp=1)
dev.off()

# pdf('tex/B_SCATTERplot.pdf')
# smoothScatter(log2(bmms[, 1]), log2(bpms[, 1]), xlab=expression(log[2] * "(MM) values"), ylab=expression(log[2] * "(PM) values"), asp=1)
# dev.off()

pdf('tex/A_INTENSITY_DISTRIBplot.pdf')
boxplot(A_batch, names=1:length(sampleNames(A_batch)), xlab="Chip A", ylab=expression(log[2] * "(intensity)"), main="Hybridization intensities distributions across 20 subjects on the A chip")
dev.off()

# pdf('tex/B_INTENSITY_DISTRIBplot.pdf')
# boxplot(B_batch, names=1:length(sampleNames(B_batch)), xlab="Chip B", ylab=expression(log[2] * "(intensity)"), main="Hybridization intensities distributions across 20 subjects on the B chip")
# dev.off()

A_batch_RMA = rma(A_batch)
# B_batch_RMA = rma(B_batch)

png('tex/A_batch_RMA.png', width=1920, height=1200)
par(mfrow = c(4, 5))
MAplot(A_batch_RMA)
dev.off()

# png('tex/B_batch_RMA.png', width=1920, height=1200)
# par(mfrow = c(4, 5))
# MAplot(B_batch_RMA)
# dev.off()