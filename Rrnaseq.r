### R code from vignette source 'Rrnaseq.Rnw'

###################################################
### code chunk number 1: Rrnaseq.Rnw:428-429 (eval = FALSE)
###################################################
## download.file("http://biocluster.ucr.edu/~tgirke/HTML_Presentations/Manuals/Workshop_Dec_12_16_2013/Rrnaseq.zip", "Rrnaseq.zip")


###################################################
### code chunk number 2: Rrnaseq.Rnw:431-433
###################################################
targets <- read.delim("./data/targets.txt")
targets


###################################################
### code chunk number 3: Rrnaseq.Rnw:446-453
###################################################
source("http://bioconductor.org/biocLite.R")
biocLite("QuasR")

library(QuasR)
targets <- read.delim("data/targets.txt")
write.table(targets[,1:2], "data/QuasR_samples.txt", row.names=FALSE, quote=FALSE, sep="\t")
sampleFile <- "./data/QuasR_samples.txt"
genomeFile <- "./data/tair10chr.fasta"
results <- "./results" # defines location where to write results
cl <- makeCluster(1) # defines number of CPU cores to use 


###################################################
### code chunk number 4: Rrnaseq.Rnw:457-461
###################################################
proj <- qAlign(sampleFile, genome=genomeFile, maxHits=1, splicedAlignment=FALSE, alignmentsDir=results, 
               clObj=cl, cacheDir=results)
	       # Note: splicedAlignment should be set to TRUE when the reads are >=50nt long  
(alignstats <- alignmentStats(proj)) # Alignment summary report


###################################################
### code chunk number 5: Rrnaseq.Rnw:473-476 (eval = FALSE)
###################################################
## library(Rsubread); library(Rsamtools)
## dir.create("results") # Note: all output data will be written to directory 'results'
## buildindex(basename="./results/tair10chr.fasta", reference="./data/tair10chr.fasta") # Build indexed reference genome


###################################################
### code chunk number 6: Rrnaseq.Rnw:480-489 (eval = FALSE)
###################################################
## targets <- read.delim("./data/targets.txt") # Import experiment design information
## input <- paste("./data/", targets$FileName, sep="")
## output <- paste("./results/", targets$FileName, ".sam", sep="")
## reference <- "./results/tair10chr.fasta"
## for(i in seq(along=targets$FileName)) {
##         align(index=reference, readfile1=input[i], output_file=output[i], nthreads=8, indels=1, TH1=2)
##         asBam(file=output[i], destination=gsub(".sam", "", output[i]), overwrite=TRUE, indexDestination=TRUE)
##         unlink(output[i])
## }


###################################################
### code chunk number 7: Rrnaseq.Rnw:501-504 (eval = FALSE)
###################################################
## library(modules) # Skip this and next line if you are not using IIGB's biocluster
## moduleload("bowtie2/2.1.0"); moduleload("tophat/2.0.8b") # loads bowtie2/tophat2 from module system
## system("bowtie2-build ./data/tair10chr.fasta ./data/tair10chr.fasta")


###################################################
### code chunk number 8: Rrnaseq.Rnw:508-525 (eval = FALSE)
###################################################
## library(Rsamtools)
## dir.create("results") # Note: all output data will be written to directory 'results'
## input <- input <- paste("./data/", targets$FileName, sep="")
## output <- paste("./results/", targets$FileName, sep="")
## reference <- "./data/tair10chr.fasta"
## for(i in seq(along=input)) {
##         unlink(paste(output[i], ".tophat", sep=""), force=TRUE, recursive=TRUE)
##         tophat_command <- paste("tophat -p 4 -g 1 --segment-length 15 -i 30 -I 3000 -o ", output[i], ".tophat ", reference, " ", input[i], sep="")
##                 # -G: supply GFF with transcript model info (preferred!)
## 		# -g: ignore all alginments with >g matches
## 		# -p: number of threads to use for alignment step
## 		# -i/-I: min/max intron lengths
## 		# --segment-length: length of split reads (25 is default)
##         system(tophat_command)
## 	sortBam(file=paste(output[i], ".tophat/accepted_hits.bam", sep=""), destination=paste(output[i], ".tophat/accepted_hits", sep=""))
##         indexBam(paste(output[i], ".tophat/accepted_hits.bam", sep=""))
## }


###################################################
### code chunk number 9: Rrnaseq.Rnw:559-566
###################################################
library(ShortRead); library(Rsamtools)
Nreads <- countLines(dirPath="./data", pattern=".fastq$")/4
bfl <- BamFileList(paste0("./results/", targets$FileName, ".bam"), yieldSize=50000, index=character())
Nalign <- countBam(bfl)
(read_statsDF <- data.frame(FileName=names(Nreads), Nreads=Nreads, Nalign=Nalign$records, 
                            Perc_Aligned=Nalign$records/Nreads*100))
write.table(read_statsDF, "results/read_statsDF.xls", row.names=FALSE, quote=FALSE, sep="\t")


###################################################
### code chunk number 10: Rrnaseq.Rnw:576-581 (eval = FALSE)
###################################################
## qQCReport(proj, pdfFilename="results/qc_report.pdf")
## source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/fastqQuality.R")
## myfiles <- paste0("data/", targets$FileName); names(myfiles) <- targets$SampleName
## fqlist <- seeFastq(fastq=myfiles, batchsize=50000, klength=8)
## pdf("results/fastqReport.pdf", height=18, width=4*length(myfiles)); seeFastqPlot(fqlist); dev.off()


###################################################
### code chunk number 11: Rrnaseq.Rnw:593-601
###################################################
library(rtracklayer); library(GenomicRanges); library(Rsamtools)
gff <- import.gff("./data/TAIR10_GFF3_trunc.gff", asRangedData=FALSE)
seqlengths(gff) <- end(ranges(gff[which(elementMetadata(gff)[,"type"]=="chromosome"),]))
subgene_index <- which(elementMetadata(gff)[,"type"] == "exon")
gffsub <- gff[subgene_index,] # Returns only gene ranges
gffsub[1:4, c(2,5)]
ids <- gsub("Parent=|\\..*", "", elementMetadata(gffsub)$group)
gffsub <- split(gffsub, ids) # Coerce to GRangesList


###################################################
### code chunk number 12: Rrnaseq.Rnw:611-619
###################################################
library(GenomicFeatures)
txdb <- makeTranscriptDbFromGFF(file="data/TAIR10_GFF3_trunc.gff",
	format="gff3",
	dataSource="TAIR",
	species="Arabidopsis thaliana")
saveDb(txdb, file="./data/TAIR10.sqlite")
txdb <- loadDb("./data/TAIR10.sqlite") 
eByg <- exonsBy(txdb, by="gene")


###################################################
### code chunk number 13: Rrnaseq.Rnw:629-642
###################################################
samples <- as.character(targets$FileName)
samplespath <- paste("./results/", samples, ".bam", sep="")
names(samplespath) <- samples
countDF <- data.frame(row.names=names(eByg))
for(i in samplespath) {
        aligns <- readGAlignmentsFromBam(i) # Substitute next two lines with this one.
        counts <- countOverlaps(eByg, aligns, ignore.strand=TRUE)
        countDF <- cbind(countDF, counts)
}
colnames(countDF) <- samples
countDF[1:4,]
write.table(countDF, "./results/countDF", quote=FALSE, sep="\t", col.names = NA)
countDF <- read.table("./results/countDF")


###################################################
### code chunk number 14: Rrnaseq.Rnw:652-658
###################################################
library(GenomicRanges)
##bfl <- BamFileList(samplespath, yieldSize=50000, index=character())
##countDF2 <- summarizeOverlaps(eByg, bfl, mode="Union", ignore.strand=TRUE)
##countDF2 <- assays(countDF2)$counts
##colnames(countDF2) <- samples
##countDF2[1:4,]


###################################################
### code chunk number 15: Rrnaseq.Rnw:668-671
###################################################
countDF3 <- qCount(proj, txdb, reportLevel="gene", orientation="any") 
countDF3[1:4,]
write.table(countDF3, "results/countDFgene.xls", col.names=NA, quote=FALSE, sep="\t")


###################################################
### code chunk number 16: Rrnaseq.Rnw:681-690
###################################################
returnRPKM <- function(counts, gffsub) {
        geneLengthsInKB <- sum(width(reduce(gffsub)))/1000 # Length of exon union per gene in kbp
        millionsMapped <- sum(counts)/1e+06 # Factor for converting to million of mapped reads.
        rpm <- counts/millionsMapped # RPK: reads per kilobase of exon model.
        rpkm <- rpm/geneLengthsInKB # RPKM: reads per kilobase of exon model per million mapped reads.
        return(rpkm)
}
countDFrpkm <- apply(countDF, 2, function(x) returnRPKM(counts=x, gffsub=eByg))
countDFrpkm[1:4,]


###################################################
### code chunk number 17: Rrnaseq.Rnw:694-695
###################################################
rpkmDFgene <- t(t(countDF3[,-1]/countDF3[,1] * 1000)/colSums(countDF3[,-1]) *1e6)


###################################################
### code chunk number 18: sampletree
###################################################
library(ape) 
d <- cor(countDFrpkm, method="spearman")
hc <- hclust(dist(1-d))
plot.phylo(as.phylo(hc), type="p", edge.col=4, edge.width=3, show.node.label=TRUE, no.margin=TRUE)


###################################################
### code chunk number 19: Exercise 1 (eval = FALSE)
###################################################
## ## 1. Read alignment
## library(QuasR)
## sampleFile <- "./data/QuasR_samples.txt"; genomeFile <- "./data/tair10chr.fasta"; results <- "./results"; cl <- makeCluster(1)
## proj <- qAlign(sampleFile, genome=genomeFile, maxHits=1, splicedAlignment=FALSE, alignmentsDir=results, clObj=cl, cacheDir=results)
## (alignstats <- alignmentStats(proj)) # Alignment summary report
## 
## ## 2. Sense and antisense read counting
## library(GenomicFeatures)
## txdb <- loadDb("./data/TAIR10.sqlite") 
## eByg <- exonsBy(txdb, by="gene")
## senseDF <- qCount(proj, txdb, reportLevel="gene", orientation="same") 
## senseDF[1:4,]
## antisDF <- qCount(proj, txdb, reportLevel="gene", orientation="opposite") 
## antisDF[1:4,]
## 
## ## 3. Identify genes with much more antisense read counts
## index <- which(rowSums(antisDF[,-1]/senseDF[,-1] >=3) >= 2)
## senseDF[index,]
## antisDF[index,]
## 
## ## 4. Plot result with ggbio
## library(ggbio); library(Rsamtools)
## genes(txdb)[names(index)] # Returns positional information for candidates
## TRLa <- readGAlignmentsFromBam("./results/SRR064154.fastq.bam", use.names=TRUE, param=ScanBamParam(which=GRanges("Chr5", IRanges(44500, 47300)))) 
## TRLb <- readGAlignmentsFromBam("./results/SRR064155.fastq.bam", use.names=TRUE, param=ScanBamParam(which=GRanges("Chr5", IRanges(44500, 47300))))
## p1 <- autoplot(TRLa, geom = "rect", aes(color = strand, fill = strand))
## p2 <- autoplot(TRLb, geom = "rect", aes(color = strand, fill = strand))
## p3 <- autoplot(txdb, which=GRanges("Chr5", IRanges(44500, 47300)), names.expr = "gene_id")
## tracks(TRLa=p1, TRLb=p2, Transcripts=p3, heights = c(0.3, 0.3, 0.4)) + ylab("")


###################################################
### code chunk number 20: Rrnaseq.Rnw:767-770
###################################################
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/colAg.R")
countDFrpkm_mean <- colAg(myMA=countDFrpkm, group=c(1,1,2,2), myfct=mean)
countDFrpkm_mean[1:4,]


###################################################
### code chunk number 21: Rrnaseq.Rnw:773-779
###################################################
countDFrpkm_mean <- cbind(countDFrpkm_mean, log2ratio=log2(countDFrpkm_mean[,2]/countDFrpkm_mean[,1]))
countDFrpkm_mean <- countDFrpkm_mean[is.finite(countDFrpkm_mean[,3]), ]
degs2fold <- countDFrpkm_mean[countDFrpkm_mean[,3] >= 1 | countDFrpkm_mean[,3] <= -1,]
degs2fold[1:4,]
write.table(degs2fold, "./results/degs2fold.xls", quote=FALSE, sep="\t", col.names = NA)
degs2fold <- read.table("./results/degs2fold.xls")


###################################################
### code chunk number 22: Rrnaseq.Rnw:789-801
###################################################
library(DESeq)
countDF <- read.table("./results/countDF")
conds <- targets$Factor
cds <- newCountDataSet(countDF, conds) # Creates object of class CountDataSet derived from eSet class
counts(cds)[1:4, ] # CountDataSet has similar accessor methods as eSet class.
cds <- estimateSizeFactors(cds) # Estimates library size factors from count data. Alternatively, one can provide here the true library sizes with sizeFactors(cds) <- c(..., ...)
cds <- estimateDispersions(cds) # Estimates the variance within replicates
res <- nbinomTest(cds, "AP3", "TRL") # Calls DEGs with nbinomTest
res <- na.omit(res)
res2fold <- res[res$log2FoldChange >= 1 | res$log2FoldChange <= -1,]
res2foldpadj <- res2fold[res2fold$padj <= 0.05, ]
res2foldpadj[1:4,1:8]


###################################################
### code chunk number 23: Rrnaseq.Rnw:811-821
###################################################
library(edgeR)
countDF <- read.table("./results/countDF")
y <- DGEList(counts=countDF, group=conds) # Constructs DGEList object
y <- estimateCommonDisp(y) # Estimates common dispersion
y <- estimateTagwiseDisp(y) # Estimates tagwise dispersion
et <- exactTest(y, pair=c("AP3", "TRL")) # Computes exact test for the negative binomial distribution.
topTags(et, n=4)
edge <- as.data.frame(topTags(et, n=50000)) 
edge2fold <- edge[edge$logFC >= 1 | edge$logFC <= -1,]
edge2foldpadj <- edge2fold[edge2fold$FDR <= 0.01, ]


###################################################
### code chunk number 24: Rrnaseq.Rnw:831-850
###################################################
library(edgeR)
countDF <- read.table("./results/countDF")
y <- DGEList(counts=countDF, group=conds) # Constructs DGEList object
## Filtering and normalization
keep <- rowSums(cpm(y)>1) >= 2; y <- y[keep, ]
y <- calcNormFactors(y)
design <- model.matrix(~0+group, data=y$samples); colnames(design) <- levels(y$samples$group) # Design matrix
## Estimate dispersion
y <- estimateGLMCommonDisp(y, design, verbose=TRUE) # Estimates common dispersions
y <- estimateGLMTrendedDisp(y, design) # Estimates trended dispersions
y <- estimateGLMTagwiseDisp(y, design) # Estimates tagwise dispersions 
## Fit the negative binomial GLM for each tag
fit <- glmFit(y, design) # Returns an object of class DGEGLM
contrasts <- makeContrasts(contrasts="AP3-TRL", levels=design) # Contrast matrix is optional
lrt <- glmLRT(fit, contrast=contrasts[,1]) # Takes DGEGLM object and carries out the likelihood ratio test. 
edgeglm <- as.data.frame(topTags(lrt, n=length(rownames(y))))
## Filter on fold change and FDR
edgeglm2fold <- edgeglm[edgeglm$logFC >= 1 | edgeglm$logFC <= -1,]
edgeglm2foldpadj <- edgeglm2fold[edgeglm2fold$FDR <= 0.01, ]


###################################################
### code chunk number 25: Rrnaseq.Rnw:861-866
###################################################
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/overLapper.R")
setlist <- list(edgeRexact=rownames(edge2foldpadj), edgeRglm=rownames(edgeglm2foldpadj), DESeq=as.character(res2foldpadj[,1]), RPKM=rownames(degs2fold))
OLlist <- overLapper(setlist=setlist, sep="_", type="vennsets")
counts <- sapply(OLlist$Venn_List, length)
vennPlot(counts=counts, mymain="DEG Comparison")


###################################################
### code chunk number 26: Rrnaseq.Rnw:880-886
###################################################
library(lattice); library(gplots)
y <- countDFrpkm[rownames(edgeglm2foldpadj)[1:20],]
colnames(y) <- targets$Factor
y <- t(scale(t(as.matrix(y))))
y <- y[order(y[,1]),]
levelplot(t(y), height=0.2, col.regions=colorpanel(40, "darkblue", "yellow", "white"), main="Expression Values (DEG Filter: FDR 1%, FC > 2)", colorkey=list(space="top"), xlab="", ylab="Gene ID")


###################################################
### code chunk number 27: Rrnaseq.Rnw:900-909
###################################################
library(GOstats); library(GO.db); library(ath1121501.db)
geneUniverse <- rownames(countDF) 
geneSample <- res2foldpadj[,1]
params <- new("GOHyperGParams", geneIds = geneSample, universeGeneIds = geneUniverse, 
        annotation="ath1121501", ontology = "MF", pvalueCutoff = 0.5,
	conditional = FALSE, testDirection = "over")
hgOver <- hyperGTest(params)
summary(hgOver)[1:4,]
htmlReport(hgOver, file = "results/MyhyperGresult.html")


###################################################
### code chunk number 28: Rrnaseq.Rnw:950-957
###################################################
library(ggbio)
AP3 <- readGAlignmentsFromBam("./results/SRR064154.fastq.bam", use.names=TRUE, param=ScanBamParam(which=GRanges("Chr1", IRanges(49457, 51457)))) 
TRL <- readGAlignmentsFromBam("./results/SRR064166.fastq.bam", use.names=TRUE, param=ScanBamParam(which=GRanges("Chr1", IRanges(49457, 51457))))
p1 <- autoplot(AP3, geom = "rect", aes(color = strand, fill = strand))
p2 <- autoplot(TRL, geom = "rect", aes(color = strand, fill = strand))
p3 <- autoplot(txdb, which=GRanges("Chr1", IRanges(49457, 51457)), names.expr = "gene_id")
tracks(AP3=p1, TRL=p2, Transcripts=p3, heights = c(0.3, 0.3, 0.4)) + ylab("")


###################################################
### code chunk number 29: Exercise 2 (eval = FALSE)
###################################################
## ## 1. List of upregulated DEGs
## setlistup <- list(edgeRexact=rownames(edge2foldpadj[edge2foldpadj$logFC >= 1, ]), 
## 		  edgeRglm=rownames(edgeglm2foldpadj[edgeglm2foldpadj$logFC >= 1, ]), 
## 		  DESeq=res2foldpadj[res2foldpadj$log2FoldChange <= -1, 1], 
## 	          RPKM=rownames(degs2fold[degs2fold$log2ratio <= -1,]))
## ## 2. List of downregulated DEGs
## setlistdown <- list(edgeRexact=rownames(edge2foldpadj[edge2foldpadj$logFC <= -1, ]), 
## 		  edgeRglm=rownames(edgeglm2foldpadj[edgeglm2foldpadj$logFC <= -1, ]), 
## 		  DESeq=res2foldpadj[res2foldpadj$log2FoldChange >= 1, 1], 
## 	          RPKM=rownames(degs2fold[degs2fold$log2ratio >= 1,]))
## 
## ## 3. Plot up and downregulated DEGs in one venn diagram 
## source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/overLapper.R") 
## OLlistup <- overLapper(setlist=setlistup, sep="_", type="vennsets")
## OLlistdown <- overLapper(setlist=setlistdown, sep="_", type="vennsets")
## counts <- list(sapply(OLlistup$Venn_List, length), sapply(OLlistdown$Venn_List, length))
## vennPlot(counts=counts, ccol=c("red", "blue"), colmode=2, mysub="Top: DEG UP; Bottom: DEG Down", yoffset=c(0.3, -0.2)) 


###################################################
### code chunk number 30: Rrnaseq.Rnw:1001-1014
###################################################
source("data/Fct/gffexonDEXSeq.R")
gffexonDEXSeq <- exons2DEXSeq(gff=gff)
ids <- as.character(elementMetadata(gffexonDEXSeq)[, "ids"])
countDFdex <- data.frame(row.names=ids)
for(i in samplespath) {
        aligns <- readGAlignmentsFromBam(i) # Substitute next two lines with this one.
        counts <- countOverlaps(gffexonDEXSeq, aligns)
        countDFdex <- cbind(countDFdex, counts)
}
colnames(countDFdex) <- samples
countDFdex[1:4,1:2]
write.table(countDFdex, "./results/countDFdex", quote=FALSE, sep="\t", col.names = NA)
countDFdex <- read.table("./results/countDFdex")


###################################################
### code chunk number 31: Rrnaseq.Rnw:1023-1042
###################################################
## library(DEXSeq)
## samples <- as.character(targets$Factor); names(samples) <- targets$FileName 
## countDFdex[is.na(countDFdex)] <- 0
## Construct ExonCountSet from scratch
## Commented by: Shiva Prasad Gaddameedi
## newExonCountSet2 is depricated and getting failed to use DEXSeqDataSet in place of newExonCountSet2
## exset <- newExonCountSet2(countDF=countDFdex) # fData(exset)[1:4,]
## exset <- DEXSeqDataSet(countDFdex, samples, design,featureID,ids)
## Performs normalization
## exset <- estimateSizeFactors(exset) 
## Evaluate variance of the data by estimating dispersion using Cox-Reid (CR) likelihood estimation
## exset <- estimateDispersions(exset) 
## Fits dispersion-mean relation to the individual CR dispersion values
## exset <- fitDispersionFunction(exset) 
## Performs Chi-squared test on each exon and Benjmini-Hochberg p-value adjustment for mutliple testing
## exset <- testForDEU(exset) 
## Estimates fold changes of exons
## exset <- estimatelog2FoldChanges(exset) 
## Obtain results in data frame
## deuDF <- DEUresultTable(exset) 
## Count number of genes with differential exon usage
## table(tapply(deuDF$padjust < 0.01, geneIDs(exset), any)) 




###################################################
### code chunk number 32: dexseq1
###################################################
## plotDEXSeq(exset, "Parent=AT1G01100", displayTranscripts=TRUE, expression=TRUE, legend=TRUE)  
## Generate many plots and write them to results directory
## mygeneIDs <- unique(as.character(na.omit(deuDF[deuDF$geneID %in% unique(deuDF$geneID),])[,"geneID"]))
## DEXSeqHTML(exset, geneIDs=mygeneIDs, path="results", file="DEU.html") 


###################################################
### code chunk number 33: Rrnaseq.Rnw:1066-1067
###################################################
sessionInfo()
