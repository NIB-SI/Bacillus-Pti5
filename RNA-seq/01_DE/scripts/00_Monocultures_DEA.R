# Monocultures experiment differential expresseion analysis


##############################
############# ROOTS

# combine all counts in to an expression matrix
files <- list.files(path="../output/mapping/Monocultures/", pattern=".ReadsPerGene.out.tab",  full.names= TRUE)

geneNames <- read.table(files[1], header = FALSE, sep = "\t", skip= 4)[,1]

counts1 <- NULL
colNs <- NULL
counti <- NULL

for(i in 1:length(files)) {
  counti <- read.table(files[i], header = FALSE, sep = "\t", skip= 4)[,4] #  row 4 are reverse stranded reads
  colNi <- unlist(strsplit(unlist(strsplit(unlist(strsplit(files[i], ".", fixed = TRUE))[3], "/", fixed = TRUE))[4],"_"))[2] # change if needed to capture sample name
  
  counts1 <- cbind(counts1, counti)
  colNs <- cbind(colNs , colNi) 
}

colNs

colnames(counts1) <- as.vector(colNs)
rownames(counts1) <- geneNames

head(counts1) # check if it makes sense!

# based on centrifuge results as well as first pass of DE analysis, sample BS472 has to be discarded
(colnames(counts1)=="BS472")[22]

counts <- counts1[,-22]

#export table with all counts
write.table(counts, file="../output/Monocultures/StGenome_mergedGFF_STAR_counts_Stranded.txt", sep="\t", row.names=TRUE)


library ("limma")
library("edgeR")
library("stringr")

colnames(counts)

# read phenodata from analytes.txt/phenodata.txt (samples used for RNA-Seq)
phenodata <- read.table("../input/phenodata.txt",  row.names=1, header = TRUE)
dim(phenodata)
head(phenodata)

## selecting columns for ROOT samples!!
phenodata <- phenodata[1:23, ]
counts <- counts[, 1:23]

#check if samples are in same order; if TRUE it's OK :)
all(rownames(phenodata)==colnames(counts))

levels(phenodata$treatment)
# assign groups to samples
group <- as.factor(as.vector(phenodata$treatment))

x <- counts

## Create a DGEList object for limma statistical analysis and specify library size i.e. number of reads sequenced per sample
# default library sizes (sum of mapped reads)!!
y <- DGEList(counts=x, group=group)


#####
#define color palette for the samples in graphs
col= c(rep("grey",4), rep("sienna1",4), rep("red1",4), rep("grey20",4), rep("sienna4",4), rep("red4",3))

#check density plot
pdf("../other/root_density_plot_before_lowexpr_filter.pdf", onefile=TRUE, family= "Helvetica")
nsamples <- ncol(x)

lcpm <- log(as.matrix(x),10)
plot(density(lcpm), col=col[1], lwd=2, ylim=c(0,0.4), las=2, main="", xlab="")
title(main="A. BEFORE REMOVAL", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
  
}
legend("topright", colnames(lcpm), text.col=col, cex = 0.75, bty="n")
dev.off()

## Filter genes based on expression -- keep genes (rows) that have over 100 counts in at least 4 samples == one group  ### had to filter more stringently here!
# filtering is necessary for voom method to work properly
keep.exprs <- rowSums(y$counts>100)>=4


## Normalization of dataset for different library sizes
y1 <- y[keep.exprs, , keep.lib.sizes=TRUE]
y1 <- calcNormFactors(y1)

## Plot QC plots using different functions e.g.:

pdf("../other/root_log10rawcounts_boxplot.pdf", onefile=TRUE, family= "Helvetica")
boxplot(log(y$counts+1,10), las=2, ylab="log10(counts)", col=col)
dev.off()

pdf("../other/root_log10filteredcounts_boxplot.pdf", onefile=TRUE, family= "Helvetica")
boxplot(log(y1$counts+1,10), las=2, ylab="log10(counts)", col=col)
dev.off()

#density plots before and after removing low expressed genes
pdf("../other/root_norm_counts_raw&filtered_densityplots.pdf", onefile=TRUE, family= "Helvetica")
opar <- par()
par(mfrow=c(1,2), cex = 0.6)
nsamples <- ncol(x)

lcpm <- log(as.matrix(x),10)
plot(density(lcpm), col=col[1], lwd=2, ylim=c(0,0.4), las=2, main="", xlab="")
title(main="A. BEFORE REMOVAL", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
#legend("topright", colnames(lcpm), text.col=col, bty="n")

lcpm <- log(as.matrix(y1),10)
plot(density(lcpm), col=col[1], lwd=2, ylim=c(0,1), las=2, main="", xlab="")
title(main="B. AFTER REMOVAL", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", colnames(lcpm), text.col=col, bty="n")
par(opar)
dev.off()


#MDS (PCA-like) graph
pdf("../other/root_norm_counts_filtered_MDS.pdf", onefile=TRUE, family= "Helvetica")
plotMDS(y1, labels=colnames(y1), col = col, cex = 0.6)
dev.off()

pdf("../other/root_norm_cf_lcpm_MDS.pdf", onefile=TRUE, family= "Helvetica")
lcpm <- log(as.matrix(y1),10)
plotMDS(lcpm, labels=colnames(lcpm), col = col, cex = 0.6)
dev.off()

###########
## limma-voom protocol
# Create design matrix
design <- model.matrix(~0+group)
colnames(design) <- levels(group)

# limma voom fit for filtered RNA-seq dataset (y1)
pdf("../other/root_voom_mean-variance_trend.pdf", onefile=TRUE, family= "Helvetica")
v <- voom(y1,design,plot=TRUE)
dev.off()
fit <- lmFit(v, design)

## Define contrasts i.e. comparisons between groups:
### !!! ###
contrastMatrix = makeContrasts("roots_inoculated_PS216_2h_A-roots_noninoculated_A",
                               "roots_inoculated_PS218_2h_A-roots_noninoculated_A",
                               "roots_inoculated_PS216_26h_B-roots_noninoculated_B",
                               "roots_inoculated_PS218_26h_B-roots_noninoculated_B",
                               "(roots_inoculated_PS216_2h_A+roots_inoculated_PS216_26h_B)-(roots_noninoculated_A+roots_noninoculated_B)",
                               "(roots_inoculated_PS218_2h_A+roots_inoculated_PS218_26h_B)-(roots_noninoculated_A+roots_noninoculated_B)",
                               levels=design)

fit2 = contrasts.fit(fit, contrastMatrix)
## Check DEG in contrasts (adj p-val cutoff 0.05, |logFC| > 0) using alternative to eBayes (just as control)
tfit <- treat(fit2)
tfit
logFCcut <- 1
dt <- decideTests(tfit, lfc=logFCcut, adjust.method = "BH") #defaults adjust.method = "BH", p.value = 0.05, lfc=0
summary(dt)
colnames(dt)

## eBayes statistics calculation
fit2 <- eBayes(fit2)
pdf("../other/root_SIGMA_vs_A_plot.pdf", onefile=TRUE, family= "Helvetica")
plotSA(fit2)
dev.off()

## make results table
results1 <- topTable(fit2, coef=1, number=1000000, sort.by="none")
results2 <- topTable(fit2, coef=2, number=1000000, sort.by="none")
results3 <- topTable(fit2, coef=3, number=1000000, sort.by="none")
results4 <- topTable(fit2, coef=4, number=1000000, sort.by="none")
results5 <- topTable(fit2, coef=5, number=1000000, sort.by="none")
results6 <- topTable(fit2, coef=6, number=1000000, sort.by="none")

#make one expression matrix with logFCs and adj.P.vals
results <- cbind(results1[,1], results1[,5],
                 results2[,1], results2[,5],
                 results3[,1], results3[,5],
                 results4[,1], results4[,5],
                 results5[,1], results5[,5],
                 results6[,1], results6[,5])

colnames(results) <- c(paste(colnames(contrastMatrix)[1], " ", colnames(results1[1])),
                       paste(colnames(contrastMatrix)[1], " ", colnames(results1[5])),
                       paste(colnames(contrastMatrix)[2], " ", colnames(results2[1])),
                       paste(colnames(contrastMatrix)[2], " ", colnames(results2[5])),
                       paste(colnames(contrastMatrix)[3], " ", colnames(results3[1])),
                       paste(colnames(contrastMatrix)[3], " ", colnames(results3[5])),
                       paste(colnames(contrastMatrix)[4], " ", colnames(results4[1])),
                       paste(colnames(contrastMatrix)[4], " ", colnames(results4[5])),
                       paste(colnames(contrastMatrix)[5], " ", colnames(results5[1])),
                       paste(colnames(contrastMatrix)[5], " ", colnames(results5[5])),
                       paste(colnames(contrastMatrix)[6], " ", colnames(results6[1])),
                       paste(colnames(contrastMatrix)[6], " ", colnames(results6[5])))

rownames(results) <- rownames(results1)
# add raw expression data
length(rownames(y1))==length(results[,1])
all(rownames(y1) == results[,1]) # if FALSE, have to do merge, not cbind!


pdf("../other/root_TMMnormLOGcpm_boxplot.pdf", onefile=TRUE, family= "Helvetica")
boxplot(cpm(y1, log=TRUE, prior.count=0.5), las=2, ylab="TMM normalized log_cpm values with prior.counts 0.5", col=col)
dev.off()

results.raw <- merge(results, y1$counts, by.x="row.names", by.y="row.names", all.x= TRUE, all.y= FALSE, sort= FALSE)
head(results.raw)
colnames(results.raw)[1] <- c("GeneID")
write.table(results.raw, file="../output/Monocultures/root_logFC_padj_rawReads.txt", sep="\t", quote=TRUE, row.names=FALSE)

#instead of raw counts here I export TMM normalized cpm values with prior.counts 0.5 
results.tmm <- merge(results, cpm(y1, log=TRUE, prior.count=0.5), by.x="row.names", by.y="row.names", all.x= TRUE, all.y= FALSE, sort= FALSE)
head(results.tmm)
colnames(results.tmm)[1] <- c("GeneID")
write.table(results.tmm, file="../output/Monocultures/root_logFC_padj_TMMcpm.txt", sep="\t", quote=TRUE, row.names=FALSE)

#unfiltered TMM-normalized results output
y_GSEA <- calcNormFactors(y)
y_GSEA <- cpm(y_GSEA, log=TRUE, prior.count=0.5)
y_GSEA <- cbind(rownames(y_GSEA), rep("NA", length(rownames(y_GSEA))), y_GSEA)
colnames(y_GSEA)[1:2] <- c("NAME","DESCRIPTION")

head(y_GSEA)  
write.table(y_GSEA, file="../output/Monocultures/root_TMM_GSEA_noDESC.txt", sep="\t", quote=FALSE, row.names=FALSE)




################################################
## LIMMA VENN DIAGRAMS FOR SELECTED CONTRASTS ##
################################################

## Define contrasts i.e. comparisons between groups:
### !!! ###

# modified vennDiagram function without box and circle borders (I further modified this one https://gist.github.com/mevers/9c846e6293db44dd37695c46b8f2b6a2#file-my-venndiagram-r)

my.vennDiagram <- function (object, include = "both", names = NULL, 
                            mar = rep(1,4), cex = c(1.5, 1, 0.7), 
                            lwd = 1, circle.col = NULL, counts.col = NULL,
                            show.include = NULL, ...)
{
  include <- as.character(include)
  LenInc <- min(length(include), 2)
  if (is(object, "VennCounts")) {
    include <- include[1]
    LenInc <- 1
  }
  else {
    if (LenInc > 1)
      z2 <- vennCounts(object, include = include[2])[,
                                                     "Counts"]
    object <- vennCounts(object, include = include[1])
  }
  z <- object[, "Counts"]
  nsets <- ncol(object) - 1
  if (nsets > 5)
    stop("Can't plot Venn diagram for more than 5 sets")
  VennZone <- object[, 1:nsets, drop = FALSE]
  VennZone <- apply(VennZone, 1, function(x) paste(x, sep = "",
                                                   collapse = ""))
  names(z) <- VennZone
  if (length(include) == 2)
    names(z2) <- VennZone
  if (is.null(names))
    names <- colnames(object)[1:nsets]
  FILL.COL <- TRUE
  if (is.null(circle.col)) {
    circle.col <- par("col")
    FILL.COL <- FALSE
  }
  if (length(circle.col) < nsets)
    circle.col <- rep(circle.col, length.out = nsets)
  if (is.null(counts.col))
    counts.col <- par("col")
  if (length(counts.col) < LenInc)
    counts.col <- rep(counts.col, length.out = LenInc)
  if (is.null(show.include))
    show.include <- as.logical(LenInc - 1)
  old.par <- par()$mar
  on.exit(par(mar = old.par))
  par(mar = mar)
  if (nsets <= 3) {
    plot(x = 0, y = 0, type = "n", xlim = c(-4, 4), ylim = c(-4,
                                                             4), xlab = "", ylab = "", axes = FALSE, ...)
    theta <- 2 * pi * (0:360)/360
    xcentres <- switch(nsets, 0, c(-1, 1), c(-1, 1, 0))
    ycentres <- switch(nsets, 0, c(0, 0), c(1, 1, -2)/sqrt(3))
    r <- 1.5
    xtext <- switch(nsets, -1.2, c(-1.2, 1.2), c(-1.2, 1.2,
                                                 0))
    ytext <- switch(nsets, 1.8, c(1.8, 1.8), c(2.4, 2.4,
                                               -3))
    for (circle in 1:nsets) {
      if (!FILL.COL)
        lines(xcentres[circle] + r * cos(theta), ycentres[circle] +
                r * sin(theta), lwd = lwd, col = circle.col[circle])
      if (FILL.COL) {
        RGB <- col2rgb(circle.col[circle])/255
        ALPHA <- 0.06
        RGB.ALP <- rgb(RGB[1, 1], RGB[2, 1], RGB[3, 1],
                       alpha = ALPHA)
        polygon(xcentres[circle] + r * cos(theta), ycentres[circle] +
                  r * sin(theta), border = NA,
                lwd = lwd, col = RGB.ALP)
      }
      text(xtext[circle], ytext[circle], names[circle],
           cex = cex)
    }
    #switch(nsets, rect(-3, -2.5, 3, 2.5), rect(-3, -2.5,
    #    3, 2.5), rect(-3, -3.5, 3, 3.3))
    showCounts <- switch(nsets, function(counts, cex, adj,
                                         col, leg) {
      text(2.3, -2.1, sprintf(""), cex = cex, col = col,
           adj = adj)
      text(0, 0, counts[2], cex = cex, col = col, adj = adj)
      if (show.include) text(-2.3, -2.1, leg, cex = cex,
                             col = col, adj = adj)
    }, function(counts, cex, adj, col, leg) {
      text(2.3, -2.1, sprintf(""), cex = cex, col = col,
           adj = adj)
      text(1.5, 0.1, counts[2], cex = cex, col = col, adj = adj)
      text(-1.5, 0.1, counts[3], cex = cex, col = col,
           adj = adj)
      text(0, 0.1, counts[4], cex = cex, col = col, adj = adj)
      if (show.include) text(-2.3, -2.1, leg, cex = cex,
                             col = col, adj = adj)
    }, function(counts, cex, adj, col, leg) {
      text(2.5, -3, sprintf(""), cex = cex, col = col, adj = adj)
      text(0, -1.7, counts[2], cex = cex, col = col, adj = adj)
      text(1.5, 1, counts[3], cex = cex, col = col, adj = adj)
      text(0.75, -0.35, counts[4], cex = cex, col = col,
           adj = adj)
      text(-1.5, 1, counts[5], cex = cex, col = col, adj = adj)
      text(-0.75, -0.35, counts[6], cex = cex, col = col,
           adj = adj)
      text(0, 0.9, counts[7], cex = cex, col = col, adj = adj)
      text(0, 0, counts[8], cex = cex, col = col, adj = adj)
      if (show.include) text(-2.5, -3, leg, cex = cex,
                             col = col, adj = adj)
    })
    if (LenInc == 1)
      adj <- c(0.5, 0.5)
    else adj <- c(0.5, 0)
    print(z)
    showCounts(counts = z, cex = cex[1], adj = adj, col = counts.col[1],
               leg = include[1])
    if (LenInc == 2)
      showCounts(counts = z2, cex = cex[1], adj = c(0.5,
                                                    1), col = counts.col[2], leg = include[2])
    return(invisible())
  }
  plot(c(-20, 420), c(-20, 420), type = "n", axes = FALSE,
       ylab = "", xlab = "", ...)
  relocate_elp <- function(e, alpha, x, y) {
    phi <- (alpha/180) * pi
    xr <- e[, 1] * cos(phi) + e[, 2] * sin(phi)
    yr <- -e[, 1] * sin(phi) + e[, 2] * cos(phi)
    xr <- x + xr
    yr <- y + yr
    cbind(xr, yr)
  }
  if (4 == nsets) {
    #rect(-20, -20, 420, 400)
    elps <- cbind(162 * cos(seq(0, 2 * pi, len = 1000)),
                  108 * sin(seq(0, 2 * pi, len = 1000)))
    if (!FILL.COL) {
      polygon(relocate_elp(elps, 45, 130, 170), border = NA,
              lwd = lwd)
      polygon(relocate_elp(elps, 45, 200, 200), border = NA,
              lwd = lwd)
      polygon(relocate_elp(elps, 135, 200, 200), border = NA,
              lwd = lwd)
      polygon(relocate_elp(elps, 135, 270, 170), border = NA,
              lwd = lwd)
    }
    if (FILL.COL) {
      RGB <- col2rgb(circle.col)/255
      ALPHA <- 0.06
      RGB.ALP1 <- rgb(RGB[1, 1], RGB[2, 1], RGB[3, 1],
                      alpha = ALPHA)
      RGB.ALP2 <- rgb(RGB[1, 2], RGB[2, 2], RGB[3, 2],
                      alpha = ALPHA)
      RGB.ALP3 <- rgb(RGB[1, 3], RGB[2, 3], RGB[3, 3],
                      alpha = ALPHA)
      RGB.ALP4 <- rgb(RGB[1, 4], RGB[2, 4], RGB[3, 4],
                      alpha = ALPHA)
      polygon(relocate_elp(elps, 45, 130, 170), border = NA,
              lwd = lwd, col = RGB.ALP1)
      polygon(relocate_elp(elps, 45, 200, 200), border = NA,
              lwd = lwd, col = RGB.ALP2)
      polygon(relocate_elp(elps, 135, 200, 200), border = NA,
              lwd = lwd, col = RGB.ALP3)
      polygon(relocate_elp(elps, 135, 270, 170), border = NA,
              lwd = lwd, col = RGB.ALP4)
    }
    text(35, 315, names[1], cex = cex[1])
    text(138, 350, names[2], cex = cex[1])
    text(262, 347, names[3], cex = cex[1])
    text(365, 315, names[4], cex = cex[1])
    text(35, 250, z["1000"], cex = cex[2], col = counts.col[1],
    )
    text(140, 315, z["0100"], cex = cex[2], col = counts.col[1])
    text(260, 315, z["0010"], cex = cex[2], col = counts.col[1])
    text(365, 250, z["0001"], cex = cex[2], col = counts.col[1])
    text(90, 282, z["1100"], cex = cex[3], col = counts.col[1])
    text(95, 110, z["1010"], cex = cex[2], col = counts.col[1])
    text(200, 52, z["1001"], cex = cex[3], col = counts.col[1])
    text(200, 292, z["0110"], cex = cex[2], col = counts.col[1])
    text(300, 110, z["0101"], cex = cex[2], col = counts.col[1])
    text(310, 282, z["0011"], cex = cex[3], col = counts.col[1])
    text(130, 230, z["1110"], cex = cex[2], col = counts.col[1])
    text(245, 81, z["1101"], cex = cex[3], col = counts.col[1])
    text(155, 81, z["1011"], cex = cex[3], col = counts.col[1])
    text(270, 230, z["0111"], cex = cex[2], col = counts.col[1])
    text(200, 152, z["1111"], cex = cex[2], col = counts.col[1])
    text(400, 15, sprintf(""), cex = cex[1], col = counts.col[1])
    if (length(include) == 2) {
      text(35, 238, z2["1000"], cex = cex[2], col = counts.col[2])
      text(140, 304, z2["0100"], cex = cex[2], col = counts.col[2])
      text(260, 304, z2["0010"], cex = cex[2], col = counts.col[2])
      text(365, 238, z2["0001"], cex = cex[2], col = counts.col[2])
      text(90, 274, z2["1100"], cex = cex[3], col = counts.col[2])
      text(95, 100, z2["1010"], cex = cex[2], col = counts.col[2])
      text(200, 43, z2["1001"], cex = cex[3], col = counts.col[2])
      text(200, 280, z2["0110"], cex = cex[2], col = counts.col[2])
      text(300, 100, z2["0101"], cex = cex[2], col = counts.col[2])
      text(310, 274, z2["0011"], cex = cex[3], col = counts.col[2])
      text(130, 219, z2["1110"], cex = cex[2], col = counts.col[2])
      text(245, 71, z2["1101"], cex = cex[3], col = counts.col[2])
      text(155, 72, z2["1011"], cex = cex[3], col = counts.col[2])
      text(270, 219, z2["0111"], cex = cex[2], col = counts.col[2])
      text(200, 140, z2["1111"], cex = cex[2], col = counts.col[2])
      text(400, -2, sprintf(""), cex = cex[1], col = counts.col[2])
      if (show.include) {
        text(10, 15, include[1], cex = cex[1], col = counts.col[1])
        text(10, -2, include[2], cex = cex[1], col = counts.col[2])
      }
    }
    return(invisible())
  }
  #rect(-20, -30, 430, 430)
  elps <- cbind(150 * cos(seq(0, 2 * pi, len = 1000)), 60 *
                  sin(seq(0, 2 * pi, len = 1000)))
  if (!FILL.COL) {
    polygon(relocate_elp(elps, 90, 200, 250), border = NA,
            lwd = lwd)
    polygon(relocate_elp(elps, 162, 250, 220), border = NA,
            lwd = lwd)
    polygon(relocate_elp(elps, 234, 250, 150), border = NA,
            lwd = lwd)
    polygon(relocate_elp(elps, 306, 180, 125), border = NA,
            lwd = lwd)
    polygon(relocate_elp(elps, 378, 145, 200), border = NA,
            lwd = lwd)
  }
  if (FILL.COL) {
    RGB <- col2rgb(circle.col)/255
    ALPHA <- 0.06
    RGB.ALP1 <- rgb(RGB[1, 1], RGB[2, 1], RGB[3, 1], alpha = ALPHA)
    RGB.ALP2 <- rgb(RGB[1, 2], RGB[2, 2], RGB[3, 2], alpha = ALPHA)
    RGB.ALP3 <- rgb(RGB[1, 3], RGB[2, 3], RGB[3, 3], alpha = ALPHA)
    RGB.ALP4 <- rgb(RGB[1, 4], RGB[2, 4], RGB[3, 4], alpha = ALPHA)
    RGB.ALP5 <- rgb(RGB[1, 5], RGB[2, 5], RGB[3, 5], alpha = ALPHA)
    polygon(relocate_elp(elps, 90, 200, 250), border = NA,
            lwd = lwd, col = RGB.ALP1)
    polygon(relocate_elp(elps, 162, 250, 220), border = NA,
            lwd = lwd, col = RGB.ALP2)
    polygon(relocate_elp(elps, 234, 250, 150), border = NA,
            lwd = lwd, col = RGB.ALP3)
    polygon(relocate_elp(elps, 306, 180, 125), border = NA,
            lwd = lwd, col = RGB.ALP4)
    polygon(relocate_elp(elps, 378, 145, 200), border = NA,
            lwd = lwd, col = RGB.ALP5)
  }
  text(50, 285, names[1], cex = cex[1])
  text(200, 415, names[2], cex = cex[1])
  text(350, 305, names[3], cex = cex[1])
  text(350, 20, names[4], cex = cex[1])
  text(100, -10, names[5], cex = cex[1])
  text(61, 231, z["10000"], cex = cex[2], col = counts.col[1])
  text(200, 332, z["01000"], cex = cex[2], col = counts.col[1])
  text(321, 248, z["00100"], cex = cex[2], col = counts.col[1])
  text(290, 84, z["00010"], cex = cex[2], col = counts.col[1])
  text(132, 72, z["00001"], cex = cex[2], col = counts.col[1])
  text(146, 253, z["11000"], cex = cex[3], col = counts.col[1])
  text(123, 191, z["10100"], cex = cex[3], col = counts.col[1])
  text(275, 155, z["10010"], cex = cex[3], col = counts.col[1])
  text(137, 149, z["10001"], cex = cex[3], col = counts.col[1])
  text(243, 271, z["01100"], cex = cex[3], col = counts.col[1])
  text(175, 270, z["01010"], cex = cex[3], col = counts.col[1])
  text(187, 120, z["01001"], cex = cex[3], col = counts.col[1])
  text(286, 193, z["00110"], cex = cex[3], col = counts.col[1])
  text(267, 238, z["00101"], cex = cex[3], col = counts.col[1])
  text(228, 108, z["00011"], cex = cex[3], col = counts.col[1])
  text(148, 213, z["11100"], cex = cex[3], col = counts.col[1])
  text(159, 255, z["11010"], cex = cex[3], col = counts.col[1])
  text(171, 144, z["11001"], cex = cex[3], col = counts.col[1])
  text(281, 178, z["10110"], cex = cex[3], col = counts.col[1])
  text(143, 166, z["10101"], cex = cex[3], col = counts.col[1])
  text(252, 148, z["10011"], cex = cex[3], col = counts.col[1])
  text(205, 258, z["01110"], cex = cex[3], col = counts.col[1])
  text(254, 248, z["01101"], cex = cex[3], col = counts.col[1])
  text(211, 121, z["01011"], cex = cex[3], col = counts.col[1])
  text(267, 214, z["00111"], cex = cex[3], col = counts.col[1])
  text(170, 234, z["11110"], cex = cex[3], col = counts.col[1])
  text(158, 172, z["11101"], cex = cex[3], col = counts.col[1])
  text(212, 142, z["11011"], cex = cex[3], col = counts.col[1])
  text(263, 183, z["10111"], cex = cex[3], col = counts.col[1])
  text(239, 235, z["01111"], cex = cex[3], col = counts.col[1])
  text(204, 193, z["11111"], cex = cex[2], col = counts.col[1])
  text(400, 7, sprintf(""), cex = cex[1], col = counts.col[1])
  if (length(include) == 2) {
    text(61, 220, z2["10000"], cex = cex[2], col = counts.col[2])
    text(200, 321, z2["01000"], cex = cex[2], col = counts.col[2])
    text(321, 237, z2["00100"], cex = cex[2], col = counts.col[2])
    text(290, 73, z2["00010"], cex = cex[2], col = counts.col[2])
    text(132, 61, z2["00001"], cex = cex[2], col = counts.col[2])
    text(146, 244, z2["11000"], cex = cex[3], col = counts.col[2])
    text(123, 180, z2["10100"], cex = cex[3], col = counts.col[2])
    text(275, 144, z2["10010"], cex = cex[3], col = counts.col[2])
    text(137, 143, z2["10001"], cex = cex[3], col = counts.col[2])
    text(243, 260, z2["01100"], cex = cex[3], col = counts.col[2])
    text(175, 259, z2["01010"], cex = cex[3], col = counts.col[2])
    text(187, 110, z2["01001"], cex = cex[3], col = counts.col[2])
    text(286, 186, z2["00110"], cex = cex[3], col = counts.col[2])
    text(267, 230, z2["00101"], cex = cex[3], col = counts.col[2])
    text(228, 97, z2["00011"], cex = cex[3], col = counts.col[2])
    text(148, 203, z2["11100"], cex = cex[3], col = counts.col[2])
    text(159, 249, z2["11010"], cex = cex[3], col = counts.col[2])
    text(171, 137, z2["11001"], cex = cex[3], col = counts.col[2])
    text(281, 171, z2["10110"], cex = cex[3], col = counts.col[2])
    text(143, 155, z2["10101"], cex = cex[3], col = counts.col[2])
    text(252, 137, z2["10011"], cex = cex[3], col = counts.col[2])
    text(205, 247, z2["01110"], cex = cex[3], col = counts.col[2])
    text(254, 242, z2["01101"], cex = cex[3], col = counts.col[2])
    text(211, 112, z2["01011"], cex = cex[3], col = counts.col[2])
    text(267, 207, z2["00111"], cex = cex[3], col = counts.col[2])
    text(170, 223, z2["11110"], cex = cex[3], col = counts.col[2])
    text(158, 162, z2["11101"], cex = cex[3], col = counts.col[2])
    text(212, 133, z2["11011"], cex = cex[3], col = counts.col[2])
    text(263, 172, z2["10111"], cex = cex[3], col = counts.col[2])
    text(239, 228, z2["01111"], cex = cex[3], col = counts.col[2])
    text(204, 182, z2["11111"], cex = cex[2], col = counts.col[2])
    text(400, -10, sprintf(""), cex = cex[1], col = counts.col[2])
    if (show.include) {
      text(10, 7, include[1], cex = cex[1], col = counts.col[1])
      text(10, -10, include[2], cex = cex[1], col = counts.col[2])
    }
  }
  invisible()
}





contrastMatrix = makeContrasts("roots_inoculated_PS216_2h_A-roots_noninoculated_A",
                               "roots_inoculated_PS218_2h_A-roots_noninoculated_A",
                               "roots_inoculated_PS216_26h_B-roots_noninoculated_B",
                               "roots_inoculated_PS218_26h_B-roots_noninoculated_B",
                               "(roots_inoculated_PS216_2h_A+roots_inoculated_PS216_26h_B)-(roots_noninoculated_A+roots_noninoculated_B)",
                               "(roots_inoculated_PS218_2h_A+roots_inoculated_PS218_26h_B)-(roots_noninoculated_A+roots_noninoculated_B)",
                               levels=design)

fit2 = contrasts.fit(fit, contrastMatrix)
fit2 <- eBayes(fit2)

dt <- decideTests(fit2, lfc = 1, adjust.method = "BH") #defaults adjust.method = "BH", p.value = 0.05, lfc=1
summary(dt)
colnames(dt)

## change dt colnames
colnames(dt) <- c("PS216 2h", "PS218 2h", "PS216 26h", "PS218 26h", "PS216 bothTPs", "PS218 bothTPs")

# single strains, 2 time points
vennDiagram(dt[, 1:4], include=c("up", "down"), counts.col=c("red", "blue"), circle.col = c("orange", "blue", "red", "magenta"), cex=c(1,0.8,0.7))
my.vennDiagram(dt[, 1:4], include=c("up", "down"), counts.col=c("red", "blue"), circle.col = c("orange", "blue", "red", "magenta"), cex=c(1,0.8,0.7))

pdf(file = "../output/Monocultures/limma_venns_Exp1_roots.pdf", family = "Helvetica-Narrow")
my.vennDiagram(dt[, 1:4], include=c("up", "down"), counts.col=c("red", "blue"), circle.col = c("orange", "blue", "red", "magenta"), cex=c(1,0.8,0.7))
my.vennDiagram(dt[, 1:2], include=c("up", "down"), counts.col=c("red", "blue"), circle.col = c("orange", "blue"), cex=c(1,0.8,0.7))
my.vennDiagram(dt[, 3:4], include=c("up", "down"), counts.col=c("red", "blue"), circle.col = c("red", "magenta"), cex=c(1,0.8,0.7))
dev.off()


























##############################
############# LEAVES

# combine all counts in to an expression matrix
files <- list.files(path="../output/mapping/Monocultures/", pattern=".ReadsPerGene.out.tab",  full.names= TRUE) 

geneNames <- read.table(files[1], header = FALSE, sep = "\t", skip= 4)[,1]

counts1 <- NULL
colNs <- NULL
counti <- NULL

for(i in 1:length(files)) {
  counti <- read.table(files[i], header = FALSE, sep = "\t", skip= 4)[,4] #  row 4 are reverse stranded reads
  colNi <- unlist(strsplit(unlist(strsplit(unlist(strsplit(files[i], ".", fixed = TRUE))[3], "/", fixed = TRUE))[4],"_"))[2]  # change if needed to capture sample name
  
  counts1 <- cbind(counts1, counti)
  colNs <- cbind(colNs , colNi) 
}

colNs

colnames(counts1) <- as.vector(colNs)
rownames(counts1) <- geneNames

head(counts1) # check if it makes sense!

# based on centrifuge results as well as first pass of DE analysis, sample BS472 has to be discarded
(colnames(counts1)=="BS472")[22]

counts <- counts1[,-22]

library ("limma")
library("edgeR")
library("stringr")

colnames(counts)

# read phenodata from analytes.txt/phenodata.txt (samples used for RNA-Seq)
phenodata <- read.table("../input/phenodata.txt",  row.names=1, header = TRUE)
dim(phenodata)
head(phenodata)

## selecting columns for LEAF samples!!
phenodata <- phenodata[24:35, ]
counts <- counts[, 24:35]


#check if samples are in same order; if TRUE it's OK :)
all(rownames(phenodata)==colnames(counts))


# assign groups to samples
group <- as.factor(as.vector(phenodata$treatment))
levels(group) 

x <- counts

## Create a DGEList object for limma statistical analysis and specify library size i.e. number of reads sequenced per sample
# default library sizes (sum of mapped reads)!!
y <- DGEList(counts=x, group=group)


#####
#define color palette for the samples in graphs
col= c(rep("grey",4), rep("green1",4), rep("green4",4))

#check density plot
pdf("../output/Monocultures/leaf_density_plot_before_lowexpr_filter.pdf", onefile=TRUE, family= "Helvetica")
nsamples <- ncol(x)

lcpm <- log(as.matrix(x),10)
plot(density(lcpm), col=col[1], lwd=2, ylim=c(0,0.4), las=2, main="", xlab="")
title(main="A. BEFORE REMOVAL", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
  
}
legend("topright", colnames(lcpm), text.col=col, cex = 0.75, bty="n")
dev.off()

## Filter genes based on expression -- keep genes (rows) that have over 50 counts in at least 4 samples == one group
# filtering is necessary for voom method to work properly
keep.exprs <- rowSums(y$counts>50)>=4


## Normalization of dataset for different library sizes
y1 <- y[keep.exprs, , keep.lib.sizes=TRUE]
y1 <- calcNormFactors(y1)

## Plot QC plots using different functions e.g.:

pdf("../output/Monocultures/leaf_log10rawcounts_boxplot.pdf", onefile=TRUE, family= "Helvetica")
boxplot(log(y$counts+1,10), las=2, ylab="log10(counts)", col=col)
dev.off()

pdf("../output/Monocultures/leaf_log10filteredcounts_boxplot.pdf", onefile=TRUE, family= "Helvetica")
boxplot(log(y1$counts+1,10), las=2, ylab="log10(counts)", col=col)
dev.off()

#density plots before and after removing low expressed genes
pdf("../output/Monocultures/leaf_norm_counts_raw&filtered_densityplots.pdf", onefile=TRUE, family= "Helvetica")
opar <- par()
par(mfrow=c(1,2), cex = 0.6)
nsamples <- ncol(x)

lcpm <- log(as.matrix(x),10)
plot(density(lcpm), col=col[1], lwd=2, ylim=c(0,0.4), las=2, main="", xlab="")
title(main="A. BEFORE REMOVAL", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
#legend("topright", colnames(lcpm), text.col=col, bty="n")

lcpm <- log(as.matrix(y1),10)
plot(density(lcpm), col=col[1], lwd=2, ylim=c(0,1), las=2, main="", xlab="")
title(main="B. AFTER REMOVAL", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", colnames(lcpm), text.col=col, bty="n")
par(opar)
dev.off()


#MDS (PCA-like) graph
pdf("../output/Monocultures/leaf_norm_counts_filtered_MDS.pdf", onefile=TRUE, family= "Helvetica")
plotMDS(y1, labels=colnames(y1), col = col, cex = 0.6)
dev.off()

pdf("../output/Monocultures/leaf_norm_cf_lcpm_MDS.pdf", onefile=TRUE, family= "Helvetica")
lcpm <- log(as.matrix(y1),10)
plotMDS(lcpm, labels=colnames(lcpm), col = col, cex = 0.6)
dev.off()

###########
## limma-voom protocol
# Create design matrix
design <- model.matrix(~0+group)
colnames(design) <- levels(group)

# limma voom fit for filtered RNA-seq dataset (y1)
pdf("../output/Monocultures/leaf_voom_mean-variance_trend.pdf", onefile=TRUE, family= "Helvetica")
v <- voom(y1,design,plot=TRUE)
dev.off()
fit <- lmFit(v, design)

## Define contrasts i.e. comparisons between groups:
### !!! ###
contrastMatrix = makeContrasts("leaves_inoculated_PS216_26h-leaves_noninoculated",
                               "leaves_inoculated_PS218_26h-leaves_noninoculated",
                               levels=design)

fit2 = contrasts.fit(fit, contrastMatrix)
## Check DEG in contrasts (adj p-val cutoff 0.05, |logFC| > 0) using alternative to eBayes (just as control)
tfit <- treat(fit2)
tfit
logFCcut <- 1
dt <- decideTests(tfit, lfc=logFCcut, adjust.method = "BH") #defaults adjust.method = "BH", p.value = 0.05, lfc=0
summary(dt)
colnames(dt)

## eBayes statistics calculation
fit2 <- eBayes(fit2)
pdf("../output/Monocultures/leaf_SIGMA_vs_A_plot.pdf", onefile=TRUE, family= "Helvetica")
plotSA(fit2)
dev.off()

## make results table
results1 <- topTable(fit2, coef=1, number=1000000, sort.by="none")
results2 <- topTable(fit2, coef=2, number=1000000, sort.by="none")

#make one expression matrix with logFCs and adj.P.vals
results <- cbind(results1[,1], results1[,5],
                 results2[,1], results2[,5])

colnames(results) <- c(paste(colnames(contrastMatrix)[1], " ", colnames(results1[1])),
                       paste(colnames(contrastMatrix)[1], " ", colnames(results1[5])),
                       paste(colnames(contrastMatrix)[2], " ", colnames(results2[1])),
                       paste(colnames(contrastMatrix)[2], " ", colnames(results2[5])))

rownames(results) <- rownames(results1)
# add raw expression data
length(rownames(y1))==length(results[,1])
all(rownames(y1) == results[,1]) # if FALSE, have to do merge, not cbind!


pdf("../output/Monocultures/leaf_TMMnormLOGcpm_boxplot.pdf", onefile=TRUE, family= "Helvetica")
boxplot(cpm(y1, log=TRUE, prior.count=0.5), las=2, ylab="TMM normalized log_cpm values with prior.counts 0.5", col=col)
dev.off()

results.raw <- merge(results, y1$counts, by.x="row.names", by.y="row.names", all.x= TRUE, all.y= FALSE, sort= FALSE)
head(results.raw)
colnames(results.raw)[1] <- c("GeneID")
write.table(results.raw, file="../output/Monocultures/leaf_logFC_padj_rawReads.txt", sep="\t", quote=TRUE, row.names=FALSE)

#instead of raw counts here I export TMM normalized cpm values with prior.counts 0.5 
results.tmm <- merge(results, cpm(y1, log=TRUE, prior.count=0.5), by.x="row.names", by.y="row.names", all.x= TRUE, all.y= FALSE, sort= FALSE)
head(results.tmm)
colnames(results.tmm)[1] <- c("GeneID")
write.table(results.tmm, file="../output/Monocultures/leaf_logFC_padj_TMMcpm.txt", sep="\t", quote=TRUE, row.names=FALSE)

#unfiltered TMM-normalized results output
y_GSEA <- calcNormFactors(y)
y_GSEA <- cpm(y_GSEA, log=TRUE, prior.count=0.5)
y_GSEA <- cbind(rownames(y_GSEA), rep("NA", length(rownames(y_GSEA))), y_GSEA)
colnames(y_GSEA)[1:2] <- c("NAME","DESCRIPTION")

head(y_GSEA)  
write.table(y_GSEA, file="../output/Monocultures/leaf_TMM_GSEA_noDESC.txt", sep="\t", quote=FALSE, row.names=FALSE)





################################################
## LIMMA VENN DIAGRAMS FOR SELECTED CONTRASTS ##
################################################

# find my.vennDiagram above in the ROOT section


contrastMatrix = makeContrasts("leaves_inoculated_PS216_26h-leaves_noninoculated",
                               "leaves_inoculated_PS218_26h-leaves_noninoculated",
                               levels=design)

fit2 = contrasts.fit(fit, contrastMatrix)
fit2 <- eBayes(fit2)

dt <- decideTests(fit2, lfc = 1, adjust.method = "BH") #defaults adjust.method = "BH", p.value = 0.05, lfc=1
summary(dt)
colnames(dt)

## change dt colnames
colnames(dt) <- c("PS216 26h", "PS218 26h")

# single strains
vennDiagram(dt, include=c("up", "down"), counts.col=c("red", "blue"), circle.col = c("brown", "blue4"), cex=c(1,0.8,0.7))
my.vennDiagram(dt, include=c("up", "down"), counts.col=c("red", "blue"), circle.col = c("brown", "magenta4"), cex=c(1,0.8,0.7))

pdf(file = "../output/Monocultures/limma_venns_Exp1_leaves.pdf", family = "Helvetica-Narrow")
my.vennDiagram(dt, include=c("up", "down"), counts.col=c("red", "blue"), circle.col = c("brown", "magenta4"), cex=c(1,0.8,0.7))
dev.off()




sessionInfo()
# R version 4.1.1 (2021-08-10)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19045)
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=Slovenian_Slovenia.1250  LC_CTYPE=Slovenian_Slovenia.1250    LC_MONETARY=Slovenian_Slovenia.1250
# [4] LC_NUMERIC=C                        LC_TIME=Slovenian_Slovenia.1250    
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] stringr_1.5.0 edgeR_3.36.0  limma_3.50.3 
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.9        compiler_4.1.1    later_1.3.0       urlchecker_1.0.1  prettyunits_1.1.1 profvis_0.3.7     remotes_2.4.2    
# [8] tools_4.1.1       digest_0.6.31     pkgbuild_1.4.0    pkgload_1.3.2     lattice_0.20-45   memoise_2.0.1     lifecycle_1.0.3  
# [15] rlang_1.0.6       shiny_1.7.4       cli_3.4.1         rstudioapi_0.14   fastmap_1.1.0     fs_1.5.2          htmlwidgets_1.6.0
# [22] devtools_2.4.5    grid_4.1.1        locfit_1.5-9.6    glue_1.6.2        R6_2.5.1          processx_3.8.0    sessioninfo_1.2.2
# [29] callr_3.7.3       purrr_0.3.5       magrittr_2.0.3    ps_1.7.2          promises_1.2.0.1  ellipsis_0.3.2    htmltools_0.5.4  
# [36] usethis_2.1.6     mime_0.12         xtable_1.8-4      httpuv_1.6.7      stringi_1.7.8     miniUI_0.1.1.1    cachem_1.0.6     
# [43] crayon_1.5.2 