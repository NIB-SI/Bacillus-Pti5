# properly merge DEG results of the two experiments and add v4 and v6 gene annotations
# Notes:
# - mapping was still done to merged v4 genome annotations GFF3 (Petek et al. 2020)
# - only uniquely mapped reads are counted (thus loosing some identical tandemly repeated genes e.g. PR1b)


#setup R

setwd("O:/DEJAVNOSTI/OMIKE/pISA-Projects/_p_Endophytes/_I_Bsubtilis/_S_NGS_mixed_cultures_exp2/_A_03_exp2-STAR_limma/scripts")
getwd()

#import the tables

exp1_shoots <- read.delim("../../../_S_NGS_monocultures/_A_RNAseq-STAR_limma/output/leaf_logFC_padj_TMMcpm.txt")
exp1_roots  <- read.delim("../../../_S_NGS_monocultures/_A_RNAseq-STAR_limma/output/root_logFC_padj_TMMcpm.txt")
exp2        <- read.delim("../output/root_logFC_padj_TMMcpm.txt")

#merge exp1 shooot and root tables
exp1 <- merge(exp1_shoots, exp1_roots, by = "GeneID", all = TRUE)
#reorder columns so that DE stats are together, remove unnecessary columns
colnames(exp1)
colnames(exp1[, c(1:5, 18:25, 6:17, 30:52)])

exp1 <- exp1[, c(1:5, 18:25, 6:17, 30:52)]


#merge exp1 and exp2
exp12 <- merge(exp1, exp2, by = "GeneID", all = TRUE)
#reorder columns
colnames(exp12)
colnames(exp12[, c(1:13, 49:58, 14:48, 59:78)])

exp12 <- exp12[, c(1:13, 49:58, 14:48, 59:78)]

#stats
dim(exp1)
dim(exp2)
dim(exp12)

#add gene functional annotations 
#translate v4 to v6 IDs using UniTato translation table 
v4v6translations <- read.delim("../input/Phureja_v4-v6.1_translations_dedupl.txt", header= TRUE, skip = 8)
v4v6translations <- v4v6translations[, c(1,2,10,3)]
mapman <- read.delim("../input/Phureja_v4-v6.1_annotations-extended.txt", header= TRUE)

fun.annot2 <- merge(v4v6translations, mapman[, c(1, 21:23)], by.x="v6.1_geneID", by.y="geneID", all.x= TRUE, all.y= TRUE, sort= FALSE)
results_matrix.addcpms.annot <- merge(fun.annot2, exp12, by.x="v4_geneID", by.y="GeneID", all.x= FALSE, all.y= TRUE, sort= FALSE)
dim(results_matrix.addcpms.annot)
### GRRRR, there were rows with duplicated geneid in the 
#fixed the dirty way in Excel by sorting and deleting duplicated rows


# get the subset of UniTatoIDs that has no v6.1 IDs
subsetNA <- results_matrix.addcpms.annot[is.na(results_matrix.addcpms.annot$Unitato_MapMan4_BINCODE), ]
subsetNA.annot <- merge(subsetNA[, -c(2:7)], fun.annot2, by.x="v4_geneID", by.y="v6.1_geneID", all.x= TRUE, all.y= FALSE, sort= FALSE)
subsetNA.annot1 <- cbind(subsetNA.annot[, 1], subsetNA.annot[, 79:84], subsetNA.annot[, 2:78] )
colnames(subsetNA.annot1)[1:2] <- c("v4_geneID", "v6.1_geneID")
colnames(subsetNA.annot1)[1:10]

subsetNotNA <- results_matrix.addcpms.annot[!is.na(results_matrix.addcpms.annot$Unitato_MapMan4_BINCODE), ]
colnames(subsetNotNA)[1:10]

#sanity check
dim(subsetNA.annot1)
dim(subsetNotNA)
dim(subsetNA.annot1)[1]+dim(subsetNotNA)[1] == dim(results_matrix.addcpms.annot)[1]

summary(colnames(subsetNotNA) == colnames(subsetNA.annot1))

results_matrix.addcpms.annot1 <- rbind(subsetNotNA, subsetNA.annot1)

colnames(results_matrix.addcpms.annot1)

write.table(results_matrix.addcpms.annot1, file = "../output/Exp1&2_DE_withMapManAnnot.tab", sep="\t", quote=TRUE, row.names=FALSE)


#sessioninfo
sessionInfo()

# R version 4.3.1 (2023-06-16 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19045)
# 
# Matrix products: default
# 
# 
# locale:
#   [1] LC_COLLATE=Slovenian_Slovenia.utf8  LC_CTYPE=Slovenian_Slovenia.utf8    LC_MONETARY=Slovenian_Slovenia.utf8 LC_NUMERIC=C                       
# [5] LC_TIME=Slovenian_Slovenia.utf8    
# 
# time zone: Europe/Ljubljana
# tzcode source: internal
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# loaded via a namespace (and not attached):
#   [1] compiler_4.3.1    fastmap_1.2.0     cli_3.6.3         htmltools_0.5.8.1 tools_4.3.1       rstudioapi_0.17.1 yaml_2.3.10       rmarkdown_2.29    knitr_1.49       
# [10] xfun_0.49         digest_0.6.37     rlang_1.1.4       evaluate_1.0.1
