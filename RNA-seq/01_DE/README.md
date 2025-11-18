# Differential expression analysis

## Inputs

* experimental design file (phenodata.txt)
* mapped read counts files (.ReadsPerGene.out.tab from the STAR mapping stage)

## Outputs

* 
* 

## Essential parts of the analysis scripts

[edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)

```
x <- counts
y <- DGEList(counts=x, group=group)
keep.exprs <- rowSums(y$counts>100)>=4
y1 <- y[keep.exprs, , keep.lib.sizes=TRUE]
y1 <- calcNormFactors(y1)
```

[limma](https://www.bioconductor.org/packages/release/bioc/html/limma.html)

```
design <- model.matrix(~0+group)
colnames(design) <- levels(group)
v <- voom(y1,design,plot=TRUE)
fit <- lmFit(v, design)
contrastMatrix = makeContrasts("treatment1-control",
                               "treatment2-control",
                               ...,
                               levels=design)  # replace treatment and control by corresponding groups
fit2 = contrasts.fit(fit, contrastMatrix)
fit2 <- eBayes(fit2)
resultsX <- topTable(fit2, coef=X, number=1000000, sort.by="none") # replace X by exact coef number
```

