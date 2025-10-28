# RNA-seq analysis pipeline

## QC
```
fastqc -t $ncores ./data/*.fq.gz -o ../output/fastqc
multiqc ../output/ -o ../output/fastqc
```

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

## Trimming

PHRED basecall quality score > 20

[TrimGalore](https://github.com/FelixKrueger/TrimGalore)
```
for i in $(ls ./data | grep fq.gz | sed s/_[12].fq.gz// | sort -u) 
do 
  trim_galore --cores $ncores --quality 20 --paired --trim-n --gzip -o ../output/trim_galore ./data/${i}_1.fq.gz ./data/${i}_2.fq.gz
done
```

## Taxonomic classification

[Centrifuge](https://github.com/DaehwanKimLab/centrifuge)

[Pavian](https://github.com/fbreitwieser/pavian)

## Mapping

[STAR](https://github.com/alexdobin/STAR)

```
STAR --genomeLoad LoadAndExit --genomeDir ../output/STAR_index

for i in $(ls ../output/trim_galore/ | grep fq.gz | sed s/_[12]_val_[12].fq.gz// | sort -u)
do  STAR \
  --genomeDir ../output/STAR_index \
  --runThreadN $ncores \
  --quantMode GeneCounts \
  --readFilesIn  ../output/trim_galore/${i}_1_val_1.fq.gz ../output/trim_galore/${i}_2_val_2.fq.gz \
  --outFileNamePrefix ../output/unique_map/_$i. \
  --outFilterMultimapNmax 1 \
  --outFilterMismatchNoverReadLmax 0.02 \
  --outSAMtype None \
  --quantTranscriptomeBan Singleend \
  --outFilterType BySJout \
  --alignSJoverhangMin 10 \
  --alignSJDBoverhangMin 1 \
  --alignIntronMin 20 \
  --alignIntronMax 10000 \
  --alignMatesGapMax 10000 \
  --readFilesCommand pigz -cd -p $ncores
done
```          


## DE

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

## GSEA

[GSEA](https://www.gsea-msigdb.org/gsea/index.jsp)

e.g. of .cls file:
```
12 3 1
# uninoculated PS216 PS218
uninoculated uninoculated uninoculated uninoculated PS216 PS216 PS216 PS216 PS218 PS218 PS218 PS218
```

[biokit](https://github.com/martingarridorc/biokit)

  ```gseaFromStats```
