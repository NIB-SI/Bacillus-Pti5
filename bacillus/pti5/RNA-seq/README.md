# RNA-seq


## Trimming

[TrimGalore](https://github.com/FelixKrueger/TrimGalore)

PHRED basecall quality score > 20

## QC

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

## Taxonomic classification

[Centrifuge](https://github.com/DaehwanKimLab/centrifuge)

[Pavian](https://github.com/fbreitwieser/pavian)

## Mapping

[STAR](https://github.com/alexdobin/STAR)

additional alignment parameters:
```
--outFilterMultimapNmax 1

--outFilterMismatchNoverReadLmax 0.02
--quantTranscriptomeBan Singleend
--outFilterType BySJout
--alignSJoverhangMin 10
--alignSJDBoverhangMin 1
--alignIntronMin 20
--alignIntronMax 10000
--alignMatesGapMax 10000
```

## DE

[edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)

[limma](https://www.bioconductor.org/packages/release/bioc/html/limma.html)

## GSEA

[GSEA](https://www.gsea-msigdb.org/gsea/index.jsp)

[biokit](https://github.com/martingarridorc/biokit)
