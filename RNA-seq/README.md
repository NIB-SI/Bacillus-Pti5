# RNA-seq analysis pipeline

## input files

- genome fasta
- general feature format (gff3)
- gene set file (gmt)
- fastq


## pipeline

### QC

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

### Trimming

PHRED basecall quality score > 20

[TrimGalore](https://github.com/FelixKrueger/TrimGalore)


### Taxonomic classification

[Centrifuge](https://github.com/DaehwanKimLab/centrifuge)

[Pavian](https://github.com/fbreitwieser/pavian)

### Mapping

[STAR](https://github.com/alexdobin/STAR)

     
### DE

[edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)


[limma](https://www.bioconductor.org/packages/release/bioc/html/limma.html)


### GSEA

[GSEA](https://www.gsea-msigdb.org/gsea/index.jsp)

[biokit](https://github.com/martingarridorc/biokit)

### Network analysis
[SKM](https://skm.nib.si/)

[SKM tools](https://github.com/NIB-SI/skm-tools)

[DiNAR](https://github.com/NIB-SI/DiNAR)

[Cytoscape](https://cytoscape.org/)

