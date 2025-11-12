# RNA-seq analysis pipeline

## Input files

- reference potato genome file (fasta)
- reference potato genome annotation file (gff3)
- Illumina paired-end RNA-seq reads (fastq)
- gene set file (gmt)
- phenotype labels file (cls)

## Steps and tools used

### Quality control (QC)

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

### Trimming

PHRED basecall quality score > 20

[TrimGalore](https://github.com/FelixKrueger/TrimGalore)


### Taxonomic classification of reads

[Centrifuge](https://github.com/DaehwanKimLab/centrifuge)

[Pavian](https://github.com/fbreitwieser/pavian)

### Read mapping

[STAR](https://github.com/alexdobin/STAR)
     
### Differential expression (DE)

[edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)

[limma](https://www.bioconductor.org/packages/release/bioc/html/limma.html)

### Gene set enrichment analysis (GSEA)

[GSEA](https://www.gsea-msigdb.org/gsea/index.jsp)

[biokit](https://github.com/martingarridorc/biokit)

### Network analysis

[SKM](https://skm.nib.si/)

[SKM tools](https://github.com/NIB-SI/skm-tools)

[DiNAR](https://github.com/NIB-SI/DiNAR)

[Cytoscape](https://cytoscape.org/)

