# Quality control and read mapping 

## 

## Input files

* reference potato genome file (fasta)
* reference potato genome annotation file (gff3)
* Illumina paired-end RNA-seq reads (fastq)



## Steps and tools used





### Taxonomic classification

* [Centrifuge](https://github.com/DaehwanKimLab/centrifuge)
* [Pavian](https://github.com/fbreitwieser/pavian)



### Quality control (QC)

* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [MultiQC](https://github.com/MultiQC/MultiQC)



### Trimming

PHRED basecall quality score > 20

* [TrimGalore](https://github.com/FelixKrueger/TrimGalore)



### Mapping

counting only uniquely mapped paired reads

* [STAR](https://github.com/alexdobin/STAR)
