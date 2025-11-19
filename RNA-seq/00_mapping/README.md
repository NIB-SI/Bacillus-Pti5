# Quality control and read mapping



## Input files

* reference potato genome file (fasta)
* reference potato genome annotation file (gff3)
* Illumina paired-end RNA-seq reads (fastq)



## Most relevant output files



### Quality control (QC)

Mutliqc reports including stats of all qc steps:

* ./output/Monocultures/multiqc/multiqc\_report.html
* ./output/Interactions/multiqc/multiqc\_report.html



### STAR mapping

In directories ./output/Monocultures/STAR/ and ./output/Interactions/STAR/



* \*.ReadsPerGene.out.tab - raw counts tables
* \*.Log.final.out - mapping statistics
* \*.Log.out -







## Steps and tools used





### Taxonomic classification

* [Centrifuge](https://github.com/DaehwanKimLab/centrifuge)
* [Pavian](https://github.com/fbreitwieser/pavian)



### QC

* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [MultiQC](https://github.com/MultiQC/MultiQC)



### Trimming

PHRED basecall quality score > 20

* [TrimGalore](https://github.com/FelixKrueger/TrimGalore)



### Mapping

counting only uniquely mapped paired reads

* [STAR](https://github.com/alexdobin/STAR)
