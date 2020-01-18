# intron_order

## Prerequisite

### R packages, better use Rstudio as I build the web viewer using shiny

packages within R

install.packages(c("shiny","networkD3","readr","DT","dplyr","igraph","dbscan","stringr","gtools") )


package within Bioconductor

install.packages("BiocManager")

BiocManager::install("Sushi")


### If user prepare to call intron splicing order from their own BAM files, need install JAVA
ORACLE JDK8


## Steps

### 1. Align FASTQ reads with splice wise aligner. 
STAR, minimap2, et.al

'samtools index <Bam file>'

### 2. Calculated intron splicing order pairs
'java -jar java/isoLarge.jar  anno/hg19_gencode_from_ucsc_nothick_nocds.bed  <bam_file> <output_file>'

put the output file under data/iso_3rd/

### 3. Build intron splicing order graph and matrix
Edit the run.R to change the working dir

Source the below R script in Rstudio.

run.R

### Suggestions and comments are welcome:  limeng@picb.ac.cn




## For non human genome

### gene_id and transcript id file are downloaded from ENSEMBL bioMart

### script to replace CDS position into TSS and TES in bed, with example:
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$1"\t"$2"\t"$9"\t"$10"\t"$11"\t"$12}' hg19_gencode_from_ucsc.bed >

hg19_gencode_from_ucsc_nothick_nocds.bed

### Convert bed file into introns, with example:
sh code/convert_bed_to_introns.sh hg19_gencode_from_ucsc.bed > hg19_gencode_intron_from_ucsc.bed



