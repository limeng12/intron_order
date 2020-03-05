
## Prerequisite

### R and R packages

packages within R
```
install.packages(c("readr","Rcpp","DT","dplyr","igraph","dbscan","stringr","gtools") )
```

### If need call intron splicing order from users' own BAM files, then need install JAVA (JRE or JDK)
ORACLE JDK8/JRE8

## Quick started
```
git clone https://github.com/limeng12/intron_order.git
```
Source the below file in Rstudio
```
intron_order/run.R
```

## Steps

### 1. Align FASTQ reads with splice-wise aligner. 
STAR, minimap2, et.al, then index the bam file
```
samtools index <Bam file>
```

### 2. Calculated intron splicing order pairs using the custome java program
```
java -jar java/isoLarge.jar  anno/hg19_gencode_from_ucsc_nothick_nocds.bed <bam_file> <output_file> <optional INT e.g. 20>
```
The last parameter is the minium length of nucleotides aligned in intron side of intron-exon junction

Please put the output file under `data/`, since the R code will treat data/ as directory of intron splicing order pairs files. 

### 3. Build intron splicing order graph and matrix
If users are not working with Rstudio, will need to edit the run.R to change the working dir to `intron_order`

Source the below R script in Rstudio.
```
intron_order/run.R
```
### Suggestions and comments are more than welcome:  limeng@picb.ac.cn or use issues



## For non-human genome

### Prepare intron splicing order pairs calling annotation

Java program used the HTSJDK to process BAM files, so need replace CDS position into TSS and TES in bed. 

The hg19_gencode_from_ucsc.bed can be easily downloaded using UCSC table browser directly

Script with an example:
```
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$1"\t"$2"\t"$9"\t"$10"\t"$11"\t"$12}' hg19_gencode_from_ucsc.bed > hg19_gencode_from_ucsc_nothick_nocds.bed
```


### Prepare transcription ID and gene symbol information.
This can be easily got from ENSEMBL BioMart server

Column names are:
```
gene_id,trans_id,gene_symbol
```
