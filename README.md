## Intron Splicing Order
http://intron-splicing-order.online:3838/iso/

## Prerequisite

### R and R packages

packages within R
```
install.packages(c("readr","Rcpp","dplyr","igraph","dbscan","stringr","gtools","rstudioapi","gridExtra") )
```
```
BiocManager::install("lpsymphony")
```

### Calling intron splicing order from users' own BAM files, then need install JAVA (JRE or JDK)
Oracle JDK8/JRE8

## Steps

### 1. Aligning FASTQ reads with splice-wise aligner. 
STAR, minimap2, et.al, then index the bam file
```
samtools index <Bam file>
```

### 2. Calculating intron splicing order pairs using the custome java program
```
java -jar java/isoLarge.jar -i anno/hg19_gencode_from_ucsc_nothick_nocds.bed -ibam <bam_file> -o <output_file> -t <optional INT e.g. 90>
```
The last parameter is the minium length of nucleotides aligned in intron side of intron-exon junction

Please put the output file under `data/`, since the R code will treat data/ as directory of intron splicing order pairs files. 

Output format

| Column    | Meaning |
| ----------|:-------:|
| Column #1 |     Transcript id                                             |
| Column #2 |     Intron 1, the coordinate of relatively slower spliced intron       |
| Column #3 |     Intron 2, the coordinate of relatively faster spliced introns (also include detected junctions)      |
| Column #4 |     Strand       |
| Column #5 |     Deprecated |
| Column #6 |     Read count supports this intron splicing order pair (intron 1 spliced after intron 2) |
| Column #7 |     Read count supports both two introns were spliced |

### 3. Calculating most likely intron splicing orders
If users are not working with Rstudio, then will need to edit the run.R to change the working dir to `intron_order`

Source the below R script in Rstudio.
```
intron_order/code/run_human.R
```


## For other genomes that are not included here

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

### Suggestions and comments are welcome:  limeng49631@aliyun.com or use issues


