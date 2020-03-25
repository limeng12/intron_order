## Intron Splicing Order
http://intron-splicing-order.online:3838/iso/

### Code Ocean
This repository is compatibility with Code Ocean to facilitate reproducibility

To reproduce the result in Code Ocean:
1. Git clone this repository into a `capsule` in Code Ocean
2. Add the prerequisite R packages below in the environment setting in capsule.
3. Set `code/run_all.R` as file to run
4. Click the `Reproducible Run`

### Intermediate results
This repository also includes intermediate results and data for this work. 

The result and annotations for each organism are stored seperatelly in each folder `human`,`zebrafish`,`fly`,`pombe`,`yeast`. 

 `all_MLO_Plot.pdf` stores intron splicing read count graph order by most likely order.
 
 `all_MLO_Table.pdf` stores intron splicing frequency matrix order by most likely order.
 
 `*no_thick.bed` stores transcripts annotations.
 
 `gene_id_tran*` stores gene_id gene symbol and transcript id.
 
 `best_order.tsv` stores the calculated most likely order for this organism. 
 
 `iso_pair` stores calculated intron splicing order pairs. 
 
## Prerequisite

### R and R packages

packages within R
```
install.packages(c("readr","Rcpp","dplyr","igraph","dbscan","stringr","gtools","rstudioapi","gridExtra") )
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
java -jar java/isoLarge.jar  anno/hg19_gencode_from_ucsc_nothick_nocds.bed <bam_file> <output_file> <optional INT e.g. 20>
```
The last parameter is the minium length of nucleotides aligned in intron side of intron-exon junction

Please put the output file under `data/`, since the R code will treat data/ as directory of intron splicing order pairs files. 

### 3. Calculating most likely intron splicing orders
If users are not working with Rstudio, then will need to edit the run.R to change the working dir to `intron_order`

Source the below R script in Rstudio.
```
intron_order/code/run_human.R
```
### Suggestions and comments are more than welcome:  limeng@picb.ac.cn or use issues



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
