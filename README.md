
## Prerequisite

### R packages, better use Rstudio as I build the web viewer using shiny

packages within R
```
install.packages(c("shiny","networkD3","readr","DT","dplyr","igraph","dbscan","stringr","gtools") )
```

package within Bioconductor
```
install.packages("BiocManager")

BiocManager::install("Sushi")
```

### Call intron splicing order from their own BAM files, need install JAVA
ORACLE JDK8/JRE8

## Quick started
```
git clone https://github.com/limeng12/intron_order.git
```
Source the below file
```
intron_order/run.R
```

## Steps

### 1. Align FASTQ reads with splice-wise aligner. 
STAR, minimap2, et.al, then index the bam file
```
samtools index <Bam file>
```

### 2. Calculated intron splicing order pairs using a custome java program
```
java -jar java/isoLarge.jar  anno/hg19_gencode_from_ucsc_nothick_nocds.bed <bam_file> <output_file> <optional INT e.g. 20>
```
The last parameter is the minium length of nucleotides aligned in intron

Please put the output file under `data/iso_3rd/`, since the R code will treat data/iso_3rd/ as directory of intron splicing order pairs files. 
### 3. Build intron splicing order graph and matrix
If users are not working with Rstudio, will need to edit the run.R to change the working dir to `intron_order`

Source the below R script in Rstudio.
```
intron_order/run.R
```
### Suggestions and comments are welcome:  limeng@picb.ac.cn



## For non human genome

### Prepare intron splicing order pairs calling annotation

Java program used the HTSJDK to process BAM files, so need replace CDS position into TSS and TES in bed. 

The hg19_gencode_from_ucsc.bed download from UCSC table browser directly

Script with an example:
```
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$1"\t"$2"\t"$9"\t"$10"\t"$11"\t"$12}' hg19_gencode_from_ucsc.bed > hg19_gencode_from_ucsc_nothick_nocds.bed
```

`intron_order/run.R` need intron annotation to annotate the intron splicing order pairs from JAVA program


### Prepare intron number annotation
Convert transcript bed file into intron bed file, with an example from human:

The hg19_gencode_from_ucsc.bed download from UCSC table browser directly
```
sh code/convert_bed_to_introns.sh hg19_gencode_from_ucsc.bed > hg19_gencode_intron_from_ucsc.bed
```

### Prepare transcription ID and gene symbol information.
This can be got from ENSEMBL BioMart server

Column names are:
```
gene_id,trans_id,gene_symbol,trans_start,trans_end,strand,chr,gene_start,gene_end
```

### Prepare for Sushi transcription structure. 
An example of human:

The `gencode.v29lift37.annotation.gtf` file can be downloaded from ENSEMBL or GENCODE 

```
code/prepare_trans_stuc_plot.R
```

