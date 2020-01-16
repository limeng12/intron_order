# intron_order

## Align FASTQ reads with splice wise aligner. 
STAR, minimap2, et.al
'samtools index <Bam file>'


## Calculated intron splicing order pairs
'java -jar isoLarge.jar  anno/hg19_gencode_from_ucsc_nothick_nocds.bed  <bam_file> <output_file>'

put the output file under data/iso_3rd/


## Build intron splicing order graph and matrix
Edit the run.R to change the working dir
Run the below R script in Rstudio
run.R


