#STAR --runThreadN 4 --runMode genomeGenerate \
#--genomeDir ~/Documents/projects/iso/result/star_index_chr21 \
#--genomeFastaFiles ~/Documents/projects/iso/result/chr21/chr21.fa \
#--sjdbGTFfile ~/Documents/projects/iso/result/chr21/gencode.v29lift37.annotation_chr21.gtf \
#--sjdbOverhang 149 \
#--genomeSAindexNbases 11
source ~/.profile

#######################################################simulation########################################
minimap2 -t 4 -ax splice ~/Documents/projects/iso/result/chr21/chr21.fa \
~/Documents/projects/iso/result/human_read_given_order_simulation_long_super_long.fasta.gz |
samtools view -@ 4 -S -b | samtools sort -@ 4 -o ./result/human_simulate_super_long.bam

samtools index ./result/human_simulate_super_long.bam

java -jar ./code/run_sh/isoLarge.jar anno/hg19_gencode_from_ucsc_nothick_nocds_chr21.bed  ./result/human_simulate_super_long.bam \
./result/iso_simulation_super_long/human_simulate_super_long_iso_large_unique_intron.tsv 


####################################################### ########################################
minimap2 -t 4 -ax splice ~/Documents/projects/iso/result/chr21/chr21.fa ~/Documents/projects/iso/result/human_read_given_order_simulation_long.fasta.gz |
samtools view -@ 4 -S -b | samtools sort -@ 4 -o  ./result/human_simulate_long.bam

samtools index ./result/human_simulate_long.bam

java -jar ./code/run_sh/isoLarge.jar anno/hg19_gencode_from_ucsc_nothick_nocds_chr21.bed  ./result/human_simulate_long.bam \
./result/iso_simulation_long/human_simulate_long_iso_large_unique_intron.tsv 

# 20 false true a3a5_bias_all_pombe.tsv
####################################################### ########################################
#ulimit -n 10000

#rm -r ./result/human_simulate_pair_STARtmp
#rm ./result/human_simulate_pair*.bam

STAR --runThreadN 4 \
--genomeDir ./result/star_index_chr21 \
--readFilesIn  ./result/human_read_given_order_simulation_pair_1.fasta.gz ./result/human_read_given_order_simulation_pair_2.fasta.gz \
--outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 \
--outFilterMismatchNoverLmax 0.04 --alignIntronMin 10 --readFilesCommand 'gunzip -c' \
--alignIntronMax 1000000 --alignMatesGapMax 1000000 \
--outSAMattributes NH HI NM MD AS nM jM jI XS --outFileNamePrefix \
./result/human_simulate_pair --outSAMtype BAM Unsorted 

samtools sort -@ 3 -o ./result/human_simulate_pair_sorted.bam ./result/human_simulate_pairAligned.out.bam

samtools index ./result/human_simulate_pair_sorted.bam

java -jar ./code/run_sh/isoLarge.jar anno/hg19_gencode_from_ucsc_nothick_nocds_chr21.bed ./result/human_simulate_pair_sorted.bam \
./result/iso_simulation_pair/human_simulate_pair_iso_large_unique_intron.tsv 
#20 false true a3a5_bias_all_pombe.tsv


######################################################################################################################################
