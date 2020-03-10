
source ~/.profile

#######################################################simulation########################################
minimap2 -t 4 -ax splice fasta/Schizosaccharomyces_pombe.ASM294v2.fa ~/Documents/projects/iso/pombe/result/pombe_read_given_order_simulation_super_long.fasta.gz |
samtools view -@ 4 -S -b | samtools sort -@ 4 -o  ./result/iso_simulation_super_long.bam

samtools index ./result/iso_simulation_super_long.bam

java -jar ../code/run_sh/isoLarge.jar Schizosaccharomyces_pombe.ASM294v2.43.chr_nothick.bed  ./result/iso_simulation_super_long.bam ./iso_simulation_super_long/pombe_simulate_super_long_iso_large_unique_intron.tsv 


minimap2 -t 4 -ax splice fasta/Schizosaccharomyces_pombe.ASM294v2.fa ~/Documents/projects/iso/pombe/result/pombe_read_given_order_simulation_long.fasta.gz |
samtools view -@ 4 -S -b | samtools sort -@ 4 -o  ./result/pombe_simulate_long.bam

samtools index ./result/pombe_simulate_long.bam

java -jar ../code/run_sh/isoLarge.jar Schizosaccharomyces_pombe.ASM294v2.43.chr_nothick.bed  ./result/pombe_simulate_long.bam ./iso_simulation_long/pombe_simulate_long_iso_large_unique_intron.tsv 


rm -r ./result/pombe_simulate_pair*_STARtmp

STAR --runThreadN 4 --quantMode GeneCounts \
--genomeDir ./star_index \
--readFilesIn  ./result/pombe_read_given_order_simulation_pair_1.fasta.gz  ./result/pombe_read_given_order_simulation_pair_2.fasta.gz \
--chimJunctionOverhangMin 20 --chimSegmentMin 20 --chimSegmentReadGapMax 3 \
--outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 \
--outFilterMismatchNoverLmax 0.04 --alignIntronMin 10 --readFilesCommand 'gunzip -c' \
--alignIntronMax 1000000 --alignMatesGapMax 1000000 \
--outSAMattributes NH HI NM MD AS nM jM jI XS --outFileNamePrefix \
./result/pombe_simulate_pair --outSAMtype BAM Unsorted 

samtools sort -@ 3 -o ./result/pombe_simulate_pair_sorted.bam ./result/pombe_simulate_pairAligned.out.bam
#./result/pombe_simulate_pair --outSAMtype BAM SortedByCoordinate 


samtools index ./result/pombe_simulate_pair_sorted.bam

java -jar ../code/run_sh/isoLarge.jar Schizosaccharomyces_pombe.ASM294v2.43.chr_nothick.bed ./result/pombe_simulate_pair_sorted.bam ./iso_simulation_pair/pombe_simulate_pair_iso_large_unique_intron.tsv 


###########################################################################################################################################################
