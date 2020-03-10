#simulate RNA-seq read given intron splicing order.
library(BSgenome);
library(BSgenome.Hsapiens.UCSC.hg19);
library(polyester);
library(Biostrings)
library(Rsamtools)
library(stringr)
library(stringi)

## incorrected predict transcript usually because the pre-mRNA is long, and read doesn't cover or can't cover this length. 
##

setwd("/Users/mengli/Documents/projects/iso/pombe");
#cat Schizosaccharomyces_pombe.ASM294v2.dna.chromosome.I.fa >> Schizosaccharomyces_pombe.ASM294v2.fa
#cat Schizosaccharomyces_pombe.ASM294v2.dna.chromosome.II.fa >> Schizosaccharomyces_pombe.ASM294v2.fa
#cat Schizosaccharomyces_pombe.ASM294v2.dna.chromosome.III.fa >> Schizosaccharomyces_pombe.ASM294v2.fa
#cat Schizosaccharomyces_pombe.ASM294v2.dna.chromosome.MT.fa >> Schizosaccharomyces_pombe.ASM294v2.fa
#samtools faidx Schizosaccharomyces_pombe.ASM294v2.fa
source("../code/generate_fragments_me.R")


number_of_frag_per_premrna<-50;
number_of_frag_per_premrna_var<-20

pombe_fa<-FaFile("fasta/Schizosaccharomyces_pombe.ASM294v2.fa");

get_seq_region<-function(x){
  if(x==""){
    return("")
  }
  
  chr<-str_c(sapply(str_split(x,":"),"[",1) );
  if(!str_detect(chr,"chr")){
    chr<-str_c("chr",chr);
    
    #chr<-str_sub(chr,4);
  }
  
  start<-as.numeric( sapply(str_split(x,":"),"[",2) );
  
  end<-as.numeric( sapply(str_split(x,":"),"[",3) );
  
  strand<-sapply(str_split(x,":"),"[",4) 
  
  gr <- GRanges(
    seqnames = c(chr),
    ranges = IRanges( start, end ),
    strand = strand
  );
  
  #a<-getSeq(BSgenome.Hsapiens.UCSC.hg19,chr,start,end,strand="+");
  
  a<-as.character( Rsamtools::getSeq(pombe_fa, gr) );
  a
  #as.character(a);
  
}

#xx<-"I:1:100:+";
#xx<-"chrI:1:10:+";
#get_seq_region(xx);


concat_one_premrna<-function(t_exons,t_introns){
  pre<-""
  for( j in 1:length(t_exons)){
    pre<- str_c(pre,get_seq_region ( t_exons[j] ), get_seq_region(t_introns[j]) );
  }
  pre
}


### get the sequence of introns and exons for each transcript
ucsc_annotation_bed<-read.table("Schizosaccharomyces_pombe.ASM294v2.43.chr_nothick.bed",
                                as.is = TRUE,sep = "\t")
trans_exons_introns<-list();

for (row_num in 1:nrow(ucsc_annotation_bed)) {
  trans_id<-ucsc_annotation_bed[row_num,4]
  
  trans_start<-ucsc_annotation_bed[row_num,2];
  
  trans_end<-ucsc_annotation_bed[row_num,3];
  
  trans_id<-ucsc_annotation_bed[row_num,4];
  
  chr<-ucsc_annotation_bed[row_num,1];
  
  strand<-ucsc_annotation_bed[row_num,6];
  
  exon_len<-as.numeric(as.character(unlist(strsplit(ucsc_annotation_bed[row_num,11], split = ","))) );
  
  if(length(exon_len)<=2){
    next;
  }
  
  exon_pos<-as.numeric(as.character(unlist(strsplit(ucsc_annotation_bed[row_num,12], split = ","))) );
  
  
  exons_start<-trans_start+exon_pos+1;
  exons_end<-exons_start+exon_len-1;
  
  
  introns_start<-(exons_end+1)[-1*length(exon_len)];
  introns_end<-(exons_start-1)[-1];
  
  
  exons<-str_c(chr,":",exons_start,":",exons_end,":",strand);
  introns<-str_c(chr,":",introns_start,":",introns_end,":",strand);
  
  if(strand=="-"){
    exons<-rev(exons);
    introns<-rev(introns);
  }
  
  introns<-c(introns,"");
  
  
  t_given_order<-sample(length(introns):1);

  
  trans_exons_introns[[trans_id]]<-list(exons=exons,
                                              introns=introns,
                                              trans_id=trans_id,
                                              orders=t_given_order ) ;
}
##

print(paste0("total number of multi-introns containing transcripts: ", length(trans_exons_introns) ) );


# bb<-concat_one_premrna(exons,introns)
# dd<-get_seq_region("chr6:348842-350821")
#given_order<-c(3,1,4,7,5,6,2,4);



#### simulated reads
#t<-1;
unlink(paste0("result/pombe_read_given_order_simulation_super_long.fasta.gz") );

unlink(paste0("result/pombe_read_given_order_simulation_long.fasta.gz") );

unlink(paste0("result/pombe_read_given_order_simulation_pair_1.fasta.gz") );

unlink(paste0("result/pombe_read_given_order_simulation_pair_2.fasta.gz") );

#  trans_one<-trans_exons_introns[[2]]

trans_order<-list();

for(trans_one in trans_exons_introns){
  
  ##build pre-mRNA sequence set for each transcript
  pre_mrna_set<-c();
  
  tt_introns<-trans_one$introns;
  
  tt_exons<-trans_one$exons;
  
  tt_trans_id<-trans_one$trans_id;
    
  given_order<-(trans_one$orders);
  
  #trans_order<-rbind(tt_trans_id,given_order);
  
  trans_order[[tt_trans_id]]<-given_order;
  
  for(i in 1:length(given_order) ){
    #tt_introns<-introns;
    tt_introns[given_order[i]]<-""
    #for(j in 1:i){
    #  tt_introns[j]<-"";
    #}
    pre_mrna_set<-c(pre_mrna_set,concat_one_premrna(tt_exons,tt_introns)  );
    
  }
  
  names(pre_mrna_set)<-str_c("XXXX",stri_rand_strings(length(pre_mrna_set), 100, pattern = "[A-Za-z0-9]"));
    
  for(i in 1:length(pre_mrna_set)){
    
    ############################################################
    pre_mrna_set_dna<-DNAStringSet(rep(pre_mrna_set[i],
                                       max(rpois(1, number_of_frag_per_premrna)  ,1) ));
    
    #####single end super long 3rd generation sequencing reads
    fragments_super_long<-generate_fragments_me(pre_mrna_set_dna, fraglen = 6000, fragsd = 2000, distr = "normal", readlen = 5000)
    
    #paired_reads<-get_reads(fragments,250,paired = TRUE);
    
    t_read_len<-max((rnorm(1,5000,1600)),1000);
    single_reads<-get_reads(fragments_super_long,t_read_len,paired = FALSE);
    
    writeXStringSet(single_reads, filepath = paste0("result/pombe_read_given_order_simulation_super_long.fasta.gz"), 
                    format = "fasta",append=TRUE,compress = TRUE);
    
    ############################################################
    pre_mrna_set_dna<-DNAStringSet(rep(pre_mrna_set[i],
                                       max(rpois(1, number_of_frag_per_premrna*3)  ,1) ));
    
    #####single end 3rd generation sequencing reads
    fragments_long<-generate_fragments_me(pre_mrna_set_dna, fraglen = 6000, fragsd = 2000, distr = "normal", readlen = 800)
   
    #paired_reads<-get_reads(fragments,250,paired = TRUE);
    
    t_read_len<-max((rnorm(1,800,250) ), 150);
    single_reads<-get_reads(fragments_long,t_read_len,paired = FALSE);
    
    writeXStringSet(single_reads, filepath = paste0("result/pombe_read_given_order_simulation_long.fasta.gz"), 
                    format = "fasta",append=TRUE,compress = TRUE);
    
    ############################################################
    #####pair end reads
    #pre_mrna_set_dna<-DNAStringSet(rep(pre_mrna_set[i],
    #                                   rpois(1, number_of_frag_per_premrna*median(width(fragments_super_long))/300)  ) );
    
    pre_mrna_set_dna<-DNAStringSet(rep(pre_mrna_set[i],
                                       max(rpois(1, number_of_frag_per_premrna)  ,1) ));
    
    
    fragments<-generate_fragments_me(pre_mrna_set_dna, fraglen = 600, distr = "normal", fragsd = 200, readlen = 300)
    
    pair_reads<-get_reads(fragments,150,paired = TRUE);
    
    
    left_reads<-pair_reads[1:length(pair_reads) %% 2==1];
    
    right_reads<-pair_reads[1:length(pair_reads) %% 2==0];
    
    
    writeXStringSet(left_reads, filepath = paste0("result/pombe_read_given_order_simulation_pair_1.fasta.gz"), 
                    format = "fasta",append=TRUE,compress = TRUE);
    
    writeXStringSet(right_reads, filepath = paste0("result/pombe_read_given_order_simulation_pair_2.fasta.gz"), 
                    format = "fasta",append=TRUE,compress = TRUE);
  
  }
  
}


save(trans_order,file="result/simulation_orders.Rd");
## one nanopore ~1-3 million reads
## 5,557,643 fragments are simulated
#  fragments_one<-fragments[1:2]
#  pair_reads<-get_reads(fragments_one,rnorm(1,800,15),paired = TRUE);


  