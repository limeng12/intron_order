#simulate RNA-seq read given intron splicing order.
library(BSgenome);
library(BSgenome.Hsapiens.UCSC.hg19);
library(polyester);
library(Biostrings)
library(Rsamtools)
library(stringr)
library(stringi)
library(dplyr)

## incorrected predict transcript usually because the pre-mRNA is long, and read doesn't cover or can't cover this length. 
## 

setwd("/Users/mengli/Documents/projects/iso/");
source("code/generate_fragments_me.R");

number_of_frag_per_premrna<-200;
number_of_frag_per_premrna_sd<-100;
max_number_premrna<-5000;
t_max_frag_len<-30000;

#max_frag_len=t_max_frag_len



chr21_fa<-FaFile("./result/chr21/chr21.fa");

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
  
  a<-as.character( Rsamtools::getSeq(chr21_fa, gr) );
  a
  #as.character(a);
  
}

#xx<-"I:1:100:+";
#xx<-"chr21:1:10:+";
#get_seq_region(xx);


concat_one_premrna<-function(t_exons,t_introns){
  pre<-""
  for( j in 1:length(t_exons)){
    pre<- str_c(pre,get_seq_region ( t_exons[j] ), get_seq_region(t_introns[j]) );
  }
  pre
}


trans_gene_id_map<-read.table("anno/hg19_ensembl_gene_id_trans_id_map.tsv",as.is=TRUE,sep="\t")[,1:2]
colnames(trans_gene_id_map)<-c("gene_id","trans_id")
#for each gene give a symbol, if a 
gene_id=unique(trans_gene_id_map[,"gene_id"])
symbol=sample( c("i","n"),length(unique(trans_gene_id_map[,"gene_id"] ) ) ,replace = TRUE)
               
gene_symbol_tbl<-data.frame(gene_id=gene_id,symbol=symbol );

trans_gene_id_map<-inner_join(trans_gene_id_map,gene_symbol_tbl,"gene_id"="gene_id") 


### get the sequence of introns and exons for each transcript
ucsc_annotation_bed<-read.table("anno/hg19_gencode_from_ucsc_nothick_nocds.bed",
                                as.is = TRUE,sep = "\t");
colnames(ucsc_annotation_bed)<-c("chr","start","end","trans_id","score","strand","thick_start","thick_end",
"a_","exon_count","exon_len","exon_start");

ucsc_annotation_bed[,"trans_id"]<-sapply(str_split(ucsc_annotation_bed[,"trans_id"],"\\." ),"[",1 )                                
    
ucsc_annotation_bed<-inner_join(ucsc_annotation_bed,trans_gene_id_map,"trans_id"="trans_id");            

ucsc_annotation_bed<-ucsc_annotation_bed[str_detect(ucsc_annotation_bed$chr,"chr21"),];

print(paste0("Total unique gene number: ", length(unique(ucsc_annotation_bed$gene_id) ) ) );
#579 genes, 1415 transcripts

#keep one transcript per gene

ucsc_annotation_bed<-ucsc_annotation_bed[order(ucsc_annotation_bed$exon_count,decreasing = TRUE),];
ucsc_annotation_bed<-ucsc_annotation_bed[!duplicated(ucsc_annotation_bed[,"gene_id"]),]


trans_exons_introns<-list();

for (row_num in 1:nrow(ucsc_annotation_bed)) {
  trans_id<-ucsc_annotation_bed[row_num,"trans_id"]
  
  trans_start<-ucsc_annotation_bed[row_num,"start"];
  
  trans_end<- ucsc_annotation_bed[row_num,"end"];
  
  trans_id<-ucsc_annotation_bed[row_num,"trans_id"];
  
  chr<-ucsc_annotation_bed[row_num,"chr"];
  
  strand<-ucsc_annotation_bed[row_num,"strand"];
  
  exon_len<-as.numeric(as.character(unlist(strsplit(ucsc_annotation_bed[row_num,"exon_len"], split = ","))) );
  
  if(chr!="chr21"){
    next;
  }
  
  if(length(exon_len)<=2){
    next;
  }
  
  exon_pos<-as.numeric(as.character(unlist(strsplit(ucsc_annotation_bed[row_num,"exon_start"], split = ","))) );
  
  
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
  
  symbol<-ucsc_annotation_bed[row_num,"symbol"];

  t_order<-c()
  if(symbol=="i"){
  	t_order<-1:length(introns)
  }else{
    t_order<-length(introns):1
  }
  
  t_order<-sample(t_order);
  
  trans_exons_introns[[trans_id]]<-list(exons=exons,
                                              introns=introns,
                                              trans_id=trans_id,
                                              orders=t_order ) ;
}
##

rm(ucsc_annotation_bed)
print(paste0("Total number of multi-introns containing transcripts: ", length(trans_exons_introns) ) );


# bb<-concat_one_premrna(exons,introns)
# dd<-get_seq_region("chr6:348842-350821")
#given_order<-c(3,1,4,7,5,6,2,4);



#### simulated reads
#t<-1;
unlink(paste0("result/human_read_given_order_simulation_long.fasta.gz") );

unlink(paste0("result/human_read_given_order_simulation_long_super_long.fasta.gz") );

unlink(paste0("result/human_read_given_order_simulation_pair_1.fasta.gz") );

unlink(paste0("result/human_read_given_order_simulation_pair_2.fasta.gz") );

#  trans_one<-trans_exons_introns[[2]]

trans_order<-list();

for(trans_one_index in 1:length(trans_exons_introns) ){
  
  trans_one<-trans_exons_introns[[trans_one_index]]
  ##build pre-mRNA sequence set for each transcript
  pre_mrna_set<-c();
  
  tt_introns<-trans_one$introns;
  
  tt_exons<-trans_one$exons;
  
  tt_trans_id<-trans_one$trans_id;
    
  #given_order<-sample(trans_one$orders);
  print(paste0(trans_one_index,":",tt_trans_id,":",length(tt_exons)) )
  
  given_order<-(trans_one$orders);
  
  #trans_order<-rbind(tt_trans_id,given_order);
  
  trans_order[[tt_trans_id]]<-given_order;
  
  for(i in 1:length(given_order) ){
    #tt_introns<-introns;
    tt_introns[given_order[i]]<-""
    #for(j in 1:i){
    #  tt_introns[j]<-"";
    #}
    pre_mrna_set<-c(pre_mrna_set,concat_one_premrna(tt_exons, tt_introns)  );
    
  }
  
  names(pre_mrna_set)<-str_c("XXXX",stri_rand_strings(length(pre_mrna_set), 100, pattern = "[A-Za-z0-9]"));
  
  for(i in 1:length(pre_mrna_set)){
    
    ######################################################SUPER LONG READ##########################################
    pre_mrna_set_dna<-DNAStringSet(rep(pre_mrna_set[i],
                                       min(rpois(1, number_of_frag_per_premrna),max_number_premrna )  ) );
    
    #####single end super long 3rd generation sequencing reads, SD for third should be large
    fragments<-generate_fragments_me(pre_mrna_set_dna, fraglen = 15000, 
                                     fragsd = 5000, distr = "normal", readlen = 12000,
                                     max_frag_len=t_max_frag_len)
    #  plot(hist(rnorm(1000,15000,5000)))
    #  plot(hist(rnorm(1000,12000,4000)))
    #paired_reads<-get_reads(fragments,250,paired = TRUE);
    
    t_read_len<-max( (rnorm(1,12000,4000) ), 2000);
    single_reads<-get_reads(fragments,t_read_len,paired = FALSE);
    
    writeXStringSet(single_reads, 
                    filepath = paste0("result/human_read_given_order_simulation_long_super_long.fasta.gz"), 
                    format = "fasta",append=TRUE,compress = TRUE);
    
    ##########################################################LONG READ###############################################
    rm(pre_mrna_set_dna,fragments,single_reads)
    
    pre_mrna_set_dna<-DNAStringSet(rep(pre_mrna_set[i], 
                                       min(rpois(1, number_of_frag_per_premrna),max_number_premrna )*12000/800  ) );
    
    #####single end 3rd generation sequencing reads
    fragments<-generate_fragments_me(pre_mrna_set_dna, fraglen = 15000, 
                                     fragsd = 5000, distr = "normal", readlen = 800,max_frag_len=t_max_frag_len)
   
    #paired_reads<-get_reads(fragments,250,paired = TRUE);
    
    t_read_len<-max((rnorm(1,800,250) ), 200);
    single_reads<-get_reads(fragments,t_read_len,paired = FALSE);
    
    writeXStringSet(single_reads, filepath = paste0("result/human_read_given_order_simulation_long.fasta.gz"), 
                    format = "fasta",append=TRUE,compress = TRUE);
    
    
    ############################################################SHORT PAIR END################################
    rm(pre_mrna_set_dna,fragments,single_reads)
    
    pre_mrna_set_dna<-DNAStringSet(rep(pre_mrna_set[i], 
                                       min(rpois(1, number_of_frag_per_premrna), max_number_premrna )*12000/300 ) );
    
    #####pair end reads
    fragments<-generate_fragments_me(pre_mrna_set_dna, fraglen = 15000,
                                     distr = "normal", fragsd = 5000, readlen = 300,max_frag_len=t_max_frag_len);
    
    pair_reads<-get_reads(fragments,150,paired = TRUE);
    
    
    left_reads<-pair_reads[1:length(pair_reads) %% 2==1];
    
    right_reads<-pair_reads[1:length(pair_reads) %% 2==0];
    
    
    writeXStringSet(left_reads, filepath = paste0("result/human_read_given_order_simulation_pair_1.fasta.gz"), 
                    format = "fasta",append=TRUE,compress = TRUE);
    
    writeXStringSet(right_reads, filepath = paste0("result/human_read_given_order_simulation_pair_2.fasta.gz"), 
                    format = "fasta",append=TRUE,compress = TRUE);
  
    
    rm(pre_mrna_set_dna,fragments,pair_reads,left_reads,right_reads)
  }
}


save(trans_order,file="result/simulation_orders_human.Rd");

## one nanopore ~1-3 million reads
## 5,557,643 fragments are simulated
#  fragments_one<-fragments[1:2]
#  pair_reads<-get_reads(fragments_one,rnorm(1,800,15),paired = TRUE);


  