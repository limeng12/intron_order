setwd("/Users/mengli/Documents/projects/iso");
library(stringr)
library(ggplot2)

get_length_diff<-function(region1,region2){
  as.numeric(sapply(strsplit(region1,":|-"),"[",2) )->region1_start;
  as.numeric(sapply(strsplit(region1,":|-"),"[",3) )->region1_end;
  
  as.numeric(sapply(strsplit(region2,":|-"),"[",2) )->region2_start;
  as.numeric(sapply(strsplit(region2,":|-"),"[",3) )->region2_end;
  
  len_diff<-abs(region1_end-region1_start-(region2_end-region2_start)  );
  
  len_diff;
  #region1_start<region2_start;
}


se1<-read.table(paste0("anno/fromGTF.SE.txt"),header = TRUE);
se1$smallES<-pmin(se1$upstreamES+1,se1$downstreamES+1);
se1$smallEE<-pmin(se1$upstreamEE,se1$downstreamEE);

se1$largeES<-pmax(se1$upstreamES+1,se1$downstreamES+1);
se1$largeEE<-pmax(se1$upstreamEE,se1$downstreamEE);


up_int<-str_c(se1$chr,":",se1$smallEE+1,"-",se1$exonStart_0base+1-1);

down_int<-str_c(se1$chr,":",se1$exonEnd+1,"-",se1$largeES-1);

all_alt_se_introns<-c(up_int,down_int);



se1<-read.table(paste0("anno/fromGTF.A5SS.txt"),header = TRUE);

short_int<-1:nrow(se1);

#up_int<-str_c(se1$chr,se1$smallEE+1,"-",exonStart_0base);

short_int[se1$strand=="-"]<-str_c(se1[se1$strand=="-",]$chr,":",se1[se1$strand=="-",]$flankingEE+1,"-",se1[se1$strand=="-",]$longExonStart_0base);

short_int[se1$strand=="+"]<-str_c(se1[se1$strand=="+",]$chr,":",se1[se1$strand=="+",]$longExonEnd+1,"-",se1[se1$strand=="+",]$flankingES+1-1);


long_int<-1:nrow(se1);

#up_int<-str_c(se1$chr,se1$smallEE+1,"-",exonStart_0base);

long_int[se1$strand=="-"]<-str_c(se1[se1$strand=="-",]$chr,":",se1[se1$strand=="-",]$flankingEE+1,"-",se1[se1$strand=="-",]$shortES);

long_int[se1$strand=="+"]<-str_c(se1[se1$strand=="+",]$chr,":",se1[se1$strand=="+",]$shortEE+1,"-",se1[se1$strand=="+",]$flankingES+1-1);



all_alt_a5_introns<-c(long_int,short_int);

len_diff_a5<-get_length_diff(long_int,short_int);

p_len_diff_a5<-ggplot()+geom_histogram(aes(x=len_diff_a5) )#+scale_x_log10();




se1<-read.table(paste0("anno/fromGTF.A3SS.txt"),header = TRUE);

short_int<-1:nrow(se1);

#up_int<-str_c(se1$chr,se1$smallEE+1,"-",exonStart_0base);

short_int[se1$strand=="+"]<-str_c(se1[se1$strand=="+",]$chr,":",se1[se1$strand=="+",]$flankingEE+1,"-",se1[se1$strand=="+",]$longExonStart_0base+1-1);

short_int[se1$strand=="-"]<-str_c(se1[se1$strand=="-",]$chr,":",se1[se1$strand=="-",]$longExonEnd+1,"-",se1[se1$strand=="-",]$flankingES+1-1);


long_int<-1:nrow(se1);

long_int[se1$strand=="+"]<-str_c(se1[se1$strand=="+",]$chr,":",se1[se1$strand=="+",]$flankingEE+1,"-",se1[se1$strand=="+",]$shortES+1-1);

long_int[se1$strand=="-"]<-str_c(se1[se1$strand=="-",]$chr,":",se1[se1$strand=="-",]$shortEE+1,"-",se1[se1$strand=="-",]$flankingES+1-1);


#as_region_neg<-(se1[se1$strand=="-",]$flankingEE+1-(se1[se1$strand=="-",]$shortES+1-1) )

all_alt_a3_introns<-c(long_int,short_int);


len_diff_a3<-get_length_diff(long_int,short_int);


se1<-read.table(paste0("anno/fromGTF.RI.txt"),header = TRUE);


#short_int<-1:nrow(se1);
se1$smallES<-pmin(se1$upstreamES+1,se1$downstreamES+1);
se1$smallEE<-pmin(se1$upstreamEE,se1$downstreamEE);

se1$largeES<-pmax(se1$upstreamES+1,se1$downstreamES+1);
se1$largeEE<-pmax(se1$upstreamEE,se1$downstreamEE);

short_int<-str_c(se1$chr,":",se1$smallEE+1,"-",se1$largeES-1);


all_alt_ri_introns<-short_int;






p_len_diff_a3<-ggplot()+geom_histogram(aes(x=len_diff_a3) )#+scale_x_log10();

source("code/multiplot.R")
pdf("result/A5SS_A3SS_alternative_region_length.pdf",width = 8,height = 4);
#print(p_len_diff_a3);
#print(p_len_diff_a5);
multiplot(p_len_diff_a5,p_len_diff_a3,cols = 2)

dev.off();



