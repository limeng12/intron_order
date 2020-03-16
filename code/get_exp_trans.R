library(stringr)

setwd("/Users/mengli/Documents/projects/iso");

trans_exp1<-read.table("anno/k562_hg19_ENCFF360MVV.tsv",header = TRUE,as.is = TRUE,sep="\t")

trans_exp2<-read.table("anno/k562_hg19_ENCFF928EIW.tsv",header = TRUE,as.is = TRUE,sep="\t")


trans_exp1<-trans_exp1[trans_exp1$TPM>1,];
trans_exp2<-trans_exp2[trans_exp2$TPM>1,];


trans_exp1$transcript_id<-sapply(str_split(trans_exp1$transcript_id,"\\."),"[",1 );
trans_exp2$transcript_id<-sapply(str_split(trans_exp2$transcript_id,"\\."),"[",1 );

cat( unique(c(trans_exp1$transcript_id,trans_exp2$transcript_id)), 
    file="anno/k562_exp_trans_id.tsv",sep="\n");


