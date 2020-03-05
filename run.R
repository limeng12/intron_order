library(rstudioapi)
this.dir<-dirname(rstudioapi::getActiveDocumentContext()$path) 
setwd(this.dir);


files_all<-list.files("data/iso_3rd/",full.names =TRUE,pattern = "*tsv");

gene_trans_id_tbl<-"./anno/hg19_ensembl_gene_id_trans_id_map.tsv";

bed_anno<-"anno/hg19_gencode_from_ucsc_nothick_nocds.bed"

most_likeli_order_output_path<-"./result/best_order.tsv"



read_count_threshold<-0;

source("code/necessary_code.R");


#################################################################################################################

##extracte intron postion from bed file
intron_pos_mat_fr<-get_intron_from_bed(bed_anno);

##build intron splicing order pairs
iso_final<-build_iso_object2(files_all,intron_pos_mat_fr,read_count_threshold=read_count_threshold);

## summary intron splicing order pairs
iso_summary<-get_iso_summary(iso_final,intron_pos_mat_fr);


######################################build intron splicing matrix, graph and most likely order#########################

t_alpha<-0.1;

##calculated intron splicing order adjacent matrix
t_igraph_list<-get_adj2(iso_final,iso_summary,0.95);

##caculated intron splicing unit
t_igraph_list<-get_members(t_igraph_list,t_alpha);


gene_trans_id_map<-read.table(gene_trans_id_tbl,header = FALSE,as.is = TRUE,sep = "\t")[,1:3];
colnames(gene_trans_id_map)<-c("gene_id","trans_id","gene_symbol");

##calculate most likely order
t_igraph_list<-cal_mlp(t_igraph_list, most_likeli_order_output_path,
                       t_alpha, read_count_threshold,gene_trans_id_map);

##
##t_igraph_list<-add_gene_symbol_intron_pos(t_igraph_list,gene_trans_id_map,intron_pos_mat_fr);


#save( t_igraph_list, file="result/t_igraph_list.Rd",version = 2);



