library(shiny);
library(networkD3);
library(Sushi);
library(readr);
library(DT);
library(dplyr);
library(igraph);
library(dbscan);
library(stringr);
library(gtools);
library(ggraph)
options(scipen=999);

source("code/utils.R");
source("code/get_intron_from_bed.R");
source("code/build_iso_object2.R",echo=TRUE);
source("code/get_iso_summary.R",echo=TRUE);

source("code/draw_3d.R");
source("code/mlp3.R");
source("code/get_adj2.R");
source("code/get_members.R");
source("code/cal_mlp_graph.R");
source("code/add_gene_symbol_intron_pos.R");



get_mlo_pipe<-function(t_bed_anno,
                       files_all,
                       t_most_likeli_order_output_path,
                       gene_trans_id_tbl,
                       read_count_threshold=0,
                       t_trans_exp_file=""){
  
  ##extracte intron postion from bed file
  intron_pos_mat_fr<-get_intron_from_bed(t_bed_anno);
  
  ##build intron splicing order pairs
  iso_final<-build_iso_object2(files_all,intron_pos_mat_fr,trans_exp_file=t_trans_exp_file);
  
  ## summary intron splicing order pairs
  iso_summary<-get_iso_summary(iso_final,intron_pos_mat_fr);
  
  ######################################build intron splicing matrix, graph and most likely order#########################
  
  
  t_alpha<-0.1;
  
  ##calculated intron splicing order adjacent matrix
  t_igraph_list<-get_adj2(iso_final,iso_summary,0.95,read_count_threshold);
  
  ##caculated intron splicing unit
  t_igraph_list<-get_members(t_igraph_list,t_alpha);
  
  
  gene_trans_id_map<-read.table(gene_trans_id_tbl,header = FALSE,as.is = TRUE,sep = "\t")[,1:3];
  colnames(gene_trans_id_map)<-c("gene_id","trans_id","gene_symbol");
  
  ##calculate most likely order
  t_igraph_list<-cal_mlp(t_igraph_list, t_most_likeli_order_output_path,
                         t_alpha,gene_trans_id_map);
  
  ##
  t_igraph_list
}


