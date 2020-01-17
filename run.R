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

#setwd("some_where_in_your_computer/intron_order");
setwd("/Users/mengli/Documents/projects/iso/intron_order");

#this.dir <- dirname(parent.frame(1)$ofile);


##############################prepare pairwise intron orders################################################################
source("code/build_iso_object.R",echo=TRUE)

#files_all<-list.files("data/iso_3rd/",full.names =TRUE,pattern = "*unique_intron.tsv");
files_all<-list.files("data/iso_3rd/",full.names =TRUE);

label<-"human"

is_large=TRUE;

gene_trans_id_tbl<-"./anno/gene_id_trans_id.tsv";

ucsc_intron_anno<-"./anno/hg19_gencode_intron_from_ucsc.bed";

if_uniq="uniq_intron";
if(is_large){
  if_uniq="not_uniq_intron";
}

t_result_path<-paste0("result/all_iso_data_",label,"_",if_uniq,".Rd") ;

build_iso_object(files_all,gene_trans_id_tbl,ucsc_intron_anno,is_large,t_result_path);



######################################build matrix graph MLO################################################################

isoform_num_produce<-100
read_count_threshold<-0
draw_and_save_graph<-FALSE # not used 
return_graph<-TRUE;
t_alpha<-0.1;

load(t_result_path);


source("code/draw.R");
#source("code/mlp.R");
source("code/draw_3d.R");
source("code/mlp3.R");
source("code/get_adj.R");
source("code/get_members.R");
source("code/cal_mlp_graph.R");

iso_slow_sumary<-iso_slow_sumary[order((iso_slow_sumary$edge_count)/
                                         (iso_slow_sumary[,"int_count"]^2+0.1),decreasing = TRUE),];

t_igraph_list<-get_adj("result/iso_test_unique_order_by_count_graph_output.pdf",
                       iso_final,iso_slow_sumary,
                       isoform_num_produce,read_count_threshold,
                       draw_and_save_graph,return_graph);

t_igraph_list<-get_members(t_igraph_list,t_alpha);

t_igraph_list<-cal_mlp(t_igraph_list,"./result/best_order.tsv",t_alpha,read_count_threshold);

t_igraph_list<-draw_3d(t_igraph_list,"./result/html/",t_alpha,FALSE);

save( t_igraph_list, file="result/t_igraph_list.Rd",version = 2);


######################################################shiny#################################################################

load("result/t_igraph_list.Rd");

load("anno/gene_trans_id_map.Rd");

load("anno/sushi_trans_file.Rd");


source("code/shiny_web.R",echo = TRUE);

