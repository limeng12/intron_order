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

options(scipen=999);

####set directory ####
if (!requireNamespace("rstudioapi", quietly = TRUE) )
  install.packages("rstudioapi")

library(rstudioapi)
this.dir<-dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(this.dir);


####key parameter####
## number of isoform output, number of isoform wish to check, 
## I sorted the output by read count support, so if produce too many isoforms, the last few may not accurate
isoform_num_produce<-100


                                      ################below are scripts################
#################################################################################################################
##############################prepare pairwise intron orders###########################################################
source("code/build_iso_object.R",echo=TRUE)

files_all<-list.files("data/iso_3rd/",full.names =TRUE);

gene_trans_id_tbl<-"./anno/hg19_ensembl_gene_id_trans_id_map.tsv";

ucsc_intron_anno<-"./anno/hg19_gencode_intron_from_ucsc.bed";

label<-"human"

##deprecated, used when consider each intron one time 
is_large=TRUE;

if_uniq="uniq_intron";
if(is_large){
  if_uniq="not_uniq_intron";
}

###tmp file to store intron splicing order pairs
t_result_path<-paste0("result/all_iso_data_",label,"_",if_uniq,".Rd") ;


build_iso_object(files_all,gene_trans_id_tbl,ucsc_intron_anno,is_large,t_result_path, TRUE);
####If the transcript id contain . in nature, then may following script instead
##build_iso_object(files_all,gene_trans_id_tbl,ucsc_intron_anno,is_large,t_result_path, TRUE);


######################################build intron splicing matrix, graph and most likely order#########################

## the read count threahold
read_count_threshold<-0

## suppose that the read count support intron 1 spliced before intron 2 is 10
## but intron 2 spliced before intron 1 is 0, then add t_alpha on both side to void 0 in likelihood multiplication
t_alpha<-0.1;

draw_and_save_graph<-FALSE;
return_graph<-TRUE;

load(t_result_path);


source("code/draw.R");
source("code/draw_3d.R");
source("code/mlp3.R");
source("code/get_adj.R");
source("code/get_members.R");
source("code/cal_mlp_graph.R");

iso_slow_sumary<-iso_slow_sumary[order((iso_slow_sumary$edge_count)/
                                         (iso_slow_sumary[,"int_count"]^2+0.1),decreasing = TRUE),];

##calculated intron splicing order adjacent matrix
t_igraph_list<-get_adj("result/iso_test_unique_order_by_count_graph_output.pdf",
                       iso_final,iso_slow_sumary,
                       isoform_num_produce,read_count_threshold,
                       draw_and_save_graph,return_graph);

##caculated intron splicing unit
t_igraph_list<-get_members(t_igraph_list,t_alpha);
##most likely order
t_igraph_list<-cal_mlp(t_igraph_list,"./result/best_order.tsv",t_alpha,read_count_threshold);

##draw intron splicing order graph
t_igraph_list<-draw_3d(t_igraph_list, paste0(getwd(),"./result/html/"),t_alpha,TRUE);

##output intron splicing order adjacent matrix
for( i in 1:length(t_igraph_list) ){

  write.table(t_igraph_list[[i]]$adjacency_matrix, file=str_c("./result/adj_matrix/",names(t_igraph_list)[i],".tsv" ),
            sep="\t",col.names = FALSE,row.names = FALSE  );
}


save( t_igraph_list, file="result/t_igraph_list.Rd",version = 2);


######################################################shiny###########################################################

load("result/t_igraph_list.Rd");

gene_trans_id_map<-read.table("anno/hg19_ensembl_gene_id_trans_id_map.tsv",
                              header = FALSE,as.is = TRUE,sep = "\t");
colnames(gene_trans_id_map)<-c("gene_id","trans_id","gene_symbol","trans_start","trans_end","strand",
                               "chr","gene_start","gene_end");


load("anno/sushi_trans_file.Rd");

source("code/shiny_web.R",echo = TRUE);

