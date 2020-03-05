library(rstudioapi)
this.dir<-dirname(rstudioapi::getActiveDocumentContext()$path) 
setwd(this.dir);


files_all<-list.files("data/",full.names =TRUE,pattern = "*tsv");

gene_trans_id_tbl<-"./anno/hg19_ensembl_gene_id_trans_id_map.tsv";

bed_anno<-"anno/hg19_gencode_from_ucsc_nothick_nocds.bed"

most_likeli_order_output_path<-"./result/best_order.tsv"


read_count_threshold<-0;

source("code/necessary_code.R");

#################################################################################################################

t_igraph_list<-get_mlo_pipe(bed_anno,files_all,most_likeli_order_output_path,
                            gene_trans_id_tbl,read_count_threshold);


#save( t_igraph_list, file="result/t_igraph_list.Rd",version = 2);

