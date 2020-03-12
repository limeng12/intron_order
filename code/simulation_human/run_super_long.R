setwd("/Users/mengli/Documents/projects/iso/");


files_all<-list.files("result/iso_simulation_super_long/",
                      full.names=TRUE,pattern = "*unique_intron.tsv");

gene_trans_id_tbl<-"anno/hg19_ensembl_gene_id_trans_id_map.tsv";

bed_anno<-"anno/hg19_gencode_from_ucsc.bed";

most_likeli_order_output_path<-"./result/best_order/best_order_simulation_super_long.tsv"


read_count_threshold<-0

read_coverage_threshold<-0

source("code/necessary_code.R");


#t_igraph_list<-get_mlo_pipe(bed_anno,files_all,most_likeli_order_output_path,
#                            gene_trans_id_tbl,read_count_threshold,"",read_coverage_threshold)

t_igraph_list<-get_mlo_pipe(t_bed_anno=bed_anno,
                            t_files_all=files_all,
                            t_most_likeli_order_output_path=most_likeli_order_output_path,
                            t_gene_trans_id_tbl=gene_trans_id_tbl,
                            t_read_count_threshold=read_count_threshold,
                            t_trans_exp_file="",
                            t_read_cov_threshold=read_coverage_threshold,
                            t_trim_trans_id_by_dot=TRUE);