
files_all<-list.files("./data/fly/iso_pair/",full.names =TRUE,pattern = "*unique_intron.tsv");

gene_trans_id_tbl<-"./data/fly/gene_id_trans_id_fly.tsv";

bed_anno<-"./data/fly/dm6_ensembl_no_thick.bed"

most_likeli_order_output_path<-"./results/best_order_fly.tsv"


read_count_threshold<-0

source("code/necessary_code.R");

#t_igraph_list<-get_mlo_pipe(bed_anno,files_all,most_likeli_order_output_path,
#                            gene_trans_id_tbl,read_count_threshold,"fly/s2_exp_trans_id.tsv");

t_igraph_list<-get_mlo_pipe(t_bed_anno=bed_anno,
                            t_files_all=files_all,
                            t_most_likeli_order_output_path=most_likeli_order_output_path,
                            t_gene_trans_id_tbl=gene_trans_id_tbl,
                            t_read_count_threshold=read_count_threshold,
                            t_trans_exp_file="./data/fly/s2_exp_trans_id.tsv",
                            t_read_cov_threshold=0.95,
                            t_trim_trans_id_by_dot=TRUE);


