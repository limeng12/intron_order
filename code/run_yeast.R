

files_all<-list.files("./data/yeast/iso_pair/",full.names =TRUE,pattern = "*unique_intron.tsv");

gene_trans_id_tbl<-"./data/yeast/gene_id_trans_id_yeast.txt";

bed_anno<-"./data/yeast/sacCer3_no_thick.bed"

most_likeli_order_output_path<-"./results/best_order_yeast.tsv"


read_count_threshold<-0

source("code/necessary_code.R");

t_igraph_list<-get_mlo_pipe(t_bed_anno=bed_anno,
                            t_files_all=files_all,
                            t_most_likeli_order_output_path=most_likeli_order_output_path,
                            t_gene_trans_id_tbl=gene_trans_id_tbl,
                            t_read_count_threshold=read_count_threshold,
                            t_trans_exp_file="",
                            t_read_cov_threshold=0.95,
                            t_trim_trans_id_by_dot=FALSE);

#save( t_igraph_list, file="pombe/t_igraph_list.Rd",version = 2);

