

files_all<-list.files("./data/pombe/iso_pair/",full.names =TRUE,pattern = "*unique_intron.tsv");

gene_trans_id_tbl<-"./data/pombe/pombe_ensembl_gene_id_trans_id_map.tsv";

bed_anno<-"./data/pombe/Schizosaccharomyces_pombe.ASM294v2.43.chr_nothick.bed"

most_likeli_order_output_path<-"./results/best_order_pombe.tsv"


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


