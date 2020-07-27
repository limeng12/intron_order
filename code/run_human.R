

files_all_iso_pairs<-list.files("./data/",full.names =TRUE,pattern = "*_intron.tsv");

gene_trans_id_tbl<-"./data/hg19_ensembl_gene_id_trans_id_map.tsv";

bed_anno<-"./data/hg19_gencode_from_ucsc_nothick_nocds.bed"

most_likeli_order_output_path<-"./results/best_order_human.tsv"

read_count_threshold<-0;

source("code/necessary_code.R");
#################################################################################################################


t_igraph_list<-get_mlo_pipe(t_bed_anno=bed_anno,
                            t_files_all=files_all_iso_pairs,
                            t_most_likeli_order_output_path=most_likeli_order_output_path,
                            t_gene_trans_id_tbl=gene_trans_id_tbl,
                            t_read_count_threshold=read_count_threshold,
                            t_trans_exp_file="./data/k562_exp_trans_id.tsv",
                            t_read_cov_threshold=0.95,
                            t_trim_trans_id_by_dot=TRUE,
                            t_reatined_intron_psi_file="rmats_human/RI.MATS.JC.txt");

#dir.create("./results/adj_matrix",showWarnings = FALSE)

#for( i in 1:length(t_igraph_list) ){
  
#  write.table(t_igraph_list[[i]]$adjacency_matrix, file=str_c("./results/adj_matrix/",names(t_igraph_list)[i],".tsv" ),
#              sep="\t",col.names = FALSE,row.names = FALSE  );
#}


#save( t_igraph_list, file="result/t_igraph_list.Rd",version = 2);

