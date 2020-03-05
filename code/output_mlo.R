


output_mlo<-function(t_igraph_list,t_gene_trans_id_map,t_intron_pos_mat_fr,output_file,t_trim_trans_id_by_dot){
  
  if(!is.null(t_gene_trans_id_map) ){
    t_gene_trans_id_map_v<-t_gene_trans_id_map[,"gene_symbol"];
    names(t_gene_trans_id_map_v)<-t_gene_trans_id_map[,"trans_id"];
    #best_orde_fr[,"gene_symbol"]<-t_gene_trans_id_map_v[best_orde_fr$transcript_id];
  }
  
  
  ##add gene symbol intron postion information
  for(g in 1:length(t_igraph_list) ){
    
    intron_pos_mat_fr_one<-t_intron_pos_mat_fr[t_intron_pos_mat_fr$trans_id==t_igraph_list[[g]]$trans_id,];
    
    intron_pos_mat_fr_one<-intron_pos_mat_fr_one[order(intron_pos_mat_fr_one$intron_order,decreasing=FALSE),];
    
    
    intron_pos_index_fr_one<-data.frame(
      pos=str_c(intron_pos_mat_fr_one[,"chr"],":",intron_pos_mat_fr_one[,"start"],"-",intron_pos_mat_fr_one[,"end"]),
      index=intron_pos_mat_fr_one[,"intron_order"],
      stringsAsFactors=FALSE
    );
    
    t_igraph_list[[g]]$intron_pos_index_fr<-intron_pos_index_fr_one
    
    
    
    if(t_trim_trans_id_by_dot){
      t_igraph_list[[g]]$trans_id<-sapply(strsplit(t_igraph_list[[g]]$trans_id,"\\."),"[",1);
    }
    
    if( !is.null(t_gene_trans_id_map) ) {
      t_igraph_list[[g]]$gene_symbol<-t_gene_trans_id_map_v[t_igraph_list[[g]]$trans_id];
    }    
    
  
  }
  
  
  
  #best_orde_fr[,"spearman_rho"]<-spearman_rho;
  #best_orde_fr[,"spearman_rho_abs"]<-spearman_rho_abs;
  #best_orde_fr[,"spearman_p_value"]<-spearman_p_value;
  unlink(output_file)
  ###output to file
  cat(c("gene_symbol$transcript_id$relative_likelihood$best_order$number_of_orders_have_same_prob$percent_coverage_order_pair$disorder_p_value$spearman_rho$spearman_rho_abs$spearman_p_value\n"),
      file=str_c(output_file ) );
  
  for(i in 1:length(t_igraph_list) ){
    cat(
      paste0(t_igraph_list[[i]]$gene_symbol,"$",
             t_igraph_list[[i]]$trans_id,"$",
             (t_igraph_list[[i]]$chi_stat),"$", 
             paste0(t_igraph_list[[i]]$best_order,collapse=","),"$",
             (t_igraph_list[[i]]$number_of_maximum_order),"$",
             (t_igraph_list[[i]]$percent_coverage_pair),"$",
             (t_igraph_list[[i]]$disorder_p_value),"$",
             
             t_igraph_list[[i]]$spearman_rho,"$",
             t_igraph_list[[i]]$spearman_rho_abs,"$",
             t_igraph_list[[i]]$spearman_p_value,
             "\n"),
      append=TRUE, file=str_c(output_file ) );
    
  }
  
  
  t_igraph_list
}


