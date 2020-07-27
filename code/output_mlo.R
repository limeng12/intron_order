
output_mlo<-function(t_igraph_list,t_gene_trans_id_map,t_intron_pos_mat_fr,output_file,t_trim_trans_id_by_dot){
  
  if(!is.null(t_gene_trans_id_map) ){
    t_gene_trans_id_map_v<-t_gene_trans_id_map[,"gene_symbol"];
    names(t_gene_trans_id_map_v)<-t_gene_trans_id_map[,"trans_id"];
    #best_orde_fr[,"gene_symbol"]<-t_gene_trans_id_map_v[best_orde_fr$transcript_id];
  }
  
  cat("\n");
  
  print("Summarise and output to file:");
  
  pb = txtProgressBar(min = 1, max = length(t_igraph_list) , initial = 0, width=100, style=3) 
  
  ##add gene symbol intron postion information
  for(g in 1:length(t_igraph_list) ){
    
    setTxtProgressBar(pb,g);
    
    if(t_trim_trans_id_by_dot){
      t_igraph_list[[g]]$trans_id<-sapply(strsplit(t_igraph_list[[g]]$trans_id,"\\."),"[",1);
      names(t_igraph_list)[g]<- t_igraph_list[[g]]$trans_id
      
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
  cat(c("gene_symbol$transcript_id$spearman_p_value_one_side_less$best_order$number_of_orders_have_same_prob$percent_coverage_order_pair$normalized_entropy$spearman_rho$spearman_rho_abs$relative_likelihood\n"),
      file=str_c(output_file ) );
  
  for(i in 1:length(t_igraph_list) ){
    cat(
      paste0(t_igraph_list[[i]]$gene_symbol,"$",
             t_igraph_list[[i]]$trans_id,"$",
             (t_igraph_list[[i]]$spearman_p_value_less),"$", 
             paste0(t_igraph_list[[i]]$best_order,collapse=","),"$",
             (t_igraph_list[[i]]$number_of_maximum_order),"$",
             (t_igraph_list[[i]]$percent_coverage_pair),"$",
             (t_igraph_list[[i]]$entropy_sum_normalized),"$",
             
             t_igraph_list[[i]]$spearman_rho,"$",
             t_igraph_list[[i]]$spearman_rho_abs,"$",
             t_igraph_list[[i]]$chi_stat,
             "\n"),
      append=TRUE, file=str_c(output_file ) );
    
  }
  
  
  t_igraph_list
}


