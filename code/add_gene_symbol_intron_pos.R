


add_gene_symbol_intron_pos<-function(t_igraph_list,t_gene_trans_id_map,t_intron_pos_mat_fr){
  
  for(g in 1:length(t_igraph_list) ){
    
    #t_igraph_list[[g]]$gene_symbol<-subset(t_gene_trans_id_map,trans_id==t_igraph_list[[g]]$trans_id,gene_symbol);
    
    
    intron_pos_mat_fr_one<-t_intron_pos_mat_fr[t_intron_pos_mat_fr$trans_id=="ENST00000237247",];
    intron_pos_mat_fr_one<-intron_pos_mat_fr_one[order(intron_pos_mat_fr_one$intron_order,decreasing=FALSE),];
    
    
    intron_pos_index_fr_one<-data.frame(
      pos=str_c(intron_pos_mat_fr_one[,"chr"],":",intron_pos_mat_fr_one[,"start"],"-",intron_pos_mat_fr_one[,"end"]),
      index=intron_pos_mat_fr_one[,"intron_order"],
      stringsAsFactors=FALSE
    );
    
    t_igraph_list[[g]]$intron_pos_index_fr<-intron_pos_index_fr_one
    
  }
  
  
  t_igraph_list
}


