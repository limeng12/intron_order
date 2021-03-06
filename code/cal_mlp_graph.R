library(stringr);
options(scipen=999);

cal_mlp<-function(t_igraph_list,t_alpha=0){
  
  cat("\n");
  print("Calculating most likely order: ")
  print( length(t_igraph_list) )
  pb = txtProgressBar(min = 0, max = length(t_igraph_list), initial = 0, width=100, style=3) 
  #setTxtProgressBar(pb,i);
  
  for(i in 1:length(t_igraph_list) ){
    
    setTxtProgressBar(pb,i);
    
    
    adj_mat<-t_igraph_list[[i]]$adjacency_matrix;
    
    if(i%%100==1){
      #print(paste0(i,":",names(t_igraph_list)[i]) );
    }
    best_order_ls<-find_path_global(adj_mat, t_alpha);
    
    t_igraph_list[[i]]$best_order<-best_order_ls$best_order;
    
    t_igraph_list[[i]]$p_value<-best_order_ls$p_value;
    
    t_igraph_list[[i]]$chi_stat<-best_order_ls$chi_stat;
    #t_igraph_list[[i]]$chi_stat<-best_order_ls$normalized_relative_likelihood;
    
    
    t_igraph_list[[i]]$permut_p<-best_order_ls$permut_p;
    
    t_igraph_list[[i]]$entropy<-best_order_ls$entropy;
      
    t_igraph_list[[i]]$number_of_maximum_order<-best_order_ls$number_of_maximum_order;
    
    t_igraph_list[[i]]$disorder_p_value<-best_order_ls$disorder_p_value;
    
    
    best_oders<-t_igraph_list[[i]]$best_order
    #best_oders<-as.numeric( str_split(t_igraph_list[[i]]$best_order,",")[[1]] );
    
    if(length(best_oders)==1){
      t_igraph_list[[i]]$spearman_p_value<-NA
      t_igraph_list[[i]]$spearman_rho<-NA
      t_igraph_list[[i]]$spearman_rho_abs<-NA
    }else{
    
      #t_igraph_list[[i]]$spearman_p_value<- cor.test(best_oders,1:length(best_oders) ,method="spearman")$p.value
      
      t_igraph_list[[i]]$spearman_p_value_greater<- cor.test(best_oders,1:length(best_oders),
                                                     alternative="greater",method="spearman")$p.value
      
      t_igraph_list[[i]]$spearman_p_value_less<- cor.test(best_oders,1:length(best_oders),
                                                     alternative="less",method="spearman")$p.value
      
      t_igraph_list[[i]]$spearman_rho<-cor(best_oders,1:length(best_oders),method="spearman") 
      t_igraph_list[[i]]$spearman_rho_abs<-abs(cor(best_oders,1:length(best_oders),method="spearman" ) )  
      
    }
  }
  
  return(t_igraph_list);
}



  #unlink(output_file);
  
  #cat(c("gene_symbol$transcript_id$p_value_log$best_order$bayesian_factor$relative_likelihood\n"),
  #    file=str_c(output_file ) );
  #cat(c("gene_symbol$transcript_id$p_value_log$best_order$number_of_orders_have_same_prob$percent_coverage_order_pair\n"),
  #    file=str_c(output_file ) );
    # for(ii in 1:nrow(adj_mat)){
    #   for(jj in 1:ncol(adj_mat)){
    #     if(adj_mat[ii,jj]+adj_mat[jj,ii]<t_read_count_threshold){
    #       
    #       adj_mat[ii,jj]<-adj_mat[jj,ii]<-0;
    #     }
    #     
    #   }
    # }