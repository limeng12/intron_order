library(stringr);
options(scipen=999);

cal_mlp<-function(t_igraph_list,output_file,t_alpha,t_read_count_threshold=0){
  
  unlink(output_file);
  
  #cat(c("gene_symbol$transcript_id$p_value_log$best_order$bayesian_factor$relative_likelihood\n"),
  #    file=str_c(output_file ) );
  #cat(c("gene_symbol$transcript_id$p_value_log$best_order$number_of_orders_have_same_prob$percent_coverage_order_pair\n"),
  #    file=str_c(output_file ) );
  
  for(i in 1:length(t_igraph_list) ){
    
    adj_mat<-t_igraph_list[[i]]$adjacency_matrix;
    
    for(ii in 1:nrow(adj_mat)){
      for(jj in 1:ncol(adj_mat)){
        if(adj_mat[ii,jj]+adj_mat[jj,ii]<t_read_count_threshold){
          
          adj_mat[ii,jj]<-adj_mat[jj,ii]<-0;
        }
        
      }
    }
    
    
    print(paste0(i,":",names(t_igraph_list)[i]) );
    
    best_order_ls<-find_path_global(adj_mat, t_alpha);
    
    t_igraph_list[[i]]$best_order<-best_order_ls$best_order;
    
    t_igraph_list[[i]]$p_value<-best_order_ls$p_value;
    
    t_igraph_list[[i]]$number_of_maximum_order<-best_order_ls$number_of_maximum_order;
    
    t_igraph_list[[i]]$disorder_p_value<-best_order_ls$disorder_p_value;
    
    # cat(
    #   
    #   str_c(t_igraph_list[[i]]$gene_symbol,"$",t_igraph_list[[i]]$trans_id,"$",(best_order_ls$p_value),"$", 
    #         paste0(best_order_ls$best_order,collapse=","),"$",best_order_ls$bs,"$",
    #         best_order_ls$rela_likeli_v,"\n"),
    #   append=TRUE, file=str_c(output_file ) );
    
  }
  
  # output_file<-"result/best_order/best_order_simulation_pair.tsv"
  
  cat(c("gene_symbol$transcript_id$p_value_log$best_order$number_of_orders_have_same_prob$percent_coverage_order_pair$disorder_p_value\n"),
      file=str_c(output_file ) );
  
  for(i in 1:length(t_igraph_list) ){
    cat(
      paste0(t_igraph_list[[i]]$gene_symbol,"$",
            t_igraph_list[[i]]$trans_id,"$",
            (t_igraph_list[[i]]$p_value),"$", 
            paste0(t_igraph_list[[i]]$best_order,collapse=","),"$",
            (t_igraph_list[[i]]$number_of_maximum_order),"$",
            (t_igraph_list[[i]]$percent_coverage_pair),"$",
            (t_igraph_list[[i]]$disorder_p_value),"\n"),
      append=TRUE, file=str_c(output_file ) );
   
  }
  
  # output_file<-"./fly/result/best_order.tsv"
  
  ##add spearman correlation, meta analysis p-values
  best_orde_fr<-as.data.frame( read_delim(output_file,col_names=TRUE,delim="$",
                                          col_types="ccdcidddddd")  );

  spearman_rho<-c();
  spearman_p_value<-c();
  spearman_rho_abs<-c()
  
  for(i in 1:nrow(best_orde_fr) ){
    
    best_oders<-as.numeric( str_split(best_orde_fr[i,"best_order"],",")[[1]] );
    
    if(length(best_oders)<3){
      spearman_p_value<-c(spearman_p_value,NA);
      spearman_rho<-c(spearman_rho,NA);
      spearman_rho_abs<-c(spearman_rho_abs,NA);
      next;
    }
    
    spearman_p_value<-c(spearman_p_value,
                        cor.test(best_oders,1:length(best_oders) ,method="spearman")$p.value  );
    spearman_rho<-c(spearman_rho,cor(best_oders,1:length(best_oders)) );
    
    spearman_rho_abs<-c(spearman_rho_abs,abs(cor(best_oders,1:length(best_oders),method="spearman" ) )  );
    
  }
  
  
  best_orde_fr[,"spearman_rho"]<-spearman_rho;
  best_orde_fr[,"spearman_rho_abs"]<-spearman_rho_abs;
  best_orde_fr[,"spearman_p_value"]<-spearman_p_value;
  
  
  meta_p_spearman<-pchisq(-2*sum( (log(spearman_p_value[spearman_p_value>0]) ) ,na.rm = TRUE),
                          df=length(spearman_p_value), lower.tail = FALSE,log.p = TRUE);
  

  best_orde_fr[,"p_value"]<- exp(best_orde_fr$p_value_log);
  
  p_values<-best_orde_fr[,"p_value_log"];
  
  ###fisher's method combine P-values
  chi_stat<- ( -2*sum( (p_values) ,na.rm = TRUE) );
  
  meta_p<-pchisq(chi_stat,df=length(p_values)*2,lower.tail = FALSE,log.p = TRUE);
  
  
  cat(str_c("####p-value of in oder splicing (meta,log): ",format(meta_p),
            "      ###p-value of in oder splicing (meta spearman,log): ",format(meta_p_spearman,scientific=TRUE),"\n"),
      file=output_file);
  
  best_orde_fr<-best_orde_fr[order(best_orde_fr$p_value_log,decreasing=TRUE),];
  
  write.table(best_orde_fr,file=output_file, 
              col.names = TRUE,row.names = FALSE, sep="$" ,append = TRUE,quote = FALSE);
  
  return(t_igraph_list);
}

