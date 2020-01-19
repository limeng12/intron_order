library(stringr);
options(scipen=999);


cal_mlp<-function(t_igraph_list,output_file,t_alpha,t_read_count_threshold=0){
  
  #unlink(str_c("result/best_order/best_order.tsv" ) );
  unlink(output_file);
  
  #cat(c("gene_symbol\ttranscript_id\tp_value_log\tbest_order\tbayesian_factor\trelative_likelihood\n") ,file=str_c(output_file ) );
  cat(c("gene_symbol\ttranscript_id\tp_value_log\tbest_order\n"),
      file=str_c(output_file ) );
  #  i<-1;
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
    
    #best_order_ls<-del_find_path_iter(adj_mat,t_alpha);
    best_order_ls<-find_path_global(adj_mat, t_alpha);
    
    
    t_igraph_list[[i]]$best_order<-best_order_ls$best_order;
    
    t_igraph_list[[i]]$p_value<-best_order_ls$p_value;
    
    
    #print( best_order_ls$best_order );
    
    
    cat(
      
      str_c(t_igraph_list[[i]]$gene_symbol,"\t",t_igraph_list[[i]]$trans_id,"\t",(best_order_ls$p_value),"\t", 
            paste0(best_order_ls$best_order,collapse=","),"\t",best_order_ls$bs,"\t",
            best_order_ls$rela_likeli_v,"\n"),
      append=TRUE, file=str_c(output_file ) );
    
    
    #write.table(best_order,file=str_c("result/best_order/",names(t_igraph_list)[i],".tsv" ),
    #            sep="\t",col.names = FALSE,row.names = FALSE  );
    
  }
  
  

  
  #setwd("/Users/mengli/Documents/projects/iso");
  
  best_orde_fr<-read.table(output_file,header=TRUE,as.is=TRUE,sep="\t" ) 
  
  
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
                        cor.test(best_oders,1:length(best_oders) )$p.value  );
    spearman_rho<-c(spearman_rho,cor(best_oders,1:length(best_oders)) );
    
    spearman_rho_abs<-c(spearman_rho_abs,abs(cor(best_oders,1:length(best_oders) ) )  );
    
  }
  
  
  
  best_orde_fr[,"spearman_rho"]<-spearman_rho;
  
  best_orde_fr[,"spearman_rho_abs"]<-spearman_rho_abs;
  
  best_orde_fr[,"spearman_p_value"]<-spearman_p_value;
  
  
  meta_p_spearman<-pchisq(-2*sum( (log(spearman_p_value[spearman_p_value>0]) ) ,na.rm = TRUE),
                          df=length(spearman_p_value), lower.tail = FALSE,log.p = TRUE);
  
  
  ######best_orde_fr$p.value<- log(best_orde_fr$p.value);
  
  best_orde_fr[,"p_value"]<- exp(best_orde_fr$p_value_log);
  
  p_values<-best_orde_fr[,"p_value_log"];
  
  ###fisher's method combine P-values
  chi_stat<- ( -2*sum( (p_values) ) );
  
  meta_p<-pchisq(chi_stat,df=length(p_values),lower.tail = FALSE,log.p = TRUE);
  
  
  cat(str_c("####P-value of in oder splicing (meta,log): ",meta_p,
            "      ###P-value of in oder splicing (meta spearman,log): ",meta_p_spearman,"\n"),
      file=output_file);
  
  
  best_orde_fr<-best_orde_fr[order(best_orde_fr$p_value_log,decreasing=TRUE),];
  
  write.table(best_orde_fr,file=output_file, 
              col.names = TRUE,row.names = FALSE, sep="\t" ,append = TRUE,quote = FALSE);
  
  
  return(t_igraph_list);
}


