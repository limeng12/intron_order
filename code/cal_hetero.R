
cal_hetero<-function(t_igraph_list, t_alpha){
  
  cat("\n");
  print("calculate intron splicing order matrix heterogeneity:");
  pb = txtProgressBar(min = 0, max = length(t_igraph_list), initial = 0, width=100, style=3) 
  
  
  for(i in 1:length(t_igraph_list) ){
    
    setTxtProgressBar(pb,i);
    
    
    t_adj_mat<-t_igraph_list[[i]]$adjacency_matrix;
    
    t_best_orders<-t_igraph_list[[i]]$best_order;
    
    t_adj_mat_order<-t_adj_mat[t_best_orders,t_best_orders];
    
    
    entropy_sum<-0
    n<-0
    t_adj_mat_li<-t_adj_mat_order;
    for(ii in 1:nrow(t_adj_mat_li)){
      for(jj in 1:ncol(t_adj_mat_li)){
        if(ii<jj){
          p<-(t_adj_mat_li[ii,jj]+t_alpha) / (t_adj_mat_li[ii,jj]+t_adj_mat_li[jj,ii]+2*t_alpha);
          entropy_sum<-entropy_sum + (-1*p*log2(p) );
          n<-n+1
          #t_adj_mat_li[i,j]<-p*log(p);
          #t_adj_mat_li[j,i]<-1-p;
        }
        

      }
      
    }
    
    t_entropy_sum_normalized <- entropy_sum/n;
    #print(t_entropy_sum_normalized)
    t_igraph_list[[i]]$entropy_sum_normalized<-t_entropy_sum_normalized
    
  }
  
  return(t_igraph_list);
  
}
