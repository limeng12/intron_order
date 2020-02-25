library(gtools)
library(stringr)
library(stringi)
library(ggplot2)
library(Rcpp)
#   scater_compare_real_mat();

#   scater_compare();

#hill climbing
#hill

#dynamic
#source("code/mlp8.R");
#dynamic CPP
sourceCpp("code/mlp9.cpp",cleanupCacheDir=FALSE);

#source("code/mlp2.R")
sourceCpp("code/mlp2.cpp",cleanupCacheDir=FALSE);

##Cpp permutation
sourceCpp("code/mlp3.cpp",cleanupCacheDir=FALSE);

#  adj_mat<-t_igraph_list[[1]]$adjacency_matrix;
#  ini_intron_order<-colnames(adj_mat)
#  t_adj_mat<-adj_mat
#  t_alpha_v=0.01

#source("code/mlp2.R")

disorder_p<-function(t_adj_mat,t_best_order,t_alpha_v){
  
  prob_bino<-0;
  for(i in 1:(length(t_best_order) -1) ){
    for(j in (i+1):length(t_best_order)){
      A_i_j<-t_adj_mat[t_best_order[i],t_best_order[j]];
      A_j_i<-t_adj_mat[t_best_order[j],t_best_order[i]];
      
      phi<-(A_i_j+t_alpha_v)/(A_i_j+A_j_i+2*t_alpha_v)
      prob_bino<-prob_bino+dbinom(A_i_j,(A_i_j+A_j_i),phi,log=TRUE);
      
    }
  }
  
  
  prob_bino_null<-0;
  for(i in 1:(length(t_best_order) -1) ){
    for(j in (i+1):length(t_best_order)){
      A_i_j<-t_adj_mat[t_best_order[i],t_best_order[j]];
      A_j_i<-t_adj_mat[t_best_order[j],t_best_order[i]];
      
      phi<-0.5
      prob_bino_null<-prob_bino_null+dbinom(A_i_j,(A_i_j+A_j_i),phi,log=TRUE);
      
    }
  }
  
  df_t<-0;
  for(i in 1:nrow(t_adj_mat)){
    for(j in 1:ncol(t_adj_mat)){
      if((i>j) && (t_adj_mat[i,j]+t_adj_mat[j,i]!=0 ) ){
        df_t<-df_t+1;
      }
      
    }
  }
  
  chi_stat<- -2* ((prob_bino_null-prob_bino) );
  
  p_t<-pchisq(chi_stat,df= df_t ,lower.tail = FALSE, log.p=TRUE);
  
  p_t;
}

#can't be too small, influcence too much to the MLE

find_path_global<-function(t_adj_mat, t_alpha_v=0.05){
  
  path_list<-NULL;
  
  ## 80% of transcrpts' intron number < 8;
  if(nrow(t_adj_mat) < 12){
    
    #path_list<-find_best_order_full2(t_adj_mat,colnames(t_adj_mat),t_alpha_v);
    path_list<-find_best_order_full2_c(t_adj_mat,t_alpha_v);
    
  #}else if(nrow(t_adj_mat) < 15) {
    
    #path_list<-del_find_path_iter(t_adj_mat, t_alpha_v);
  #  path_list<-hill(t_adj_mat,t_alpha_v);
    
  }else if(nrow(t_adj_mat) < 20){
    #path_list<-find_opti_dynam(t_adj_mat,t_alpha_v);
    path_list<-find_opti_dynam_r_cpp(t_adj_mat,t_alpha_v);
    
  }else{
    #path_list<-del_find_path_iter(t_adj_mat, t_alpha_v);
    #path_list<-hill(t_adj_mat,t_alpha_v);
    path_list<-hill_c(t_adj_mat,t_alpha_v);
    
  }
  
  #best_score<-path_list$best_score;
  
  best_order<-path_list$best_order;
  best_score<-calp2(t_adj_mat+t_alpha_v,best_order,t_alpha_v) ;
  
  
  number_of_maximum_order<-path_list$number_of_maximum_order
  
  m_ini_intron_order<-colnames(t_adj_mat)[order(as.numeric(colnames(t_adj_mat)))];
  
  avg_li<-calp2(t_adj_mat+t_alpha_v,m_ini_intron_order,t_alpha_v) ;
  
  
  #p-value of disorder
  
  disorder_p_value<-disorder_p(t_adj_mat,best_order,t_alpha_v);
  
  
  #chi_stat<- -2* log(exp(avg_li-best_score) );
  
  chi_stat<- -2* ((avg_li-best_score) );
  
  ##p-value of splicing is in order
  df_t<-0;
  for(i in 1:nrow(t_adj_mat)){
    for(j in 1:ncol(t_adj_mat)){
      if((i>j) && (t_adj_mat[i,j]+t_adj_mat[j,i]!=0 ) ){
      
      if(
          ( which( m_ini_intron_order==i ) < which( m_ini_intron_order==j ) ) !=
          ( which( best_order==i ) < which( best_order==j ) )
          
      ) {
        df_t<-df_t+1;
        
      }
        
      }
      
    }
  }
  
  #pt_t<-log(0.5);
  
  #if(chi_stat!=0){
    p_t<-pchisq(chi_stat,df= df_t ,lower.tail = FALSE, log.p=TRUE);
  #}
  
  if(is.infinite(p_t)){
    p_t<- -10000;
  }
  
  if(!is.null(number_of_maximum_order)){
    print(str_c("Best path number: ", number_of_maximum_order ) );
  }
  
  list(p_value=p_t,
       best_order=best_order,bs=chi_stat*(-1)/2,
       rela_likeli_v=rela_likeli(0,avg_li,df_t,best_score), 
       number_of_maximum_order=number_of_maximum_order,
       disorder_p_value=disorder_p_value);
  
}


rela_likeli<-function(k_1,l_1, k_2,l_2){
  
  aic_1<-2*k_1-2*(l_1);
  
  aic_2<-2*k_2-2*(l_2);
  
  a<-exp( (aic_1-aic_2)/2 )
  
  a
}


calp2<-function(t_adj_mat, t_intron_order, t_alpha_v){
  
  p_all<-0;
  
  #for(i in 1:(length(t_intron_order)-1) ){
  order_comb<-combn(t_intron_order,2);
  
  
  #calculate the prob
  for(i in 1:  ncol(order_comb)){
    
    order_first<-order_comb[1,i];
    order_next<-order_comb[2,i];
    
    ##skipped orders that no read count support
    if( (t_adj_mat[order_first,order_next]==t_alpha_v) && (t_adj_mat[order_next,order_first]==t_alpha_v) ){
      next;
    }
    
    p_one<-log( t_adj_mat[order_first,order_next]/
                     (t_adj_mat[order_first,order_next]+t_adj_mat[order_next,order_first] ) );
    #print(p_one)
    
    p_all<-p_all+p_one;
    
  }
  #  print(p_all)
  p_all
}


find_best_order_full2<-function(adj_mat,ini_intron_order,t_alpha_v,plot_dis=FALSE,output_id="",
                                suff_index_start=1,suff_index_end=length(ini_intron_order) ){
  
  #best_intron_order<-m_ini_intron_order<-ini_intron_order;
  
  #allPerms_mat<-allPerms(ini_intron_order);
  #int_num<-length(ini_intron_order)
  t_adj_mat<-adj_mat+t_alpha_v;
  
  
  int_num<-suff_index_end-suff_index_start+1;
  
  allPerms_mat<-permutations(n = int_num, r = int_num, v = ini_intron_order[suff_index_start:suff_index_end]);
  
  li<-rep(NA,nrow(allPerms_mat) );
  
  for( i in 1:nrow(allPerms_mat) ){
    if(i %% 10000 ==1){
      #   print(i);
    }
    
    t_ini_intron_order<-ini_intron_order;
    
    t_ini_intron_order[suff_index_start:suff_index_end]<-as.character(allPerms_mat[i,])
    
    li[i]<-calp2(t_adj_mat, t_ini_intron_order,t_alpha_v);
    
  }
  
  
  if(plot_dis){
    a<-stri_rand_strings(1,100)
    
    print(paste0("result/mlv_dis/",output_id,".jpeg"))
    
    jpeg(paste0("result/mlv_dis/",output_id,".jpeg"),width = 600, height=500 );
    li_order<-sort(li);
    
    #qplot(1:length(li_order),li_order,xlab="sorted intron order",ylab="Log likelihood value");
    plot(1:length(li_order),li_order,xlab="sorted intron order",
         ylab="Log likelihood value",pch=4);# type="s",cex=0.5
    dev.off();
    
  }
  
  best_intron_order<-
    as.character(allPerms_mat[which.max(li),]);
  
  number_of_maximum_order<-sum(li==max(li));
  
  #print(str_c("Best path number: ", number_of_maximum_order ) );
  
  
  list(best_score=max(li),best_order=best_intron_order,number_of_maximum_order=number_of_maximum_order);
  
  #return(best_intron_order);
  
}



permute_find_opti<-function(t_adj_mat,init_order, t_alpha_v, only_best=TRUE ){
  
  best_order_list<-list()
  best_order_list[[1]]<-init_order;
  
  best_order_v<-calp2(t_adj_mat, init_order, t_alpha_v);
  
  #calp2(t_adj_mat, t_intron_order, t_alpha_v)
  
  for(i in 1:length(init_order)){
    
    for(j in 1:length(init_order)){
      
      for(k in 1:length(best_order_list)){
        if(i<j){
          t_intron_order<-best_order_list[[k]];
          
          t<-t_intron_order[i];
          t_intron_order[i]<-t_intron_order[j];
          
          t_intron_order[j]<-t;
          #t_intron_order[i]<-t_intron_order[j];
          
          t_best_order_v<-calp2(t_adj_mat, t_intron_order, t_alpha_v);
          
          if(t_best_order_v>best_order_v){
            best_order_list[[k]]<-t_intron_order
            
            best_order_v<-t_best_order_v
          }
          
          if((!only_best) && (t_best_order_v==best_order_v) && (length(best_order_list)<100)){
            best_order_list[[length(best_order_list)+1]]<-t_intron_order

          }
        }
      }
    }
    
    
  }
  
  return(best_order_list)
}

# 
# test<-function( t_alpha_v){
#   
#   t<-2
#   
#   adj_mat <- matrix(nrow = t,ncol = t);
#   
#   colnames(adj_mat)<-1:t
#   rownames(adj_mat)<-1:t
#   
#   
#   for( i in 1:nrow(adj_mat)){
#     for( j in 1:ncol(adj_mat)){
#       if(i!=j){
#         #p<-runif(1)
#         adj_mat[i,j]<-rpois(1,20);
#         adj_mat[j,i]<-rpois(1,15);
#         
#         adj_mat[i,j]<-adj_mat[j,i]<-1;
#         
#       }
#       if(i==j){
#         adj_mat[i,j]<-0;
#       }
#       
#     }
#   }
#   
#   
#   t_adj_mat2<-adj_mat[1:t,1:t, drop=FALSE];
#   
#   #t_adj_mat3<-t_adj_mat2;
# 
#   
#   path1<-find_path_global(t_adj_mat2,t_alpha_v);
#   
#   disorder_p_v<-disorder_p(t_adj_mat2,path1$best_order,t_alpha_v);
#   
#   
#   #print("Test MLO algorithm")
#   return(path1) ;
# }
# 
# 
# test(0.1);

