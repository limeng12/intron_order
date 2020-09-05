library(gtools)
library(stringr)
library(stringi)
#library(ggplot2)
library(Rcpp)

#dynamic programming
#source("code/mlp8.R");
#sourceCpp("code/mlp9.cpp",cleanupCacheDir=FALSE);
sourceCpp("code/mlp15.cpp",cleanupCacheDir=FALSE);

#hill climbing
#source("code/mlp2.R")
sourceCpp("code/mlp2.cpp",cleanupCacheDir=FALSE);

##Cpp permutation
sourceCpp("code/mlp3.cpp",cleanupCacheDir=FALSE);

#linear programming
source("code/mlp16.R");


#can't be too small, influcence too much to the MLE

find_path_global<-function(t_adj_mat, t_alpha_v=0.1,is_verbose=FALSE){
  
  path_list<-NULL;
  
  ## 80% of transcrpts' intron number < 8;
  if(nrow(t_adj_mat) < 12){
    
    #path_list<-find_best_order_full2(t_adj_mat,colnames(t_adj_mat),t_alpha_v);
    path_list<-find_best_order_full2_c(t_adj_mat,t_alpha_v);
  }else{
  #}else if(nrow(t_adj_mat) < 13){
    #path_list<-find_opti_dynam(t_adj_mat,t_alpha_v);
    #path_list<-find_opti_dynam_r_cpp(t_adj_mat,t_alpha_v);
  #  path_list<-find_opti_dynam_r_cpp_bit(t_adj_mat,t_alpha_v);
    
  #}else if(nrow(t_adj_mat) < 56){
    #path_list<-del_find_path_iter(t_adj_mat, t_alpha_v);
    #path_list<-hill(t_adj_mat,t_alpha_v);
    #path_list<-hill_c(t_adj_mat,t_alpha_v);
    path_list<-lp_kenemy(t_adj_mat,t_alpha_v);
  }
  #}else{
  #  path_list<-hill_c(t_adj_mat,t_alpha_v);
    
  #}
  
  #best_score<-path_list$best_score;
  
  best_order<-path_list$best_order;
  best_score<-calp2(t_adj_mat,best_order,t_alpha_v) ;
  permut_p<-path_list$permut_p;
  entropy_one<-path_list$entropy;
  
  if(!is.na(entropy_one)){
    entropy_one<-entropy_one/log2( factorial(length(best_order) ) );
  }
  
  worst_score<-calp2(t_adj_mat,rev(best_order),t_alpha_v) ;
  
  
  #print(entropy_one);
  
  number_of_maximum_order<-path_list$number_of_maximum_order
  
  m_ini_intron_order<-colnames(t_adj_mat)[order(as.numeric(colnames(t_adj_mat)))];
  
  in_order_spliced_v<-calp2(t_adj_mat,m_ini_intron_order,t_alpha_v) ;
  
  
  #p-value of disorder
  #disorder_p_value<-disorder_p(t_adj_mat,best_order,t_alpha_v);
  disorder_p_value<-NA
  
  #chi_stat<- -2* log(exp(avg_li-best_score) );
  
  chi_stat<- -2* ((in_order_spliced_v-best_score) );
  
  normalized_relative_likelihood<-(exp(in_order_spliced_v)-0.99*exp(worst_score))/
    (exp(best_score)-0.98*exp(worst_score));
    
  normalized_relative_likelihood<-log(normalized_relative_likelihood);
  
  relative_likelihood<- ((in_order_spliced_v-best_score) );
  
  
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
  
  if(is_verbose&&(!is.null(number_of_maximum_order)) ){
    print(str_c("Best path number: ", number_of_maximum_order ) );
  }
  
  if((!is.na(permut_p))&&(permut_p<=0)){
    permut_p=1/1000000;
  }
    
  list(
    best_order=best_order,
    number_of_maximum_order=number_of_maximum_order,
    entropy=entropy_one,
    normalized_relative_likelihood=normalized_relative_likelihood,
    chi_stat=relative_likelihood,
    
    p_value=p_t,
    permut_p= log(permut_p),
       
    bs=chi_stat*(-1)/2,
    rela_likeli_v=NA, 
    disorder_p_value=disorder_p_value);
  
}


calp2<-function(t_adj_mat, t_intron_order, t_alpha_v){
  
  t_adj_mat<-t_adj_mat+t_alpha_v
  
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
  #print(p_all)
  p_all
}

  
  #best_intron_order<-m_ini_intron_order<-ini_intron_order;
  #allPerms_mat<-allPerms(ini_intron_order);
  #int_num<-length(ini_intron_order)
 

find_best_order_full2<-function(adj_mat,ini_intron_order,t_alpha_v,plot_dis=FALSE,output_id="",
                                suff_index_start=1,suff_index_end=length(ini_intron_order) ){
 
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
  #return(best_intron_order);

  
  list(best_score=max(li),best_order=best_intron_order,number_of_maximum_order=number_of_maximum_order);
  
  
}

