library(gtools)
library(stringr)
##
#  adj_mat<-t_igraph_list[[1]]$adjacency_matrix;
#  ini_intron_order<-colnames(adj_mat)
#  t_adj_mat<-adj_mat
#  t_alpha_v=0.01

#source("code/mlp2.R")

#can't be too small, influcence too much to the MLE



find_path_global<-function(t_adj_mat, t_alpha_v=0.05){
  
  path_list<-NULL;
  
  ## 80% of transcrpts' intron number < 8;
  if(nrow(t_adj_mat) < 8){
    
    path_list<-find_best_order_full2(t_adj_mat,colnames(t_adj_mat),t_alpha_v);
    
  }else{
    
    path_list<-del_find_path_iter(t_adj_mat,t_alpha_v);
    
  }
  
  best_score<-path_list$best_score
  
  best_order<-path_list$best_order;
  
  
  m_ini_intron_order<-colnames(t_adj_mat)[order(as.numeric(colnames(t_adj_mat)))];
  avg_li<-calp2(t_adj_mat+t_alpha_v,m_ini_intron_order,t_alpha_v) ;
  
  
  #chi_stat<- -2* log(exp(avg_li-best_score) );
  
  chi_stat<- -2* ((avg_li-best_score) );
  
  ##p-value of splicing is in order
  df_t<-0;
  for(i in 1:nrow(t_adj_mat)){
    for(j in 1:ncol(t_adj_mat)){
      if((i>j) && (t_adj_mat[i,j]+t_adj_mat[j,i]!=0 ) ){
        df_t<-df_t+1;
      }
      
    }
  }
  
  p_t<-pchisq(chi_stat,df= df_t ,lower.tail = FALSE, log.p=TRUE);
  
  
  if(is.infinite(p_t)){
    p_t<- -10000;
  }
  
  list(p_value=p_t,
       best_order=best_order,bs=chi_stat*(-1)/2,rela_likeli_v=rela_likeli(0,avg_li,df_t,best_score) );
  
}


rela_likeli<-function(k_1,l_1, k_2,l_2){
  
  aic_1<-2*k_1-2*(l_1);
  
  aic_2<-2*k_2-2*(l_2);
  
  a<-exp( (aic_1-aic_2)/2 )
  
  a
}





del_find_path_iter<-function(t_adj_mat2, t_alpha_v=0.05){
  # t_adj_mat<-t_adj_mat[1:11,1:11]
  
  t_adj_mat<-t_adj_mat2
  
  ##likehood matrix
  t_adj_mat_li<-t_adj_mat;
  
  colnames(t_adj_mat_li)<-colnames(t_adj_mat);
  rownames(t_adj_mat_li)<-rownames(t_adj_mat);
  
  
  for(i in 1:nrow(t_adj_mat_li)){
    for(j in 1:ncol(t_adj_mat_li)){
      if(i==j){
        t_adj_mat_li[i,j]<-t_adj_mat_li[j,i]<-0;
      }
      
      if(t_adj_mat_li[i,j]+t_adj_mat_li[j,i]==0){
        next;
      }
      
      if(i>j){
        p<- (t_adj_mat_li[i,j]+t_alpha_v)/(t_adj_mat_li[i,j]+t_adj_mat_li[j,i]+2*t_alpha_v);
        
        t_adj_mat_li[i,j]<-p
        t_adj_mat_li[j,i]<-1-p
      }
      
    }
    
  }
  
  star_grow<-c();
  
  #m_best_order_list<-list();
  max_order_count<-5000;
  m_best_order_list<-list();
  
  ###consier second-best rowSums
  max_idex_consider_num<-as.integer( nrow(t_adj_mat_li)/10 );
  
  if(max_idex_consider_num<3){
    max_idex_consider_num<- 3;
  }
  
  #if(nrow(t_adj_mat_li))
  #del_iter(t_adj_mat_li,star_grow);
  #[[str_c(sample(1:100000000,1),
  #                       sample(1:100000000,1) ,"_")]] <-
    
  best_order_list<-del_iter(t_adj_mat_li, star_grow, max_idex_consider_num, m_best_order_list, max_order_count);
  
  
  li_c<-c();
  
  t_adj_mat<-t_adj_mat2+t_alpha_v;
  
  
  for(i in 1:length(best_order_list ) ){
    
    li_c<-c(li_c,calp2(t_adj_mat,best_order_list[[i]],t_alpha_v) );
    
  }
  
  
  max_index<-which(max(li_c)==li_c);
  
  print(str_c("Best path number: ", length(max_index) ) );
  
  max_index<-which.max(li_c);
  
  
  best_order<-best_order_list[[max_index]];
  
  best_score<-calp2(t_adj_mat,best_order_list[[max_index]],t_alpha_v);
  
  t_igraph_list<-best_order_list[[max_index]];
  
  
  list(best_score=best_score,best_order=best_order);
 
  #c( p_t, best_order );
  
}




del_iter<-function(t_adj_mat_li, star_grow, t_max_idex_consider_num, best_order_list, max_order_count=1000){
  #max_idex_consider_num<-4;
  
  
  if(length(best_order_list)> max_order_count){
    return(best_order_list);
  }
  
  
  if( nrow(t_adj_mat_li)==0 ){
    best_order_list[[str_c(sample(1:100000000,1), sample(1:100000000,1) ,"_") ]]<-star_grow 
    
    return(best_order_list);
    
  }
  
  #while(  (nrow(t_adj_mat_li)>0 ) ){
  
  #the sum likehood  of each node spliced before all others
  row_sum_li<-rowSums(t_adj_mat_li);
  
  max_index<-which(max(row_sum_li)==row_sum_li);
  
  
  row_sum_li_t<-row_sum_li;
  
  while(length(max_index)<t_max_idex_consider_num  ){
    
    row_sum_li_t[max(row_sum_li_t)==row_sum_li_t]<-0;
    
    if( sum(row_sum_li_t)==0 ){
      break;
    }
    
    max_index<-c(max_index, which(max(row_sum_li_t)==row_sum_li_t) );
    #max_index<-which(max(row_sum_li)=row_sum_li);
    
  }
  
  
  #print(length(max_index) );
  
  star_grow_back<-star_grow;
  
  t_adj_mat_li_back<-t_adj_mat_li;
  
  for( i in 1:length(max_index) ){
    
    t_adj_mat_li<-t_adj_mat_li_back;
    
    best_start<-rownames(t_adj_mat_li)[ max_index[i] ];
    
    #if(is.na(best_start) )
    #{print("bug");
    #  View(max_index)
    #  View(t_adj_mat_li);
    #  stop();
    # }

    star_grow_new<-c(star_grow_back, best_start);
    #t_adj_mat_li[best_start,]<-;
    #t_adj_mat_li[,best_start]<-;
    
    
    t_adj_mat_li<-t_adj_mat_li[setdiff(rownames(t_adj_mat_li),best_start), , drop=FALSE];
    
    
    ##if ncol(t_adj_mat_li) is 0, then only one element in t_adj_mat_li_back, thus max_index's length is 1
    if(ncol(t_adj_mat_li)>0){
      t_adj_mat_li<-t_adj_mat_li[,setdiff(colnames(t_adj_mat_li),best_start),drop=FALSE]
    }else{
      best_order_list[[str_c(sample(1:100000000,1), sample(1:100000000,1) ,"_") ]]<-star_grow_new;
      
      return(best_order_list);
    }
      
    
    
    if(nrow(t_adj_mat_li)>0){
      #best_order_list[[length(best_order_list)+1]]<-
      #  del_iter(t_adj_mat_li,star_grow_new);

      tmp_list<-list();
      
      tmp_list<-del_iter(t_adj_mat_li,star_grow_new,
               t_max_idex_consider_num,best_order_list,max_order_count)
      
      best_order_list<-unique( c(best_order_list, tmp_list) );
      
    }else{
      best_order_list[[str_c(sample(1:100000000,1), sample(1:100000000,1) ,"_") ]]<-star_grow_new 
      return(best_order_list);
    }
    
    
  }
  
  return(best_order_list);
  
}



calp2<-function(t_adj_mat,t_intron_order,t_alpha_v){
  
  #t_adj_mat<-t_adj_mat+t_alpha_v;
  
  p_all<-0;
  
  #combn(t_intron_order,2)
  
  #for(i in 1:(length(t_intron_order)-1) ){
  order_comb<-combn(t_intron_order,2);
  
  ###find all the order pairs that consistent of the input order(t_intron_order)
  #order_comb<-matrix(ncol=0,nrow = 2);
  #for(i in 1:(length(t_intron_order)-1) ){
  
  #  one_pair<-as.character(combn(t_intron_order[(i+1):length(t_intron_order)],1));
  
  #  order_comb<-cbind(order_comb,matrix(c(rep(i,length(one_pair)),one_pair),byrow = TRUE,nrow = 2 )   );
  #}
  
  
  #calculate the prob
  for(i in 1:  ncol(order_comb)){
    
    order_first<-order_comb[1,i];
    order_next<-order_comb[2,i];
    
    ##skipped orders that no read count support
    if( (t_adj_mat[order_first,order_next]==t_alpha_v) && (t_adj_mat[order_next,order_first]==t_alpha_v) ){
      next;
    }
    
    p_one<-log(    t_adj_mat[order_first,order_next]/
                     (t_adj_mat[order_first,order_next]+t_adj_mat[order_next,order_first] ) );
    #print(p_one)
    
    p_all<-p_all+p_one;
    
  }
  #  print(p_all)
  p_all
}




find_best_order_full2<-function(adj_mat,ini_intron_order,t_alpha_v,
                                suff_index_start=1,suff_index_end=length(ini_intron_order)){
  
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
  
  best_intron_order<-
    as.character(allPerms_mat[which.max(li),]);
  
  
  list(best_score=max(li),best_order=best_intron_order);
  
  #return(best_intron_order);
  
}



scater_compare<-function(){
  dyn<-c()
  
  opti<-c();
  
  for( i in 1:1000){
    a<-compare(0.1);
    
    dyn<-c(dyn,a[1])
    
    opti<-c(opti,a[2])
    
  }
  
  pdf("result/compare_dyn_opti.pdf")
  scatter.smooth(dyn,opti,xlab="likelihood value by dynamic",
                 ylab="likelihood value by permutation",
                 main="random shuffer 1000 times"
                 );
  
  
  dev.off();
}


#  t_alpha_v<-0.1
compare<-function( t_alpha_v){
  
  t<-200
  adj_mat <- matrix(nrow = t,ncol = t);
  
  colnames(adj_mat)<-1:t
  rownames(adj_mat)<-1:t
  
  
  for( i in 1:nrow(adj_mat)){
    for( j in 1:ncol(adj_mat)){
      if(i!=j){
        #p<-runif(1)
        adj_mat[i,j]<-rpois(1,20);
        adj_mat[j,i]<-rpois(1,15);
      }
      if(i==j){
        adj_mat[i,j]<-0;
      }
      
    }
  }
  
  
  t_adj_mat2<-adj_mat[1:t,1:t, drop=FALSE];
  
  t_adj_mat3<-t_adj_mat2;
  
  for(i in 1:nrow(t_adj_mat3) ){
    for(j in 1:ncol(t_adj_mat3) ){
      
      if(i>j){
        
        p_minus<-t_adj_mat3[i,j]-t_adj_mat3[j,i];
        
        p_minus_rev<-t_adj_mat3[j,i]-t_adj_mat3[i,j];
        
        t_adj_mat3[i,j]<-p_minus;
        
        t_adj_mat3[j,i]<-p_minus_rev;
      }
    }
  }
  
  path1<-del_find_path_iter(t_adj_mat2,t_alpha_v);
  
  #path3<-del_find_path_iter2(t_adj_mat2,t_alpha_v);
  
  
  path2<-find_best_order_full2(t_adj_mat2,colnames(t_adj_mat2),t_alpha_v);
  
  
  print(path1$best_order);
  dyna<-calp2(t_adj_mat2,path1$best_order,t_alpha_v);
  
  
  #print(path3$best_order)
  #calp2(t_adj_mat2,path3$best_order,0.05);
  
  
  print(path2$best_order)
  opti<-calp2(t_adj_mat2,path2$best_order,t_alpha_v);
  
  return( c( dyna,opti) );
  
}


