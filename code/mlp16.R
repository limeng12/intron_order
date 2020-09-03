library(lpSolve)
library(ggm)
library(gtools)
#library(Rglpk)


lp_kenemy<-function(t_adj_mat,t_alpha_v,verbose_hill=FALSE){
  
  
  t_adj_mat<-t_adj_mat+t_alpha_v;
  
  t_adj_mat_li<-t_adj_mat;
  for(i in 1:nrow(t_adj_mat_li)){
    for(j in 1:ncol(t_adj_mat_li)){
      if(i>j){
        p<-t_adj_mat_li[i,j]/(t_adj_mat_li[i,j]+t_adj_mat_li[j,i])
        
        t_adj_mat_li[i,j]<-log(p)
        t_adj_mat_li[j,i]<-log(1-p)
      }
      
      if(i==j){
        t_adj_mat_li[i,j]<-0;
      }
      
    }
    
  }
  
  n<-nrow(t_adj_mat_li);
  
  f.obj<-as.vector(t(t_adj_mat_li));
  
  constraints1<-matrix(nrow=n*(n-1)/2,ncol = n*n)
  t<-1
  #one_constraint<-rep(0,n*n);
  for(i in 1:nrow(t_adj_mat_li) ){
    for(j in 1:ncol(t_adj_mat_li)){
      if(i>j){
        #one_constraint<-rep(0,n*n);
        
        constraints1[t, (i-1)*n+j ]<-1;
        constraints1[t, (j-1)*n+i ]<-1;
        
        #constraints1[t,]<-one_constraint;
        t<-t+1        
      }
      
    }
    
  }
  
  
  
  f.dir<-rep("==", nrow(constraints1) );
  
  
  all_combns<-combn(n, 3)
  
  constraints2<-matrix(0, nrow=ncol(all_combns)*6,ncol = n*n);
  t<-1
  
  
  for(i in 1:ncol(all_combns)){
    all_per<-permutations(3,3,all_combns[,i])
    for(j in 1:nrow(all_per) ){
      #one_constraint<-rep(0,n*n);
      
      constraints2[t, (all_per[j,1]-1)*n+all_per[j,2] ]<-1
      constraints2[t, (all_per[j,2]-1)*n+all_per[j,3] ]<-1
      constraints2[t, (all_per[j,3]-1)*n+all_per[j,1] ]<-1
      
      
      #constraints2[t,]<-one_constraint;
      t<-t+1
      
    }
    
  }
  
  f.con<-rbind(constraints1,constraints2);
  
  f.dir<-c(f.dir, rep(">=", nrow(constraints2) ) );
  
  rm(constraints1,constraints2);
  
  
  f.rhs<-rep(1,length(f.dir) );
  
  #bounds<-list(lower=list(ind=1:(n*n),val=rep(0,n*n)),upper=list(ind=1:(n*n),val=rep(1,n*n) ));
  print("Got all coefficients")
  max <- TRUE
  a<-lp ("max", f.obj, f.con, f.dir, f.rhs, all.int=TRUE, all.bin=TRUE);
  
  #a<-Rglpk_solve_LP(obj=f.obj, mat=f.con, dir=f.dir, rhs=f.rhs, max = max,bounds=bounds, types="B",
  #                  control = list("verbose" =TRUE, "canonicalize_status" = FALSE) )
  
  #, tm_limit=100000
  best_order<-topOrder(matrix(a$solution,byrow = TRUE,nrow = n) );
  
  permut_p<-NA
  entropy<-NA
  
  list(best_order=best_order, permut_p=permut_p, entropy=entropy);
  
}



get_test_matrix<-function(){
  
  t<-12
  t_alpha_v<-0.1
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
  t_adj_mat<-t_adj_mat2
  
  
  lp_kenemy(t_adj_mat,0.01);
  find_opti_dynam_r_cpp_bit(t_adj_mat, t_alpha_v);
  
  
}

