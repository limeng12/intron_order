library(gtools);
library(stringr);
library(ggplot2);
library(tictoc)

#library(RAppArmor)

setwd("/Users/mengli/Documents/projects/iso");
#setrlimit(2)

#hill
#source("code/mlp2.R")
#permutation
source("code/mlp3.R")
#dynamic
source("code/mlp8.R")
#dynamic CPP
#sourceCpp("code/mlp9.cpp");


# scater_compare()
scater_compare_simulate<-function(){
  dyn<-c();
  
  opti<-c();
  
  hill_v<-c();
  
  dyn_cpp<-c();
  
  for( i in 1:1000){
    
    if(i %% 10 ==0){
      print(paste0("cur= ",i) )
    }
    
    a<-compare(0.1);
    
    dyn<-c(dyn,a[1]);
    
    opti<-c(opti,a[2]);
    
    hill_v<-c(hill_v,a[3]);
    
    dyn_cpp<-c(dyn_cpp,a[4]);
  }
  
  pdf("result/compare_hill_permut_simulation.pdf", width=8, height=8);
  
  plot(hill_v,opti,xlab="Likelihood values by hill-climbing",
       ylab="Likelihood values by permutation"
       ,main="Simulation data (n=1,000)"
  );
  
  # plot(dyn,opti,xlab="Likelihood values by dynamic programming ",
  #      ylab="Likelihood values by permutation"
  #      ,main="Simulation data (n=1,000)"
  # );
  
  plot(dyn_cpp,opti,xlab="Likelihood values by dynamic programming CPP",
       ylab="Likelihood values by permutation"
       ,main="Simulation data (n=1,000)"
  );
  
  
  dev.off();
}


#  t_alpha_v<-0.1
compare<-function( t_alpha_v){
  
  t<-sample(6:11,1)
  #t<-8
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
  
  
  #path1<-find_opti_dynam(t_adj_mat2,t_alpha_v);
  
  #print(paste0("permutation: ",toc()) )
  path2<-find_best_order_full2_c(t_adj_mat2,t_alpha_v);
  
  #print(paste0("hill climbing: ",toc()) )
  path3<-hill_c(t_adj_mat2,t_alpha_v);
  
  #print(path1$best_order);
  path4<-find_opti_dynam_r_cpp(t_adj_mat2,t_alpha_v);
  
  
  #dyna<-calp2(t_adj_mat2,path1$best_order,t_alpha_v);
  
  opti<-calp2(t_adj_mat2,path2$best_order,t_alpha_v);
  
  hill_v<-calp2(t_adj_mat2,path3$best_order,t_alpha_v);
  
  dyn_cpp<-calp2(t_adj_mat2,path4$best_order,t_alpha_v);

  
  return( c( 0,opti,hill_v,dyn_cpp) );
  
}


#   scater_compare_real_mat()
scater_compare_real_mat<-function(){
  dyn<-c()
  
  dyn_cpp<-c()
  
  opti<-c();
  
  hill_v<-c();
  
  mat_files<-list.files("./result/adj_matrix/");
  
  n<-1;
  
  for( i in 1:length(mat_files) ){
    
    mat_file<-mat_files[i];
    
    mat<-read.table(paste0("./result/adj_matrix/", mat_file),header = FALSE );
    if(nrow(mat)<6 || nrow(mat)>11 ){
      next;
    }
    
    if(n>1000){
      break;
    }
    
    print(paste0("current matrix file: ", n));
    
    n<-n+1
    
    a<-compare_real_mat(mat_file,0.1);
    
    dyn<-c(dyn,a[1]);
    
    opti<-c(opti,a[2]);
    
    hill_v<-c(hill_v,a[3]);
    
    dyn_cpp<-c(dyn_cpp,a[4])
    
  }
  
  
  pdf("result/compare_hill_permut_human_mat.pdf", width=8, height=8);
  
  #par(mfrow=c(1,2))
  
  #plot(dyn,opti,xlab="Likelihood values by algorithm developed here",
  #               ylab="Likelihood values by permutation"
  #,main="Compare 100 times based on random generated matrix"
  #);
  
  plot(hill_v,opti,xlab="Likelihood values by hill-climbing ",
       ylab="Likelihood values by permutation"
       ,main=paste0("Real data (n=",n,")" )
  );
  
  # plot(dyn,opti,xlab="Likelihood values by dynamic programming ",
  #      ylab="Likelihood values by permutation"
  #      ,main=paste0("Real data (n=",n,")" )
  # );
  
  plot(dyn_cpp,opti,xlab="Likelihood values by dynamic programming CPP",
       ylab="Likelihood values by permutation"
       ,main=paste0("Real data (n=",n,")" )
  );
  
  dev.off();
}



compare_real_mat<-function(mat_file, t_alpha_v=0.1){
  
  #t<-7
  
  mat<-read.table(paste0("./result/adj_matrix/", mat_file),header = FALSE );
  
  m_output_id<-str_split(mat_file,"\\.")[[1]][1]
            
  colnames(mat)<-1:ncol(mat)
  
  rownames(mat)<-1:nrow(mat)
  
  t_adj_mat2<-as.matrix(mat);
  
  #print(paste0("dynamic: ",toc()) )
  
  
  
  #path1<-find_opti_dynam(t_adj_mat2,t_alpha_v);
  
  #print(paste0("permutation: ",toc()) );
  path2<-find_best_order_full2_c(t_adj_mat2,t_alpha_v);
  
  #print(paste0("hill climbing: ",toc()) );
  path3<-hill_c(t_adj_mat2,t_alpha_v);
  
  path4<-find_opti_dynam_r_cpp(t_adj_mat2,t_alpha_v);
  

  #dyna<-calp2(t_adj_mat2,path1$best_order,t_alpha_v);
  
  opti<-calp2(t_adj_mat2,path2$best_order,t_alpha_v);
  
  hill_v<-calp2(t_adj_mat2,path3$best_order,t_alpha_v);
  
  dyna_cpp<-calp2(t_adj_mat2,path4$best_order,t_alpha_v);
  
  return( c( 0,opti,hill_v,dyna_cpp) );
  
}



##generate simulated read count matrix with random read count
##try to prove the hill and permutation are almost same  
scater_compare_simulate();

##fill read count matrix with human real dataset, intron number <8 and >6
scater_compare_real_mat();

