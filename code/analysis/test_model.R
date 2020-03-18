library(ggplot2)
library(Cairo)
#setwd("/Users/mengli/Documents/projects/iso");
source("code/utils.R")
source("code/mlp3.R")

generate_read_count_mat<-function(t_total_number_intron=50,random_noise=TRUE){
  
  t<-t_total_number_intron;
  t_alpha_v<-0.1;
  
  given_orders<-sample(1:t);
  
  t_adj_mat2<-matrix(nrow=t, ncol=t);
  t_adj_mat2<-apply(t_adj_mat2,c(1,2),function(x){return(0)});
  colnames(t_adj_mat2)<-1:t
  rownames(t_adj_mat2)<-1:t
  
  for(i in 1:(length(given_orders)-1) ){
    for(j in (i+1):length(given_orders)){
      if(i!=j){
        t_adj_mat2[given_orders[i],given_orders[j]]<-max(1,rpois(1,15) );
      }
      
    }
  }
  
  if(random_noise){
    for(i in 1:nrow(t_adj_mat2) ){
      for(j in 1:ncol(t_adj_mat2)){
        
        if(i!=j){
          t_adj_mat2[i,j]<-t_adj_mat2[i,j]+rpois(1,3);
          
        }
      }
      
    }
    
  }
  
  #path3<-hill(t_adj_mat2,t_alpha_v);
  
  #hill_v<-calp2(t_adj_mat2,path3$best_order,t_alpha_v);
  
  #a<-cor(given_orders,as.numeric( path3$best_order) )
  list(mat=t_adj_mat2,given_orders=given_orders)
  #a
}


test_test_large_intron_given_order_times<-function(t_sim_matrix_number){
  
  ##average number of intron in each transcript
  total_number_intron<-7;
  min_total_number_intron<-2;
  #total_number_intron<-3;
  #min_total_number_intron<-2;
  
  
  #sim_times<-50;
  sim_matrix_number<-t_sim_matrix_number;
  t_alpha_v<-0.1;
  
  
  #######generate read count matrix
  all_read_mat<-vector("list", sim_matrix_number);
  
  for( k in 1:sim_matrix_number){
    
    gen_mat_obj<-generate_read_count_mat( max(rpois(1,total_number_intron), min_total_number_intron),FALSE );
    
    all_read_mat[[k]]$mat<-gen_mat_obj$mat;
    all_read_mat[[k]]$given_orders<-gen_mat_obj$given_orders;
    
  }
  
  
  all_read_mat_backup<-all_read_mat
  ######impute the adjacent matrix into 0 to about 100%, 80%, 65%, 50%, 
  t_intron_pair_cov_threshod_v<-c(100,99,95,90,80,70,60);
  
  test_cors_mat<-matrix(nrow=0,ncol=2)
  
  for(t_intron_pair_cov_threshod in t_intron_pair_cov_threshod_v){
    all_read_mat<-all_read_mat_backup
    
    #  t_intron_pair_cov_threshod<-80
    ###  random erase percent of intron splicing order pairs
    for(k in 1:length(all_read_mat) ){
      m_mat<-all_read_mat[[k]]$mat
      
      mat_index_pos<-0
      back_up_value<-0;
      
      while( cal_intron_pair_cov(m_mat) > t_intron_pair_cov_threshod/100 ){
        
        mat_index_pos<-sample(1:(nrow(m_mat)*nrow(m_mat) ) ,1);
        
        back_up_value<-m_mat[ mat_index_pos ];
        m_mat[ mat_index_pos ]<-0;
        
      }
      
      m_mat[ mat_index_pos ]<-back_up_value;
      
      all_read_mat[[k]]$mat<-m_mat;
    }
    
    
    ######
    ##calculate correlation between given orders and calculated orders
    test_cors<-c();
    
    for( k in 1:length(all_read_mat)){
      print(paste0("Curent kept > %: ", t_intron_pair_cov_threshod,"   mat: ", k));
      
     # print(paste0("Curent kept %: ", k));
      
      #path3<-hill(all_read_mat[[k]]$mat,t_alpha_v);
      path3<-find_path_global(all_read_mat[[k]]$mat,t_alpha_v);
      #path3<-del_find_path_iter(all_read_mat[[k]]$mat, t_alpha_v);
      
      #hill_v<-calp2(t_adj_mat2,path3$best_order,t_alpha_v);
      
      a<-cor(all_read_mat[[k]]$given_orders,as.numeric( path3$best_order) ,method="spearman");
      
      test_cors<-c(test_cors,a );
      
    }
    
    test_cors_mat<-rbind(test_cors_mat,cbind(rep( paste0("\u2265",t_intron_pair_cov_threshod,"%") ,
                                              length(test_cors)),test_cors))
    
  }
  
  test_cors_mat_fr<-data.frame(label=test_cors_mat[,1], 
                               cors=as.numeric(as.character(test_cors_mat[,2]) ) ,
                               stringsAsFactors = FALSE);
  
  test_cors_mat_fr$label <- factor(test_cors_mat_fr$label,
                                        levels = c(paste0("\u2265",t_intron_pair_cov_threshod_v,"%")  ) );
  
  CairoPDF("result/test_model_large_intron_hill_erase_percent.pdf", width=8, height=4);
  
  p<-ggplot(test_cors_mat_fr)+geom_violin(aes(x="",y=cors, fill=label))+theme_minimal()+
    ylim(-1.1,1.1)+ylab("")+xlab("")+
    facet_grid(cols=vars(label), scales = "free", space = "free")+theme(legend.position="none")+
    ggtitle(paste0("Correlation between simulated orders and calculated orders: n=1,000") )+
    theme(legend.position="none",text = element_text(size=8));
  
  # ,length(all_read_mat)
  
  print(p);
  dev.off();
  
}

#p<-p+

##generate read count matrix with given orders (add some poisson noise)
##try to prove as long as the matrix is correct, the model result will be correct, even for large number of introns
#test_test_large_intron_given_order_times(1000);

