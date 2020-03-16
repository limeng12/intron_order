
library(dbscan);


get_members<-function(t_igraph_list,t_alpha){
  
  for( i in 1:length(t_igraph_list) ){
    mermbers<-get_members_matrix(t_igraph_list[[i]]$adjacency_matrix,t_alpha);
    
    t_igraph_list[[i]]$members<-mermbers;
    
  }
  
  return(t_igraph_list);
}


get_members_matrix<-function(tt_adj_mat,t_alpha_v=0.05){
  
  if( nrow(tt_adj_mat) <=2 ){
    
    members<-rep(0,nrow(tt_adj_mat))
    names(members)<-colnames(tt_adj_mat);
    
    
    for(i in 1:length(members) ){
      
      if(members[i]==0){
        members[i]<-"not_in_unit(0)";
      }else{
        members[i]<-str_c("unit_",members[i]);
      }
      
      
    }
    
    return(members)
  }
  
  
  t_adj_mat_p<-tt_adj_mat+t_alpha_v;
  
  for( i in 1:nrow(t_adj_mat_p)){
    
    for(j in 1:ncol(t_adj_mat_p)){
      
      if(i==j){
        t_adj_mat_p[i,j]<-0.5;
        next;
      }
      if( (t_adj_mat_p[i,j]==0) && (t_adj_mat_p[j,i]==0) ){
        next;
      }
      
      if(i>j){
        
        t_adj_mat_p[i,j]<-(t_adj_mat_p[i,j])/(t_adj_mat_p[i,j]+t_adj_mat_p[j,i]);
        
        t_adj_mat_p[j,i]<-1-t_adj_mat_p[i,j]
      }
    }
    
  }
  
  
  #members<-dbscan(t_adj_mat_p,2)$cluster;
  #dis_mat<-(2-cor(t_adj_mat_p,method = "pearson",use = "pairwise.complete.obs") )
  #dis_mat[is.na(dis_mat)]<-3;
  #dist_obj<-as.dist( dis_mat );
  #a<-optics(dist_obj,minPts = 2);
  #sqrt( sum( (t_adj_mat_p[,1]-t_adj_mat_p[,2])^2 ) )
  #spearman corrlation of 0.6 as threshold
  #eps_t<-a$eps/2;
  
  
  a<-optics(t_adj_mat_p,minPts = 2);
  
  eps_t<-(2-0.8);
  
  if(!is.na(a$eps_cl) ){
    eps_t<-a$eps_cl
  }
  
  extractDBSCAN(a,eps_cl=eps_t )->b;
  
  members<-b$cluster;
  
  names(members)<-colnames(t_adj_mat_p);
  
  
  for(i in 1:length(members) ){
    
    if(members[i]==0){
      members[i]<-"not_in_unit(0)";
    }else{
      members[i]<-str_c("unit_",members[i]);
    }
    
  }
  
  members;
  
}



#   tt_adj_mat<-t_igraph_list[[1]]$adjacency_matrix;
#   t_alpha_v=0.001
#library(RColorBrewer)
#library(PReMiuM);
#cols1<-brewer.pal(12, name = 'Set3') ;

#cols2<-brewer.pal(8, name = 'Dark2') ;

#cols3<-brewer.pal(12, name = 'Paired') 

#cols<-c(cols1,cols2,cols3);

