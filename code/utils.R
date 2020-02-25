
if_smaller_f<-function(region1,region2){
  sapply(strsplit(region1,":|-"),"[",2)->region1_start;
  sapply(strsplit(region2,":|-"),"[",2)->region2_start;
  
  region1_start<region2_start;
}


strReverse <- function(x)
  sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")

cal_intron_pair_cov<-function(t_mat){
  count_of_nonzero_ele<-0;
  
  for(i in 1:nrow(t_mat)){
    for(j in 1:ncol(t_mat)){
      if( (t_mat[i,j]+t_mat[j,i]>0) && (i>j) ){
        count_of_nonzero_ele<-count_of_nonzero_ele+1;
      }
      
    }
    
  }
  
  #number of introns
  n<-nrow(t_mat)
  percent_coverage_pair<-count_of_nonzero_ele/( n*(n-1)/2+0.0 );
  
  percent_coverage_pair
}
