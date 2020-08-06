//[[Rcpp::plugins(cpp11)]]

//C++ implementation of dynamic
#include <map>
#include <vector>
#include <set>
#include <math.h>       /* log */
#include <unordered_map>
#include <algorithm>
#include <iterator>
#include <Rcpp.h>


double calp2_c(double **t_read_count_mat_li,std::vector<int>& order_arr){
  
  double p_sum=0;
  int dim=order_arr.size();
  
  for(int i=0;i<dim-1;i++){
    for(int j=i+1;j<dim;j++){
      p_sum=p_sum+t_read_count_mat_li[order_arr[i]][order_arr[j]];
    }
    
  }
  
  return p_sum;
}


double calp2_p_v(double **t_read_count_mat_li, int t_dim, int sim_times=1000000){
  
  std::vector<int> init_order_p_v(t_dim);
  
  for(int i = 0; i <t_dim; i++)
    init_order_p_v[i] = i;
  
  
  
  double p_v=0;
  
  std::sort( init_order_p_v.begin(),init_order_p_v.end() ); 
  
  float in_order_li=calp2_c(t_read_count_mat_li,init_order_p_v);
  
  
  double sum_of_less_than_in_order=0;
  double sim_sum=0;
  
  for(int i=0;i< sim_times;i++){
    std::random_shuffle ( init_order_p_v.begin(), init_order_p_v.end() );
    
    double tmp_li=calp2_c(t_read_count_mat_li,init_order_p_v);
    
    sim_sum+= std::exp(tmp_li);
    
    if(tmp_li<=in_order_li ){
      sum_of_less_than_in_order+= std::exp(tmp_li);
    }
    
  }
  //number_of_less_than_in_order
  
  double permut_p=sum_of_less_than_in_order/sim_sum;
  
  return permut_p;
}




std::vector<int> hill_iter_c(double **t_read_count_mat_li,std::vector<int> init_order, std::set<int> full_order){
  
  double best_score=INT_MIN;
  std::vector<int> best_order;
  
  std::set<int> candidates = full_order;
  std::set<int>::iterator it;
  std::vector<int>::iterator it_v;
  
  //candidates.erase(init_order.begin(),init_order.end());
  for(it_v=init_order.begin();it_v!=init_order.end();it_v++){
    candidates.erase(*it_v);
  }
  
  
  for(it=candidates.begin();it!=candidates.end();it++){
    
    for( int j=0;j<(init_order.size()+1);j++ ){
      std::vector<int> tmp_order=init_order;
      tmp_order.insert(tmp_order.begin()+j,*it);
      
      double tmp_score=calp2_c(t_read_count_mat_li,tmp_order);
      
      if(tmp_score>best_score){
        
        best_score=tmp_score;
        best_order=tmp_order;
      }
      
    }
  }
  
  return best_order;
}


// [[Rcpp::export]]
Rcpp::List hill_c(Rcpp::NumericMatrix t_read_count_mat,double t_alpha_v=0.1){
  int dim = t_read_count_mat.nrow();
  
  double **read_count_mat_li=new double *[dim];
  for(int i = 0; i <dim; i++)
    read_count_mat_li[i] = new double[dim];
  
  
  for(int i=0;i<dim;i++){
    for(int j=0;j<dim;j++){
      if(i==j){
        read_count_mat_li[i][j]=read_count_mat_li[j][i]=0;
        continue;
      }
      
      if( (t_read_count_mat(i,j)+t_read_count_mat(j,i) )==0 ) {
        read_count_mat_li[i][j]=read_count_mat_li[j][i]=std::log(0.5);
        continue;
      }
      
      if(i>j){
        double p= (t_read_count_mat(i,j)+t_alpha_v)/(t_read_count_mat(i,j)+t_read_count_mat(j,i)+2*t_alpha_v);
        
        read_count_mat_li[i][j]=std::log(p);
        read_count_mat_li[j][i]=std::log(1-p);
      }
      
    }
  }
  
  //std::vector<int> li;
  //li.reserve(std::);
  //std::cout<<"calculate log prob"<<std::endl;
  
  
  //https://www.geeksforgeeks.org/all-permutations-of-an-array-using-stl-in-c/
  std::vector<int> init_order(2);
  std::set<int> full_order;
  double max_li=INT_MIN;
  for(int i = 0; i <dim; i++)
    full_order.insert(i);
  
  
  for(int i=0;i<dim;i++){
    for(int j=0;j<dim;j++){
      
      if( (read_count_mat_li[i][j]>max_li) && (i!=j) ){
        max_li=read_count_mat_li[i][j];
        init_order[0]=i;
        init_order[1]=j;
      }
      
    }
  }
  
  //std::cout << "init" << std::endl;
  
  while (init_order.size() < dim ) {
    init_order=hill_iter_c(read_count_mat_li,init_order,full_order);
    
  }
  //std::cout << "iteration over" << std::endl;
  
  //full_order=init_order;
  max_li=calp2_c(read_count_mat_li,init_order);
  
  
  for(int i = 0; i <dim; i++)
    init_order[i] = init_order[i]+1;
  
  double permut_p=calp2_p_v(read_count_mat_li,dim);
  //double permut_p=0;
  
  Rcpp::List L = Rcpp::List::create(Rcpp::Named("best_order") = init_order ,
                                    Rcpp::_["entropy"]=NA_REAL,
                                    Rcpp::_["best_score"] = max_li,
                                    Rcpp::_["number_of_maximum_order"]=NA_INTEGER,
                                    Rcpp::_["permut_p"]=permut_p
                                    );
  
  return L; 
}




