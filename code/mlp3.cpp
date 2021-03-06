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

inline float calp2_c(float **t_read_count_mat_li,std::vector<int>& order_arr, int t_dim){
  
  float p_sum=0;
  //int dim=order_arr.size();
  
  for(int i=0;i<t_dim-1;i++){
    for(int j=i+1;j<t_dim;j++){
      p_sum=p_sum+t_read_count_mat_li[order_arr[i]][order_arr[j]];
    }
    
  }
  
  return p_sum;
}


// [[Rcpp::export]]
Rcpp::List find_best_order_full2_c(Rcpp::NumericMatrix t_read_count_mat,double t_alpha_v=0.1){
  int dim = t_read_count_mat.nrow();
  
  float **read_count_mat_li=new float *[dim];
  for(int i = 0; i <dim; i++)
    read_count_mat_li[i] = new float[dim];
  
  
  for(int i=0;i<dim;i++){
    for(int j=0;j<dim;j++){
      if(i==j){
        read_count_mat_li[i][j]=read_count_mat_li[j][i]=0;
        continue;
      }
      
      if( (t_read_count_mat(i,j)+t_read_count_mat(j,i) )==0 ){
        read_count_mat_li[i][j]=read_count_mat_li[j][i]=std::log(0.5);
        continue;
      }
      
      if(i>j){
        float p= (t_read_count_mat(i,j)+t_alpha_v)/(t_read_count_mat(i,j)+t_read_count_mat(j,i)+2*t_alpha_v);
        
        read_count_mat_li[i][j]=std::log(p);
        read_count_mat_li[j][i]=std::log(1-p);
      }
      
    }
  }
  
  std::vector<float> li(std::pow(2,dim));
  //li.reserve(std::);
  
  //https://www.geeksforgeeks.org/all-permutations-of-an-array-using-stl-in-c/
  std::vector<int> init_order(dim);
  std::vector<int> best_order(dim);
  float max_li=INT_MIN;
  
  for(int i = 0; i <dim; i++){
    init_order[i] = i;
  }
  
  float in_order_li=calp2_c(read_count_mat_li,init_order,dim);
  
  std::sort( init_order.begin(),init_order.end() ); 
  
  // Find all possible permutations 
  //std::cout << "Possible permutations are:\n"; 
  do {
    float tmp_li=calp2_c(read_count_mat_li,init_order,dim);
    li.push_back(tmp_li);
    
    if(tmp_li>max_li){
      max_li=tmp_li;
      std::copy( std::begin(init_order),std::end(init_order),best_order.begin() ) ;
      
    }
    
    //display(a, n); 
  } while ( std::next_permutation(init_order.begin(), init_order.end()) );
  
  
  int number_of_maximum_order=0;
  //std::count(li.begin(), li.end(), max_li);
  
  double sum_of_less_than_in_order=0;
  
  double entropy=0;
  
  double prob_sum=0;
  
  
  std::vector<float>::iterator li_d;
  
  for(li_d=li.begin();li_d!=li.end();li_d++){
    prob_sum+=std::exp(*li_d);
  }
  
  for(li_d=li.begin();li_d!=li.end();li_d++){
    if(std::abs(*li_d-max_li)<=std::numeric_limits<double>::epsilon() ){
      number_of_maximum_order++;
    }
    
    //std::cout << "probability:" << *li_d << std::endl; 
    double prob_one=std::exp(*li_d)/prob_sum;
    
    double entropy_one=prob_one*std::log2(prob_one);
    
    //std::cout << "entropy_one:" << entropy_one << std::endl; 
    
    entropy+=entropy_one;
    
    //in_order_li
    if(*li_d<=in_order_li ){
      //sum_of_less_than_in_order+= std::exp(*li_d);
      sum_of_less_than_in_order+= prob_one;
      
    }
    
  }
  
  double permut_p=sum_of_less_than_in_order;
  
  //double permut_p=sum_of_less_than_in_order/prob_sum;
  //double permut_p=0;
  
  for(int i = 0; i <dim; i++){
    best_order[i] = best_order[i]+1;
  }
  
  Rcpp::List L = Rcpp::List::create(Rcpp::Named("best_order") = best_order,
                                    Rcpp::_["entropy"]=-1*entropy,
                                    Rcpp::_["best_score"] = max_li,
                                    Rcpp::_["number_of_maximum_order"]=number_of_maximum_order,
                                    Rcpp::_["permut_p"]=permut_p
                                    );
  
 return L; 
}



