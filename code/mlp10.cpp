// [[Rcpp::plugins(cpp11)]]


// this algorithm can be speed up using bit operation on set
//C++ implementation of dynamic
#include <map>
#include <vector>
#include <set>
#include <bitset>
#include <math.h>       /* log */
#include <unordered_map>
#include <Rcpp.h>

//using namespace Rcpp;
//using namespace std;



double calp2_c(float **t_read_count_mat_li,std::vector<int>& order_arr){
  
  double p_sum=0;
  int dim=order_arr.size();
  
  for(int i=0;i<dim-1;i++){
    for(int j=i+1;j<dim;j++){
      p_sum=p_sum+t_read_count_mat_li[order_arr[i]][order_arr[j]];
    }
    
  }
  
  return p_sum;
}


double calp2_p_v(float **t_read_count_mat_li, int t_dim, int sim_times=1000000){
  
  std::vector<int> init_order_p_v(t_dim);
  
  for(int i = 0; i <t_dim; i++)
    init_order_p_v[i] = i;
  
  
  
  //double p_v=0;
  
  std::sort( init_order_p_v.begin(),init_order_p_v.end() ); 
  
  float in_order_li=calp2_c(t_read_count_mat_li,init_order_p_v);
  
  
  double number_of_less_than_in_order=0;
  
  for(int i=0;i< sim_times;i++){
    std::random_shuffle ( init_order_p_v.begin(), init_order_p_v.end() );
    
    double tmp_li=calp2_c(t_read_count_mat_li,init_order_p_v);
    
    if(tmp_li<=in_order_li ){
      number_of_less_than_in_order++;
    }
    
  }
  //number_of_less_than_in_order
  
  double permut_p=number_of_less_than_in_order/sim_times;
  
  return permut_p;
}




struct CompareSet_bit{
  
  bool operator()(const std::bitset<32>& lhs, const std::bitset<32>& rhs) const
  {
    //if(lhs.size()!=rhs.size()){
    //  return lhs.size() < rhs.size();
    //}
    
    if(lhs.count()!=rhs.count()){
      return lhs.count()<rhs.count();
    }
    
    return lhs.to_ulong()<rhs.to_ulong();
    
  }
};

inline float get_log_sum_c_bit(int t_sink, std::bitset<32> t_parents,float **mat_li, int t_dim ){
  
  float t_log_sum=0;
  
  //std::set<int>::iterator it;
  //for (it = t_parents.begin(); it != t_parents.end(); ++it){
  //for (std::size_t i = 0; i < t_parents.size(); ++i) {
  for (std::size_t i = 0; i < t_dim; ++i) {
  
    if( t_parents.test(i)  && (mat_li[i][t_sink]!=0)) {
        //t_log_sum = t_log_sum + std::log(mat_li[i][t_sink]);
        t_log_sum = t_log_sum + (mat_li[i][t_sink]);
        
      
    }
  }
  
  return t_log_sum;
}


int combinations(int set, int at, int r, int n,std::set<std::bitset<32>,CompareSet_bit>& subsets) {
  //std::cout<<'tx';
  
  //std::set<std::bitset<32>,CompareSet_bit > subsets;
  // Return early if there are more elements left to select than what is available.
  int elementsLeftToPick = n - at;
  if (elementsLeftToPick < r) {return 0;}
  //std::cout<<'tx';
  // We selected 'r' elements so we found a valid subset!
  if (r == 0) {
    subsets.insert(std::bitset<32>(set));
    //std::cout<<'tx';
    
  } else {
    for (int i = at; i < n; i++) {
      // Try including this element
      //std::cout<<'tx';
      
      set ^= (1 << i);
      
      combinations(set, i + 1, r - 1, n,subsets);
      
      //subsets.insert(t_set.begin(),t_set.end());
      
      // Backtrack and try the instance where we did not include this element
      set ^= (1 << i);
    }
  }
}


std::set<std::bitset<32>,CompareSet_bit > get_all_combinations(int n) {
  std::set<std::bitset<32>,CompareSet_bit > subsets;
  
  for(int r=1; r<=n;r++){
    combinations(0, 0, r, n,subsets);
    //subsets.insert( t_set.begin(),t_set.end() );
  }
  return subsets;
}


std::vector<double> find_opti_dynam_c_bit(float **read_count_mat_li,double t_alpha_v, int dim){
  
  //double read_count_mat_li[dim][dim];
  
  //for(int i=0;i<dim;i++){
  //for(int j=0;j<dim;j++) {
  //std::cout << read_count_mat_li[i][j] << '\n';
  
  //}
  //}
  
  //std::vector<int> v_set(dim);
  
  //for(int i=0;i<v_set.size();i++){
  //  v_set[i]=i;
  //}
  
  std::set<std::bitset<32>,CompareSet_bit > v_set_power = get_all_combinations(dim);
  
  std::unordered_map<std::bitset<32>, float> scores;
  std::unordered_map<std::bitset<32>, int> sinks;
  
  scores.reserve(pow(2,dim));
  sinks.reserve(pow(2,dim));
  
  std::set<std::bitset<32>,CompareSet_bit >::iterator it;
  
  for (it = v_set_power.begin(); it != v_set_power.end(); ++it){
    
    std::bitset<32> one_sub=*it;
    if(one_sub.size()==0){
      continue;
    }
    
    //std::string w_set="";
    std::bitset<32> w_set=one_sub;
    //std::for_each(one_sub.begin(), one_sub.end(), [&](const int &piece){
    //  w_set =w_set+":"+ std::to_string(piece); });
    //for (std::size_t i = 0; i < bs.size(); ++i) {
    
    
    scores[w_set]=0.0;
    sinks[w_set]= -1;
    
    if(one_sub.count()>2){
      std::set<int>::iterator iit;
      
      //for (iit = one_sub.begin(); iit != one_sub.end(); ++iit){
      for (std::size_t i = 0; i < dim; ++i) {
        if(!one_sub.test(i)){
          continue;
        }
        //upvars<-cset_difference(i,j);
        //std::bitset<32> upvars;
        //upvars=one_sub;
        //upvars.erase(*iit);
        one_sub.reset(i);
        
        //std::string upvars_str="";
        //std::for_each(upvars.begin(), upvars.end(), [&](const int &piece){
        //  upvars_str =upvars_str+":"+ std::to_string(piece); });
        
        float skore=scores[one_sub];
        skore=skore+get_log_sum_c_bit(i, one_sub, read_count_mat_li,dim);
        
        if( (skore > scores[w_set]) || (sinks[w_set]== -1)  ){
          scores[w_set]=skore;
          sinks[w_set]=i;
          
        }
        one_sub.set(i);
        
        
      }
      
    }else if(one_sub.count()==2){
      int index_one=0;
      int index_two=0;
      
      bool a=true;
      //for (std::size_t i = 0; i < one_sub.size(); ++i) {
      for (std::size_t i = 0; i < dim; ++i) {
          
        if(a&one_sub.test(i)){
          index_one=i;
          a=false;
        }
        
        if((!a)&one_sub.test(i)){
          
          index_two=i;
        }
      }
      
      
      if(read_count_mat_li[index_one][index_two]>read_count_mat_li[index_two][index_one]){
        //scores[w_set]=std::log(read_count_mat_li[index_one][index_two]);
        scores[w_set]=(read_count_mat_li[index_one][index_two]);
        
        sinks[w_set]=index_two;
      }else{
        //scores[w_set]=std::log(read_count_mat_li[index_two][index_one]);
        scores[w_set]=(read_count_mat_li[index_two][index_one]);
        
        sinks[w_set]=index_one;
      }
      
    }else{
      scores[w_set]=0;
      //sinks[w_set]=*(one_sub.begin());
      //for (std::size_t i = 0; i < one_sub.size(); ++i) {
      for (std::size_t i = 0; i < dim; ++i) {
      
        if(one_sub.test(i)){
          sinks[w_set]=i;
        }
      }
      
    }
  }
  
  //std::string w_set="";
  //std::for_each(v_set.begin(), v_set.end(), [&](const int &piece){
  //  w_set =w_set+":"+ std::to_string(piece); });
  
  
  std::vector<double> ord(dim+1);
  std::bitset<32> left(pow(2,(dim))-1 );
  //for(int i=0;i<dim;i++){
  //  left.insert(i);
  //}
  
  std::bitset<32> full_set(pow(2,(dim))-1 );
  
  float best_socre=scores[full_set];
  
  
  for(int i=dim-1; i>=0;i-- ){
    //std::string upvars_str="";
    //std::for_each(left.begin(), left.end(), [&](const int &piece){
    //  upvars_str=upvars_str+":"+ std::to_string(piece); });
    
    //upvars_str=left;
    //ord[i]=sinks[upvars_str] ;
    ord[i]=sinks[left] ;
    
    //left<- setdiff(left,ord[i]);
    //left.erase(ord[i]);
    left.reset(ord[i]);
  }
  
  
  ord[dim]=best_socre;
  
  return ord;
}




// [[Rcpp::export]]
Rcpp::List find_opti_dynam_r_cpp_bit(Rcpp::NumericMatrix t_read_count_mat,double t_alpha_v=0.1){
  int dim = t_read_count_mat.nrow();
  
  //int dim=5;
  double **m_read_count_mat=new double *[dim];
  
  for(int i = 0; i <dim; i++)
    m_read_count_mat[i] = new double[dim];
  
  
  for(int i=0;i<dim;i++){
    for(int j=0;j<dim;j++){
      m_read_count_mat[i][j]= t_read_count_mat(i,j) ;
    }
  }
  
  
  float **m_read_count_mat_li=new float *[dim];
  for(int i = 0; i <dim; i++)
    m_read_count_mat_li[i] = new float[dim];
  
  
  for(int i=0;i<dim;i++){
    for(int j=0;j<dim;j++){
      if(i==j){
        m_read_count_mat_li[i][j]=m_read_count_mat_li[j][i]=0;
        continue;
      }
      
      if( (m_read_count_mat[i][j]+m_read_count_mat[j][i])==0){
        m_read_count_mat_li[i][j]=m_read_count_mat_li[j][i]=std::log(0.5);
        continue;
      }
      
      if(i>j){
        float p= (m_read_count_mat[i][j]+t_alpha_v)/(m_read_count_mat[i][j]+m_read_count_mat[j][i]+2*t_alpha_v);
        
        m_read_count_mat_li[i][j]=std::log(p);
        m_read_count_mat_li[j][i]=std::log(1-p);
      }
      
    }
  }
  
  
  std::vector<double> ord=find_opti_dynam_c_bit(m_read_count_mat_li,t_alpha_v,dim);
  Rcpp::NumericVector out(dim );
  
  for(int i=0;i<dim;i++){
    out[i]=ord[i]+1;
  }
  
  
  //std::cout<<dim;
  //double permut_p=calp2_p_v(m_read_count_mat_li,dim);
  double permut_p=0;
  
  Rcpp::List L = Rcpp::List::create(Rcpp::Named("best_order") = out ,
                                    Rcpp::_["best_score"] = ord[dim],
                                    Rcpp::_["number_of_maximum_order"]=NA_INTEGER,
                                    Rcpp::_["permut_p"]=permut_p
                                    );
  
  return L;
}



