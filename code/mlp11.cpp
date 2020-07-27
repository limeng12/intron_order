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

//bigger 32 sotre child, lower 32 store parents
std::unordered_map<int, std::unordered_map<std::bitset<32>, float> > build_parents_tbl(int t_dim,float **mat_li){
  
  //std::unordered_map<std::bitset<64>, float> parents;
  
  std::unordered_map<int, std::unordered_map<std::bitset<32>, float> > parents;
  parents.reserve(t_dim);
  for (std::size_t i = 0; i < t_dim; ++i) {
    (parents[i]).reserve(std::pow(2,t_dim));
  }
  //parents.reserve(t_dim*pow(2,t_dim));

  std::vector<std::bitset<32> > sets;
  
  for (std::size_t i = 0; i < t_dim; ++i) {
    //i is the child
    
    sets.clear();
    sets.reserve(pow(2,t_dim));
    //std::list<std::bitset<32> >::iterator it;
    
    //std::bitset<64> one_set_64;one_set_64.set(32+i);
    
    for (std::size_t j = 0; j < t_dim; ++j) {
      if(j==i){
        continue;
      }
      
      //for (it = sets.begin(); it != sets.end(); ++it) {
      int cur_size_vec=sets.size();
      for(std::size_t t=0;t<cur_size_vec;t++){
        
        std::bitset<32> t_one_set_32( sets.at(t) );
        if(t_one_set_32.test(j)){
          continue;
        }
        
        float score=(parents[i])[t_one_set_32];
        
        t_one_set_32.set(j);
        //sets.insert( t_one_set_32 );
        sets.push_back(t_one_set_32);
        
        (parents[i])[t_one_set_32]=score+mat_li[j][i];
        
        //std::bitset<64> t_one_set_64;t_one_set_64.set(32+i);
        //t_one_set_64=t_one_set_64.to_ulong()+t_one_set_32.reset(j).to_ulong();
        //std::bitset<64> t_one_set_64_addj(t_one_set_64);t_one_set_64_addj.set(j);
        
        //parents[t_one_set_64_addj]=parents[t_one_set_64]+mat_li[j][i];
      }
      std::bitset<32> t_one_set_32;
      
      t_one_set_32.set(j);
      //sets.insert(one_set_32);
      sets.push_back(t_one_set_32);
      (parents[i])[t_one_set_32]=mat_li[j][i];
      
      //std::bitset<64> t_one_set_64;t_one_set_64.set(32+i);
      //t_one_set_64.set(j);
      
    }
    
  }
  
  std::cout<<"build parent table"<<std::endl;
  return parents;
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


std::vector<double> find_opti_dynam_c_bit(double **read_count_mat,double t_alpha_v, int dim){
  
  //double read_count_mat_li[dim][dim];
  
  float **read_count_mat_li=new float *[dim];
  
  for(int i = 0; i <dim; i++)
    read_count_mat_li[i] = new float[dim];
  
  
  for(int i=0;i<dim;i++){
    for(int j=0;j<dim;j++){
      if(i==j){
        read_count_mat_li[i][j]=read_count_mat_li[j][i]=0;
        continue;
      }
      
      if( (read_count_mat[i][j]+read_count_mat[j][i])==0){
        read_count_mat_li[i][j]=read_count_mat_li[j][i]=std::log(0.5);
        continue;
      }
      
      if(i>j){
        float p= (read_count_mat[i][j]+t_alpha_v)/(read_count_mat[i][j]+read_count_mat[j][i]+2*t_alpha_v);
        
        read_count_mat_li[i][j]=std::log(p);
        read_count_mat_li[j][i]=std::log(1-p);
      }
      
    }
  }
  
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
  
  std::unordered_map<int, std::unordered_map<std::bitset<32>, float> > parents_score=build_parents_tbl(dim,read_count_mat_li);
  
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
        //skore=skore+get_log_sum_c_bit(i, one_sub, read_count_mat_li,dim);
        //std::bitset<64> one_child_parents;one_child_parents.set(i+32);
        //one_child_parents=one_child_parents.to_ulong()+one_sub.to_ulong();
        
        skore=skore+(parents_score[i])[one_sub];
        
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
  
  std::vector<double> ord=find_opti_dynam_c_bit(m_read_count_mat,t_alpha_v,dim);
  Rcpp::NumericVector out(dim );
  
  for(int i=0;i<dim;i++){
    out[i]=ord[i]+1;
  }
  //std::cout<<dim;
  
  Rcpp::List L = Rcpp::List::create(Rcpp::Named("best_order") = out ,
                                    Rcpp::_["best_score"] = ord[dim],
                                                               Rcpp::_["number_of_maximum_order"]=NA_INTEGER);
  
  return L;
}



