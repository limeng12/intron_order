// [[Rcpp::plugins(cpp11)]]

//C++ implementation of dynamic
#include <map>
#include <vector>
#include <set>
#include <math.h>       /* log */
#include <unordered_map>
#include <Rcpp.h>

//using namespace Rcpp;
//using namespace std;

struct CompareSet{
  
  bool operator()(const std::set<int>& lhs, const std::set<int>& rhs) const
  {
    if(lhs.size()!=rhs.size()){
      return lhs.size() < rhs.size();
    }
    
    return lhs<rhs;
    
  }
};

inline double get_log_sum_c(int t_sink, std::set<int>& t_parents,double **mat_li ){
  
  double t_log_sum=0;
  
  std::set<int>::iterator it;
  for (it = t_parents.begin(); it != t_parents.end(); ++it){
    //if(mat_li[*it][t_sink]!=0)
      t_log_sum=t_log_sum + std::log(mat_li[*it][t_sink]);
    
  }
  
  return t_log_sum;
}

// https://www.geeksforgeeks.org/find-distinct-subsets-given-set/
std::vector<std::set<int> >  get_power_set(std::vector<int> arr) {
  int n=arr.size();
  
  std::vector<std::set<int> > list;
  list.reserve(std::pow(2,n));

  
  /* Run counter i from 000..0 to 111..1*/
  for (int i = 0; i < (int) pow(2, n); i++)
  {
    std::set<int> subset;
    
    // consider each element in the set
    for (int j = 0; j < n; j++)
    {
      // Check if jth bit in the i is set. If the bit
      // is set, we consider jth element from set
      if ((i & (1 << j)) != 0)
        //subset += to_string(arr[j]) + "|";
        subset.insert(arr[j]);
    }
    
    // if subset is encountered for the first time
    // If we use set<string>, we can directly insert
    //if (find(list.begin(), list.end(), subset) == list.end())
    list.push_back(subset);
  }
  return list;
}


std::vector<double> find_opti_dynam_c(double **read_count_mat,double t_alpha_v, int dim){
  
  //double read_count_mat_li[dim][dim];
  
  double **read_count_mat_li=new double *[dim];
  double **read_count_mat_li_log=new double *[dim];
  
  for(int i = 0; i <dim; i++){
    read_count_mat_li[i] = new double[dim];
    read_count_mat_li_log[i] = new double[dim];
  }
  for(int i=0;i<dim;i++){
    for(int j=0;j<dim;j++){
      if(i==j){
        read_count_mat_li[i][j]=read_count_mat_li[j][i]=0;
        read_count_mat_li_log[i][j]=read_count_mat_li_log[j][i]=0;
        
        continue;
      }
      
      if( (read_count_mat[i][j]+read_count_mat[j][i])==0){
        read_count_mat_li[i][j]=read_count_mat_li[j][i]=0.5;
        read_count_mat_li_log[i][j]=read_count_mat_li_log[j][i]=std::log(0.5);
        
        continue;
      }
      
      if(i>j){
        double p= (read_count_mat[i][j]+t_alpha_v)/(read_count_mat[i][j]+read_count_mat[j][i]+2*t_alpha_v);
        
        read_count_mat_li[i][j]=p;
        read_count_mat_li[j][i]=1-p;
        
        read_count_mat_li_log[i][j]=std::log(p);
        read_count_mat_li_log[j][i]=std::log(1-p);
      }
      
    }
  }
  
  //for(int i=0;i<dim;i++){
    //for(int j=0;j<dim;j++) {
      //std::cout << read_count_mat_li[i][j] << '\n';
      
    //}
  //}
  
  std::vector<int> v_set(dim);
  
  for(int i=0;i<v_set.size();i++){
    v_set[i]=i;
  }
  
  std::vector<std::set<int> > v_set_power = get_power_set(v_set);
  
  //std::sort(v_set_power.begin(), v_set_power.end(), CompareSet() );
  //v_set_power.erase( unique( v_set_power.begin(), v_set_power.end() ), v_set_power.end() );
  
  std::set<std::set<int>,CompareSet> v_set_power_set( v_set_power.begin(), v_set_power.end() );
  //vec.assign( s.begin(), s.end() );
  
  
  std::unordered_map<std::string, double> scores;
  std::unordered_map<std::string, int> sinks;
  
  scores.reserve(pow(2,dim));
  sinks.reserve(pow(2,dim));
  
  std::set<std::set<int> >::iterator it;
  
  for (it = v_set_power_set.begin(); it != v_set_power_set.end(); ++it){
    
    std::set<int> one_sub=*it;
    if(one_sub.size()==0){
      continue;
    }
    
    std::string w_set="";
    
    std::for_each(one_sub.begin(), one_sub.end(), [&](const int &piece){
      w_set =w_set+":"+ std::to_string(piece); });
    
    
    scores[w_set]=0.0;
    sinks[w_set]= -1;
    
    if(one_sub.size()>2){
      std::set<int>::iterator iit;
      
      for (iit = one_sub.begin(); iit != one_sub.end(); ++iit){
        //upvars<-cset_difference(i,j);
        //std::set<int> upvars;
        //upvars=one_sub;
        //upvars.erase(*iit);
        
        //double log_sum=0;
        double skore=0;
        
        std::string upvars_str="";
        std::for_each(one_sub.begin(), one_sub.end(), [&](const int &piece){
          if(piece!=*iit){
            upvars_str =upvars_str+":"+ std::to_string(piece);
            skore=skore+(read_count_mat_li_log[piece][*iit]);
          }
          });
        
        skore=skore+scores[upvars_str];
        //double skore=scores[upvars_str]+log_sum;
        //skore=skore+get_log_sum_c(*iit, upvars, read_count_mat_li);
        
        if( (sinks[w_set]== -1) || (skore > scores[w_set])){
          scores[w_set]=skore;
          sinks[w_set]=*iit;
          
        }
        
      }
      
    }else if(one_sub.size()==2){
      
      int index_one=*(one_sub.begin());
      
      int index_two=*(one_sub.rbegin());
      
      if(read_count_mat_li_log[index_one][index_two]>read_count_mat_li_log[index_two][index_one]){
        scores[w_set]=(read_count_mat_li_log[index_one][index_two]);
        sinks[w_set]=index_two;
      }else{
        scores[w_set]=(read_count_mat_li_log[index_two][index_one]);
        sinks[w_set]=index_one;
      }
      
    }else{
      scores[w_set]=0;
      sinks[w_set]=*(one_sub.begin());
    }
    
  }
  
  
  std::string w_set="";
  std::for_each(v_set.begin(), v_set.end(), [&](const int &piece){
    w_set =w_set+":"+ std::to_string(piece); });
  
  
  std::vector<double> ord(dim+1);
  std::set<int> left;
  for(int i=0;i<dim;i++){
    left.insert(i);
  }
  
  double best_socre=scores[w_set];
  
  for(int i=dim-1; i>=0;i-- ){
    std::string upvars_str="";
    std::for_each(left.begin(), left.end(), [&](const int &piece){
      upvars_str=upvars_str+":"+ std::to_string(piece); });
    
    ord[i]=sinks[upvars_str] ;
    //left<- setdiff(left,ord[i]);
    left.erase(ord[i]);
  }
  
  ord[dim]=best_socre;
  
  return ord;
}


// [[Rcpp::export]]
Rcpp::List find_opti_dynam_r_cpp(Rcpp::NumericMatrix t_read_count_mat,double t_alpha_v=0.1){
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
  
  std::vector<double> ord=find_opti_dynam_c(m_read_count_mat,t_alpha_v,dim);
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



