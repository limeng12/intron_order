// [[Rcpp::plugins(cpp11)]]


// this algorithm can be speed up using bit operation on set
//C++ implementation of dynamic
#include <cstdint>
#include <map>
#include <vector>
#include <set>
#include <bitset>
#include <math.h>       /* log */
#include <unordered_map>
#include <Rcpp.h>

//using namespace Rcpp;
//using namespace std;

//https://stackoverflow.com/questions/47981/how-do-you-set-clear-and-toggle-a-single-bit
/* a=target variable, b=bit number to act upon 0-n */
#define BIT_SET(a,b) ((a) |= (1ULL<<(b)))
#define BIT_CLEAR(a,b) ((a) &= ~(1ULL<<(b)))
#define BIT_FLIP(a,b) ((a) ^= (1ULL<<(b)))
#define BIT_CHECK(a,b) (!!((a) & (1ULL<<(b))))        // '!!' to make sure this returns 0 or 1


//https://stackoverflow.com/questions/109023/how-to-count-the-number-of-set-bits-in-a-32-bit-integer
inline int numberOfSetBits(uint32_t i)
{
  // Java: use int, and use >>> instead of >>
  // C or C++: use uint32_t
  i = i - ((i >> 1) & 0x55555555);
  i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
  return (((i + (i >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}


struct CompareSet_bit_int {
  
  bool operator()(const uint32_t lhs, const uint32_t rhs) const
  {
    //if(lhs.size()!=rhs.size()){
    //  return lhs.size() < rhs.size();
    //}
    
    if(numberOfSetBits(lhs)!=numberOfSetBits(rhs) ) {
      return numberOfSetBits(lhs)<numberOfSetBits(rhs);
    }
    
    return lhs<rhs;
    
  }
};

bool compare_two_set(const uint32_t lhs, const uint32_t rhs) 
{
  //if(lhs.size()!=rhs.size()){
  //  return lhs.size() < rhs.size();
  //}
  
  if(numberOfSetBits(lhs)!=numberOfSetBits(rhs) ) {
    return numberOfSetBits(lhs)<numberOfSetBits(rhs);
  }
  
  return lhs<rhs;
  
}


double calp2_c(float **t_read_count_mat_li,std::vector<int>& order_arr) {
  
  double p_sum=0;
  int dim=order_arr.size();
  
  for(int i=0; i<dim-1; i++) {
    for(int j=i+1; j<dim; j++) {
      p_sum=p_sum+t_read_count_mat_li[order_arr[i]][order_arr[j]];
    }
    
  }
  
  return p_sum;
}


double calp2_p_v(float **t_read_count_mat_li, int t_dim, int sim_times=1000000) {
  
  std::vector<int> init_order_p_v(t_dim);
  
  for(int i = 0; i <t_dim; i++)
    init_order_p_v[i] = i;
  
  
  
  //double p_v=0;
  
  std::sort( init_order_p_v.begin(),init_order_p_v.end() );
  
  float in_order_li=calp2_c(t_read_count_mat_li,init_order_p_v);
  
  
  double sum_of_less_than_in_order=0;
  
  double sim_sum=0;
  
  
  for(int i=0; i< sim_times; i++) {
    std::random_shuffle ( init_order_p_v.begin(), init_order_p_v.end() );
    
    double tmp_li=calp2_c(t_read_count_mat_li,init_order_p_v);
    
    sim_sum+= std::exp(tmp_li);
    
    if(tmp_li<=in_order_li ) {
      sum_of_less_than_in_order+=std::exp(tmp_li);
    }
    
  }
  //number_of_less_than_in_order
  
  double permut_p=sum_of_less_than_in_order/sim_sum;
  
  return permut_p;
}


inline float get_log_sum_c_bit(int t_sink, uint32_t t_parents,float **mat_li, int t_dim ) {
  
  float t_log_sum=0;
  
  //std::set<int>::iterator it;
  //for (it = t_parents.begin(); it != t_parents.end(); ++it){
  //for (std::size_t i = 0; i < t_parents.size(); ++i) {
  for (std::size_t i = 0; i < t_dim; ++i) {
    
    //if( t_parents.test(i)  && (mat_li[i][t_sink]!=0)) {
    //BIT_CHECK
    if( BIT_CHECK(t_parents,i)  && (mat_li[i][t_sink]!=0)) {
      
      //t_log_sum = t_log_sum + std::log(mat_li[i][t_sink]);
      t_log_sum = t_log_sum + (mat_li[i][t_sink]);
      
      
    }
  }
  
  return t_log_sum;
}


int combinations(uint32_t set, int at, int r, int n, std::vector<uint32_t>& subsets) {
  //std::cout<<'tx';
  
  //std::set<std::bitset<32>,CompareSet_bit > subsets;
  // Return early if there are more elements left to select than what is available.
  int elementsLeftToPick = n - at;
  if (elementsLeftToPick < r) {
    return 0;
  }
  //std::cout<<'tx';
  // We selected 'r' elements so we found a valid subset!
  if (r == 0) {
    subsets.push_back(set);
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

int binomial_coe(int n,int r) {
  if(r==0 || r==n)
    return 1;
  return binomial_coe(n-1,r)+binomial_coe(n-1,r-1);
}

std::vector<uint32_t> get_all_combinations(int n) {
  std::vector<uint32_t> subsets(pow(2,n));
  
  for(int r=1; r<=n; r++) {
    combinations(0, 0, r, n,subsets);
    //subsets.insert( t_set.begin(),t_set.end() );
  }
  return subsets;
}


std::vector<double> find_opti_dynam_c_bit(float **read_count_mat_li,double t_alpha_v, int dim) {
  
  
  std::vector<uint32_t> v_set_power_vec=get_all_combinations(dim);
  
  std::sort(v_set_power_vec.begin(), v_set_power_vec.end(), compare_two_set);
  v_set_power_vec.erase(std::unique(v_set_power_vec.begin(), v_set_power_vec.end(), compare_two_set) );
  
  //std::set<uint32_t,CompareSet_bit_int > v_set_power = get_all_combinations(dim)
  //std::vector<uint32_t> v_set_power_vec(v_set_power.begin(),v_set_power.end());
  std::vector<uint32_t>::iterator it;
  //v_set_power.clear();
  
  
  //std::cout<<"Got all combinations"<< std::endl;
  
  std::vector<float> scores(pow(2,dim),INT_MIN);
  std::vector<unsigned char> sinks(pow(2,dim),-1);
  
  //std::map<std::bitset<32>, int, CompareSet_bit> sinks;
  
  //std::cout<<"Allocate arrays"<< std::endl;
  
  
  for (it = v_set_power_vec.begin(); it != v_set_power_vec.end(); ++it) {
    
    uint32_t one_sub=*it;
    
    uint32_t number_bit= numberOfSetBits(one_sub);
    
    if(number_bit==0) {
      continue;
    }
    
    
    //uint32_t w_set=one_sub;
    
    if(number_bit>2) {
      //std::set<int>::iterator iit;
      
      //for (iit = one_sub.begin(); iit != one_sub.end(); ++iit){
      for (std::size_t i = 0; i < dim; ++i) {
        //if(! one_sub.test(i)) {
        
        
        if(!BIT_CHECK(one_sub,i)) {
          continue;
        }
        
        BIT_CLEAR(one_sub,i);
        
        float skore=scores[one_sub];
        skore=skore+get_log_sum_c_bit(i, one_sub, read_count_mat_li,dim);
        
        BIT_SET(one_sub,i);
        
        
        if( (skore > scores[one_sub] )   ) {
          
          scores[one_sub]=skore;
          sinks[one_sub]=i;
          
        }
        //one_sub.set(i);
        
      }
      
    } else if(number_bit==2) {
      int index_one=0;
      int index_two=0;
      
      bool a=true;
      //for (std::size_t i = 0; i < one_sub.size(); ++i) {
      for (std::size_t i = 0; i < dim; ++i) {
        
        //if(a&one_sub.test(i)) {
        if(a&&BIT_CHECK(one_sub,i) ) {
          
          index_one=i;
          a=false;
        }
        
        if((!a)&&BIT_CHECK(one_sub,i) ) {
          
          index_two=i;
        }
      }
      
      
      if(read_count_mat_li[index_one][index_two]>read_count_mat_li[index_two][index_one]) {
        //scores[w_set]=std::log(read_count_mat_li[index_one][index_two]);
        scores[one_sub]=(read_count_mat_li[index_one][index_two]);
        
        sinks[one_sub]=index_two;
      } else {
        //scores[w_set]=std::log(read_count_mat_li[index_two][index_one]);
        scores[one_sub]=(read_count_mat_li[index_two][index_one]);
        
        sinks[one_sub]=index_one;
      }
      
    } else {
      scores[one_sub]=0;
      //sinks[w_set]=*(one_sub.begin());
      //for (std::size_t i = 0; i < one_sub.size(); ++i) {
      for (std::size_t i = 0; i < dim; ++i) {
        
        if( BIT_CHECK(one_sub,i) ) {
          sinks[one_sub]=i;
        }
      }
      
      
    }
    
  }
  
  
  std::vector<int> ord(dim+1);
  uint32_t left=(pow(2,(dim))-1 );
  
  
  uint32_t full_set=(pow(2,(dim))-1 );
  
  float best_socre=scores[full_set];
  
  
  for(int i=dim-1; i>=0; i-- ) {
    
    
    ord[i]=sinks[left] ;
    
    BIT_CLEAR(left,ord[i] );
  }
  
  
  ord[dim]=best_socre;
  
  
  std::vector<double> ord_double;
  for(int i=0; i<ord.size(); i++) {
    ord_double.push_back(ord[i]);
  }
  
  return ord_double;
}




// [[Rcpp::export]]
Rcpp::List find_opti_dynam_r_cpp_bit(Rcpp::NumericMatrix t_read_count_mat,double t_alpha_v=0.1) {
  int dim = t_read_count_mat.nrow();
  
  //int dim=5;
  double **m_read_count_mat=new double *[dim];
  
  for(int i = 0; i <dim; i++)
    m_read_count_mat[i] = new double[dim];
  
  
  for(int i=0; i<dim; i++) {
    for(int j=0; j<dim; j++) {
      m_read_count_mat[i][j]= t_read_count_mat(i,j) ;
    }
  }
  
  
  float **m_read_count_mat_li=new float *[dim];
  for(int i = 0; i <dim; i++)
    m_read_count_mat_li[i] = new float[dim];
  
  
  for(int i=0; i<dim; i++) {
    for(int j=0; j<dim; j++) {
      if(i==j) {
        m_read_count_mat_li[i][j]=m_read_count_mat_li[j][i]=0;
        continue;
      }
      
      if( (m_read_count_mat[i][j]+m_read_count_mat[j][i])==0) {
        m_read_count_mat_li[i][j]=m_read_count_mat_li[j][i]=std::log(0.5);
        continue;
      }
      
      if(i>j) {
        float p= (m_read_count_mat[i][j]+t_alpha_v)/(m_read_count_mat[i][j]+m_read_count_mat[j][i]+2*t_alpha_v);
        
        m_read_count_mat_li[i][j]=std::log(p);
        m_read_count_mat_li[j][i]=std::log(1-p);
      }
      
    }
  }
  
  
  std::vector<double> ord=find_opti_dynam_c_bit(m_read_count_mat_li,t_alpha_v,dim);
  Rcpp::NumericVector out(dim );
  
  for(int i=0; i<dim; i++) {
    out[i]=ord[i]+1;
  }
  
  
  //std::cout<<dim;
  double permut_p=calp2_p_v(m_read_count_mat_li,dim);
  //double permut_p=0;
  
  Rcpp::List L = Rcpp::List::create(Rcpp::Named("best_order") = out ,
                                    Rcpp::_["entropy"]=NA_REAL,
                                    Rcpp::_["best_score"] = ord[dim],
                                                               Rcpp::_["number_of_maximum_order"]=NA_INTEGER,
                                                               Rcpp::_["permut_p"]=permut_p
  );
  
  return L;
}



