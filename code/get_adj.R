#setwd("/Users/mengli/Documents/projects/iso");
library(igraph);
library(networkD3);
library(stringi);

source("code/map.R");
source("code/to_gephi.R");

source("code/triangle_star_shapes.R");
source("code/utils.R")

get_adj<-function(output_path,t_iso_final,t_iso_slow_sumary,number_transcript_calculated,
                  save_file=FALSE,return_graph_list=FALSE){
  
  
  igraph_list<-list();
  
  #post_evidence_count_all<-matrix(nrow = 0,ncol=4);
  
  for(g in 1:min(nrow(t_iso_slow_sumary),number_transcript_calculated ) ){
    
    iso_1<-unique(t_iso_final[t_iso_final[,"id"] %in% t_iso_slow_sumary[g,"id"],
                              c("gene_symbol","id","gencode_intron_o_first",
                                "gencode_intron_o_next","strand","read_count","first","nexti") ]);
    
    iso_1[,"gencode_intron_o_first"]<-as.character(iso_1[,"gencode_intron_o_first"]);
    iso_1[,"gencode_intron_o_next"]<-as.character(iso_1[,"gencode_intron_o_next"]);
    
    all_ints<-unique( c(iso_1[,"gencode_intron_o_first"],iso_1[,"gencode_intron_o_next"]) );
    
    ##all the intron numbers
    intron_number_position_index<-str_count(all_ints[1],"_")+1
    
    t_gene_symbol<-as.character(iso_1[1,"gene_symbol"]);
    t_trans_id<-iso_1[1,"id"];
    if(is.na(t_gene_symbol) ){
      t_gene_symbol<-""
    }
    
    print(str_c( g,":",t_gene_symbol,":",t_trans_id) );
    
    
    ##total number of introns
    t_total_intron_count<- t_iso_slow_sumary[g,"int_count"];
    
    ##make intron splicing order matrix
    adjacency_matrix<-matrix(0,nrow=(t_total_intron_count),ncol=(t_total_intron_count) );
    rownames(adjacency_matrix)<-(1:t_total_intron_count);
    colnames(adjacency_matrix)<-(1:t_total_intron_count);
    
    ##
    gencode_intron_o_first_number<-
      as.character(sapply(str_split((iso_1[,"gencode_intron_o_first"]),"_"),"[", intron_number_position_index) );
    
    gencode_intron_o_next_number<-
      as.character(sapply(str_split((iso_1[,"gencode_intron_o_next"]),"_"),"[", intron_number_position_index) );
    
    
    for(i in 1:nrow(iso_1)){
        adjacency_matrix[ gencode_intron_o_first_number[i] , gencode_intron_o_next_number[i] ]<-(iso_1[i,"read_count"]);
    }
    
    ##add node size and color
    node_size<-rep(10, t_total_intron_count)
    names(node_size)<-1:t_total_intron_count; 
    
    node_color=rep("green", t_total_intron_count );
    names(node_color)<-1:t_total_intron_count;
    
    
    percent_coverage_pair<-cal_intron_pair_cov(adjacency_matrix);
    
    
    igraph_list[[str_c(t_trans_id) ]]<-
      list(
        gene_symbol=t_gene_symbol,
        trans_id=t_trans_id,
        adjacency_matrix=adjacency_matrix,
        percent_coverage_pair=percent_coverage_pair,
        #coverage_pair=count_of_nonzero_ele,
        intron_pair_count=(t_iso_slow_sumary[g,"intron_pair_count"]),
        #g2=g2,
        node_size=node_size,
        node_color=node_color
      );
    
  }
  
  
  if(return_graph_list){
    return(igraph_list)
  }
  
  
}


