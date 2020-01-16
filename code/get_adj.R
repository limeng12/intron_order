#setwd("/Users/mengli/Documents/projects/iso");
library(igraph);
library(networkD3);
library(stringi)
strReverse <- function(x)
  sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")

source("code/map.R");
source("code/to_gephi.R");

source("code/triangle_star_shapes.R");


get_adj<-function(output_path,t_iso_final,t_iso_slow_sumary,number_transcript_calculated,read_count_t,
                  save_file=FALSE,return_graph_list=FALSE){
  
  
  igraph_list<-list();
  
  #post_evidence_count_all<-matrix(nrow = 0,ncol=4);
  
  for(g in 1:min(nrow(t_iso_slow_sumary),number_transcript_calculated ) ){
    
    iso_1<-unique(t_iso_final[t_iso_final[,"id"] %in% t_iso_slow_sumary[g,"id"],
                              c("gene_symbol","id","gencode_intron_o_first",
                                "gencode_intron_o_next","strand","read_count","first","nexti") ]);
    
    iso_1<-iso_1[iso_1[,"read_count"]>=read_count_t,];
    
    
    iso_1[,"gencode_intron_o_first"]<-as.character(iso_1[,"gencode_intron_o_first"]);
    iso_1[,"gencode_intron_o_next"]<-as.character(iso_1[,"gencode_intron_o_next"]);
    
    
    all_ints<-unique(c(iso_1[,"gencode_intron_o_first"],iso_1[,"gencode_intron_o_next"]))
    
    intron_number_position_index<-str_count(all_ints[1],"_")+1
    
    
    t_gene_symbol<-as.character(iso_1[1,"gene_symbol"]);
    t_trans_id<-iso_1[1,"id"];
    if(is.na(t_gene_symbol) ){
      t_gene_symbol<-""
    }
    
    print(str_c( g,":",t_gene_symbol,":",t_trans_id) );
    
    
    
    t_total_intron_count<- t_iso_slow_sumary[g,"int_count"];
    
    
    
    adjacency_matrix<-matrix(0,nrow=(t_total_intron_count),ncol=(t_total_intron_count) );
    rownames(adjacency_matrix)<-(1:t_total_intron_count);
    colnames(adjacency_matrix)<-(1:t_total_intron_count);
    
    
    for(i in 1:nrow(iso_1)){
        gencode_intron_o_first_number<-
          as.character(sapply(str_split((iso_1[i,"gencode_intron_o_first"]),"_"),"[", intron_number_position_index) );
        
        gencode_intron_o_next_number<-
          as.character(sapply(str_split((iso_1[i,"gencode_intron_o_next"]),"_"),"[", intron_number_position_index) );
        
        adjacency_matrix[ gencode_intron_o_first_number , gencode_intron_o_next_number ]<-(iso_1[i,"read_count"]);
    }
      
    
    
    g2<-graph_from_adjacency_matrix(adjacency_matrix,mode="directed",weighted=TRUE);
      
    post_evidence_count_t_v<-rep(0,length(V(g2)));
      
    for( i in 1:length(V(g2)) ){
      #print(i);
      this_v_order<-as.numeric( names(V(g2)[i]) );
        
      neighbors<-neighbors(g2,  V(g2)[i],mode="in");
        
      #neighbors_orders<-as.numeric( sapply(strsplit(names(neighbors),"_"),"[",2) );
      neighbors_orders<-as.numeric( names(neighbors) )
        
      post_evidence_count<-sum(neighbors_orders>this_v_order,na.rm = TRUE);
      post_evidence_count_t_v[i]<-post_evidence_count

        
    }
      
    
    node_size=map(post_evidence_count_t_v,c(8.2,50.6) );
    node_size[is.na(node_size)]=min(node_size,na.rm = TRUE);
      
    names(node_size)<-V(g2)$name;# sapply(str_split(V(g2)$name,":"),"[",1);
      
    node_color=rep("green",length(V(g2)) );
    names(node_color)<-V(g2)$name;# sapply(str_split(V(g2)$name,":"),"[",1 );
      
      
    igraph_list[[str_c(t_trans_id) ]]<-
      list(
        gene_symbol=t_gene_symbol,
        trans_id=t_trans_id,
        adjacency_matrix=adjacency_matrix,
        edge_count=(iso_slow_sumary[g,"edge_count"]),
        g2=g2,
        node_size=node_size,
        node_color=node_color
      );
    
  }
  
  
  if(return_graph_list){
    return(igraph_list)
  }
  
  
}


