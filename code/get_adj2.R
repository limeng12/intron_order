#setwd("/Users/mengli/Documents/projects/iso");
library(igraph);
#library(networkD3);
library(stringi);

source("code/map.R");
source("code/to_gephi.R");

source("code/triangle_star_shapes.R");
source("code/utils.R")

# t_iso_final<-iso_final
# t_iso_slow_sumary<-iso_summary
# intron_pair_cov_threshold<-0.95

get_adj2<-function(t_iso_final,t_iso_slow_sumary,intron_pair_cov_threshold,t_read_count_threshold){
  
  isoform_num_produce<-sum(t_iso_slow_sumary[,"percent_intron_pair_coverage"]>intron_pair_cov_threshold);
  
  print(paste0("Number of multi introns transcripts that contain intron pairs >",
               sprintf("%1.2f%%", intron_pair_cov_threshold*100),"= ", isoform_num_produce ) );
  
  
  igraph_list<-list();
  
  #t_iso_final<-f
  #t_iso_slow_sumary<-g
  
  isoform_num_produce<-sum(t_iso_slow_sumary[,"percent_intron_pair_coverage"] > intron_pair_cov_threshold);
  
  #only seleted well detected transcripts in iso_slow_summary object
  t_iso_slow_sumary<-t_iso_slow_sumary[t_iso_slow_sumary[,"percent_intron_pair_coverage"] > intron_pair_cov_threshold,];
  
  #only seleted well detected transcripts in iso object
  t_iso_final<-t_iso_final[t_iso_final[,"id"] %in% t_iso_slow_sumary[,"trans_id"],];
  
  
  #t_iso_final<-t_iso_final[order(t_iso_final[,"id"]),];
  #row_num_group_index<-group_split(data.frame(id=t_iso_final[,"id"], row_num=1:nrow(t_iso_final)  ),id,keep = TRUE);
  
  #row_num_group_index<-group_keys(data.frame(id=t_iso_final[,"id"], row_num=1:nrow(t_iso_final) ) );
  #row_num_group_index<-split(data.frame(id=t_iso_final[,"id"], row_num=1:nrow(t_iso_final)  ),id );
  
  cat("\n");
  
  print("Build adjacent matrix: ");
  pb = txtProgressBar(min = 0, max = nrow(t_iso_slow_sumary ), initial = 0, width=100, style=3) ;
  
  
  for(g in 1:nrow(t_iso_slow_sumary)  ){
    
    setTxtProgressBar(pb ,g);
    
    #t_trans_id<-iso_1[1,"id"];
    t_trans_id<-t_iso_slow_sumary[g,"trans_id"];
    #t_trans_id<-as.character( as.data.frame(row_num_group_index[[g]]) [1,"id"]);
    
    #if( (t_iso_slow_sumary[g,"percent_intron_pair_coverage"] < intron_pair_cov_threshold) ){
    #  next;
    #}
    iso_1<-unique(t_iso_final[t_iso_final[,"id"] %in% t_iso_slow_sumary[g,"trans_id"],])
    
    #iso_1<-unique(t_iso_final[ as.data.frame(row_num_group_index[[g]])[,"row_num"], ] );
    
    
    t_gene_symbol<-"";
    
    if(is.na(t_gene_symbol) ){
      t_gene_symbol<-""
    }
    
    #print(str_c( g,":",t_gene_symbol,":",t_trans_id) );
    
    
    ##total number of introns
    t_total_intron_count<- t_iso_slow_sumary[g,"intron_count"];
    
    ##make intron splicing order matrix
    adjacency_matrix<-matrix(0,nrow=(t_total_intron_count),ncol=(t_total_intron_count) );
    rownames(adjacency_matrix)<-(1:t_total_intron_count);
    colnames(adjacency_matrix)<-(1:t_total_intron_count);
    
    jc_pair_matrix<-matrix(0,nrow=(t_total_intron_count),ncol=(t_total_intron_count) );
    rownames(jc_pair_matrix)<-(1:t_total_intron_count);
    colnames(jc_pair_matrix)<-(1:t_total_intron_count);
    
    ##
    #gencode_intron_o_first_number<-
    #  as.character(sapply(str_split((iso_1[,"gencode_intron_o_first"]),"_"),"[", intron_number_position_index) );
    
    #gencode_intron_o_next_number<-
    #  as.character(sapply(str_split((iso_1[,"gencode_intron_o_next"]),"_"),"[", intron_number_position_index) );
    
    gencode_intron_o_first_number<-iso_1[,"intron_order_first"];
    gencode_intron_o_next_number<-iso_1[,"intron_order_next"];
    
    for(i in 1:nrow(iso_1)){
      adjacency_matrix[ gencode_intron_o_first_number[i], gencode_intron_o_next_number[i] ]<-(iso_1[i,"read_count"]);
      
      jc_pair_matrix[ gencode_intron_o_first_number[i], gencode_intron_o_next_number[i] ]<-(iso_1[i,"read_count_jc"]);
    }
    
    jc_pair_matrix[is.na(jc_pair_matrix)]<-0;
    
    
    for(ii in 1:nrow(adjacency_matrix)){
      for(jj in 1:ncol(adjacency_matrix)){
        if(adjacency_matrix[ii,jj]+adjacency_matrix[jj,ii]<t_read_count_threshold){
          
          adjacency_matrix[ii,jj]<-adjacency_matrix[jj,ii]<-0;
        }
        
      }
    }
    
    percent_coverage_pair<-cal_intron_pair_cov(adjacency_matrix);
    
    
    if(percent_coverage_pair < intron_pair_cov_threshold ){
      next;
    }
    
    ##add node size and color
    node_size<-rep(10, t_total_intron_count);
    names(node_size)<-1:t_total_intron_count; 
    
    node_color=rep("green", t_total_intron_count );
    names(node_color)<-1:t_total_intron_count;
    
    igraph_list[[str_c(t_trans_id) ]]<-
      list(
        gene_symbol=t_gene_symbol,
        trans_id=t_trans_id,
        strand=iso_1[1,"strand"],
        adjacency_matrix=adjacency_matrix,
        jc_pair_matrix=jc_pair_matrix,
        percent_coverage_pair=percent_coverage_pair,
        #coverage_pair=count_of_nonzero_ele,
        intron_pair_count=(t_iso_slow_sumary[g,"intron_pair_count"]),
        #g2=g2,
        node_size=node_size,
        #index_pos=intron_pos_index_fr$pos,
        #index_pos_fr=intron_pos_index_fr,
        
        node_color=node_color
      );
    
  }
  
  #print(str_c( "Total number of transcripts have adj matrix: ",length(igraph_list)) );
  
  #if(return_graph_list){
  return(igraph_list)
  #}
  
  
}


