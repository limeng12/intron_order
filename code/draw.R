#setwd("/Users/mengli/Documents/projects/iso");
library(igraph);
library(networkD3);
library(stringi)
strReverse <- function(x)
  sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")

source("code/map.R");
source("code/to_gephi.R");

source("code/triangle_star_shapes.R");


draw<-function(output_path,t_iso_final,t_iso_slow_sumary,count,read_count_t,save_file=FALSE,return_graph_list=FALSE){
  
  igraph_list<-list();
  
  post_evidence_count_all<-matrix(nrow = 0,ncol=4);
  
  #pdf(paste0("result/iso_test_unique_",is_large,"_",label,".pdf"),width = 40,height = 40);
  
  pdf(paste0(output_path),width = 40,height = 40);
  
  for(g in 1:min(nrow(t_iso_slow_sumary),count ) ){
    
    iso_1<-unique(t_iso_final[t_iso_final[,"id"] %in% t_iso_slow_sumary[g,"id"],
                              c("gene_symbol","id","gencode_intron_o_first",
                                "gencode_intron_o_next","strand","read_count","first","nexti") ]);
    
    iso_1<-iso_1[iso_1[,"read_count"]>=read_count_t,];
    
    #adjacency_matrix<-dcast(iso_1, gencode_intron_o_first~gencode_intron_o_next,length);
    #iso_1<-iso_1[str_detect(iso_1[,"gencode_intron_o_first"], iso_slow_sumary[g,"id"] ) & 
    #               str_detect(iso_1[,"gencode_intron_o_first"], iso_slow_sumary[g,"id"] ) ,]
    if(nrow(iso_1)<2){
      next;
    }
    
    t_gene_symbol<-as.character(iso_1[1,"gene_symbol"]);
    if(is.na(t_gene_symbol)){
      t_gene_symbol=""
    }
    
    print( str_c(g,"  :  ",t_gene_symbol,":",iso_1[1,"id"]) );
    
    t_trans_id<-iso_1[1,"id"];
    
    
    iso_1[,"gencode_intron_o_first"]<-as.character(iso_1[,"gencode_intron_o_first"]);
    iso_1[,"gencode_intron_o_next"]<-as.character(iso_1[,"gencode_intron_o_next"]);
    
    all_ints<-unique(c(iso_1[,"gencode_intron_o_first"],iso_1[,"gencode_intron_o_next"]))
    
    intron_number_index<-str_count(all_ints[1],"_")+1
    
    
    t_total_intron_count<-gencode_intron_o_frame_intron_count[t_trans_id];
    
    
    #adjacency_matrix<-matrix(0,nrow=length(all_ints),ncol=length(all_ints) );
    #rownames(adjacency_matrix)<-(all_ints);
    #colnames(adjacency_matrix)<-(all_ints);
    
    adjacency_matrix<-matrix(0,nrow=(t_total_intron_count),ncol=(t_total_intron_count) );
    rownames(adjacency_matrix)<-(1:t_total_intron_count);
    colnames(adjacency_matrix)<-(1:t_total_intron_count);
    
    
    is_alt_map_se<-c()
    is_alt_map_a5<-c()
    is_alt_map_a3<-c()
    is_alt_map_ri<-c()
    
    #if(cosi_raw_group_vector_format)
    
    cosi_map<-list();
    
    posp_introns_map<-c();
    
    
    for(i in 1:nrow(iso_1)){
      
      gencode_intron_o_first_number<-
        as.character(sapply(str_split((iso_1[i,"gencode_intron_o_first"]),"_"),"[", intron_number_index) );
      
      gencode_intron_o_next_number<-
        as.character(sapply(str_split((iso_1[i,"gencode_intron_o_next"]),"_"),"[", intron_number_index) );
      
      # next doesn't spliced, # first spliced
      # adjacency_matrix: row directed to col
      # 
      #adjacency_matrix[ iso_1[i,"gencode_intron_o_first"] , iso_1[i,"gencode_intron_o_next"] ]<-1
      #adjacency_matrix[ iso_1[i,"gencode_intron_o_next"] , iso_1[i,"gencode_intron_o_first"] ]<-(iso_1[i,"read_count"])
      
      adjacency_matrix[ gencode_intron_o_first_number , gencode_intron_o_next_number ]<-(iso_1[i,"read_count"]);
      
      #adjacency_matrix[ iso_1[i,"gencode_intron_o_first"] , iso_1[i,"gencode_intron_o_next"] ]<-(iso_1[i,"read_count"])
      
      
      
      
      if(iso_1[i,"first"] %in% all_alt_se_introns){
        is_alt_map_se<-c(is_alt_map_se,iso_1[i,"gencode_intron_o_first"] )
      }
      if(iso_1[i,"first"] %in% all_alt_a5_introns){
        is_alt_map_a5<-c(is_alt_map_a5,iso_1[i,"gencode_intron_o_first"] )
      }
      if(iso_1[i,"first"] %in% all_alt_a3_introns){
        
        is_alt_map_a3<-c(is_alt_map_a3,iso_1[i,"gencode_intron_o_first"] )
      }
      if(iso_1[i,"first"] %in% all_alt_ri_introns){
        is_alt_map_ri<-c(is_alt_map_ri,iso_1[i,"gencode_intron_o_first"] )
      }
      
      if(iso_1[i,"nexti"] %in% all_alt_se_introns){
        is_alt_map_se<-c(is_alt_map_se,iso_1[i,"gencode_intron_o_next"] )
      }
      if(iso_1[i,"nexti"] %in% all_alt_a5_introns){
        is_alt_map_a5<-c(is_alt_map_a5,iso_1[i,"gencode_intron_o_next"] )
      }
      if(iso_1[i,"nexti"] %in% all_alt_a3_introns){
        
        is_alt_map_a3<-c(is_alt_map_a3,iso_1[i,"gencode_intron_o_next"] )
      }
      if(iso_1[i,"nexti"] %in% all_alt_ri_introns){
        is_alt_map_ri<-c(is_alt_map_ri,iso_1[i,"gencode_intron_o_next"] )
      }
      
      if(iso_1[i,"nexti"] %in% names(cosi_raw_group_vector_format) ){
        cosi_map[[  iso_1[i,"gencode_intron_o_next"] ]]<-cosi_raw_group_vector_format[iso_1[i,"nexti"] ];
      }
      
      if(iso_1[i,"first"] %in% names(cosi_raw_group_vector_format) ){
        cosi_map[[  iso_1[i,"gencode_intron_o_first"]  ]]<-cosi_raw_group_vector_format[iso_1[i,"first"] ];
        
      }
      
      if(iso_1[i,"nexti"] %in% posp_introns ){
        #cosi_map[[  iso_1[i,"gencode_intron_o_next"] ]]<-cosi_raw_group_vector_format[iso_1[i,"nexti"] ];
        posp_introns_map<-c(posp_introns_map,iso_1[i,"gencode_intron_o_next"])
      }
      
      if(iso_1[i,"first"] %in% posp_introns ){
        #cosi_map[[  iso_1[i,"gencode_intron_o_first"]  ]]<-cosi_raw_group_vector_format[iso_1[i,"first"] ];
        posp_introns_map<-c(posp_introns_map,iso_1[i,"gencode_intron_o_first"])
        
      }
      
    }
    
    
    #is_alt_map_se<-str_sub(is_alt_map_se,17+7);
    #is_alt_map_a5<-str_sub(is_alt_map_a5,17+7);
    #is_alt_map_a3<-str_sub(is_alt_map_a3,17+7);
    #is_alt_map_ri<-str_sub(is_alt_map_ri,17+7);
    
    #names(cosi_map)<-str_sub(names(cosi_map),17+7);
    #posp_introns_map<-str_sub(posp_introns_map,17+7);
    
    
    #rownames(adjacency_matrix)<-str_sub(all_ints,17+7);
    #colnames(adjacency_matrix)<-str_sub(all_ints,17+7);
    
    
    is_alt_map_se<-sapply(str_split((is_alt_map_se),"_"),"[",intron_number_index) ;
    is_alt_map_a5<-sapply(str_split((is_alt_map_a5),"_"),"[",intron_number_index) ;
    is_alt_map_a3<-sapply(str_split((is_alt_map_a3),"_"),"[",intron_number_index) ;
    is_alt_map_ri<-sapply(str_split((is_alt_map_ri),"_"),"[",intron_number_index) ;
    
    names(cosi_map)<-sapply(str_split((names(cosi_map)),"_"),"[",intron_number_index) ;
    posp_introns_map<-sapply(str_split((posp_introns_map),"_"),"[",intron_number_index) ;
    
    
    #rownames(adjacency_matrix)<-sapply(str_split((all_ints),"_"),"[",intron_number_index) ;
    #colnames(adjacency_matrix)<-sapply(str_split((all_ints),"_"),"[",intron_number_index) ;
    
    
    
    g2<-graph_from_adjacency_matrix(adjacency_matrix,mode="directed",weighted=TRUE);
    
    #  degree(g2,mode="in);
    col <- rep("green",length(V(g2)))
    shape <- rep("square",length(V(g2)))
    
    post_evidence_count_t_v<-rep(0,length(V(g2)));
    #tkplot(g2)
    #neighbors(g2,  V(g2)[1],mode="in")->a
    
    for( i in 1:length(V(g2)) ){
      #this_v_order<-as.numeric( sapply(strsplit(names(V(g2)[i]),"_"),"[",2) );
      #shape[i]<-"circle";
      this_v_order<-as.numeric( names(V(g2)[i]) )
      #t_is_alt<-true
      if(  names(V(g2)[i]) %in% is_alt_map_se ){
        col[i]="yellow";
      }
      if(  names(V(g2)[i]) %in% is_alt_map_a5 ){
        #col[i]="pink";
        col[i]="cyan";
      }
      
      if(  names(V(g2)[i]) %in% is_alt_map_a3 ){
        #col[i]="red";
        #col[i]="blue";
        col[i]="wheat"
      }
      if(  names(V(g2)[i]) %in% is_alt_map_ri ){
        col[i]="white";
      }
      
      neighbors<-neighbors(g2,  V(g2)[i],mode="in");
      
      #neighbors_orders<-as.numeric( sapply(strsplit(names(neighbors),"_"),"[",2) );
      neighbors_orders<-as.numeric( names(neighbors) )
      
      post_evidence_count<-sum(neighbors_orders>this_v_order,na.rm = TRUE);
      post_evidence_count_t_v[i]<-post_evidence_count
      if((post_evidence_count>0) ){
        shape[i]<-"circle"
      }
      
      #post_evidence_count_all<-c(post_evidence_count_all,post_evidence_count)
      #post_evidence_count_all<-rbind(post_evidence_count_all,
      #                               c(t_gene_symbol,t_trans_id,names(V(g2))[i],post_evidence_count) );
      
    }
    
    #
    
    new_name<-c();
    
    vertex.frame.color_t<-rep("white",length(V(g2)));
    vertex.label.color_t<-rep("black",length(V(g2)));
    
    for( i in 1:length(V(g2)) ){
      t_new_name<-names(V(g2)[i]);
      
      
      if(names(V(g2)[i]) %in% posp_introns_map ){
        #vertex.frame.color_t[i]<-"black"
        vertex.frame.color_t[i]<-"white"
        
        vertex.label.color_t[i]<-"black"
        shape[i]<-"triangle"
        post_evidence_count_t_v[i]<-max(post_evidence_count_t_v,na.rm = TRUE);
        
      }
      
      if(t_new_name %in% names(cosi_map) ){
        t_new_name<-paste0( names(V(g2)[i]),":" ,as.character(cosi_map[[t_new_name]] ) );
      }
      
      new_name<-c(new_name,t_new_name)
    }
    
    #ay <- set.vertex.attribute(g2, "name", new_name )
    
    V(g2)$label <- new_name
    V(g2)$name <- new_name
    
    #layout_with_graphopt,charge=0.0008
    #layout_with_fr
    #layout_with_lgl
    #layout_nicely
    #layout_with_kk
    # t_trans_id<-iso_1[1,"id"];
    t_total_intron_count<-gencode_intron_o_frame_intron_count[t_trans_id];
    
    main=str_c(t_gene_symbol,":",t_trans_id,":",iso_1[1,"strand"],":",iso_1[1,"gencode_intron_o_next"],":",
               t_total_intron_count);
    
    #############################################################################################################          
    #g3<<-g2;

    
    #print(post_evidence_count_t_v);
    #tt<<-post_evidence_count_t_v;
    
    node_size=map(post_evidence_count_t_v,c(8.2,50.6) );
    node_size[is.na(node_size)]=min(node_size,na.rm = TRUE);
    
    #print(V(g2)$name);
    
    names(node_size)<-V(g2)$name;# sapply(str_split(V(g2)$name,":"),"[",1);
    
    
    node_color=col;
    names(node_color)<-V(g2)$name;# sapply(str_split(V(g2)$name,":"),"[",1 );
    
    
    
    #links=t_net_one$links,
    #nodes=t_net_one$nodes,
    igraph_list[[str_c(t_trans_id,"_",t_gene_symbol) ]]<-
      list(
           adjacency_matrix=adjacency_matrix,
           edge_count=(iso_slow_sumary[g,"edge_count"]),
           g2=g2,
           node_size=node_size,
           node_color=node_color
      );
    
    
    plot(g2,mark.shape=0.01,layout=layout_with_graphopt(g2,niter=5000),
         #vertex.color="green",
         cex=2,
         vertex.label.cex=1.0,edge.arrow.size=1,
         #vertex.size=map(degree(g2),c(1,5)),
         vertex.size=map(post_evidence_count_t_v,c(2.4,4.5) ),
         vertex.size2=map(post_evidence_count_t_v,c(2.4,4.5) ),
         edge.width=map( E(g2)$weight, c(0.25,6) ), 
         vertex.frame.color=vertex.frame.color_t,
         vertex.color=col,vertex.shape=shape,
         vertex.label.color=vertex.label.color_t,
         edge.curved=0.1  );
    
    title(main,cex.main=5);
    
    # as_adj_edge_list(g2)
    

    #names(igraph_list)[g]<-str_c(t_gene_symbol,":",t_trans_id);
    if(save_file){
      
      write_graph(g2,file=paste0("./data/graph/",main,".graphml"), format="graphml" );
      
      saveAsGEXF(g2,file=paste0("./data/graph/",main,".gexf") );
    }
    
    
  }
  
  
  dev.off();
  
  
  
  
  
  #post_evidence_count_all_fr<-data.frame(gene=post_evidence_count_all[,1],
  #                                       trans_id=post_evidence_count_all[,2],
  #                                       intron_order=post_evidence_count_all[,3],
  #                                       post_downstream_first_spliced_intron_count=post_evidence_count_all[,4]);
  
  
  #write.table(post_evidence_count_all_fr,file=paste0("result/post_downstream_first_spliced_intron_count_unique_intron_",
  #                                                   is_large,"_",label,".tsv"),
  #            sep = "\t",quote = FALSE,col.names = TRUE,row.names = FALSE);
  
  if(return_graph_list){
    return(igraph_list)
  }
  
}

