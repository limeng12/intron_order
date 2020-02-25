library(networkD3);
source("code/get_members.R")


# i<-1

draw_3d<-function(tt_igraph_list,output_path,t_alpha_v=0.05,output_graph=FALSE){
  
    
    # i<-1;
    # tt_igraph_list<-t_igraph_list
    for( i in 1:length(tt_igraph_list) ){
      
      
      print(paste0(i,":",names(tt_igraph_list)[i]) );
      
      #mermbers <- membership(cluster_walktrap(tt_igraph_list[[i]]) );
      
      #t_net_one<-igraph_to_networkD3(tt_igraph_list[[i]],mermbers);
      
      #t_net_one$links$value<-t_net_one$links$value/(max(t_net_one$links$value))*10;
      
      t_charge<- (-2*(   tt_igraph_list[[i]]$intron_pair_count )  );
      
      if(t_charge< -80){
        t_charge<- -80
      }
      
      if(t_charge> -20){
        t_charge<- -20
      }
      print(t_charge);
      
      #t_linkDistance<-500/iso_slow_sumary[i,"int_count"];
      t_linkDistance<-(1*(   tt_igraph_list[[i]]$intron_pair_count )  );
      
      if(t_linkDistance>30){
        t_linkDistance<-30
      }
      
      if(t_linkDistance<10){
        t_linkDistance<-10
      }
      
      
      #mermbers <- membership(cluster_optimal(g2) );
      
      #mermbers <- membership(cluster_walktrap(g2) );
      #mermbers<-get_members_entropy_weight(tt_igraph_list[[i]]$adjacency_matrix);
      #mermbers<-get_members_matrix(tt_igraph_list[[i]]$adjacency_matrix,t_alpha_v);
      mermbers<-tt_igraph_list[[i]]$members;
      
      #mermbers <- membership(cluster_fast_greedy(g2) );
      
      g2<-graph_from_adjacency_matrix(tt_igraph_list[[i]]$adjacency_matrix,mode="directed",weighted=TRUE);
      
      #t_net_one<-igraph_to_networkD3(tt_igraph_list[[i]]$g2,mermbers);
      t_net_one<-igraph_to_networkD3(g2,mermbers);

      
      #edge_width=map( E(g2)$weight, c(0.25,6) );
      
      t_net_one$nodes[,"node_size"]<-tt_igraph_list[[i]]$node_size[as.character(t_net_one$nodes$name)];
      t_net_one$nodes[,"node_color"]<-tt_igraph_list[[i]]$node_color[as.character(t_net_one$nodes$name)];
      
      
      #t_net_one$links
      t_net_one$links$value<-t_net_one$links$value/(max(t_net_one$links$value))*10;
      
      
      tt_igraph_list[[i]]$one_net<-t_net_one
      
      
      if(output_graph){
        g_one_net<-forceNetwork(Links = t_net_one$links, Nodes = t_net_one$nodes,
                                Source = 'source', Target = 'target', Value='value', NodeID = 'name',
                                Group = 'group',arrows=TRUE,Nodesize='node_size',
                                opacity=1, bounded=TRUE, legend=TRUE,
                                charge = t_charge, zoom=FALSE, opacityNoHover = 0.9,
                                radiusCalculation = JS(" Math.sqrt(d.nodesize)+1")
        );
        
        saveNetwork(g_one_net,file=paste0(output_path,
                                          names(tt_igraph_list)[i] ,"_d3_network.html") );
      }
      #                           linkDistance=t_linkDistance
      
      #  linkDistance = JS("function(d) {   if(d.value*10<5){ return 5; }  return d.value * 10;}"
      
      #  View(tt_igraph_list[[i]]$nodes);
      
      # linkDistance = 25, linkDistance = 30, linkDistance = 35,
      # opacity=1,linkDistance = 1,bounded=TRUE,
      # charge = -7,zoom=FALSE,opacityNoHover = 0.9
      
      #arrows = TRUE
      #print(a);
      #b<-simpleNetwork(Links = t_net_one$links, Nodes = t_net_one$nodes,
      #                Source = 'source', Target = 'target', NodeID = 'name',
      #                Group = 'group');
      # linkDistance = 5,,bounded=TRUE,charge = -10,zoom=FALSE,arrows=TRUE,
      

    }
  
  
  return(tt_igraph_list);
  
}
#saveNetwork(b,file="result/d3_network_test.html")

