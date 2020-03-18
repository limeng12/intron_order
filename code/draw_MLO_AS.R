library(igraph);
library(stringr);
library(Sushi);
library(readr);
library(DT);
library(ggraph);
library(tidygraph);


#   tt_igraph_list<-tt_igraph_list[1:2]
#   k<-1

draw_mlo_plot<-function(tt_igraph_list,output_path){
  
  int_num<-sapply(tt_igraph_list,function(x){
   
    length(x$best_order) 
  })
  
  tt_igraph_list<-tt_igraph_list[order(int_num,decreasing = TRUE)]
  
  #pdf("result/all_MLO_Plot.pdf",width = 12,height=8);
  
  pdf(output_path,width = 12,height=8);
  
  #    k<-1145    k<-1144
  for(k in 1:min(length(tt_igraph_list), .Machine$integer.max) ){
    
    
    orders<-tt_igraph_list[[k]]$best_order;
    
    print(paste0(tt_igraph_list[[k]]$gene_symbol,":",tt_igraph_list[[k]]$trans_id,":",k) );
    
    
    ###sort intro position by order
    index_pos_fr<-tt_igraph_list[[k]]$index_pos_fr;
    
    intron_pos<-c();
    
    for(i in 1:nrow(index_pos_fr) ){
      if(orders[i] %in% index_pos_fr[,"index"]){
        intron_pos<-c(intron_pos,index_pos_fr[orders[i],"pos"]);
      }else{
        intron_pos<-c(intron_pos,NA);
      }
    }
    
    
    NodeList <- data.frame(nodes=orders, x=(0:(length(orders)-1)), 
                           y=rep(0,length(orders)),
                           labels=orders,
                           r=rep(0.02,length(orders)),
                           intron_pos=intron_pos
                           );
    
    col<-c();
    
    for(i in 1:nrow(NodeList)){
      col[i]<-"green"
      
      if( NodeList[i,"intron_pos"] %in% all_alt_se_introns){
        col[i]="yellow";NodeList[i,"labels"]<-str_c(NodeList[i,"labels"],":SE");
      }
      if( NodeList[i,"intron_pos"] %in% all_alt_a5_introns){
        col[i]="cyan";  NodeList[i,"labels"]<-str_c(NodeList[i,"labels"],":A5");
      }      
      if( NodeList[i,"intron_pos"] %in% all_alt_a3_introns){
        col[i]="wheat"; NodeList[i,"labels"]<-str_c(NodeList[i,"labels"],":A3");
      }
      if( NodeList[i,"intron_pos"] %in% all_alt_ri_introns){
        col[i]="white"; NodeList[i,"labels"]<-str_c(NodeList[i,"labels"],":RI");
      }
      
    }
    
    NodeList$col<-col;
    
    NodeList$x<-NodeList$x/max(NodeList$x)*2-1
    
    
    
    m_adj_matris<-tt_igraph_list[[k]]$adjacency_matrix;
    
    from<-c();
    to<-c();
    weight<-c();
    #from<-c(from,orders[i])
    #to<-c(to,orders[j]);
    
    #weight<-c(weight,m_adj_matris[orders[i],orders[j]]);
    for(i in 1:length(orders) ){
      for(j in 1:length(orders)){
        if(i!=j){
          from<-c(from,orders[i]);
          to<-c(to,orders[j]);
          weight<-c(weight,m_adj_matris[orders[i],orders[j]]);
        }
        
      }
    }
    
    
    EdgeList <- data.frame(from,to,weight)
    
    g3<- graph_from_data_frame(vertices = NodeList, d= EdgeList, directed = TRUE)
    
    
    par(mar=c(0,0,0,0)+.1)
    #curve_multiple(g3,start=0.1)
    #label_size=5,label.size=5
    p<-ggraph(g3, layout = 'linear') + 
      geom_edge_arc(aes(label=weight),strength = 0.5,arrow = arrow(length=unit(4,'mm')),
                    end_cap=circle(5,'mm'),alpha=0.3)+
      geom_node_point(size=8,colour=col)+geom_node_text(aes(label=labels) ,size=5)+
      theme_graph(foreground = 'steelblue', fg_text_colour = 'white',
                  border=FALSE,base_size = 5,title_size=12,base_family="Helvetica")+
      ggtitle(paste0(tt_igraph_list[[k]]$gene_symbol,"-",tt_igraph_list[[k]]$trans_id) );
    
    #geom_edge_fan(aes(alpha = stat(index)), show.legend = FALSE) + 
    #geom_node_point(aes(size = Popularity)) + 
    
    print(p);
    
    
  }
  
  dev.off();
}




