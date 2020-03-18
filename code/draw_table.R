library(reshape2)

# k<-1

draw_mol_table_plot<-function(t_igraph_list,output_path){
  
  int_num<-sapply(t_igraph_list,function(x){
    
    length(x$best_order) 
  });
  
  t_igraph_list<-t_igraph_list[order(int_num,decreasing = TRUE)];
  
  # output_path<-"result/test_fre_table.pdf";
  
  pdf(output_path,width = 12,height=10);
  
  # 
  for(k in 1:min(length(t_igraph_list), .Machine$integer.max) ){
    
    print(paste0(t_igraph_list[[k]]$gene_symbol,":",t_igraph_list[[k]]$trans_id,":",k) );
    
    
    adj<-t_igraph_list[[k]]$adjacency_matrix
    
    orders<-t_igraph_list[[k]]$best_order;
    
    t_alpha_v<-0.1;
  
    t_adj_mat_li<-adj+t_alpha_v;
    
    for(i in 1:nrow(t_adj_mat_li)){
      for(j in 1:ncol(t_adj_mat_li)){
        if(i>j){
          p<-t_adj_mat_li[i,j]/(t_adj_mat_li[i,j]+t_adj_mat_li[j,i])
          
          t_adj_mat_li[i,j]<-(p);
          t_adj_mat_li[j,i]<-(1-p);
          
        }
        if(i==j){
          t_adj_mat_li[i,j]<-0
        }
        
      }
      
    }
    

    
    t_adj_mat_li_format_tmp<-t_adj_mat_li[(orders), orders];
    
    t_adj_mat_li<-t_adj_mat_li_format_tmp;
    
    #  aaa<-matrix(c(rep(1,5),rep(2,5),rep(3,5),rep(4,5),rep(5,5)),nrow=5)
    
    colnames(t_adj_mat_li)<-1:length(orders)
    rownames(t_adj_mat_li)<-1:length(orders)
    t_adj_mat_li_melt<-melt(t_adj_mat_li);
    
    t_adj_mat_li_melt[,"fre"]<- as.factor(t_adj_mat_li_melt[,"value"]>0.5)
    levels(t_adj_mat_li_melt[,"fre"])<-c("<0.5", ">0.5")
    

    t_adj_mat_li_melt[,"value"]<-format(t_adj_mat_li_melt[,"value"],digits = 2);
    
    # trans="reverse
    p<-ggplot(t_adj_mat_li_melt ) + 
      geom_tile( aes(x = Var2, y = Var1, fill = fre), color = "black", size = 0.5,height=1,width=1)+
      geom_text( aes(x = Var2, y = Var1, label = value))+
      scale_x_continuous(expand = c(0, 0),breaks=1:length(orders),labels = orders,position="top") +
      scale_y_continuous(expand = c(0, 0),trans="reverse",breaks=1:length(orders),labels=orders)+
      theme(
            plot.title = element_text(size = rel(1.2)),
            #axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            legend.title = element_blank(),
            legend.position = "right")+
      ggtitle(paste0(t_igraph_list[[k]]$gene_symbol,"-",t_igraph_list[[k]]$trans_id) )+
      scale_fill_brewer(palette = "Set3");
      
    
    print(p);
  }
  
  dev.off();
  
}