library(shiny);
library(networkD3);
library(igraph);
library(stringr);
library(Sushi);
library(readr);
library(DT);

source("code/get_members.R");
source("code/mlp3.R");

###############################################################human##################################################

#load(file="result/t_igraph_list.Rd");

#load(file="result/gene_trans_id_map.Rd");

#load(file="anno/sushi_trans_file.Rd");

######################################################################################################################

###############################################################pombe##################################################
#setwd("/Users/mengli/Documents/projects/iso");

#load(file="pombe/t_igraph_list.Rd");

#load(file="pombe/sushi_trans_file.Rd");

######################################################################################################################

all_gene_symbol_graph<-as.character(sapply(t_igraph_list,"[[","gene_symbol") );

all_trans_id_graph<-as.character(sapply(t_igraph_list,"[[","trans_id") );
#all_trans_id_graph
#all_gene_symbol_graph


init_gene_symbol<-t_igraph_list[[2]]$gene_symbol;
#save(gene_trans_id_map,t_igraph_list,file="result/t_igraph_list.Rd");

#gene_trans_id_map<-read.table("anno/hg19_ensembl_gene_id_trans_id_map.tsv",
#                              header = FALSE,as.is = TRUE,sep = "\t");

#colnames(gene_trans_id_map)<-c("gene_id","trans_id","gene_symbol","trans_start","trans_end","strand",
#                               "chr","gene_start","gene_end");

#options(shiny.port = 7775);

#options(shiny.host = "10.10.115.115");

filter_matrix_read_count<-function(t_adj_matrix,t_read_count){
  for(i in 1:nrow(t_adj_matrix)){
    for(j in 1:ncol(t_adj_matrix)){
      if(t_adj_matrix[i,j]+t_adj_matrix[j,i]<t_read_count){
        t_adj_matrix[i,j]<-
          t_adj_matrix[j,i]<-0;
        
      }
      
    }
    
  }
  return(t_adj_matrix);
}



#load("result/all_iso_data_uniq.Rd");

#  t_igraph_list


server <- function(input, output, session) {
  
  
  #m_charge<- -30;#<-input$charge;
  #m_linkDistance<-50;#;<-input$linkDistance;
  #m_nodesize<-1;#<-input$nodesize;
  #m_fontsize<-10;#<-input$fontsize;
  #m_trans<-"ENST00000541578";
  print("init");
  
  re_gra<-function(){
  }
  
  #t_net_one<-NULL;
  #cur_graph_content<-NULL;
  
  
  observe({
    
    g<-input$gene;
    #print(g);
    
    trans_ids<-gene_trans_id_map[gene_trans_id_map$gene_symbol==g,"trans_id"];
    #print(trans_ids);
    trans_ids<-sort(intersect(trans_ids,all_trans_id_graph) );
    
    updateSelectInput(session, "transcript",
                      label = paste0("Transcript"),
                      choices = trans_ids,
                      selected = trans_ids[1] )
    
  });
  
  
  observe({
    #m_charge<<-input$charge;
    #m_linkDistance<<-input$linkDistance;
    #m_nodesize<<-input$nodesize;
    #m_fontsize<<-input$fontSize;
    
    #print(paste0("observe",m_fontsize) );
    #output$g2<-renderForceNetwork( data()$net );
    
    #cur_graph_content<-re_gra();
    
    #cur_graph_content
  });
  
  
  
  data<-eventReactive(input$update,{
    
    # ENST00000368575;
    #m_trans<<-input$transcript;
    #cur_graph_content<<-re_gra();
    
    #print(cur_graph_content)
    #global_g_content<<-cur_graph_content;
    #print(m_charge);
    
    #cur_graph_content
    
  });
  
  
  output$trans_struc<-renderImage({
    
    
    #m_trans<<-input$transcript;
    g<-input$gene;
    #print(g);
    m_trans<-input$transcript;
    
    #trans_ids<-gene_trans_id_map[gene_trans_id_map$gene_symbol==g,"trans_id"];
    
    #  g<-"SNHG1"
    
    chrom = paste0("chr", unique(gene_trans_id_map[gene_trans_id_map$gene_symbol==g,"chr"]) );
    
    #chromstart = min(gene_trans_id_map[gene_trans_id_map$gene_symbol==g,"gene_start"]);
    
    #chromend = max(gene_trans_id_map[gene_trans_id_map$gene_symbol==g,"gene_end"]);
    
    #  m_trans<-"ENST00000282074"
    
    chromm = paste0("chr", unique(gene_trans_id_map[gene_trans_id_map$trans_id==m_trans,"chr"]) );
    
    chromstart = min(gene_trans_id_map[gene_trans_id_map$trans_id==m_trans,"trans_start"]);
    
    chromend = max(gene_trans_id_map[gene_trans_id_map$trans_id==m_trans,"trans_end"]);
    
    
    outfile <- tempfile(fileext=".png");
    
    chromstart<-as.numeric(chromstart)-10;
    chromend<-as.numeric(chromend)+10;
    
    
    for_sushi_tmp<-for_sushi[(for_sushi$start>=chromstart) &
                               (for_sushi$end<=chromend)&
                               (for_sushi$chrom==chromm), ];
    
    print(paste0(chromm,":",chromstart,"-",chromend) );
    #View(for_sushi_tmp);
    
    #a<<-for_sushi_tmp
    
    #outfile<-"result/"
    png(outfile,width = 900,height=350, res=100 );
    #png(outfile);
    par( mar=c(0,0,2,0) );
    
    pg = plotGenes(for_sushi_tmp,chrom,chromstart,chromend,types=for_sushi_tmp$type,
                     bheight =0.2, plotgenetype = "box",bentline=FALSE,
                   labeloffset = .4,labeltext = TRUE,fontsize=0.5,packrow=TRUE);
    
    labelgenome(chrom, chromstart, chromend, side=3, n=3, scale="Mb");
    
    dev.off();
    
    list(src = outfile,
         alt = "No trans")
    
  })
  
  
  output$mlp<-renderText( {
    
    m_trans<-input$transcript;
    
    m_read_count<-input$readCount
    
    m_alpha<-input$alpha;
    
    m_adj_matris<-filter_matrix_read_count(t_igraph_list[[m_trans]]$adjacency_matrix,m_read_count)
    
    #t_alpha<-0.1;
    best_order_ls<-find_path_global(m_adj_matris, m_alpha);
    
    
    #t_igraph_list[[i]]$best_order<-best_order_ls$best_order;
    
    #t_igraph_list[[i]]$p_value<-best_order_ls$p_value;
    
    orders=best_order_ls$best_order;
    p_value<-best_order_ls$p_value;
    
    #print(orders)
    #orders=t_igraph_list[[m_trans]]$best_order;
    
    paste0("Most likely order: ",paste0( orders,collapse = "-->" ),
           paste0( "\n\n in order spliced P-value (log): ",format(p_value ))
    )
    
    });
  
  
  
  output$mlp2<-renderText( {
    m_trans<-input$transcript;
    
    p_value=t_igraph_list[[m_trans]]$p_value;
    
    m_alpha<-input$alpha;
    
    #m_adj_matris<-filter_matrix_read_count(t_igraph_list[[m_trans]]$adjacency_matrix,m_read_count)
    
    #t_alpha<-0.1;
    #best_order_ls<-find_path_global(m_adj_matris, m_alpha);
    
    
    #t_igraph_list[[i]]$best_order<-best_order_ls$best_order;
    
    #t_igraph_list[[i]]$p_value<-best_order_ls$p_value;
    
    
    
    paste0( "in order spliced P-value (log): ",format(p_value ))
    
    });
  
  
  output$downloadData <- downloadHandler(
  
    
  filename = function() {
    g<-input$gene;
    
    m_trans<<-input$transcript;
    
    paste(g,"_",m_trans, ".tsv", sep = "")
  },
  
  content = function(file_a) {
    
    g<-input$gene;
    
    m_trans<<-input$transcript;
    
    
    adj<-t_igraph_list[[m_trans]]$adjacency_matrix;
    
    
    m_trans<-input$transcript;
    
    p_value=t_igraph_list[[m_trans]]$p_value;
    
    p_value_str<-paste0( "in order spliced P-value (log): ",p_value) 
    
    
    
    m_trans<-input$transcript;
    
    orders<-t_igraph_list[[m_trans]]$best_order;
    
    orders_str<-paste0("Most likely order: ",paste0( orders,collapse = "-->" ) ) 
    
    cat( paste0(p_value_str,"\n",orders_str,"\n"),file=file_a,append = FALSE);
    
    write.table(adj, file=file_a, row.names = TRUE,sep="\t",col.names = TRUE,append = TRUE);
    
  } )
  
  
  output$mytable = DT::renderDataTable({
    m_trans<<-input$transcript;
    
    adj<-t_igraph_list[[m_trans]]$adjacency_matrix;
    
    adj
  })
  
  
  output$g2<-renderForceNetwork( {
    
    m_charge<<-input$charge;
    m_linkDistance<<-input$linkDistance;
    m_nodesize<<-input$nodesize;
    m_fontsize<<-input$fontSize;
    
    m_trans<<-input$transcript;
    
    m_read_count<-input$readCount
    
    m_adj_matris<-filter_matrix_read_count(t_igraph_list[[m_trans]]$adjacency_matrix,m_read_count)
    
    m_alpha<-input$alpha;
    #best_order_ls<-find_path_global(m_adj_matris, t_alpha);
    
    
    #t_igraph_list[[i]]$best_order<-best_order_ls$best_order;
    
    #t_igraph_list[[i]]$p_value<-best_order_ls$p_value;
    
    
    
    
    #t_net_one<-t_igraph_list[[m_trans]]$one_net;
    
    g2<-graph_from_adjacency_matrix(m_adj_matris,mode="directed",weighted=TRUE);
    #mermbers<-t_igraph_list[[i]]$members;
    mermbers<-get_members_matrix(m_adj_matris,m_alpha);
    
    t_net_one<-igraph_to_networkD3(g2,mermbers);
    
    
    #edge_width=map( E(g2)$weight, c(0.25,6) );
    t_net_one$nodes[,"node_size"]<-rep(15,length(V(g2) ));
    
    #t_net_one$nodes[,"node_size"]<-tt_igraph_list[[i]]$node_size[as.character(t_net_one$nodes$name)];
    #t_net_one$nodes[,"node_color"]<-tt_igraph_list[[i]]$node_color[as.character(t_net_one$nodes$name)];
    
    
    #t_net_one$links
    t_net_one$links$value<-t_net_one$links$value/(max(t_net_one$links$value))*10;
    
    
    
    #print(m_trans);
    
    #if(is.null(t_net_one)){
    #  return(cur_graph_content);
    #}
    
    #cur_trans_id<-input$transcript;
    #print(paste0("aaa",m_trans) );
    par(mar=c(0,0,0,0));
    # orders=t_igraph_list[[m_trans]]$best_order,
    # p_value=t_igraph_list[[m_trans]]$p_value,
    #orders=best_order_ls$best_order,
    #p_value=-best_order_ls$p_value,
    
    t<-list(
            net=forceNetwork(Links = t_net_one$links, Nodes = t_net_one$nodes,
                             Source = 'source', Target = 'target', Value='value', NodeID = 'name',
                             Group = 'group',arrows=TRUE,Nodesize='node_size',
                             opacity=1, bounded=TRUE, legend=TRUE,linkDistance=m_linkDistance,
                             charge = m_charge, zoom=FALSE, opacityNoHover = 0.9,
                             radiusCalculation = JS( paste0(" Math.sqrt(d.nodesize)*",m_nodesize ) ),
                             fontSize=m_fontsize)
            
    );
    #print("dfdf");
    
    
    return(t$net);
    
    
  } );
  
  
}


ui <- fluidPage(
  
  #fontSize
  #titlePanel("Intron splicing order"),
  
  textOutput("mlp"), 
  #textOutput("mlp2"),
  #hr(),
  
  
  sidebarLayout(
    sidebarPanel(
      # init_gene_symbol
      #selectInput("gene","gene symbol",choices=unique(gene_trans_id_map$gene_symbol),selected="RAB2" ),
      #textInput("gene","gene symbol",value="MRPL15" ),
      
      #selectInput("transcript","Transcript", choices=gene_trans_id_map[gene_trans_id_map$gene_symbol=="SNHG1","trans_id"],selected="ENST00000260102"),
      
      # all_trans_id_graph
      # all_gene_symbol_graph
      
      #textInput("gene","gene",value=init_gene_symbol ),
      
      selectInput("gene","gene", choices=sort(all_gene_symbol_graph),
                  selected=init_gene_symbol ),
      
      
      selectInput("transcript","Transcript", choices=sort(gene_trans_id_map[gene_trans_id_map$gene_symbol==init_gene_symbol,"trans_id"]),
                  selected=gene_trans_id_map[gene_trans_id_map$gene_symbol==init_gene_symbol,"trans_id"][1] ),
      
      
      sliderInput("readCount", "readCount ",0 , min = 0,
                  max = 50, step = 1),
      sliderInput("alpha", "alpha ",0.1 , min = 0,
                  max = 1, step = 0.01),
      #actionButton("update","update"),
      sliderInput("charge", "charge ", -30, min = -300,
                  max = 0, step = 1),
      sliderInput("linkDistance", "linkDistance ",200 , min = 0,
                  max = 300, step = 1),
      sliderInput("nodesize", "nodesize ",2 , min = 0,
                  max = 20, step = 1),
      sliderInput("fontSize", "fontSize ",20 , min = 1,
                  max = 50, step = 1),
      
      downloadButton("downloadData", "Download matrix"),
      
      width=3),
    
    mainPanel(
      #tabPanel(
      
      tabsetPanel(
        
        tabPanel("Force Network", forceNetworkOutput("g2",height="400px") ),
        
        
        #             ),
        #tabPanel("trans ",
        tabPanel("read count adjacent matrix", DT::dataTableOutput("mytable") )
        
        
        
      ),
      tabPanel("transcript structure",plotOutput("trans_struc",height="200px") )
      
      
      #tabPanel("Most likely order P-value", ),
      
      #tabPanel("browser",includeHTML("shiny/bro.html") );
      
      
      # )
      
    )
  )
  
)


shinyApp(ui = ui, server = server);

#includeHTML("shiny/bro.html")


