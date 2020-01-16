library(dplyr);
library(stringr);
library(igraph);


if_smaller_f<-function(region1,region2){
  sapply(strsplit(region1,":|-"),"[",2)->region1_start;
  sapply(strsplit(region2,":|-"),"[",2)->region2_start;
  
  region1_start<region2_start;
}

# if transcript id contain dot
# trim_trans_id_by_dot=TRUE

build_iso_object<-function(files_all,gene_trans_id_tbl,ucsc_intron_anno,is_large,
                           result_path,trim_trans_id_by_dot=TRUE){
  
  ##################################get intron splicing order from annotation##################################################
  
  gencode_intron<-read.table(ucsc_intron_anno,header = FALSE,sep="\t",as.is = TRUE);
  colnames(gencode_intron)<-c("chr","start","end","id","score","strand");
  
  
  #gencode_intron_o_frame<-data.frame(
  #  trans_id= sapply(strsplit(gencode_intron[,4],"_|\\.") ,"[",1),
  #  intron_id=  as.numeric(sapply(strsplit(gencode_intron[,4],"_|\\.") ,"[",5))+1
  #)
  gencode_intron[,"trans_id"]<-sapply(strsplit(gencode_intron[,4],"_") ,"[",1);
  
  
  if(  trim_trans_id_by_dot ){
    gencode_intron[,"trans_id"]<-sapply(strsplit(gencode_intron[,4],"_|\\.") ,"[",1);
  }
  
  gencode_intron[,"intron_order"]<-as.numeric(sapply(strsplit(gencode_intron[,4],"_|\\.") ,"[",5) )+1;
  
  gencode_intron_o_frame_summary<-gencode_intron %>% group_by(trans_id) %>% summarise(intron_count=max(intron_order) );
  
  gencode_intron_o_frame_intron_count<-gencode_intron_o_frame_summary$intron_count;
  names(gencode_intron_o_frame_intron_count)<-gencode_intron_o_frame_summary$trans_id;
  
  
  
  gencode_intron[,"max_intron"]<-gencode_intron_o_frame_intron_count[gencode_intron[,"trans_id"]];
  
  gencode_intron[gencode_intron$strand=="-","intron_order"]<-gencode_intron[gencode_intron$strand=="-","max_intron"]-
    gencode_intron[gencode_intron$strand=="-","intron_order"]+1;
  
  
  
  gencode_intron[,"gencode_intron_region"]<-str_c(gencode_intron[,"chr"],":",gencode_intron[,"start"]+1,"-",gencode_intron[,"end"]); 
  
  
  gencode_intron[,"gencode_intron_o"]<-str_c(
    gencode_intron[,"trans_id"],
    "_",
    "intron",
    "_",
    gencode_intron[,"intron_order"]
    
  );
  
  
  gencode_intron_o_frame_first<-data.frame(
    gencode_intron_o_first=gencode_intron[,"gencode_intron_o"],
    gencode_intron_region=gencode_intron[,"gencode_intron_region"],
    trans_id=gencode_intron[,"trans_id"],
    
    stringsAsFactors = FALSE
    
  )
  
  gencode_intron_o_frame_next<-data.frame(
    gencode_intron_o_next=gencode_intron[,"gencode_intron_o"],
    gencode_intron_region=gencode_intron[,"gencode_intron_region"],
    trans_id=gencode_intron[,"trans_id"],
    
    stringsAsFactors = FALSE
  )
  #names(gencode_intron_region)<-gencode_intron_o; 
  
  
  #######################################merge iso(intron splicing order) files##############################################
  
  #files_all<-list.files("/Users/mengli/Documents/projects/iso/data/iso/",full.names =TRUE,pattern = "*unique_intron.tsv");
  #gene_ids<-sapply(str_split(files_all,"\\_tab"),"[",1);
  #if(is_large){
  #files_all<-list.files("/Users/mengli/Documents/projects/iso/data/iso_large/",full.names =TRUE,pattern = "*unique_intron.tsv");
  #    files_all<-list.files("/Volumes/mengli/eric_nu/iso_large/",full.names =TRUE,pattern = "*unique_intron.tsv");
  #}
  
  #iso_raw<-read.table(paste0("data/iso/",files_all[1]),header = FALSE, sep = "\t", as.is = TRUE)
  
  iso_raw<-read.table(paste0(files_all[1]),header = FALSE, sep = "\t", as.is = TRUE)
  
  if(ncol(iso_raw)==6){
    iso_raw<-iso_raw[,c(1:4,6)];
  }
  colnames(iso_raw)<-c("id","nexti","first","strand","read_count");
  
  for(i in files_all[-1]){
    
    print(i);
    #if(file.size( paste0("data/iso/",i) )==0 || (!file.exists( paste0("data/iso/",i) ) ) ){
    if(file.size( paste0(i) )==0 || (!file.exists( paste0(i) ) ) ){
      next;
    }
    
    iso_tmp<-read.table(paste0(i),header = FALSE, sep = "\t", as.is = TRUE)
    if(ncol(iso_tmp)==6){
      iso_tmp<-iso_tmp[,c(1:4,6)];
    }
    colnames(iso_tmp)<-c("id","nexti","first","strand","read_count");
    
    
    iso_raw<-rbind(iso_raw,iso_tmp);
    
  }
  
  colnames(iso_raw)<-c("id","nexti","first","strand","read_count");
  
  
  print("group togeter iso files ");
  
  library(dplyr);
  iso<-as.data.frame(iso_raw %>% group_by(id,nexti,first,strand) %>% summarise(read_count=sum(read_count) ) );
  
  #colnames(iso)<-c("id","nexti","first","strand","read_count");
  
  
  iso<-iso[iso[,"read_count"]>0,];
  
  if(trim_trans_id_by_dot){
    iso[,"id"]<-sapply(strsplit(iso[,"id"],"\\."),"[",1);
  }
  
  gene_id_trans_id_map<-read.table(gene_trans_id_tbl,
                                   header = FALSE,as.is = TRUE,sep="\t",quote="$",comment.char ="#")[,1:5];
  colnames(gene_id_trans_id_map)<-c("gene_id","trans_id","gene_symbol","trans_start","trans_end");
  
  
  iso<-left_join(iso, gene_id_trans_id_map,by=c("id"="trans_id") );
  
  
  ##only annotated introns
  iso<-iso[(iso[,"nexti"] %in% gencode_intron[,"gencode_intron_region"]) &
             (iso[,"first"] %in% gencode_intron[,"gencode_intron_region"]),];
  # colnames(iso)<-c("id","nexti","first","strand","read_count","gene_id");
  
  
  iso_final<-c();
  
  
  
  if(is_large){
    iso_final<-inner_join(iso,gencode_intron_o_frame_next,by=c("nexti"="gencode_intron_region","id"="trans_id")) %>%
      inner_join(gencode_intron_o_frame_first,by=c("first"="gencode_intron_region","id"="trans_id") );
    
    iso_final<-iso_final[(!is.na(iso_final$nexti))&(!is.na(iso_final$first)) & 
                           (!is.na(iso_final$id)) & (!is.na(iso_final$gencode_intron_o_first)) &
                           (!is.na(iso_final$gencode_intron_o_next)),];
    
    ###only keep  iso trans id eaul to intron trans id 
    iso_final<-iso_final[(sapply(str_split(iso_final[,"gencode_intron_o_first"],"_"),"[",1 )== iso_final[,"id"] ) & 
                           (sapply(str_split(iso_final[,"gencode_intron_o_next"],"_"),"[",1 )== iso_final[,"id"] ) ,];
    
    
  }else{
    ###mark the intron number based on annotation for get coordinated based intron splicing order
    ###should be full_join, but too slow
    iso_final<-full_join(iso,gencode_intron_o_frame_next,by=c("nexti"="gencode_intron_region")) %>%
      full_join(gencode_intron_o_frame_first,by=c("first"="gencode_intron_region") );
    
    ##this method join the intron splicing order file and annotation by intron coordinates 
    ## the trans id got from intron splicing order is for a specific isoform, 
    ##so it needs to be replaced by the same coordinate of other isoform for further adj matrix
    iso_final[,"id"]<-sapply(str_split(iso_final[,"gencode_intron_o_first"],"_"),"[",1 );
    
  }
  
  
  iso_final<-iso_final[iso_final$gencode_intron_o_first!=iso_final$gencode_intron_o_next,];
  
  
  
  
  #iso_final<-iso_final[(str_sub(iso_final[,"gencode_intron_o_first"],1,15)== iso_final[,"id"] ) & 
  #                       (str_sub(iso_final[,"gencode_intron_o_next"],1,15)== iso_final[,"id"] ) ,];
  
  
  ####################################################mark intron in order or not############################################
  
  
  if_smaller<-if_smaller_f(iso_final[,"nexti"],iso_final[,"first"] );
  
  is_slow<-((iso_final[,"strand"]=="+") &if_smaller) | ((iso_final[,"strand"]=="-") & (!if_smaller) );
  
  is_fast<-((iso_final[,"strand"]=="-") &if_smaller) | ((iso_final[,"strand"]=="+") & (!if_smaller) ) ;
  
  
  is_fast<-(is_fast&(!is_slow));
  
  
  print(paste0("# of not in order introns: ",sum(is_slow) ) );
  print(paste0("# of in order introns: ",sum(is_fast) ) );
  
  ##is fast
  # iso_final[((iso_final[,"strand"]=="-") &if_smaller) | ((iso_final[,"strand"]=="+") & (!if_smaller)),"is_low"]<-FALSE;
  iso_final[is_fast,"is_low"]<-FALSE
  
  ##is slow, slow overlap fast.
  # iso_final[((iso_final[,"strand"]=="+") &if_smaller) | ((iso_final[,"strand"]=="-") & (!if_smaller)) ,"is_low"]<-TRUE;
  
  iso_final[is_slow,"is_low"]<-TRUE
  
  
  #write.table(iso_final,file="result/iso_large_final_unique.tsv",sep="\t",quote = FALSE,col.names = TRUE,row.names = FALSE);
  
  
  #iso_final[,"is_fast"]<-is_fast;
  
  
  #iso_slow_sumary<- as.data.frame(iso_final %>% group_by(id) %>% summarise(count_slow=sum(is_low,na.rm = TRUE),
  #                                                                         count_fast=sum(!is_low,na.rm = TRUE),
  #                                                                         trans_len=mean(trans_end)-mean(trans_start))  );
  
  
  
  ##################add transcript start, end position and edge count in the iso_slow_sumary object############################
  
  iso_slow_sumary2<- as.data.frame(iso_final %>% group_by(gene_symbol,id,nexti,trans_start,trans_end) %>%
                                     summarise(is_low2 = (sum(is_low,na.rm=TRUE)>0)  ) );
  
  
  print(paste0("# of not in order introns after correction ") );
  
  print(paste0("# of not in order introns: ",sum(iso_slow_sumary2$is_low2) ) );
  print(paste0("# of in order introns: ",sum(!iso_slow_sumary2$is_low2) ) );
  
  
  #cat(iso_slow_sumary2[iso_slow_sumary2$is_low2,"nexti"],sep="\n",file="result/slow_intron.txt");
  
  #cat(iso_slow_sumary2[!(iso_slow_sumary2$is_low2),"nexti"],sep="\n",file="result/fast_intron.txt");
  
  
  iso_slow_sumary<-as.data.frame(iso_slow_sumary2 %>% group_by(gene_symbol,id) %>%
                                   summarise(count_slow=sum(is_low2,na.rm = TRUE),
                                             count_fast=sum(!is_low2,na.rm = TRUE),
                                             trans_len=mean(trans_end)-mean(trans_start))  );
  
  #iso_slow_sumary[,"int_count"]<-gencode_intron_o_frame_intron_count[iso_slow_sumary[,"id"]]
  
  iso_slow_sumary[,"int_count"]<-gencode_intron_o_frame_intron_count[iso_slow_sumary[,"id"]];
  
  
  
  #number of connection 
  iso_edge_count<-as.data.frame(iso_final %>% group_by(gene_symbol,id) %>%
                                  summarise(edge_count=n()) );
  
  
  
  iso_slow_sumary<-inner_join(iso_slow_sumary,iso_edge_count,by=c("gene_symbol"="gene_symbol","id"="id") );
  
  
  #gencode_intron_o_frame_intron_count[iso_slow_sumary[1:5,"id"] ];
  
  trans_slow_intr_percent<-iso_slow_sumary$count_slow/(iso_slow_sumary$count_slow+iso_slow_sumary$count_fast);
  
  #summary(trans_slow_intr_percent);
  
  
  #iso_slow_sumary<-iso_slow_sumary[order(iso_slow_sumary$count_fast+iso_slow_sumary$count_slow,decreasing = TRUE),];
  
  #iso_slow_sumary<-iso_slow_sumary[order((iso_slow_sumary$count_slow+iso_slow_sumary$count_fast)/
  #(iso_slow_sumary$trans_len),decreasing = TRUE),];
  
  
  
  #################################################output and saving##################################################
  
  #library(reshape2);
  
  #library(igraph);
  
  #  t_iso_slow_sumary<-iso_slow_sumary
  #  t_iso_final<-iso_final
  #  label<-"test"
  #  count=1
  #  read_count_t=2
  #  g<-2
  #source("code/draw.R")
  #order by post evidence
  
  #if_uniq<-NULL;
  
  # is_large
  #if(is_large){
  #  if_uniq="not_uniq_intron"
  #}else{
  #  if_uniq="uniq_intron"
  #}
  
  if(!exists("cosi_raw_group_vector_format")){
    cosi_raw_group_vector_format<-NULL;
  }
  
  if(!exists("posp_introns")){
    posp_introns<-NULL
  }
  
  if(!exists("all_alt_se_introns")){
    all_alt_se_introns<-NULL
  }
  
  if(!exists("all_alt_a5_introns")){
    all_alt_a5_introns<-NULL
  }
  
  if(!exists("all_alt_a3_introns")){
    all_alt_a3_introns<-NULL
  }
  
  if(!exists("all_alt_ri_introns")){
    all_alt_ri_introns<-NULL
  }
  
  
  save(iso_final,iso_slow_sumary,gencode_intron_o_frame_intron_count,
       all_alt_se_introns,all_alt_a5_introns,all_alt_a3_introns,
       all_alt_ri_introns,
       cosi_raw_group_vector_format,posp_introns,
       is_large,
       file=result_path);
  
  #iso_slow_sumary$gene_symbol
  
  #pdf("result/test_igraph.pdf")
  #adjacency_matrix_t<-matrix(c(1,1,0,0),nrow = 2,byrow=TRUE);
  
  #g<-graph_from_adjacency_matrix(adjacency_matrix_t,mode="directed",weighted=TRUE);
  
  #plot(g)
  
  #dev.off();
  
}



