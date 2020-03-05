library(dplyr);
library(stringr);
library(igraph);


# if transcript id contain dot
# trim_trans_id_by_dot=TRUE

build_iso_object2<-function(files_all,intron_anno,
                           read_count_threshold =0, trans_exp_file=""){
  
  ##################################get intron splicing order from annotation##################################################
  is_large=TRUE
  
  #gencode_intron<-read.table(ucsc_intron_anno,header = FALSE,sep="\t",as.is = TRUE);
  #colnames(gencode_intron)<-c("chr","start","end","id","score","strand");
  
  
  #gencode_intron[,"trans_id"]<-sapply(strsplit(gencode_intron[,4],"_") ,"[",1);
  
  
  #if(  trim_trans_id_by_dot ){
  #  gencode_intron[,"trans_id"]<-sapply(strsplit(gencode_intron[,4],"_|\\.") ,"[",1);
  #}
  
  #gencode_intron[,"intron_order"]<-as.numeric(sapply(strsplit(gencode_intron[,4],"_|\\."),
  #                                                   "[",intron_index_in_intron_file) )+1;
  
  

  #intron_anno[,"max_intron"]<-gencode_intron_o_frame_intron_count[intron_anno[,"trans_id"]];
  
  #intron_anno[intron_anno$strand=="-","intron_order"]<-intron_anno[gencode_intron$strand=="-","max_intron"]-
  #  intron_anno[intron_anno$strand=="-","intron_order"]+1;
  
  intron_anno[,"gencode_intron_region"]<-str_c(intron_anno[,"chr"],":",
                                                  intron_anno[,"start"],"-",intron_anno[,"end"]); 
  
  
  intron_anno[,"gencode_intron_o"]<-str_c(
    intron_anno[,"trans_id"],
    "_",
    "intron",
    "_",
    intron_anno[,"intron_order"]
    
  );
  
  
  intron_o_frame_first<-data.frame(
    gencode_intron_o_first=intron_anno[,"gencode_intron_o"],
    gencode_intron_region=intron_anno[,"gencode_intron_region"],
    trans_id=intron_anno[,"trans_id"],
    intron_order_first=intron_anno[,"intron_order"],
    
    stringsAsFactors = FALSE
    
  )
  
  intron_o_frame_next<-data.frame(
    gencode_intron_o_next=intron_anno[,"gencode_intron_o"],
    gencode_intron_region=intron_anno[,"gencode_intron_region"],
    trans_id=intron_anno[,"trans_id"],
    intron_order_next=intron_anno[,"intron_order"],
    
    stringsAsFactors = FALSE
    
  )
  
  #######################################merge iso(intron splicing order) files##############################################
  ## read iso files into iso_raw
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
  
  
  ## group togeter iso files 
  iso<-as.data.frame(iso_raw %>% dplyr::group_by(id,nexti,first,strand) %>% dplyr::summarise(read_count=sum(read_count) ) );
  
  
  ##not filter by read count
  iso<-iso[iso[,"read_count"]>=read_count_threshold,];
  
  if(trim_trans_id_by_dot){
    iso[,"id"]<-sapply(strsplit(iso[,"id"],"\\."),"[",1);
  }
  
  if(trans_exp_file!=""){
    exp_trans<-readLines(trans_exp_file)
    
    iso<-iso[iso$id %in% exp_trans,];
  }
  
  ## add gene symbol, trans start, trans end information
  # gene_id_trans_id_map<-read.table(gene_trans_id_tbl,
  #                                  header = FALSE,as.is = TRUE,sep="\t",quote="$",comment.char ="#")[,1:5];
  # colnames(gene_id_trans_id_map)<-c("gene_id","trans_id","gene_symbol","trans_start","trans_end");
  # gene_id_trans_id_map$trans_start<-as.numeric(gene_id_trans_id_map$trans_start)
  # gene_id_trans_id_map$trans_end<-as.numeric(gene_id_trans_id_map$trans_end)
  # iso<-left_join(iso, gene_id_trans_id_map,by=c("id"="trans_id") );
  ##keep only annotated introns
  #  gencode_intron
  iso<-iso[(iso[,"nexti"] %in% intron_anno[,"gencode_intron_region"]) &
             (iso[,"first"] %in% intron_anno[,"gencode_intron_region"]),];
  
  
  ##mark intron orders in iso as iso_final
  iso_final<-c();
  
  if(is_large){
    iso_final<-inner_join(iso,intron_o_frame_next,by=c("nexti"="gencode_intron_region","id"="trans_id")) %>%
      inner_join(intron_o_frame_first,by=c("first"="gencode_intron_region","id"="trans_id") );
    
    iso_final<-iso_final[(!is.na(iso_final$nexti))&(!is.na(iso_final$first)) & 
                           (!is.na(iso_final$id)) & (!is.na(iso_final$gencode_intron_o_first)) &
                           (!is.na(iso_final$gencode_intron_o_next)),];
    
    ###only keep iso trans id equal to intron trans id 
    iso_final<-iso_final[(sapply(str_split(iso_final[,"gencode_intron_o_first"],"_"),"[",1 )== iso_final[,"id"] ) & 
                           (sapply(str_split(iso_final[,"gencode_intron_o_next"],"_"),"[",1 )== iso_final[,"id"] ) ,];
    
    
  }else{
    ###mark the intron number based on annotation for get coordinated based intron splicing order
    ###should be full_join, but too slow
    iso_final<-full_join(iso,intron_o_frame_next,by=c("nexti"="gencode_intron_region")) %>%
      full_join(intron_o_frame_first,by=c("first"="gencode_intron_region") );
    
    ##this method join the intron splicing order file and annotation by intron coordinates 
    ## the trans id got from intron splicing order is for a specific isoform, 
    ##so it needs to be replaced by the same coordinate of other isoform for further adj matrix
    iso_final[,"id"]<-sapply(str_split(iso_final[,"gencode_intron_o_first"],"_"),"[",1 );
    
  }
  
  ###only keep intron pairs that are different
  iso_final<-iso_final[iso_final$gencode_intron_o_first!=iso_final$gencode_intron_o_next,];
  
  print(paste0("Total intron splicing order pairs: ", nrow(iso_final)) );
  

  iso_final;
  
}


