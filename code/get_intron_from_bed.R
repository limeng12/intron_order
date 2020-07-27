library(stringr)
#gencode_intron<-read.table(ucsc_intron_anno,header = FALSE,sep="\t",as.is = TRUE);
#colnames(gencode_intron)<-c("chr","start","end","id","score","strand");
#t_bed_file<-"anno/hg19_gencode_from_ucsc.bed"


get_exon_intron_from_bed<-function(t_bed_file,
                              m_trim_trans_id_by_dot=TRUE){
  
  
  bed_anno<-read.table(t_bed_file,header = FALSE,as.is = TRUE);
  colnames(bed_anno)<-c("chr","start","end","trans_id","score","strand","CDS_start","CDS_end",
                                   "","exon_count","exon_len","exon_start")
  
  bed_anno<-bed_anno[!duplicated(bed_anno$trans_id),];
  
  exon_pos_mat<-matrix(nrow = sum(bed_anno[,"exon_count"]),ncol=2);
  exon_pos_index<-1
  
  intron_pos_mat<-matrix(nrow = sum(bed_anno[,"exon_count"]-1),ncol = 9);
  intron_pos_index<-1
  
  exon_list<-str_split(bed_anno[,"exon_len"], ",")
  exon_exon_list<-str_split(bed_anno[,"exon_start"], ",")
  
  #if(  m_trim_trans_id_by_dot ){
  #  transcript_id_list<-str_split(bed_anno[,"trans_id"],"_|\\.") ;
  #}
  #transcript_id_list<-bed_anno
  
  cat("\n");
  
  print("Get exons and introns from BED file");
  
  pb = txtProgressBar(min = 1, max = nrow(bed_anno), initial = 0, width=100, style=3) 
  
  for (row_num in 1:nrow(bed_anno) ) {
    
    setTxtProgressBar(pb,row_num);
    
    
    exon<-as.numeric( exon_list[[row_num]] );
    exon_pos<-as.numeric( exon_exon_list[[row_num]] );
    
    ##remove , after last exons
    exon<-exon[-1*length(exon)];
    exon_pos<-exon_pos[-1*length(exon_pos)];
    
    exon_count<-length(exon)
    
    if(exon_count<=2)
    {
      next
    }
    
    
    transcript_id<-bed_anno[row_num,"trans_id"];
    #if(  trim_trans_id_by_dot ){
    #  transcript_id<-transcript_id_list[[row_num]][1] ;
    #}
    #print(paste0(row_num,":",transcript_id) );
    
    transcript_start_site<-bed_anno[row_num,"start"];
    #transcript_end_site<-bed_anno[row_num,"end"];
    strand<-bed_anno[row_num,"strand"];
    chr<-bed_anno[row_num,"chr"];
    
    
    for (i in 1:(length(exon)-1)) {
      start<-transcript_start_site+(exon_pos[i])+(exon[i])+1;
      end<-(exon_pos[i+1])+transcript_start_site+1-1;
      
      if(strand=="+"){
        id<-str_c(transcript_id,"_",i)
        intron_index<-i
      }else{
        id<-str_c(transcript_id,"_",(length(exon)-i) )
        intron_index<-(length(exon)-i)
      }
      
      #intron_pos_mat<-rbind(intron_pos_mat,c(chr,start,end,id,-1,strand) );
      intron_pos_mat[intron_pos_index,]<-c(chr,start,end,id,-1,strand,transcript_id,
                                           intron_index,str_c(chr,":",start,"-",end));
      intron_pos_index<-intron_pos_index+1;
    }
    
    for (i in 1:(length(exon) )) {
      
      end<-transcript_start_site+(exon_pos[i])+(exon[i])+1-1;
      start<-(exon_pos[i])+transcript_start_site+1-1+1;
      
      
      exon_pos_mat[exon_pos_index,]<-c(transcript_id,str_c(chr,":",start,"-",end));
    
      exon_pos_index<-exon_pos_index+1;
    }
    
  }
  
  exon_pos_mat_fr<-as.data.frame(exon_pos_mat, stringsAsFactors=FALSE);
  colnames(exon_pos_mat_fr)<-c("trans_id","exon_pos")
  
  cat("\n");
  
  print(paste0("Total number of introns in the annotation: ", intron_pos_index));
  
  intron_pos_mat<-intron_pos_mat[1:(intron_pos_index-1),];
  
  intron_pos_mat_fr<-as.data.frame(intron_pos_mat,stringsAsFactors = FALSE);
  
  colnames(intron_pos_mat_fr)<-c("chr","start","end","id","score","strand","trans_id","intron_order","intron_pos");
  
  intron_pos_mat_fr$intron_order<-as.numeric(intron_pos_mat_fr$intron_order);
  intron_pos_mat_fr$start<-as.numeric(intron_pos_mat_fr$start);
  intron_pos_mat_fr$end<-as.numeric(intron_pos_mat_fr$end);
  
  list(intron=intron_pos_mat_fr, exons=exon_pos_mat_fr)
  
  
}



