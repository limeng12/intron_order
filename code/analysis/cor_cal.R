#setwd("/Users/mengli/Documents/projects/iso");

#setwd("D:\\SIBS\\Wangzefeng Lab\\project\\")
library(stringr)
library(ggplot2)
library(readr)

source("code/analysis/get_score_seq.R")

#species<-"pombe"
#species<-"fly"
#species<-"human"

#Schizosaccharomyces_pombe.ASM294v2.43.chr_nothick.bed

cal_intron_metric_cor<-function(anno_save_load,best_order_path,output_label="",rebuild_anno=FALSE,ucsc_annotation_bed_path="",fa_path=""){
  

  #best_order<-read.table(best_order_path,sep = "\t",header = TRUE,as.is = TRUE,fill=TRUE);
  
  best_order<-as.data.frame(read_delim(best_order_path,delim="$",
                                    col_names = TRUE,comment="#",col_types="ccdciddddd") )
  
  best_order<-best_order[is.na( best_order$number_of_orders_have_same_prob ) |
                               (best_order$number_of_orders_have_same_prob ==1),];
  
  
  if(rebuild_anno){
    ucsc_annotation_bed<-read.table(ucsc_annotation_bed_path,sep = "\t",header = FALSE,as.is = TRUE)
    colnames(ucsc_annotation_bed)<-c("chr","start","end","trans_id","score","strand","CDS_start","CDS_end",
                                     "","exon_count","exon_len","exon_start")
    
    #ucsc_annotation_bed<-as.matrix(ucsc_annotation_bed)
    # intron_len<-matrix(nrow = nrow(ucsc_annotation_bed))
    # intron_5tss_len<-matrix(nrow = nrow(ucsc_annotation_bed))
    # intron_3tss_len<-matrix(nrow = nrow(ucsc_annotation_bed))
    all_intron_len<-c();
    all_intron_5tss_len<-c();
    all_intron_3tss_len<-c();

    all_trans_id<-c()
    
    all_intron_5ss_score<-c()
    all_intron_3ss_score<-c()
    
    #unlink("/var/tmp/ss5_seq.txt");
    #unlink("/var/tmp/ss3_seq.txt");
    
    
    exon_list<-str_split(ucsc_annotation_bed[,"exon_len"], ",")
    exon_exon_list<-str_split(ucsc_annotation_bed[,"exon_start"], ",")
    
    
    for (row_num in 1:nrow(ucsc_annotation_bed)) {
      
      exon<-as.numeric( exon_list[[row_num]] );
      exon_pos<-as.numeric( exon_exon_list[[row_num]] );
      #exon<-as.numeric(str_split(ucsc_annotation_bed[row_num,"exon_len"], ",")[[1]]);
      #exon_pos<-as.numeric(str_split(ucsc_annotation_bed[row_num,"exon_start"],  ",")[[1]]);
      
      ##remove , after last exons
      exon<-exon[-1*length(exon)];
      exon_pos<-exon_pos[-1*length(exon_pos)];
      
      transcript_id<-ucsc_annotation_bed[row_num,"trans_id"];
      print(paste0(row_num,":",transcript_id) );
      
      transcript_start_site<-ucsc_annotation_bed[row_num,"start"];
      transcript_end_site<-ucsc_annotation_bed[row_num,"end"];
      strand<-ucsc_annotation_bed[row_num,"strand"];
      chr<-ucsc_annotation_bed[row_num,"chr"];
      exon_count<-ucsc_annotation_bed[row_num,"exon_count"]
      
      intron_5ss_seq<-c();
      intron_3ss_seq<-c();
      
      if(str_detect(chr,"GL") || str_detect(chr,"gl") ){
        next;
      }
      
      if(str_detect(chr,"random")){
        next;
      }
      #if(length(exon)==length(exon_pos))
      if(exon_count<=2)
      {
        next
      }
      
      all_trans_id<-c(all_trans_id,transcript_id)
      #intron_5ss_score_v<-c();
      #intron_3ss_score_v<-c();
        
      intron_length<-c()
      intron_5tss_distance<-c()
      intron_3tss_distance<-c()
      intron_5ss_seq<-c()
      intron_3ss_seq<-c()
      
      ##intron length
      for (i in 1:(length(exon)-1)) {
        # calculate intron length
        intron_length[i]<-(exon_pos[i+1])-(exon_pos[i])-(exon[i]);
        #intron_5ss_score<-0; intron_3ss_score<-0;
  
      }
        
      ##introns' 5ss and 3ss score;
      for (i in 1:(length(exon)-1)) {
          
        intron_3ss_region<-c();intron_5ss_region<-c();
          
        if(strand=="+"){
          pos_5ss<-transcript_start_site+(exon_pos[i])+(exon[i])+1;
          pos_3ss<-(exon_pos[i+1])+transcript_start_site+1-1;
          intron_5ss_region<-str_c(chr,":",pos_5ss-3,":",pos_5ss+5,":",strand)
          intron_3ss_region<-str_c(chr,":",pos_3ss-19,":",pos_3ss+3,":",strand)
        }
          
        if(strand=="-"){
          pos_3ss<-transcript_start_site+(exon_pos[i])+(exon[i])+1;
          pos_5ss<-(exon_pos[i+1])+transcript_start_site+1-1;
          intron_5ss_region<-str_c(chr,":",pos_5ss-5,":",pos_5ss+3,":",strand)
          intron_3ss_region<-str_c(chr,":",pos_3ss-3,":",pos_3ss+19,":",strand)
        }
        
        one5_seq<-get_seq_region(intron_5ss_region,fa_path)
        one3_seq<-get_seq_region(intron_3ss_region,fa_path)
        
        while(str_detect(one5_seq,"N") ){one5_seq<-str_replace(one5_seq,"N",sample(c("A","C","G","T"), 1) )}
        while(str_detect(one3_seq,"N") ){one3_seq<-str_replace(one3_seq,"N",sample(c("A","C","G","T"), 1) )}
        #print(one5_seq)
        #print(one3_seq)
        
        intron_5ss_seq<-c(intron_5ss_seq, one5_seq )
          
        intron_3ss_seq<-c(intron_3ss_seq, one3_seq )
          
        #intron_5ss_score<-get_score_5ss();
        #intron_3ss_score<-get_score_3ss();
        #intron_5ss_score_v[i]<-intron_5ss_score;
        #intron_3ss_score_v[i]<-intron_3ss_score;
      }
        
      ##introns' 5ss and 3ss score;
      if(ucsc_annotation_bed[row_num,"strand"]=="+"){
        for (i in 1:(length(exon)-1)) {
          # calculate intron distance to 5' tss for strain "+"
          intron_5tss_distance[i]<-as.numeric(exon_pos[i])+as.numeric(exon[i])
        }
          
        for (i in 1:(length(exon)-1)) {
          # calculate intron distance to 3' tss for strain "+"
          intron_3tss_distance[i]<-as.numeric(exon_pos[i+1])
        }
          
      }
        
      if(ucsc_annotation_bed[row_num,"strand"]=="-"){
        for (i in 1:(length(exon)-1)) {
          # calculate intron distance to 5' tss for strain "-"
          intron_5tss_distance[i]<-(transcript_end_site)-(transcript_start_site)-(exon_pos[i+1]);
        }
  
        for (i in 1:(length(exon)-1)) {
          # calculate intron distance to 3' tss for strain "-"
          intron_3tss_distance[i]<-(transcript_end_site)-(transcript_start_site)-(exon_pos[i])-(exon[i])
        }
          
        intron_length<-rev(intron_length);
        intron_5tss_distance<-rev(intron_5tss_distance)
        intron_3tss_distance<-rev(intron_3tss_distance)
      }
      
          
      #print(all_lengths)
      all_intron_len<-c(all_intron_len,paste0(intron_length,collapse = ",") )
      all_intron_5tss_len<-c(all_intron_5tss_len,paste0(intron_5tss_distance,collapse = ","))
      all_intron_3tss_len<-c(all_intron_3tss_len,paste0(intron_3tss_distance,collapse = ","))
      #all_intron_5ss_score[row_num]<-intron_5ss_score_v_c
      #all_intron_3ss_score[row_num]<-intron_3ss_score_v_c
      
      #print( paste0(get_score_5ss_file(intron_5ss_seq),collapse = ",") );
      #print( paste0(get_score_3ss_file(intron_3ss_seq),collapse = ",") );
      
      all_intron_5ss_score<-c(all_intron_5ss_score,paste0(get_score_5ss_file(intron_5ss_seq),collapse = ",") );
      all_intron_3ss_score<-c(all_intron_3ss_score,paste0(get_score_3ss_file(intron_3ss_seq),collapse = ",") );
    }
    
    #print("Got intron length, distance");
    #all_intron_5ss_score<-get_score_5ss_file(intron_5ss_seq);
    
    #all_intron_3ss_score<-get_score_3ss_file(intron_3ss_seq);
    
    
    print("Got i5ss and 3ss score");
    ucsc_annotation_intron_length<-data.frame(trans_id=all_trans_id,
                                              intron_len=all_intron_len,
                                              intron_5tss_len=all_intron_5tss_len,
                                              intron_3tss_len=all_intron_3tss_len,
                                              intron_5ss_score_v=all_intron_5ss_score,
                                              intron_3ss_score_v=all_intron_3ss_score,
                                              stringsAsFactors=FALSE);
    
    save(ucsc_annotation_intron_length,file=anno_save_load);
    
  }
  
  
  load(anno_save_load);
  
  print("Begin order calculation");

  order_intron_len<-function(x,y){
    z<-""
    if(length(x)==length(y)){
      for (i in 1:length(y)) {
        z[i]<-y[x[i]]
      }
    }
    return(z)
  }
  
  #chr<-ucsc_annotation_bed[row_num,1] 
  
  # head(intron_len,10)
  # head(intron_5tss_len,5)
  # head(intron_3tss_len,5);
  
  # head(ucsc_annotation_intron_length)
  # dim(ucsc_annotation_intron_length)
  #write.table(ucsc_annotation_intron_length,"ucsc_annotation_intron_length_5tss_3tss_distance.bed",
  #sep="\t",quote=F,row.names = F)
  #pombe_ucsc_annotation_intron_length_5tss_distance.bed
  #ucsc_annotation_intron_length<-read.table("ucsc_annotation_intron_length_5tss_3tss_distance.bed",header = T,sep="\t");
  
  #for (row_num in 1:nrow(ucsc_annotation_intron_length)) {
  # intron_length<-as.vector(ucsc_annotation_intron_length[row_num,11])
  #transcript_id<-as.vector(ucsc_annotation_intron_length[row_num,4])
  #transcript<-gsub(pattern = ".\\d_\\d+",replacement = "",x = transcript_id)
  #predict_order<-as.vector(best_order[grep(pattern=transcript, x=rownames(best_order), value=TRUE),2])
  #order<-as.numeric(as.character(unlist(strsplit(predict_order, split = ","))))
  #ordered_len<-order_intron_len(order,len)
  #}
  
  correlation_intron_length<-c()
  correlation_intron_5tss_length<-c()
  correlation_intron_3tss_length<-c()
  correlation_intron_5tss_len_length<-c()
  correlation_intron_5tss_score<-c()
  correlation_intron_3tss_score<-c();
  
  for (row_num in 1:nrow(best_order) ) {
    predict_order<-as.vector(best_order[row_num,"best_order"])
    order<-as.numeric(as.character(unlist(strsplit(predict_order, split = ","))))
    
    #transcript<-strsplit((best_order[row_num,"transcript_id"]),split="_")[[1]][1];
    
    transcript<-best_order[row_num,"transcript_id"];
    print(paste0(row_num,":",transcript) );
    
    
    target_row_num_in_genom<-which(str_detect(ucsc_annotation_intron_length[,"trans_id"],transcript) );
    
    intron_length<-ucsc_annotation_intron_length[target_row_num_in_genom,"intron_len"];
    
    intron_5tss_distance<-ucsc_annotation_intron_length[target_row_num_in_genom,"intron_5tss_len"];
    
    intron_3tss_distance<-ucsc_annotation_intron_length[target_row_num_in_genom,"intron_3tss_len"];
      
    intron_5ss_score_v<-ucsc_annotation_intron_length[target_row_num_in_genom,"intron_5ss_score_v"];
    
    intron_3ss_score_v<-ucsc_annotation_intron_length[target_row_num_in_genom,"intron_3ss_score_v"];
    
    
    len<-as.numeric(str_split(intron_length,",")[[1]]);
    distance_5tss<-as.numeric(str_split(intron_5tss_distance,  ",")[[1]])
    distance_3tss<-as.numeric(str_split(intron_3tss_distance,  ",")[[1]])
    intron_5ss_score<-as.numeric(str_split(intron_5ss_score_v,  ",")[[1]])
    intron_3ss_score<-as.numeric(str_split(intron_3ss_score_v, ",")[[1]]);
    
    if(length(order)!=length(len)){
      print(paste0("Number of order is not equal to number of intron length: ", transcript));
      correlation_intron_length[row_num]<-NA;correlation_intron_5tss_length[row_num]<-NA;
      correlation_intron_3tss_length[row_num]<-NA;correlation_intron_5tss_score[row_num]<-NA;
      correlation_intron_3tss_score[row_num]<-NA;
      next;
    }
    
    
    ordered_len<-order_intron_len(order,len)
    ordered_5tss_distance<-order_intron_len(order,distance_5tss)
    ordered_3tss_distance<-order_intron_len(order,distance_3tss)
    ordered_5tss_len_distance<-order_intron_len(order,distance_5tss+len);
    
    ordered_5tss_score<-order_intron_len(order,intron_5ss_score)
    ordered_3tss_score<-order_intron_len(order,intron_3ss_score)
    
    
    correlation_intron_length[row_num]<-cor(1:length(len),as.numeric(ordered_len),method = "spearman")
    
    correlation_intron_5tss_length[row_num]<-cor(1:length(ordered_5tss_distance),
                                                 as.numeric(ordered_5tss_distance),method = "spearman")
    
    correlation_intron_3tss_length[row_num]<-cor(1:length(ordered_3tss_distance),
                                                 as.numeric(ordered_3tss_distance),method = "spearman")
    
    correlation_intron_5tss_score[row_num]<-cor(1:length(ordered_5tss_score),
                                                 as.numeric(ordered_5tss_score),method = "spearman")
    
    correlation_intron_3tss_score[row_num]<-cor(1:length(ordered_3tss_score),
                                                 as.numeric(ordered_3tss_score),method = "spearman")
    
    
    #correlation_intron_5tss_len_length[row_num]<-cor(1:length(ordered_5tss_len_distance),
    #                                                 as.numeric(ordered_5tss_len_distance),method = "spearman")
    
    #row_num<-row_num+1;
  }
  
  
  
  best_order_cor<-cbind(as.matrix(best_order),correlation_intron_length,
                        correlation_intron_5tss_length,correlation_intron_3tss_length,
                        correlation_intron_5tss_score,correlation_intron_3tss_score)
  
  best_order_cor<-data.frame(gene_id=best_order$gene_symbol,
                             trans_id=best_order$transcript_id,
                             best_order=best_order$best_order,
                             #p_value_log=best_order$p_value_log,
                             correlation_intron_length=correlation_intron_length,
                             correlation_intron_5tss=correlation_intron_5tss_length,
                             correlation_intron_3tss=correlation_intron_3tss_length,
                             correlation_intron_5tss_score=correlation_intron_5tss_score,
                             correlation_intron_3tss_score=correlation_intron_3tss_score
                             );
  



  p_all_mat_label<-c(rep("Intron length",nrow(best_order_cor)),
                         rep("Intron distance to TSS",nrow(best_order_cor)),
                         rep("Intron 5'ss socre",nrow(best_order_cor)),
                         rep("Intron 3'ss socre",nrow(best_order_cor)))
  
  p_all_mat_v<-c(best_order_cor[,"correlation_intron_length"],
                best_order_cor[,"correlation_intron_5tss"],
                best_order_cor[,"correlation_intron_5tss_score"],
                best_order_cor[,"correlation_intron_3tss_score"])
  
  p_all_mat<-data.frame(label=p_all_mat_label,
                        value=p_all_mat_v,
                        species=rep(output_label,length(p_all_mat_v) ),
                        stringsAsFactors = FALSE)
  


}

#p_all<-ggplot(p_all_mat)+geom_violin(aes(x=p_all_mat_label,y=p_all_mat_v))+theme_minimal()+xlab("")+ylab("")

#print(p_all);

#p_all<-ggplot(p_all_mat)+geom_boxplot(aes(x=p_all_mat_label,y=p_all_mat_v))+theme_minimal()+xlab("")+ylab("")

#print(p_all);
#dev.off();

#source("code/multiplot.R");
#colnames(best_order_cor)<-c("gene_id","trans_id","p_value_log","best_order","bayesian_factor","relative_likelihood",
#                            "spearman_rho","spearman_rho_abs","spearman_p_value","p_value",
#                            "correlation_intron_length","correlation_intron_5tss","correlation_intron_3tss",
#                            "correlation_intron_5tss_score","correlation_intron_3tss_score");

#best_order_cor_ordered<-best_order_cor[order(best_order_cor[,"p_value_log"]),];


##if only two elements, then cor = -1 or +1
#intron number >2
#best_order_cor<-best_order_cor[str_count(best_order_cor[,"best_order"],",")>1,]

#pdf(output_path, width = 7,height = 6);

#sort(correlation_intron_5tss_length[correlation_intron_5tss_length>0])
# p1<-ggplot()+geom_histogram(aes(as.numeric(best_order_cor[,"correlation_intron_length"])) )+
#   xlab("Correlation of intron length with splicing order");
# 
# p2<-ggplot()+geom_histogram(aes(as.numeric(best_order_cor[,"correlation_intron_5tss"])  ) )+
#   xlab("Correlation of introns' 5 prime positions to TSS with splicing order");
#  p3<-ggplot()+geom_histogram(aes(as.numeric(best_order_cor[,"correlation_intron_3tss"])  ) )+
#    xlab("Correlation of introns' 3 prime positions to TSS with splicing order");

# p4<-ggplot()+geom_histogram(aes(as.numeric(best_order_cor[,"correlation_intron_5tss_score"])  ) )+
#   xlab("Correlation of introns' 5 prime MaxEntScore with splicing order");
# p5<-ggplot()+geom_histogram(aes(as.numeric(best_order_cor[,"correlation_intron_3tss_score"])  ) )+
#   xlab("Correlation of introns' 3 prime MaxEntScore with splicing order");
# 

#multiplot(p1,p2,p4,p5,cols=1);

