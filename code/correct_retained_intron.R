library(stringr)

#  a,"anno/rmats_human/RI.MATS.JC.txt",b
# t_intron_pos_mat_fr<-b


# t_igraph_list<-a
# t_reatined_intron_psi_file<-"arabidopsis/rmats_ara/RI.MATS.JC.txt"

# t_exon_pos_mat_fr<-e_f
# t_intron_pos_mat_fr<-i_f
# iso_final<-f
# iso_summary<-g


build_retained_intron_transcript<-function(t_intron_pos_mat_fr, t_exon_pos_mat_fr, t_reatined_intron_psi_file){
  

  #this list store the introns that in the inclusion form this retained introns
  other_introns_in_retained_intron<-list();
  
  
  t_reatined_intron_psi<-read.table(t_reatined_intron_psi_file, sep="\t", header = TRUE,quote = "$");
  
  
  cat("\n");
  
  print("Find introns in the inclusion isoform of the retained introns:");
  pb = txtProgressBar(min = 1, max = nrow(t_reatined_intron_psi), initial = 0, width=100, style=3);
  
  for( i in 1:nrow(t_reatined_intron_psi) ){
    
    setTxtProgressBar(pb,i);
    
    whole_inclusion_exon<-str_c(t_reatined_intron_psi[i,"chr"],":",
                          (t_reatined_intron_psi[i,"riExonStart_0base"]+1),"-",
                          t_reatined_intron_psi[i,"riExonEnd"]);
    
    inclusion_region<-str_c(t_reatined_intron_psi[i,"chr"],":",
                          t_reatined_intron_psi[i,"upstreamEE"]+1,"-",
                          t_reatined_intron_psi[i,"downstreamES"]);
    
    
    trans_ids_contain_retained_introns<-t_exon_pos_mat_fr[t_exon_pos_mat_fr$exon_pos==whole_inclusion_exon,"trans_id"];
    
    ##introns that are in the inclusion form of transcripts
    introns_in_transcripts<-t_intron_pos_mat_fr[t_intron_pos_mat_fr$trans_id %in% trans_ids_contain_retained_introns,"intron_pos" ];
    
    
    other_introns_in_retained_intron[[inclusion_region]]<-unique(introns_in_transcripts);
    
  }
  
  other_introns_in_retained_intron
  
}


correct_retained_intron<-function(t_igraph_list, t_reatined_intron_psi_file, t_intron_pos_mat_fr,t_exon_pos_mat_fr){
  
  other_introns_in_retained_intron<-build_retained_intron_transcript(t_intron_pos_mat_fr,
                                                                     t_exon_pos_mat_fr,t_reatined_intron_psi_file)
  
  ##re-assign the reads of retained intron based on PSI values
  if(t_reatined_intron_psi_file==""){
    return(t_igraph_list)
  }
  
  t_reatined_intron_psi<-read.table(t_reatined_intron_psi_file,sep="\t",header = TRUE,quote = "$");
  
  rep1_psi<-(sapply(str_split( t_reatined_intron_psi[,"IncLevel1"],","),"[",1) );
  rep2_psi<-(sapply(str_split( t_reatined_intron_psi[,"IncLevel1"],","),"[",2) );
  
  if(length(rep2_psi)==sum(is.na(rep2_psi))){
    rep2_psi<-rep1_psi;
  }
  
  rep1_psi[(rep1_psi=="NA")]<-0;   rep2_psi[(rep2_psi=="NA")]<-0;
  rep1_psi<-as.numeric(rep1_psi);  rep2_psi<-as.numeric(rep2_psi)

  
  t_reatined_intron_psi_vec<-(rep1_psi+rep2_psi)/2;
  #t_reatined_intron_psi_vec<-(rep1_psi)/1;
  
  
  names(t_reatined_intron_psi_vec)<-str_c(t_reatined_intron_psi$chr,":",t_reatined_intron_psi$upstreamEE+1,"-",
                                          t_reatined_intron_psi$downstreamES);
  
  t_reatined_intron_psi_vec<-t_reatined_intron_psi_vec[t_reatined_intron_psi_vec>0];
  
  cat("",file="results/test_correct_RI.txt",append=FALSE)
  
  
  cat("\n");
  print("Corrected retained introns:");
  pb = txtProgressBar(min = 1, max = length(t_igraph_list), initial = 0, width=100, style=3) 
  
  for(i in 1:length(t_igraph_list)){
    
    setTxtProgressBar(pb,i);
    
    
    trans_id<-t_igraph_list[[i]]$trans_id;
    
    one_igraph_list<-t_igraph_list[[i]];
    
    m_adj_matrix<-one_igraph_list$adjacency_matrix;
    
    m_jc_pair_matrix<-one_igraph_list$jc_pair_matrix;
    
    m_index_pos_fr<-one_igraph_list$index_pos_fr;
    
    one_igraph_list[["raw_m_adj_matrix"]]<-m_adj_matrix;
    
    ###one transcript may contain two or more retained intron events
    m_adj_matrix_old<-m_adj_matrix;
    
    for(j in 1:nrow(m_index_pos_fr)){
      
        if(m_index_pos_fr[j,"pos"] %in% names(t_reatined_intron_psi_vec) ){
        
          #print("here")
          psi_correct_value<-(1 - t_reatined_intron_psi_vec[m_index_pos_fr[j,"pos"] ]);
          
          
          #get the introns that in the inclusion form of the isoforms and test if in the transcripts
          need_to_correct_index<-m_index_pos_fr[,"pos"] %in% other_introns_in_retained_intron[[m_index_pos_fr[j,"pos"]]];
          
          if(is.na(need_to_correct_index)|| (length(need_to_correct_index)==0)){
            next;
          }
          
          m_jc_pair_count<-0;
          
          intron_index<-m_index_pos_fr[j,"index"];
          
          if(TRUE){
            cat(paste(trans_id,":",need_to_correct_index,":",psi_correct_value,"\n"),file="results/test_correct_RI.txt",append=TRUE)
            
            cat(paste(format(m_jc_pair_matrix[intron_index,need_to_correct_index],digits=5),collapse = ", "),
                file="results/test_correct_RI.txt",append=TRUE)
            cat("\n",file="results/test_correct_RI.txt",append=TRUE)
            
            cat(paste(format(m_adj_matrix_old[intron_index,need_to_correct_index],digits=5),collapse = ", "),
                file="results/test_correct_RI.txt",append=TRUE)
            cat("\n",file="results/test_correct_RI.txt",append=TRUE)
            
            
            cat(paste(format(m_adj_matrix_old[need_to_correct_index, intron_index],digits=5),collapse = ", "),
                file="results/test_correct_RI.txt",append=TRUE)
            cat("\n",file="results/test_correct_RI.txt",append=TRUE)
          }
          # row spliced before column
          corrected_value <- 
            m_adj_matrix_old[need_to_correct_index, intron_index] -
            psi_correct_value*(m_adj_matrix_old[need_to_correct_index, intron_index]+
                                 m_adj_matrix_old[intron_index,need_to_correct_index]+
                                 m_jc_pair_matrix[intron_index,need_to_correct_index]);
          
          corrected_value[corrected_value<0]<-0;
          
          m_adj_matrix[need_to_correct_index, intron_index] <-corrected_value;
          #m_adj_matrix[corrected_value<0, intron_index]<-0;

          if(TRUE){
            cat(paste(format(m_adj_matrix[need_to_correct_index, intron_index],digits=5), collapse = ", "),
                file="results/test_correct_RI.txt",append=TRUE)
            
            cat("\n",file="results/test_correct_RI.txt",append=TRUE)
            cat("\n",file="results/test_correct_RI.txt",append=TRUE)
            cat("\n",file="results/test_correct_RI.txt",append=TRUE)
          }
          #m_adj_matrix[, intron_index]<0 m_adj_matrix[, intron_index]<-0;
          
        }
      
    }
    
    one_igraph_list$adjacency_matrix<-m_adj_matrix;
    
    t_igraph_list[[i]]<-one_igraph_list;
    
  }
  
  t_igraph_list
}

#correct_retained_intron(a,"anno/rmats_human/RI.MATS.JC.txt",b);

