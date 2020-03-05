

# t_iso_final<-iso_final
# t_anno_intron<-intron_pos_mat_fr

get_iso_summary<-function(t_iso_final,t_anno_intron){
  
  
  iso_slow_sumary<-t_anno_intron %>% dplyr::group_by(trans_id) %>%
    dplyr::summarise(intron_count=max(intron_order) );
  
  iso_slow_sumary<-as.data.frame(iso_slow_sumary);
  
  
  ####get intron coverage pair
  small_intron<-apply(t_iso_final[c("nexti","first")],1,function(x){ min(x)});
  
  large_intron<-apply(t_iso_final[c("nexti","first")],1,function(x){ max(x)});
  
  t_iso_final[,"intron_pair"]<-str_c(small_intron,large_intron);
  
  
  iso_edge_count<-as.data.frame(t_iso_final %>% dplyr::group_by(id) %>%
                                  dplyr::summarise(intron_pair_count=n_distinct(intron_pair) ) );
  
  iso_slow_sumary<-inner_join(iso_slow_sumary,iso_edge_count,by=c("trans_id"="id") );
  
  #################################################output and saving##################################################
  
  iso_slow_sumary[,"percent_intron_pair_coverage"]<-(iso_slow_sumary$intron_pair_count)/
    (iso_slow_sumary[,"intron_count"]*(iso_slow_sumary[,"intron_count"]-1)/2);
  
  
  iso_slow_sumary<-iso_slow_sumary[order( iso_slow_sumary[,"percent_intron_pair_coverage"],decreasing = TRUE),];
  
  print(paste0("Number of multi introns transcripts  detected =", nrow(iso_slow_sumary) ) );
  
  
  iso_slow_sumary
  
}


