source("code/analysis/cor_cal.R",echo=TRUE);

pombe_best_order_path<-"./results/best_order_pombe.tsv";
pombe_anno_save_path<-"data/pombe/pombe_intron_len_score.Rd";
#pombe_ucsc_annotation_bed_path<-"data/pombe/Schizosaccharomyces_pombe.ASM294v2.43.chr_nothick.bed";
#pombe_fasta_path<-FaFile("data/pombe/Schizosaccharomyces_pombe.ASM294v2.fa");

fly_best_order_path<-"./results/best_order_fly.tsv";
fly_anno_save_path<-"data/fly/fly_intron_len_score.Rd";
#fly_ucsc_annotation_bed_path<-"fly/dm6_ensembl_no_thick.bed";
#fly_fasta_path<-FaFile("/Volumes/mengli/anno/anno_dm6/Drosophila_melanogaster.BDGP6.22.dna.toplevel.fa");


zebrafish_best_order_path<-"./results/best_order_zebrafish.tsv";
zebrafish_anno_save_path<-"./data/zebrafish/zebrafish_intron_len_score.Rd";
#zebrafish_ucsc_annotation_bed_path<-"zebrafish/GRCz11_ensembl_nothick.bed";
#zebrafish_fasta_path<-FaFile("zebrafish/fasta/Danio_rerio.GRCz11.dna.primary_assembly.fa");


human_best_order_path<-"./results/best_order_human.tsv";
human_anno_save_path<-"./data/human_intron_len_score.Rd";
#human_ucsc_annotation_bed_path<-"anno/hg19_gencode_from_ucsc_nothick_nocds.bed";
#human_fasta_path<-FaFile("/Volumes/mengli/anno/human_gencode/GRCh37.primary_assembly.genome.fa");


p_all_mat_pombe<-cal_intron_metric_cor(
pombe_anno_save_path,pombe_best_order_path,"pombe");

p_all_mat_fly<-cal_intron_metric_cor(
fly_anno_save_path,fly_best_order_path,"fly");

p_all_mat_zebrafish<-cal_intron_metric_cor(
zebrafish_anno_save_path,zebrafish_best_order_path,"zebrafish");

p_all_mat_human<-cal_intron_metric_cor(
human_anno_save_path,human_best_order_path,"human");



p_all_mat_all<-rbind(p_all_mat_pombe,p_all_mat_fly,p_all_mat_zebrafish,p_all_mat_human);

number_trans_human<-nrow(p_all_mat_human)/4
number_trans_zebrafish<-nrow(p_all_mat_zebrafish)/4
number_trans_fly<-nrow(p_all_mat_fly)/4
number_trans_pombe<-nrow(p_all_mat_pombe)/4


p_all_mat_all$species<-str_replace_all(p_all_mat_all$species,"human",paste0("Humans (n=",number_trans_human,")" )  )
p_all_mat_all$species<-str_replace_all(p_all_mat_all$species,"zebrafish",paste0("Zebrafish (n=",number_trans_zebrafish,")" )  )
p_all_mat_all$species<-str_replace_all(p_all_mat_all$species,"fly",paste0("Fruit fly (n=",number_trans_fly,")" )  )
p_all_mat_all$species<-str_replace_all(p_all_mat_all$species,"pombe",paste0("S. pombe (n=",number_trans_pombe,")" )  )

p_all_mat_all$species<-factor(p_all_mat_all$species,
                            levels =c(paste0("S. pombe (n=",number_trans_pombe,")" ),
                                      paste0("Fruit fly (n=",number_trans_fly,")" ) ,
                                      paste0("Zebrafish (n=",number_trans_zebrafish,")" ) ,
                                      paste0("Humans (n=",number_trans_human,")" )) )

p_all_mat_all$label<-factor(p_all_mat_all$label,levels=c("Intron length",
                                                         "Intron distance to TSS",
                                                         "Intron 5'ss socre",
                                                         "Intron 3'ss socre"));


p_all<-ggplot(p_all_mat_all)+geom_violin(aes(x="",y=value,fill=paste0(label,species)) )+theme_minimal()+xlab("")+ylab("")+
facet_grid(cols=vars(label),rows=vars(species), scales = "free", space = "free")+
  theme(legend.position="none")+ylim(-1.1,1.1)+
  scale_fill_manual(values=c( rep(hue_pal()(4)[1],4 ),
                              rep(hue_pal()(4)[2], 4 ),
                              rep(hue_pal()(4)[3], 4 ),rep(hue_pal()(4)[4], 4 )) );



p_all_s<-ggplot(p_all_mat_all,aes(x="",y=value,fill=paste0(label,species)) )+
  geom_violin( )+theme_minimal()+xlab("")+ylab("")+ geom_point(stat="summary", fun.y="mean") + 
  facet_grid(rows=vars(label),cols=vars(species), scales = "free", space = "free")+
  theme(legend.position="none")+ylim(-1.1,1.1)+
  scale_fill_manual(values=c( rep(hue_pal()(4)[1:4],1 ),
                              rep(hue_pal()(4)[1:4], 1 ),
                              rep(hue_pal()(4)[1:4], 1 ),rep(hue_pal()(4)[1:4], 1 )) );


pdf("results/cor_order_intron_length_5ss_score.pdf",width=8,height=12)

print(p_all_s);
#print(p_all);

dev.off();



