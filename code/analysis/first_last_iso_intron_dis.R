###the fre of last intron and first intron
library(ggplot2)
library(stringr)
library(scales)
library(readr)

setwd("/Users/mengli/Documents/projects/iso");
source("code/multiplot.R");


get_first_intron_spliced_fre<-function(t_best_order, output_label){
  
  t_best_order<-t_best_order[is.na( t_best_order$number_of_orders_have_same_prob ) |
                             (t_best_order$number_of_orders_have_same_prob ==1),];
  
  
  trans_order_list<-str_split(t_best_order[,"best_order"],",");
  
  pos_relative_first_intron<-c();
  pos_relative_last_intron<-c();
  
  for(i in 1:length(trans_order_list) ){
    
    one_trans<-trans_order_list[[i]];
    
    pos_relative<-(which(one_trans=="1")-1)/(length(one_trans)-1);
    
    pos_relative_first_intron<-c(pos_relative_first_intron,pos_relative);
    
    
    pos_relative<-(which(one_trans==as.character(length(one_trans) ) ) -1)/(length(one_trans)-1);
    
    pos_relative_last_intron<-c(pos_relative_last_intron,pos_relative);
    
  }
  
  
  p3<-ggplot()+geom_histogram(aes(x=pos_relative_first_intron) )+
    xlab("First intron order/intron number")+xlim(-0.1,1.1)+theme_minimal()+ylab("")
  
  
  p4<-ggplot()+geom_histogram(aes(x=pos_relative_last_intron) )+
    xlab("Last intron order/intron number")+xlim(-0.1,1.1)+theme_minimal()+ylab("")
  
  
  p5<-ggplot()+geom_violin(aes(x="",y=pos_relative_first_intron,fill="" ) )+
    xlab("First intron order/intron number")+ylim(-0.1,1.1)+theme_minimal()+ylab("")
  
  
  p6<-ggplot()+geom_violin(aes(x="",y=pos_relative_last_intron,fill="" ) )+
    xlab("Last intron order/intron number")+ylim(-0.1,1.1)+theme_minimal()+ylab("")
  
  
  #multiplot(p3,p4,cols=2)
  
  #multiplot(p5,p6,cols = 2);
  
  pos_relative_first_last_fr<-data.frame(rela_orders=c(pos_relative_first_intron,pos_relative_last_intron),
                                         label=c(rep(paste0("First intron order/intron number"),
                                                     length(pos_relative_first_intron)),
                                                 rep(paste0("Last intron order/intron number"),
                                                     length(pos_relative_last_intron)) ),
                                         species=rep(output_label, length(pos_relative_first_intron)+
                                                       length(pos_relative_last_intron))  )

}

best_order_human<-read.table("result/best_order/best_order.tsv",sep = "$",header = TRUE,as.is = TRUE,fill = TRUE );
#output_human<-"result/intron_splicing_order_first_last_dis.pdf"
pos_relative_first_last_fr_human<-get_first_intron_spliced_fre(best_order_human,"Human")


best_order_zebrafish<-read.table("zebrafish/result/best_order.tsv",sep = "$",header = TRUE,as.is = TRUE,fill = TRUE );
#output_human<-"result/intron_splicing_order_first_last_dis.pdf"
pos_relative_first_last_fr_zebrafish<-get_first_intron_spliced_fre(best_order_zebrafish,"Zebrafish")


#best_order_fly<-read.table("fly/result/best_order.tsv",sep = "\t",header = TRUE,as.is = TRUE,fill = TRUE );
best_order_fly<-as.data.frame(read_delim("fly/result/best_order.tsv",delim="$",
                                  col_names = TRUE,comment="#",col_types="ccdciddddd") )

#output_fly<-"fly/result/intron_splicing_order_first_last_dis.pdf"
pos_relative_first_last_fr_fly<-get_first_intron_spliced_fre(best_order_fly,"Fly")


best_order_pombe<-read.table("pombe/result/best_order.tsv",sep = "$",header = TRUE,as.is = TRUE,fill = TRUE );
#output_pombe<-"pombe/result/intron_splicing_order_first_last_dis.pdf"
pos_relative_first_last_fr_pombe<-get_first_intron_spliced_fre(best_order_pombe,"Pombe")



pos_relative_first_last_fr_all<-rbind(pos_relative_first_last_fr_human,
                                      pos_relative_first_last_fr_zebrafish,
                                      pos_relative_first_last_fr_fly,
                                      pos_relative_first_last_fr_pombe);



pos_relative_first_last_fr_all$species<-str_replace_all(pos_relative_first_last_fr_all$species,
                                                        "Fly",paste0("Fly (n=",nrow(pos_relative_first_last_fr_fly)/2,")" )  )
pos_relative_first_last_fr_all$species<-str_replace_all(pos_relative_first_last_fr_all$species,
                                                        "Zebrafish",
                                                        paste0("Zebrafish (n=",nrow(pos_relative_first_last_fr_zebrafish)/2,")" )  )

pos_relative_first_last_fr_all$species<-str_replace_all(pos_relative_first_last_fr_all$species,
                                                        "Human",paste0("Human (n=",nrow(pos_relative_first_last_fr_human)/2,")" )  )
pos_relative_first_last_fr_all$species<-str_replace_all(pos_relative_first_last_fr_all$species,
                                                        "Pombe",paste0("Pombe (n=",nrow(pos_relative_first_last_fr_pombe)/2,")" )  )


pos_relative_first_last_fr_all$species<-factor(pos_relative_first_last_fr_all$species,
                                               c(paste0("Pombe (n=",nrow(pos_relative_first_last_fr_pombe)/2,")" ),
                                                 paste0("Fly (n=",nrow(pos_relative_first_last_fr_fly)/2,")" ) ,
                                                 paste0("Zebrafish (n=",nrow(pos_relative_first_last_fr_zebrafish)/2,")" ),
                                                 paste0("Human (n=",nrow(pos_relative_first_last_fr_human)/2,")" )))

pos_relative_first_last_fr_all$label<-factor(pos_relative_first_last_fr_all$label,
                                             c("First intron order/intron number","Last intron order/intron number"))


p_all<-ggplot(pos_relative_first_last_fr_all)+geom_violin(aes(x="",y=rela_orders,fill=paste0(label,species ) ) )+
  xlab("")+ylim(-0.1,1.1)+theme_minimal()+ylab("")+
  facet_grid(cols=vars(label),rows=vars(species), scales = "free", space = "free")+theme(legend.position="none")+
  scale_fill_manual(values=c( rep(hue_pal()(1), 4 ), rep(hue_pal()(2)[2], 4 ) ) );
  

data_summary <- function(x) {
  mu <- mean(x)
  sigma1 <- mu-sd(x)
  sigma2 <- mu+sd(x)
  return(c(y=mu,ymin=sigma1,ymax=sigma2))
}

p_all_s<-ggplot(pos_relative_first_last_fr_all,aes(x="",y=rela_orders,fill=paste0(label,species ) ))+
  geom_violin( )+ geom_point(stat="summary", fun.y="mean") + 
  xlab("")+ylim(-0.1,1.1)+theme_minimal()+ylab("")+
  facet_grid(rows=vars(label),cols=vars(species), scales = "free", space = "free")+theme(legend.position="none")+
  scale_fill_manual(values=c( rep(hue_pal()(4)[1:4], 4 ), rep(hue_pal()(4)[1:4], 4 ) ) );
#geom_boxplot(width=0.1, color="black", alpha=0.8)+

pdf("result/intron_splicing_order_first_last_dis.pdf",width = 8,height = 6);

#print(p_first_rela)
#print(p_last_rela)
print(p_all_s);

#print(p_all);

#multiplot(p_human, p_fly,p_pombe ,cols=1)
dev.off();



