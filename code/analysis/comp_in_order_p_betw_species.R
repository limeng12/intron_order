library(stringr)
library(readr)

fisher_meta<-function(p_values){

  chi_stat<- ( -2*sum( (p_values) ,na.rm = TRUE) );
  
  meta_p<-pchisq(chi_stat,df=length(p_values),lower.tail = FALSE,log.p = TRUE);
  meta_p
}

setwd("/Users/mengli/Documents/projects/iso");
source("code/multiplot.R");
#source("code/analysis/combine_rho.R");

human_order<-read.table("result/best_order/best_order.tsv",header = TRUE,as.is = TRUE,sep = "$")
human_order<-human_order[is.na( human_order$number_of_orders_have_same_prob ) |
                           (human_order$number_of_orders_have_same_prob ==1),];

zebrafish_order<-read.table("zebrafish/result/best_order.tsv",header = TRUE,as.is = TRUE,sep = "$")
zebrafish_order<-zebrafish_order[is.na( zebrafish_order$number_of_orders_have_same_prob ) |
                           (zebrafish_order$number_of_orders_have_same_prob ==1),];

fly_order<-as.data.frame(read_delim("fly/result/best_order.tsv",delim="$",
                                  col_names = TRUE,comment="#",col_types="ccdciddddd") )
#fly_order<-read.table("fly/result/best_order.tsv",header = TRUE,as.is = TRUE,sep = "$")

fly_order<-fly_order[is.na( fly_order$number_of_orders_have_same_prob ) |
                           (fly_order$number_of_orders_have_same_prob ==1),];

pombe_order<-read.table("pombe/result/best_order.tsv",header = TRUE,as.is = TRUE,sep = "$")
pombe_order<-pombe_order[is.na( pombe_order$number_of_orders_have_same_prob ) |
                           (pombe_order$number_of_orders_have_same_prob ==1),];



human_order_p_value<-human_order$relative_likelihood;

zebrafish_order_p_value<-zebrafish_order$relative_likelihood;

fly_order_p_value<-fly_order$relative_likelihood;

pombe_order_p_value<-pombe_order$relative_likelihood;


order_p_value_fr<-data.frame(
  label=c(rep("Human", length(human_order_p_value)), 
          rep("Zebrafish",length(zebrafish_order_p_value)),
          rep("Fly",length(fly_order_p_value) ),
            rep("Pombe", length(pombe_order_p_value))),
  sample_size=str_c("n=",c(rep(length(human_order_p_value),length(human_order_p_value)),
                rep(length(zebrafish_order_p_value),length(zebrafish_order_p_value)),
                rep(length(fly_order_p_value),length(fly_order_p_value)),
                rep(length(pombe_order_p_value),length(pombe_order_p_value)))
                ),
  value=c(human_order_p_value,zebrafish_order_p_value,fly_order_p_value,pombe_order_p_value)
);


order_p_value_fr$label <- factor(order_p_value_fr$label,
                                 levels = c("Pombe","Fly","Zebrafish","Human") ) 





human_order_cor<-human_order[!is.na(human_order$spearman_rho),]
zebrafish_order_cor<-zebrafish_order[!is.na(zebrafish_order$spearman_rho),]
fly_order_cor<-fly_order[!is.na(fly_order$spearman_rho),]
pombe_order_cor<-pombe_order[!is.na(pombe_order$spearman_rho),]


order_p_value_spearman_rho_fr<-data.frame(
  label=c(rep("Human", length(human_order_cor[,"spearman_rho"])), 
          rep("Zebrafish", length(zebrafish_order_cor[,"spearman_rho"])),
          rep("Fly",length(fly_order_cor[,"spearman_rho"]) ),
          rep("Pombe", length(pombe_order_cor[,"spearman_rho"]))),
  sample_size=str_c("n=",c(rep(length(human_order_cor[,"spearman_rho"]),length(human_order_cor[,"spearman_rho"])),
                rep(length(zebrafish_order_cor[,"spearman_rho"]),length(zebrafish_order_cor[,"spearman_rho"])),
                rep(length(fly_order_cor[,"spearman_rho"]),length(fly_order_cor[,"spearman_rho"])),
                rep(length(pombe_order_cor[,"spearman_rho"]),length(pombe_order_cor[,"spearman_rho"])))
  ),
  value=c(human_order_cor[,"spearman_rho"],zebrafish_order_cor[,"spearman_rho"],
          fly_order_cor[,"spearman_rho"],pombe_order_cor[,"spearman_rho"])
);

order_p_value_spearman_rho_fr$label <- factor(order_p_value_spearman_rho_fr$label,
                                          levels = c("Pombe","Fly","Zebrafish","Human") ) 


library(ggplot2)
data_summary <- function(x) {
  x<-x[!is.na(x)]
  mu <- mean(x)
  sigma1 <- mu-sd(x)
  sigma2 <- mu+sd(x)
  print( c(y=mu,ymin=sigma1,ymax=sigma2) )
  return(c(y=mu,ymin=sigma1,ymax=sigma2))
}


p1<-ggplot(order_p_value_fr,aes(x=sample_size, y=value, fill=label))+
  geom_violin()+  geom_point(stat="summary", fun.y="mean") + 
  #geom_errorbar(stat="summary", fun.data="mean_se", fun.args = list(mult = 1.96), width=0.1) +
  
  facet_grid(cols=vars(label), scales = "free", space = "free")+ylab("Log relative likelihood")+
  theme_minimal()+theme(legend.position="none")+xlab("")+ylim(-10,0);


p3<-ggplot(order_p_value_spearman_rho_fr,aes(x=sample_size, y=value,fill=label))+
  geom_violin()+geom_point(stat="summary", fun.y="mean") + 
  facet_grid(cols=vars(label), scales = "free", space = "free")+ylim(-1.1,1.1)+
  theme_minimal()+ylab("Spearman Rho")+theme(legend.position="none")+xlab("")


library(gridExtra)

pdf("result/comp_in_order_species.pdf",width = 8, height=6)
grid.arrange(p1,p3,nrow=2,heights=c(1,1) );


dev.off();



#order_p_value_fr[,"type"]<-rep("log P-value",nrow(order_p_value_fr) );

# p2<-ggplot(order_p_value_spearman_fr)+geom_violin(aes(x=sample_size, y=value,fill=label))+
#   facet_grid(cols=vars(label), scales = "free", space = "free")+ylim(-0.1,1.1)+
#   theme_minimal()+ylab("Spearman P-values")+theme(legend.position="none")+xlab("")
# 


#order_p_value_spearman_rho_fr[,"type"]<-rep("Rho",nrow(order_p_value_spearman_rho_fr) );


# p1_p3<-ggplot(rbind(order_p_value_fr, order_p_value_spearman_rho_fr) )+
#   geom_violin(aes(x=sample_size, y=value,fill=label))+
#   facet_grid(cols=vars(label),rows=vars(type), scales = "free", space = "fixed")+
#   theme_minimal()+theme(legend.position="none")+xlab("")+ylab("")


#grid.arrange(p1,p3,nrow=2)
#print(p1_p3)
#print(p1);
#print(p3);
#print(p2);

