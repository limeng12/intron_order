################################################plot simulated result####################################
library(stringr)
library(ggplot2)
library(latex2exp)

setwd("/Users/mengli/Documents/projects/iso/pombe");
source("../code/multiplot.R");
load("result/simulation_orders.Rd");

trans_order_long<-trans_order;

#total_transcript_long<-0;

meet_cond_index<-c();
# where <4 means intron number <3,i.e skip transcript that contain only two introns
for(i in 1:length(trans_order_long) ){
  if( length(trans_order_long[[i]])<4){
    #trans_order_long[[i]]<-
    next;
  }
  meet_cond_index<-c(meet_cond_index,i);
  
  #total_transcript_long<-total_transcript_long+1;
}

trans_order_long<-trans_order_long[meet_cond_index];

total_transcript_long<-length(trans_order_long);




spearma_re_super_long<-read.table("result/best_order_simulation_super_long.tsv",sep="$",as.is = TRUE,header = TRUE);

spearma_r_vec_super_long<-c()

for( i in 1:nrow(spearma_re_super_long) ){
  
  orders<-as.numeric(str_split(spearma_re_super_long[i,"best_order"],",")[[1]] );
  
  trans_id<-spearma_re_super_long[i,"transcript_id"];
  
  #given_orders<-trans_order_long[trans_order_long[,1]==trans_id,2];
  if(is.null(trans_order_long[[trans_id]]) || is.na(trans_order_long[[trans_id]])){
    next;
  }
  given_orders_numer<-trans_order_long[[trans_id]];
  #given_orders_numer<-as.numeric(str_split(given_orders,",")[[1]] );
  given_orders_numer<-given_orders_numer[-1*which.max(given_orders_numer)];
  
  spearma_r<-cor(orders,given_orders_numer,method="spearman");
  
  spearma_r_vec_super_long<-c(spearma_r_vec_super_long,spearma_r);
}

trans_order_long_char<-sapply(trans_order_long,paste0,collapse=",");

#spearma_r_vec_mat_long<-cbind(simulation_re_long[,"gene_symbol"],simulation_re_long[,"transcript_id"],
#                              simulation_re_long[,"best_order"],
#                              trans_order_long_char[simulation_re_long[,"transcript_id"]],spearma_r_vec_long)



simulation_re_long<-read.table("result/best_order_simulation_long.tsv",sep="$",as.is = TRUE,header = TRUE);

spearma_r_vec_long<-c()

for( i in 1:nrow(simulation_re_long) ){
  
  orders<-as.numeric(str_split(simulation_re_long[i,"best_order"],",")[[1]] );
  
  trans_id<-simulation_re_long[i,"transcript_id"];
  
  #given_orders<-trans_order_long[trans_order_long[,1]==trans_id,2];
  if(is.null(trans_order_long[[trans_id]]) || is.na(trans_order_long[[trans_id]])){
    next;
  }
  given_orders_numer<-trans_order_long[[trans_id]];
  #given_orders_numer<-as.numeric(str_split(given_orders,",")[[1]] );
  given_orders_numer<-given_orders_numer[-1*which.max(given_orders_numer)];
  
  spearma_r<-cor(orders,given_orders_numer,method="spearman");
  
  spearma_r_vec_long<-c(spearma_r_vec_long,spearma_r);
}

trans_order_long_char<-sapply(trans_order_long,paste0,collapse=",");

# spearma_r_vec_mat<-cbind(simulation_re[,"gene_symbol"],simulation_re[,"transcript_id"],simulation_re[,"best_order"],
#                          trans_order_long_char[simulation_re[,"transcript_id"]],spearma_r_vec)



simulation_re_pair<-read.table("result/best_order_simulation_pair.tsv",sep="$",as.is = TRUE,header = TRUE);


spearma_r_vec_pair<-c()

for( i in 1:nrow(simulation_re_pair) ){
  orders<-as.numeric(str_split(simulation_re_pair[i,"best_order"],",")[[1]] );
  
  trans_id<-simulation_re_pair[i,"transcript_id"];
  
  #given_orders<-trans_order_long[trans_order_long[,1]==trans_id,2];
  if(is.null(trans_order_long[[trans_id]]) || is.na(trans_order_long[[trans_id]])){
    next;
  }
  
  given_orders_numer<-trans_order_long[[trans_id]];
  given_orders_numer<-given_orders_numer[-1*which.max(given_orders_numer)];
  #if( length(given_orders_numer)!= length(orders) ){
  #  next;
  #}
  
  spearma_r<-cor(orders,given_orders_numer,method="spearman")
  
  spearma_r_vec_pair<-c(spearma_r_vec_pair,spearma_r);
}


#trans_order_long_char<-sapply(trans_order_long,paste0,collapse=",");

# spearma_r_vec_pair_mat<-cbind(simulation_re_pair[,"transcript_id"],simulation_re_pair[,"best_order"],
#                               trans_order_long_char[simulation_re_pair[,"transcript_id"]],spearma_r_vec_pair)



#setwd("/Users/mengli/Documents/projects/iso/pombe");

spearma_r_vec_long_all<-c(spearma_r_vec_super_long,spearma_r_vec_long,spearma_r_vec_pair);

spearma_r_vec_long_name_all<-c(rep(paste0("super long read"," L=5,000\nRho>0.99=",sum(spearma_r_vec_super_long>0.99),""),
                                   length(spearma_r_vec_super_long)),
                               rep(paste0("long read"," L=800\nRho>0.99=",sum(spearma_r_vec_long>0.99),""),
                                   length(spearma_r_vec_long)),
                               rep(paste0("pair-end read"," L=150\nRho>0.99=",sum(spearma_r_vec_pair>0.99),""),
                                   length(spearma_r_vec_pair)) );

spearma_r_vec_long_name_all <- factor(spearma_r_vec_long_name_all,
                                      levels =c(
                                        paste0("pair-end read"," L=150\nRho>0.99=",sum(spearma_r_vec_pair>0.99),""),
                                       
                                        paste0("long read"," L=800\nRho>0.99=",sum(spearma_r_vec_long>0.99),""),
                                        paste0("super long read"," L=5,000\nRho>0.99=",sum(spearma_r_vec_super_long>0.99),"") ) )

spearma_r_vec_long_mat<-data.frame(spearma_r_vec_long_all=spearma_r_vec_long_all,
                                   spearma_r_vec_long_name_all=spearma_r_vec_long_name_all)

p1<-ggplot(spearma_r_vec_long_mat)+
  #ggtitle( TeX(paste("N=", total_transcript_long, " (intron number $\\geq$ 3) multi-intron containing transcripts simulation") ) )+
  geom_violin(aes(x="",y=spearma_r_vec_long_all,fill=spearma_r_vec_long_name_all))+
  theme_minimal()+xlab("")+ylab("")+ylim(-1.1,1.1)+
  facet_grid(cols=vars(spearma_r_vec_long_name_all), scales = "free", space = "free")+
  theme(legend.position="none",text = element_text(size=11) ,
        axis.text.y = element_text(size=13) )
  
# p3<-ggplot()+ggtitle("")+
#   geom_violin(aes(x="~5,000",y=spearma_r_vec_long))+
#   xlab(paste0( " (Rho>0.99=",sum(spearma_r_vec_long>0.99),")") )+
#   theme_minimal()+ylab("transcript number")+scale_y_log10()
# 
# p4<-ggplot()+ggtitle("")+
#   geom_violin(aes(x="~800",y=spearma_r_vec))+
#   xlab(paste0( " (Rho>0.99=",sum(spearma_r_vec>0.99),")") )+
#   theme_minimal()+ylab("transcript number")+scale_y_log10()
# 
# p5<-ggplot()+ggtitle("1300 multi-intron transcripts simulation")+
#   geom_violin(aes(x="~150",y=spearma_r_vec_pair))+
#   xlab(paste0(" (Rho>0.99=",sum(spearma_r_vec_pair>0.99),")") )+
#   theme_minimal()+ylab("transcript number")+scale_y_log10()

#multiplot(p1,p2,cols=2);
#multiplot(p5,p4,p3,pcols=3);


pdf("result/pombe_simulation_reads_cor.pdf",width = 5.5, height=4);

print(p1);

dev.off();


# p1<-ggplot()+
#   geom_violin(aes(x="Distribution of spearman rho long read",y=spearma_r_vec),fill="#E69F00")+
#   ylab("")+
#   xlab("")+ylim(-1,1)+theme_minimal();
# 
# p2<-ggplot()+
#   geom_violin(aes(x="Distribution of spearman rho pair-end short read",y=spearma_r_vec_pair),fill="#E69F00")+
#   ylab("")+
#   xlab("")+ylim(-1,1)+theme_minimal();
# 
