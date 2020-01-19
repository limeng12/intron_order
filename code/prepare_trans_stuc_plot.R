anno<-read.table("./anno/gencode.v29lift37.annotation.gtf",sep="\t",as.is=TRUE,header = FALSE);

exon_anno<-anno[anno[,3]=="exon",];

colnames(exon_anno)<-c("chr","source","type","start","end","score","strand","","attr");


exon_anno[,"trans"]<-sapply(strsplit(exon_anno[,"attr"],";"),
                            function(x){  return(paste0(x[grepl("transcript_id",x)]) ) ;  }  );

exon_anno[,"trans"]<-str_sub(exon_anno[,"trans"],16,30);

exon_anno[exon_anno$strand=="+","strand"]<- 1;
exon_anno[exon_anno$strand=="-","strand"]<- -1;

exon_anno$start<-exon_anno$start-1;

exon_anno_one<-exon_anno[exon_anno$trans=="ENST00000226105",];



for_sushi<-exon_anno[,c("chr","start","end","trans","score","strand","type")];
colnames(for_sushi)<-c("chrom","start","end","gene","score","strand","type");
save(for_sushi,file="anno/sushi_trans_file.Rd",version = 2);
