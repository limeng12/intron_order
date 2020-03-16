exon_anno<-read.table("../abs/anno/gencode.v29lift37.annotation_exon.gtf",sep="\t",as.is=TRUE,header = FALSE);
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
save(for_sushi,file="result/sushi_trans_file.Rd",version = 2);
#write.table(for_sushi,file="anno/sushi_trans_file.tsv",col.names=TRUE,row.nanes=FALSE,sep="\t");
#for_sushi<-read.table("anno/sushi_trans_file.tsv",sep="\t",as.is=TRUE,header = FALSE);

load(file="anno/sushi_trans_file.Rd");

gene_trans_id_map<-read.table("anno/hg19_ensembl_gene_id_trans_id_map.tsv",
                              header = FALSE,as.is = TRUE,sep = "\t");

colnames(gene_trans_id_map)<-c("gene_id","trans_id","gene_symbol","trans_start","trans_end","strand",
                               "chr","gene_start","gene_end");


#Sushi_data = data(package = 'Sushi')
data("Sushi_genes.bed", package = 'Sushi' );

#data("Sushi_5C.bedpe", package = 'Sushi' );
chrom = "chr15"

chromstart = 72998000
chromend = 73020000

for_sushi_tmp<-for_sushi[(for_sushi$start>=chromstart) &
                           (for_sushi$end<=chromend)&
                           (for_sushi$chrom==chrom), ];

png("result/trans_struc.png")
pg = plotGenes(for_sushi_tmp,chrom,chromstart,chromend,types=for_sushi_tmp$type,
               bheight =0.2, plotgenetype = "box",bentline=FALSE,
               labeloffset = .4,labeltext = TRUE);
labelgenome(chrom, chromstart,chromend,side=3,n=3,scale="Mb");
dev.off();
