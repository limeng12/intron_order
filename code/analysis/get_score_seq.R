#pombe_fa<-FaFile("fasta/Schizosaccharomyces_pombe.ASM294v2.fa");
library(Rsamtools)
library(GenomicRanges)
library(BSgenome);

#setwd("/Users/mengli/Documents/projects/iso");


get_score_5ss_file<-function(seq){
  
  cat(paste0(seq) , sep="\n",file="/var/tmp/ss5_seq.txt",append = FALSE);
  
  cur_dir<-getwd();
  
  #writeLines(seq,con = file("/var/tmp/test5.seq") );
  
  setwd("/Users/mengli/Documents/software/maxentscan");
  
  system(" perl score5.pl /var/tmp/ss5_seq.txt > /var/tmp/test5.score");
  
  score_tbl<-read.table("/var/tmp/test5.score",as.is = TRUE);
  
  setwd(cur_dir);
  
  score_tbl[,2];
  
}

get_score_3ss_file<-function(seq){
  
  cat(paste0(seq) , sep="\n",file="/var/tmp/ss3_seq.txt",append = FALSE);
  
  cur_dir<-getwd();
  
  #writeLines(seq,con = file("/var/tmp/test5.seq") );
  
  setwd("/Users/mengli/Documents/software/maxentscan");
  
  system(" perl score3.pl /var/tmp/ss3_seq.txt > /var/tmp/test3.score");
  
  score_tbl<-read.table("/var/tmp/test3.score",as.is = TRUE);
  
  setwd(cur_dir);
  
  score_tbl[,2];
  
}



get_score_5ss<-function(seq){
  cur_dir<-getwd()
  
  writeLines(seq,con = file("/var/tmp/test5.seq") );
  
  setwd("/Users/mengli/Documents/software/maxentscan");
  
  system(" perl score5.pl /var/tmp/test5.seq > /var/tmp/test5.score");
  
  score_tbl<-read.table("/var/tmp/test5.score",as.is = TRUE);
  
  setwd(cur_dir)
  
  score_tbl[,2];
  
}

get_score_3ss<-function(seq){
  
  writeLines(seq,con = file("/var/tmp/test3.seq") );
  
  cur_dir<-getwd()
  
  setwd("/Users/mengli/Documents/software/maxentscan");
  
  system(" perl score3.pl /var/tmp/test3.seq > /var/tmp/test3.score")
  
  score_tbl<-read.table("/var/tmp/test3.score",as.is = TRUE);
  
  setwd(cur_dir)
  
  score_tbl[,2];
  
}

get_seq_region<-function(x,t_fa){
  
  fa<-t_fa
  if(x==""){
    return("")
  }
  
  chr<-str_c(sapply(str_split(x,":"),"[",1) );
  if(!str_detect(chr,"chr")){
    chr<-str_c("chr",chr);
    
    #chr<-str_sub(chr,4);
  }
  

  
  
  start<-as.numeric( sapply(str_split(x,":"),"[",2) );
  
  end<-as.numeric( sapply(str_split(x,":"),"[",3) );
  
  strand<-sapply(str_split(x,":"),"[",4);
  
  if(str_starts(chr,"chrUn")){
    return( paste0(rep("N",(end-start+1) ) ,collapse="")   );
  }
  
  
  gr <- GRanges(
    seqnames = c(chr),
    ranges = IRanges( start, end ),
    strand = strand
  );
  #print(gr)
  
  #a<-getSeq(BSgenome.Hsapiens.UCSC.hg19,chr,start,end,strand="+");
  
  a<-as.character( Rsamtools::getSeq(fa, gr) );
  a
  #as.character(a);
  
}


test<-function(){
  
  h_fasta_path<-FaFile("/Volumes/mengli/anno/human_gencode/GRCh37.primary_assembly.genome.fa")
  
  p_fasta_path<-FaFile("/Volumes/mengli/anno/fission_yeast/Schizosaccharomyces_pombe.ASM294v2.dna.toplevel.fa")
  
  get_seq_region("chr1:1:1000:+",h_fasta_path);
  
  get_seq_region("chrI:18304:18312:+", p_fasta_path)->a
  
  get_score_5ss(a)
  
}
