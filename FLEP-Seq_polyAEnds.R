library(data.table)
library(GenomicFeatures)
library(GenomicAlignments)


PolyASite<-function(gr,perfix,exbytx,flanking=100,polyA_only=TRUE){
  start(gr) <- start(gr)-flanking
  end(gr)   <- end(gr)+flanking
  Transcript<-names(gr)
  TxStart   <-start(gr)
  TxEnd     <-end(gr)
  GeneStrand=as.character(strand(gr))
  Param     <-ScanBamParam(which=gr)
  exon      <-exbytx[Transcript]
  exon<-endoapply(exon,function(x,flanking) {
    Frist.idx  <-mcols(x)$exon_rank==1
    Last.idx   <-mcols(x)$exon_rank==length(x)
    x[Frist.idx]<-resize(x[Frist.idx],width = width(x)[Frist.idx]+flanking,fix = 'end')
    x[Last.idx]<-resize(x[Last.idx] ,width = width(x)[Last.idx]+flanking,fix = 'start')
    return(x)
  },flanking=flanking)
  
  file.bam=paste(perfix,"bam",sep = ".")
  file.adapter=paste(perfix,"adapter.result.txt",sep = ".") 
  
  # Bam
  alignment.dt  <- FetchAlignments(bam = file.bam,exbyselect = exon,Param=Param)
  # adapter information
  # adapter.result.txt is output by adapterFinder.py
  adapter.info.dt<-fread(file.adapter,header=FALSE)[,c(1,2,3,4,5,6,7,9,11,12,14,15)
  ][,setnames(.SD,c("ReadID","read_align_strand","rna_strand","read_length",
                    "primer_type","genome_align_start","genome_align_end","polyA_type",
                    "f_primer_start","f_align_end",
                    "r_primer_start","r_align_start"))
  ][primer_type %in% c("R-F","F-R")
  ][rna_strand == GeneStrand]
  
  # columns:
  # read_core_id              7fab6cbf-a9ed-4427-951f-741515ddba0b,chr1,7026,8687
  # read_align_strand         +            
  # rna_strand                -            
  # read_length               839
  # primer_type	              R-F          
  # genome_align_start        71  <- coordinate of read
  # genome_align_end          820 <- coordinate of read
  # primer_score              1
  # polyA_type                T
  # f_primer_type             R
  # f_primer_start            1   <- coordinate of adapter 
  # f_align_end               70  <- coordinate of read
  # r_primer_type             F
  # r_primer_start            1   <- coordinate of adapter
  # r_align_start             821 <- coordinate of read
  
  # f_ : 5' read
  # r_ : 3' read
  
  if((GeneStrand=="-")){
    adapter.info.dt[,c("ReadID","Chr","end","start"):=tstrsplit(ReadID,split=",")][,start:=as.integer(start)][,end:=as.integer(end)]
  }else{
    adapter.info.dt[,c("ReadID","Chr","start","end"):=tstrsplit(ReadID,split=",")][,start:=as.integer(start)][,end:=as.integer(end)]
  }
  # adapter.info.dt <- adapter.info.dt[ReadID==qname]
  adapter.info.dt <- adapter.info.dt[ReadID %in% alignment.dt$ReadID]
  adapter.info.dt[r_align_start <= genome_align_end, end:=end+(r_align_start-genome_align_end -1 )  ]
  adapter.info.dt[f_align_end >= genome_align_start, start:=start+(f_align_end-genome_align_start+1)]
  adapter.info.dt <- adapter.info.dt[alignment.dt,on="ReadID",nomatch=0L]
  
  if((GeneStrand=="-")){
    adapter.info.dt[read_align_strand=='+',pad_seq:=substr(Alignment,f_align_end+1,genome_align_start-1)]
    adapter.info.dt[read_align_strand=='-',pad_seq:=substr(Alignment,start = read_length - r_align_start +1 +1,stop = read_length - genome_align_end +1 - 1)]
    adapter.info.dt<-adapter.info.dt[read_align_strand=='+',pad_seq2:=SubsetAlignment(start=genome_align_start,extend=20 ,Alignment,CIGAR,Strand=read_align_strand)]
    adapter.info.dt<-adapter.info.dt[read_align_strand=='-',pad_seq2:=SubsetAlignment(start=genome_align_end  ,extend=-20,Alignment,CIGAR,Strand=read_align_strand)]
    adapter.info.dt[ ,`:=`(polyA=grepl("TTTTT$",pad_seq),
                           pad_seqA=lapply(regexec("^T{1,}",pad_seq2),attr,"match.length"))]
    adapter.info.dt[,pad_seqA:=as.numeric(pad_seqA)][ pad_seqA>=1,end:=end+pad_seqA]
    adapter.info.dt<-adapter.info.dt[end>=TxStart]
  }else{
    adapter.info.dt[read_align_strand=='+',pad_seq:=substr(Alignment,genome_align_end +1,r_align_start-1)]
    adapter.info.dt[read_align_strand=='-',pad_seq:=substr(Alignment,start = read_length - genome_align_start+1+1 ,stop = read_length - f_align_end +1-1)]
    adapter.info.dt<-adapter.info.dt[read_align_strand=='-',pad_seq2:=SubsetAlignment(start=genome_align_start,extend=20,Alignment,CIGAR,Strand=read_align_strand)]
    adapter.info.dt<-adapter.info.dt[read_align_strand=='+',pad_seq2:=SubsetAlignment(start=genome_align_end,extend=-20,Alignment,CIGAR,Strand=read_align_strand)]
    
    adapter.info.dt[,`:=`(polyA=grepl("^AAAAA",pad_seq),
                          pad_seqA=lapply(regexec("A{1,}$",pad_seq2),attr,"match.length"))]
    adapter.info.dt[,pad_seqA:=as.numeric(pad_seqA)][ pad_seqA>=1,end:=end-pad_seqA]
    adapter.info.dt<-adapter.info.dt[end<=TxEnd]
  }
  polyA.info<-adapter.info.dt[,.SD,.SDcols=c("ReadID","polyA")]
  if(polyA_only){
    adapter.info.dt<-adapter.info.dt[polyA==TRUE]
  }
  
  if((GeneStrand=="-")){
    End.Count<-adapter.info.dt[,.(Read=.N),by=.(end)][order(-Read)][,end:=abs(end-TxEnd)+1-flanking]  
  }else{
    End.Count<-adapter.info.dt[,.(Read=.N),by=.(end)][order(-Read)][,end:=abs(end-TxStart)+1-flanking]
    
  }
  
  setnames(End.Count,"end","Position")
  return(list(Count=End.Count[order(Position)],polyA=polyA.info))
}

SubsetAlignment<-function(start,extend,Alignment,CIGAR,Strand){
  CIGAR.list    <- strsplit(sub("^\\d+","",CIGAR),split="[0-9]+")
  feq.list      <- lapply(strsplit(CIGAR,split="[A-z]"),as.integer)
  sub.seqs      <- sapply(1:length(CIGAR.list),function(x) {
    CIGAR.splited <- rep(CIGAR.list[[x]],feq.list[[x]])
    Alignment.i   <-strsplit(Alignment[x],"")[[1]]
    Align.info    <- data.table(Alignment=Alignment.i,CIGAR=CIGAR.splited[!CIGAR.splited %in% c("D","N")]
    )
    if(Strand[x]=="-") {
      Align.info <-Align.info[,Position:=.N:1][CIGAR %in% c("M"),Matched.Position:=.N:1]
    }else{
      Align.info <- Align.info <-Align.info[,Position:=1:.N][CIGAR %in% c("M"),Matched.Position:=1:.N][]
    }
    
    Matched.start<-Align.info[Position==start[x]]$Matched.Position
    Matched.end  <-ifelse(extend>0,Matched.start+extend-1,Matched.start+extend+1)
    sub.seq      <-paste0(Align.info[Matched.Position %in% seq(Matched.start,Matched.end) ]$Alignment,collapse = "")
    return(sub.seq)
  })
  return(sub.seqs)
}

buildfst <- function(DT,prefix){
  file_fst     <-paste0(prefix,".fst")
  file_fst_cidx<-paste0(file_fst,".cidx")
  DT_cidx      <-DT[, .(to=as.numeric(.N)), by = list(Transcript)
  ][,to:=cumsum(to)
  ][,from:=c(1,head(to,-1)+1)]
  
  write.fst(DT,file_fst, 100)
  write.fst(DT_cidx,file_fst_cidx, 100)
}

FetchAlignments<-function(bam,exbyselect,Param){
  # match transcript
  Alns        <- readGappedReads(file = bam, use.names=TRUE,param = Param)
  Hits.result <- findSpliceOverlaps(query=Alns, subject=exbyselect, ignore.strand=TRUE,singleEnd=TRUE)
  Hits.idx    <- queryHits(Hits.result[mcols(Hits.result)$compatible])
  Hits        <- Alns[Hits.idx]
  Hits.seq    <- qseq(Hits)
  return(data.table(
    "ReadID"=names(Hits),
    "CIGAR" =cigar(Hits),
    "Alignment"=as.character(Hits.seq)))
}


setwd("R:/Qsync/MY_script/GitHub/XRNs-DegradomeAnalysis/")
Gff.file<-"./dataset/TAIR10_GFF3_genes_transposons.gtf"
txdb   <- makeTxDbFromGFF(file = Gff.file,format = "gtf")
tx     <- transcripts(txdb,use.names=TRUE,columns="gene_id")[c("AT4G27440.1","AT1G12090.1","AT1G15690.1")]
tx     <- split(tx,names(tx))
exbytx <- exonsBy(txdb, "tx",use.names=TRUE)[c("AT4G27440.1","AT1G12090.1","AT1G15690.1")]

stat.count  <- fread(input = "./dataset/Bam_FLEP_seq/Mapping_reads.tsv",col.names = c("Sample","MappingCount"))
Raw.WT.1    <- data.table(Transcript="raw",Position=as.numeric(0),Read=stat.count[Sample=="FLEP-seq_WT_rep1"]$MappingCount)
Raw.WT.2    <- data.table(Transcript="raw",Position=as.numeric(0),Read=stat.count[Sample=="FLEP-seq_WT_rep2"]$MappingCount)
WT.1.info   <-lapply(tx,PolyASite,perfix="./dataset/Bam_FLEP_seq/FLEP-seq_WT_rep1",exbytx=exbytx,polyA_only=TRUE)
WT.2.info   <-lapply(tx,PolyASite,perfix="./dataset/Bam_FLEP_seq/FLEP-seq_WT_rep2",exbytx=exbytx,polyA_only=TRUE)

# Count 
# the coordinates are exon + intron
WT.1        <- rbindlist(lapply(WT.1.info,function(x) x$Count),idcol = "Transcript")[,TP40M:=Read/stat.count[Sample=="FLEP-seq_WT_rep1"]$MappingCount*40*10^6][]
WT.2        <- rbindlist(lapply(WT.2.info,function(x) x$Count),idcol = "Transcript")[,TP40M:=Read/stat.count[Sample=="FLEP-seq_WT_rep2"]$MappingCount*40*10^6][]

# Convert coordinates to Genome
txStart<-rbindlist(lapply(tx,data.frame),idcol = "Transcript")[,tx_start:=ifelse(strand=="+",start,end)][,.(Transcript,seqnames,tx_start,strand)]

WT.1_genome<-WT.1[txStart,on="Transcript"
][strand=="+",Position_chr:=tx_start+Position-1
][strand=="-",Position_chr:=tx_start-Position+1
][,.SD,.SDcols=c("Transcript","seqnames","strand","Position_chr","Read","TP40M")
][,setnames(.SD,"Position_chr","Position")]

WT.2_genome<-WT.2[txStart,on="Transcript"
][strand=="+",Position_chr:=tx_start+Position-1
][strand=="-",Position_chr:=tx_start-Position+1
][,.SD,.SDcols=c("Transcript","seqnames","strand","Position_chr","Read","TP40M")
][,setnames(.SD,"Position_chr","Position")]

dir.create("FLEP-seq_polyAEnds",showWarnings = FALSE)
fwrite(x = WT.1_genome,file = "./FLEP-seq_polyAEnds/3Ends_polyARead_FLEP-seq_WT_R1.tsv",sep = "\t")
fwrite(x = WT.2_genome,file = "./FLEP-seq_polyAEnds/3Ends_polyARead_FLEP-seq_WT_R2.tsv",sep = "\t")

# polyA reads and others
polyA.info<-rbindlist( list("WT_rep1"=rbindlist(lapply(WT.1.info,function(x) x$polyA),idcol = "Transcript"),
                            "WT_rep2"=rbindlist(lapply(WT.2.info,function(x) x$polyA),idcol = "Transcript")),idcol = "Repeat")
fwrite(polyA.info,file = "./FLEP-seq_polyAEnds/Target.polyA.qnames.tsv",sep="\t")

