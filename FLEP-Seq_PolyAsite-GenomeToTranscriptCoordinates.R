# written by Bo-Han Hou
library(data.table)
library(openxlsx)
Chr2cDNA <- function(model,transcripts,positions){
  model_exon<-copy(model[Type %in% c("exon","extended")])
  PosDT<-data.table(OriOrder=1:length(positions),"Transcript"=transcripts,"Position"=positions,"Position_dup"=positions)
  setkey(PosDT     ,Transcript,Position,Position_dup)
  setkey(model_exon,Transcript,Start   ,End)
  # calculate model cDNA position
  Hits<-foverlaps(PosDT,model_exon, type="within",nomatch=NA)[order(OriOrder)]
  Hits[Strand=="+",Position_cDNA:=Position-Start+Start_cDNA]
  Hits[Strand=="-",Position_cDNA:=Start-Position+End_cDNA]
  return(Hits$Position_cDNA)
}
sort_model<-function(x,strand){
  if(strand=="+"){
    return(x[order(End,Start,decreasing=FALSE),])
  }else{
    return(x[order(Start,End,decreasing=TRUE),] )
  }
}
Get_Attribute <- function(Type,Attributes,AttName){
  switch(
    Type,
    GTF = sub(sprintf('.*%s "(.*?)".*',AttName),'\\1',Attributes),
    GFF = sub(sprintf('.*%s=(.*?);.*',AttName),'\\1',Attributes)
  )
}
Fetch_GFF3 <- function(model.file,Representative_Genes_models){
  GFF3   <-fread(model.file,header = FALSE,blank.lines.skip=TRUE,
                 colClasses = c("character","character","character","numeric","numeric","character","character","character","character"),
                 col.names = c("Ref","Source","Type","Start","End","Score","Strand","Phase","Attribute"))
  mRNA   <-GFF3[grep("mRNA|Transcript",Type,perl=TRUE)
  ][,{Temp<-unlist(strsplit(Attribute,split = ";"));
  .(Gene=sub("Parent=","",grep("Parent=",Temp,value=TRUE)),
    Transcript=sub("ID=","",grep("ID=",Temp,value=TRUE)))}][]
  feature<-GFF3[grep("exon|CDS|five_prime_UTR|three_prime_UTR",Type,perl=TRUE,ignore.case = TRUE)
  ][grep("exon",Type,ignore.case = TRUE),`:=`(Type="exon")
  ][grep("five_prime_UTR",Type,ignore.case = TRUE),Type:="UTR5"
  ][grep("three_prime_UTR",Type,ignore.case = TRUE),Type:="UTR3"
  ][,`:=`(Transcript=sub("Parent=","",grep("Parent=",unlist(strsplit(Attribute,split = ";")),value=TRUE)))]
  # Split feature when parent of feature contains multiple values
  if(feature[grep(",",Transcript),.N]){
    feature <- feature[, strsplit(Transcript, ",") ,by = names(feature)
    ][,Transcript:=NULL
    ][,setnames(.SD,"V1","Transcript")]  
  }
  
  GFF3   <-mRNA[feature,on="Transcript",nomatch=0L
  ][,sort_model(.SD,.BY),by=.(Strand)
  ][order(Ref,Gene,Transcript)
  ][,.SD,.SDcols=c("Ref","Type","Start","End","Strand","Gene","Transcript")]
  
  # If Representative gene model contains one more intron.
  GFF3 <- GFF3[Transcript %in% Representative_Genes_models]
  
  # Assign cDNA position
  GFF3.Exon<-GFF3[Type=="exon"
  ][,Len:=End-Start+1
  ][,End_cDNA:=cumsum(Len),by=.(Gene,Transcript)
  ][,Start_cDNA:=End_cDNA-Len+1
  ][,.SD,.SDcols=c("Ref","Type","Start","End","Strand","Gene","Transcript","Len","Start_cDNA","End_cDNA")]
  
  GFF3.Features <-GFF3[Type!="exon"
  ][,Len:=End-Start+1
  ][Strand=="+",Start_cDNA:=Chr2cDNA(GFF3.Exon,Transcript,Start)
  ][Strand=="-",Start_cDNA:=Chr2cDNA(GFF3.Exon,Transcript,End)
  ][,End_cDNA:=Start_cDNA+Len-1]
  model <- rbind(GFF3.Exon,GFF3.Features)
  
  model<-model[order(Ref,Gene,Transcript,Start_cDNA,End_cDNA),
  ][,`:=`(num=.N,order=1:.N),by=.(Transcript,Type)
  ][,.SD,.SDcols=c("Ref","Type","Start","End","Strand","Gene","Transcript","Start_cDNA","End_cDNA","Len","num","order")]
  
  extended.fw.exon <-model[Type=="exon"][num==order][Strand=="+"][,Type:="extended"]
  extended.rc.exon <-model[Type=="exon"][num==order][Strand=="-"][,Type:="extended"]
  extended.fw.exon[,`:=`(Start=End+1,Start_cDNA=End_cDNA+1,order=order+1)][,`:=`(End=Start+999,End_cDNA=Start_cDNA+999)][,Len:=End-Start+1]
  extended.rc.exon[,`:=`(End=Start-1,Start_cDNA=End_cDNA+1,order=order+1)][,`:=`(Start=End-999,End_cDNA=Start_cDNA+999)][,Len:=End-Start+1]
  
  model<-rbind(model,extended.fw.exon,extended.rc.exon)[order(Ref,Gene,Transcript,order)]
  
  model[Type=="exon",TxLen:=sum(Len),by=.(Transcript)]
  return(model)
}

setwd(dir = "R:/Qsync/XRNs/")
GFF.file <-"./dataset/TAIR10_GFF3_genes_transposons.gff" 
# File is from https://www.arabidopsis.org/download/index-auto.jsp?dir=%2Fdownload_files%2FGenes%2FTAIR10_genome_release%2FTAIR10_gff3
RepGene.file <-"./dataset/TAIR10_representative_gene_models.txt"
# File is from https://www.arabidopsis.org/download/index-auto.jsp?dir=%2Fdownload_files%2FGenes%2FTAIR10_genome_release%2FTAIR10_gene_lists
Representative_Genes_models <- fread(file = RepGene.file,
                                     sep="\t",
                                     header = FALSE,
                                     col.names = "ID",
                                     skip = "AT")[grep("AT\\d+",ID)]$ID
Gene2Representative<-setNames(Representative_Genes_models,sub(".\\d+$","",Representative_Genes_models))

model.dt <- Fetch_GFF3(model.file = GFF.file, Representative_Genes_models = Representative_Genes_models)
# Representative coding genes list
Coding.genes <- model.dt[Type=="CDS",unique(Transcript)]

# PolyA site datasheet from Mo etal. doi: 10.1186/s13059-021-02543-4
PolyA.dt <- as.data.table(read.xlsx(xlsxFile = "./Reference/13059_2021_2543_MOESM2_ESM.xlsx",sheet = 1,startRow = 2))
PolyA.dt <- PolyA.dt[,.SD,.SDcols=c(1,3,6)][,setnames(.SD,c("Gene","Strand","PolyAsite"))][,c(1,3,2)]

# Mo etal defined PolyA sites have 1-nt shifting between the genes from + and - strands.
# for example:
# gene is on the strand "+"
# NNNNNNNNNNNNNNNNAAANNNNNNNNNNNNNNNNNNNN Genome
# ------------------> FELP Alignment
# ------------------> FELP Alignment
# ------------------> FELP Alignment
#                   ^
#                   |
#                   ploy A site. The position is located at 3' end of alignment
#
# gene is on the strand "-"
# NNNNNNNNNNNNNNNNNNAAANNNNNNNNNNNNNNNNNN Genome
#                   <------------------ FELP Alignment
#                   <------------------ FELP Alignment
#                   <------------------ FELP Alignment
#                  ^
#                  |
#                  ploy A site. The position is "1-nt" downstream of 3' end of alignment

# plus 1 nt to poly A site of genes which are located in "-" strand.
PolyA.dt <- PolyA.dt[Strand == '-',PolyAsite:=PolyAsite+1][]
# Give Representative transcript ID
PolyA.dt <- PolyA.dt[,Transcript:=Gene2Representative[Gene]][,c(1,4,2,3)][Transcript %in% Coding.genes]
# genome position to transcript position
PolyA.dt <- PolyA.dt[,PolyAsite.cDNA:=Chr2cDNA(model = model.dt[Type %in% c("exon","extended")][Gene %in% PolyA.dt$Gene],
                                               transcripts = Transcript,
                                               positions = PolyAsite)][]
TxLen.dt<-model.dt[Type=="exon",unique(.SD),.SDcols=c("Transcript","TxLen")]

# Filter out poly A sites which are located on > 140-nt downstream region
RemainPloyA<-PolyA.dt[TxLen.dt,on="Transcript",nomatch=0L][(PolyAsite.cDNA-TxLen+1)<=140]

# Final,remaining poly A sites are located from 3'UTR to 140 bp downstream region
# 5'UTR     CDS       3'UTR 
# -----===============------|<---140nt--->|
#                     |<----------------->| select polyA sites 
FinalPloyA<-RemainPloyA[model.dt[Type %in% c("UTR3","extended")],on="Transcript",nomatch=0L
][PolyAsite.cDNA<=End_cDNA & PolyAsite.cDNA>=Start_cDNA
][,.SD,.SDcols=c("Transcript","Strand","Type","PolyAsite.cDNA")]

fwrite(FinalPloyA,file = "./dataset/PolyA-sites_TAIR10_3UTR+140bp_FLEP-Seq.tsv",sep = "\t")
