# Some functions of Guitar are modified from source code v2.10.0 by Bo-Han Hou
# Modify scripts to reduce RAM requirement & improve performance
# Name of modified function end with ".fast"


library(Guitar)
library(data.table)
library(ggplot2)
library(cowplot)
library(scales)
#### Modified GuitarPlot functions ####
process.reprot <- function(text){
  message(format(Sys.time(), "%Y/%m/%d %H:%M:%S"), " ... ", text)
}

GuitarPlotFast <- function(txGTF = NULL, txGFF = NULL, txGenomeVer = NULL, txTxdb = NULL, 
                           txGuitarTxdb = NULL, txGuitarTxdbSaveFile = NA, stBedFiles = NULL, 
                           stGRangeLists = NULL, stGroupName = NULL, stAmblguity = 5, 
                           stSampleNum = 3, stSampleModle = c("Equidistance","random"), txfiveutrMinLength = 100, 
                           txcdsMinLength = 100, txthreeutrMinLength = 100, txlongNcrnaMinLength = 100, 
                           txlncrnaOverlapmrna = FALSE, txpromoterLength = 1000, txtailLength = 1000, 
                           txAmblguity = 5, txPrimaryOnly = FALSE, txTxComponentProp = NULL, 
                           txMrnaComponentProp = NULL, txLncrnaComponentProp = NULL, 
                           mapFilterTranscript = TRUE, headOrtail = TRUE, enableCI = TRUE, 
                           pltTxType = c("tx", "mrna", "ncrna"), overlapIndex = 1, 
                           siteLengthIndex = 1, adjust = 1, CI_ResamplingTime = 1000, 
                           CI_interval = c(0.025, 0.975), miscOutFilePrefix = NA) 
{
  stSampleModle<-match.arg(stSampleModle)
  pltTxType    <-match.arg(pltTxType,several.ok=TRUE)
  
  genomeVersion2Txdb <- list(hg18 = "TxDb.Hsapiens.UCSC.hg18.knownGene", 
                             hg19 = "TxDb.Hsapiens.UCSC.hg19.knownGene", 
                             hg38 = "TxDb.Hsapiens.UCSC.hg38.knownGene", 
                             mm9 = "TxDb.Mmusculus.UCSC.mm9.knownGene", 
                             mm10 = "TxDb.Mmusculus.UCSC.mm10.knownGene")
  if (headOrtail) {
    txpromoterLength <- txpromoterLength
    txtailLength <- txtailLength
  }  else {
    txpromoterLength <- 0
    txtailLength <- 0
  }
  
  process.reprot("Start")
  process.reprot("Make/Load GuitarTxdb")
  guitarTxdb <- getGuitarTxdb.fast(txGTF = txGTF, txGFF = txGFF, 
                                   txGenomeVer = txGenomeVer, txTxdb = txTxdb, txGuitarTxdb = txGuitarTxdb, 
                                   txfiveutrMinLength = txfiveutrMinLength, txcdsMinLength = txcdsMinLength, 
                                   txthreeutrMinLength = txthreeutrMinLength, txlongNcrnaMinLength = txlongNcrnaMinLength, 
                                   txlncrnaOverlapmrna = txlncrnaOverlapmrna, txpromoterLength = txpromoterLength, 
                                   txtailLength = txtailLength, txAmblguity = txAmblguity, 
                                   txTxComponentProp = txTxComponentProp, txMrnaComponentProp = txMrnaComponentProp, 
                                   txLncrnaComponentProp = txLncrnaComponentProp, txPrimaryOnly = txPrimaryOnly, 
                                   pltTxType = pltTxType, genomeVersion2Txdb)
  
  
  if (!(is.na(txGuitarTxdbSaveFile))) {
    txGuitarTxdbSaveFile <- paste(txGuitarTxdbSaveFile,"GuitarTxdb", sep = "-")
    saveRDS(guitarTxdb, file = txGuitarTxdbSaveFile)
  }
  
  if (any(!pltTxType %in% guitarTxdb$txTypes)) {
    pltTxType.notfound <- pltTxType[!pltTxType %in% guitarTxdb$txTypes]
    process.reprot(paste(pltTxType.notfound," is not found in tx annotation"))
    stop(paste("Cannot plot distribution for", pltTxType.notfound))
  }
  process.reprot("import BED file")
  sitesGroup <- Guitar:::.getStGroup(stBedFiles = stBedFiles, 
                                     stGRangeLists = stGRangeLists, 
                                     stGroupName = stGroupName)
  
  GroupNames <- names(sitesGroup)
  sitesGroupNum <- length(sitesGroup)
  sitesPointsNormlize <- list()
  sitesPointsRelative <- list()
  pointWeight <- list()
  for (i in seq_len(sitesGroupNum)) {
    GroupName = GroupNames[[i]]
    process.reprot(paste("Sample", stSampleNum, "points for",GroupName, sep = " "))
    
    sitesPoints <- samplePoints.fast(sitesGroup[i], stSampleNum = stSampleNum, 
                                     stAmblguity = stAmblguity, pltTxType = pltTxType, 
                                     stSampleModle = stSampleModle, mapFilterTranscript = mapFilterTranscript,guitarTxdb)
    
    for (txType in pltTxType) {
      process.reprot(paste("Normalize for",txType,sep = " "))
      sitesPointsNormlize[[txType]][[GroupName]] <- normalize.fast(sitesPoints, guitarTxdb, txType, overlapIndex, siteLengthIndex)
      sitesPointsRelative[[txType]][[GroupName]] <- sitesPointsNormlize[[txType]][[GroupName]][[1]]
      pointWeight[[txType]][[GroupName]] <- sitesPointsNormlize[[txType]][[GroupName]][[2]]
    }
  }
  p.list=list()
  process.reprot("Create plot")
  for (txType in pltTxType) {
    print(paste("start figure plotting for", txType,"..."))
    title <- switch(txType,
                    tx="Distribution on Transcript",
                    mrna="Distribution on mRNA", 
                    ncrna="Distribution on ncRNA")
    densityDataframe_CI <- Guitar:::.generateDensity_CI(sitesPointsRelative[[txType]],pointWeight[[txType]], CI_ResamplingTime, adjust = adjust, enableCI = enableCI)
    p.list[[txType]] <- Guitar:::.plotDensity_CI(densityDataframe_CI, componentWidth = guitarTxdb[[txType]]$componentWidthAverage_pct, headOrtail, title, enableCI = enableCI)
  }
  process.reprot("End")
  return(p.list)
} 
getGuitarTxdb.fast <-function (txGTF = NULL, txGFF = NULL, 
                               txGenomeVer = NULL, txTxdb = NULL, 
                               txGuitarTxdb = NULL, txfiveutrMinLength = 100, 
                               txcdsMinLength = 100, txthreeutrMinLength = 100, 
                               txlongNcrnaMinLength = 100, txlncrnaOverlapmrna = FALSE, 
                               txpromoterLength = 1000, txtailLength = 1000, 
                               txAmblguity = txAmblguity, txTxComponentProp = NULL, 
                               txMrnaComponentProp = NULL, txLncrnaComponentProp = NULL, 
                               txPrimaryOnly = FALSE, pltTxType = c("tx", "mrna", "ncrna"), genomeVersion2Txdb) 
{
  if (!(is.null(txGuitarTxdb))) {
    if (!file.exists(txGuitarTxdb)) {
      stop(paste0("\"txGuitarTxdb\" is not exist, please double check!"))
    }
    process.reprot(paste("load", txGuitarTxdb, "as GuitarTxdb...", sep = " "))
    guitarTxdb <- readRDS(txGuitarTxdb)
    if(!is.null(txTxComponentProp)){
      guitarTxdb$tx$componentWidthAverage_pct      <- setNames(txTxComponentProp/sum(txTxComponentProp),c("promoter","rna","tail"))
      guitarTxdb$tx$componentStartAverage_pct      <- setNames(c(0,head(cumsum(guitarTxdb$tx$componentWidthAverage_pct),-1)),c("promoter","rna","tail"))
    }
    if(!is.null(txLncrnaComponentProp)){
      guitarTxdb$ncrna$componentWidthAverage_pct   <- setNames(txLncrnaComponentProp/sum(txLncrnaComponentProp),c("promoter","ncrna","tail"))
      guitarTxdb$ncrna$componentStartAverage_pct   <- setNames(c(0,head(cumsum(guitarTxdb$ncrna$componentWidthAverage_pct),-1)),c("promoter","ncrna","tail"))
    }
    if(!is.null(txMrnaComponentProp)){
      guitarTxdb$mrna$componentWidthAverage_pct    <- setNames(txMrnaComponentProp/sum(txMrnaComponentProp),c("promoter","utr5","cds","utr3","tail"))
      guitarTxdb$mrna$componentStartAverage_pct    <- setNames(c(0,head(cumsum(guitarTxdb$mrna$componentWidthAverage_pct),-1)),c("promoter","utr5","cds","utr3","tail"))
    }
    
    # adjust position, width and ratio of tx if length of promoter or tail in guitarTxdb are different from user give
    for (txType in unlist(guitarTxdb$txTypes)){
      txpromoterLengthOri <- guitarTxdb[[txType]]$componentWidth[1,1]
      txtailLengthOri <- guitarTxdb[[txType]]$componentWidth[1,ncol(guitarTxdb[[txType]]$componentWidth)]
      if (txpromoterLengthOri != txpromoterLength || txtailLengthOri != txtailLength ){
        tx <- unlist(guitarTxdb[[txType]]$tx)
        idx.exon.frist    <- start(guitarTxdb[[txType]]$tx@partitioning)
        idx.exon.last     <- end(guitarTxdb[[txType]]$tx@partitioning)
        plus.strand       <- strand(tx[idx.exon.last])=="+"
        
        if ( txpromoterLengthOri != txpromoterLength) {
          guitarTxdb[[txType]]$componentWidth[,1]<-txpromoterLength
          start(tx[idx.exon.frist][plus.strand]) <- start(tx[idx.exon.frist][plus.strand])+txpromoterLengthOri-txpromoterLength
          end(tx[idx.exon.last][!plus.strand])  <- end(tx[idx.exon.last][!plus.strand])-txpromoterLengthOri+txpromoterLength
          
        }
        if ( txtailLengthOri != txtailLength ) {
          guitarTxdb[[txType]]$componentWidth[,ncol(guitarTxdb[[txType]]$componentWidth)]<-txtailLength
          end(tx[idx.exon.last][plus.strand])     <- end(tx[idx.exon.last][plus.strand])-txtailLengthOri+txtailLength
          start(tx[idx.exon.frist][!plus.strand]) <- start(tx[idx.exon.frist][!plus.strand])+txtailLengthOri-txtailLength
        }
        
        guitarTxdb[[txType]]$tx         <- relist(tx,guitarTxdb[[txType]]$tx)
        guitarTxdb[[txType]]$endPoint   <- t(apply(guitarTxdb[[txType]]$componentWidth,1,cumsum))
        guitarTxdb[[txType]]$startPoint <- guitarTxdb[[txType]]$endPoint-guitarTxdb[[txType]]$componentWidth +1
        guitarTxdb[[txType]]$txLength   <- rowSums(guitarTxdb[[txType]]$componentWidth)
        
        guitarTxdb[[txType]]$componentWidthPtc     <- guitarTxdb[[txType]]$componentWidth/rowSums(guitarTxdb[[txType]]$componentWidth)
        componentWidth_ratio_avg        <- colSums(guitarTxdb[[txType]]$componentWidthPtc)
        guitarTxdb[[txType]]$componentWidthAverage <- floor(componentWidth_ratio_avg/sum(componentWidth_ratio_avg) * 500 + 0.5)
      }
    }
    
    
  }else{
    
    txdb <- Guitar:::.getTxdb(txGTF = txGTF, txGFF = txGFF, txGenomeVer = txGenomeVer, 
                              txTxdb = txTxdb, genomeVersion2Txdb)
    guitarTxdb <- makeGuitarTxdb(txdb, txfiveutrMinLength = txfiveutrMinLength, 
                                 txcdsMinLength = txcdsMinLength, txthreeutrMinLength = txthreeutrMinLength, 
                                 txlongNcrnaMinLength = txlongNcrnaMinLength, txlncrnaOverlapmrna = txlncrnaOverlapmrna, 
                                 txpromoterLength = txpromoterLength, txtailLength = txtailLength, 
                                 txAmblguity = txAmblguity, txTxComponentProp = txTxComponentProp, 
                                 txMrnaComponentProp = txMrnaComponentProp, txLncrnaComponentProp = txLncrnaComponentProp, 
                                 txPrimaryOnly = txPrimaryOnly, pltTxType = pltTxType)
  }
  return(guitarTxdb)
}

samplePoints.fast <- function(sitesGrangelists, stSampleNum = 5, stAmblguity = 5, 
                              pltTxType = c("tx", "mrna", "ncrna"), stSampleModle = "Equidistance", 
                              mapFilterTranscript = FALSE, guitarTxdb) 
{
  mkSitesPoints <- function(sitesWidth,stSampleModle,stSampleNum){
    Sample.fun <- switch(stSampleModle,
                         Equidistance=function(x, i) {
                           if (i == 1) {
                             round(x/2)
                           }    else {
                             round(seq(1, x - 1, length.out = i))
                           }
                         },
                         random=function(x, i) {
                           if (i == 1) {
                             round(x/2)
                           }
                           else {
                             sort(sample(x, i, replace = FALSE))
                           }
                         })
    # build sitesPoints matrix based on unique sitesWidth
    sitesPointsVector <- vector(mode = "integer",length(sitesWidth))
    sitesWidth.uni    <- matrix(sort(unique(sitesWidth)), ncol = 1)
    sitesWidth.idx    <- setNames(1:length(sitesWidth.uni), sitesWidth.uni)
    sitesPoints.uni   <- matrix(data = apply(sitesWidth.uni, MARGIN = 1, FUN = Sample.fun, i=stSampleNum),
                                ncol = stSampleNum, byrow = T)
    # Covert sitesWidth to sitesPointsVector
    sitesPointsVector <- t(sitesPoints.uni[sitesWidth.idx[as.character(sitesWidth)],])
    return(sitesPointsVector)
  }
  
  sitesGRangeDataframe <- list()
  stSampleNum <- 2 * stSampleNum - 1
  
  for (txType in pltTxType) {
    
    process.reprot("Step1: Convert genomic coordinates to transcript coordinates")
    mapsiteGRanges <- GRangesListmapToTranscripts.fast(sitesGrangelists[[1]], 
                                                       mapFilterTranscript, 
                                                       guitarTxdb[[txType]]$tx)
    sitesWidth <- width(mapsiteGRanges)
    process.reprot(paste("Step2: Build sitesPoints by",stSampleModle,sep=" "))
    sitesPointsVector <- mkSitesPoints(sitesWidth = sitesWidth,stSampleModle = stSampleModle,stSampleNum = stSampleNum)
    
    process.reprot("Step3: Build GRange")
    sitesGRangeDataframe[[txType]]         <- mapsiteGRanges[rep(seq_len(length(mapsiteGRanges)), each=stSampleNum),]
    start(sitesGRangeDataframe[[txType]])  <- sitesPointsVector+start(sitesGRangeDataframe[[txType]]) - 1L
    end(sitesGRangeDataframe[[txType]])    <- start(sitesGRangeDataframe[[txType]]) + 1L
    strand(sitesGRangeDataframe[[txType]]) <- "*"
    
    # use table() to count and sign them
    xHits.times     <- table(mapsiteGRanges$xHits)
    dt <- data.table(sitesLength     = sitesWidth, 
                     xHits           = mapsiteGRanges$xHits, 
                     pointsOverlapTx = as.vector(xHits.times[as.character(mapsiteGRanges$xHits)]))
    mcols(sitesGRangeDataframe[[txType]]) <- dt[rep(seq_len(.N), each=stSampleNum), ]
    
  }
  return(sitesGRangeDataframe)
}
GRangesListmapToTranscripts.fast <- function(site, mapFilterTranscript = FALSE, transcripts)
{ 
  names(site) <- 1:length(site)
  if(mapFilterTranscript) {
    xWidthes <- sum(width(site))
    names(xWidthes) <- names(site)
  }
  if(is(site,"CompressedGRangesList")){
    tx_coord <- mapToTranscripts(unlist(site), transcripts,ignore.strand=FALSE)
  }else{
    tx_coord <- mapToTranscripts(site, transcripts,ignore.strand=FALSE)
  }
  
  xHit_txHit_joint <- paste(names(tx_coord), tx_coord$transcriptsHits, sep='-')
  tx_coord_grouped <- split(tx_coord, xHit_txHit_joint)
  mapping_reduced  <- reduce(tx_coord_grouped)
  
  # some sites map to tx has multiple regions after reduce because of isoform of tx
  tx.region.num           <- elementNROWS(mapping_reduced)
  tx.region.selected      <- names(tx.region.num[tx.region.num==1])
  tx_coord_filtered       <- unlist(mapping_reduced[tx.region.selected])
  mcols(tx_coord_filtered)<- data.table(name=names(tx_coord_filtered))[,c("xHits","txHits"):=tstrsplit(name,split='-',type.convert=TRUE)][,c(2:3)]
  
  # remove hits whose length smaller than sites because of isoform of tx
  if(mapFilterTranscript) {
    tx_coord_filtered_width <- width(tx_coord_filtered)
    tx_coord_filtered <- tx_coord_filtered[tx_coord_filtered_width == xWidthes[tx_coord_filtered$xHits]]
  }
  idx <- GenomicRanges:::get_out_of_bound_index(tx_coord_filtered)
  if (length(idx) != 0) {
    tx_coord_filtered <- tx_coord_filtered[-idx]
  }
  return(tx_coord_filtered)
}

normalize.fast <- function(sitesGRanges, guitarTxdb, txType, overlapIndex, siteLengthIndex) 
{
  sitesPointsPositionNormalize <- list()
  sitesPointsPosTx        <- list()
  sitesPointsPosTx        <- end(sitesGRanges[[txType]])
  names(sitesPointsPosTx) <- seqnames(sitesGRanges[[txType]])
  startPointMat           <- guitarTxdb[[txType]]$startPoint[names(sitesPointsPosTx), ]
  startPointDiffer        <- startPointMat - sitesPointsPosTx
  
  max_which.fast<-function(x){
    x[x>=0] <- -Inf # assign -Inf to prevent from max.col() catch them
    setNames(max.col(x),rownames(x))
  }
  sitesPointsComponet <- max_which.fast(startPointDiffer)
  sitesPointsPositionComponet <- startPointDiffer[cbind(seq_along(sitesPointsComponet), sitesPointsComponet)] * -1
  
  # normalize
  sitesPointsComponetWidthAvg  <- guitarTxdb[[txType]]$componentWidthAverage_pct[sitesPointsComponet]
  sitesPointsComponetStart_pct <- guitarTxdb[[txType]]$componentStartAverage_pct[sitesPointsComponet]
  componentWidthMat            <- guitarTxdb[[txType]]$componentWidth[names(sitesPointsPosTx),]
  sitesPointsComponetWidth     <- componentWidthMat[cbind(seq_along(sitesPointsComponet), sitesPointsComponet)]
  sitesPointsPositionNormalize <- sitesPointsPositionComponet/sitesPointsComponetWidth * sitesPointsComponetWidthAvg + sitesPointsComponetStart_pct
  names(sitesPointsPositionNormalize) <- sitesGRanges[[txType]]$xHits
  
  # weight
  sitesComponet_pct       <- guitarTxdb[[txType]]$componentWidthPtc[names(sitesPointsComponet), ]
  col.levels              <- sort(unique(sitesPointsComponet))
  sitesPointsComponet_pct <- setNames(numeric(length(sitesPointsComponet)),colnames(sitesComponet_pct)[sitesPointsComponet])
  for (i in col.levels){
    sitesPointsComponet_pct[sitesPointsComponet == i] <- sitesComponet_pct[sitesPointsComponet == i,i]
  }
  
  sitesPointsWeight        <- sitesPointsComponetWidthAvg/(sitesGRanges[[txType]]$pointsOverlapTx^overlapIndex)/sitesPointsComponet_pct * (sitesGRanges[[txType]]$sitesLength^siteLengthIndex)
  names(sitesPointsWeight) <- sitesGRanges[[txType]]$xHits
  return(list(sitesPointsPositionNormalize, sitesPointsWeight))
}

#### Plot function####
getColors <- function(n)
{
  if (requireNamespace('scales')){
    scales::hue_pal()(n)
  }
}

Density.plot <- function(dt,colorset=NULL,ncol=NULL,facet.lable.show=FALSE){
  if(is.null(colorset)){
    if(is.factor(dt$ID)){
      colorset <- setNames(getColors(length(levels(dt$ID))),levels(dt$ID))
    }else{
      colorset <- setNames(getColors(length(unique(dt$ID))),unique(dt$ID))
    }
  }else{
    colorset <- colorset[names(colorset) %in% unique(dt$ID)]
  }
  p<-ggplot()+
    geom_area(data = dt,alpha=0.5,position = "identity",size =1,aes(x=x,y=density,fill=ID,color=ID)) +
    scale_fill_manual(name = "ID",values = colorset) +
    scale_colour_manual(name = "ID",values = colorset) +
    scale_y_continuous(expand = expansion(mult=c(0,0.05)))+
    labs(x=NULL)+
    facet_wrap(~group,ncol = ncol)+
    theme_half_open()+
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          legend.title = element_blank(),
          legend.position = 'top',
          legend.justification = "center")
  if(!facet.lable.show){
    p<-p+theme(strip.background = element_blank(),
            strip.text = element_blank())
  }
  return(p)
}

MakeModelStructure <- function(feature.width,height,label=c("","5'UTR","CDS","3'UTR","")){
  if(length(feature.width)!=5){
    stop("feature.width should have 5 values!")
  }
  if( !is.null(label) && length(label)!=5){
    stop("label should have 5 values!")
  }
  feature.width     <- feature.width/sum(feature.width)
  model.height      <- height*1.06
  model.height.diff <- model.height-height
  Feature.model <- data.table(
    comp  = c("upstream","5'UTR","CDS","3'UTR","downstream"),
    label = label,
    width = feature.width,
    xmin  = c(0,head(cumsum(feature.width),-1)) + 0.001,
    xmax  = cumsum(feature.width),
    ymin  = model.height,
    ymax  = model.height,
    label_pos = model.height
  )[,`:=`(mid = (xmax+xmin)/2,
          label_pos = label_pos + 12/10 * model.height.diff)
  ][grep("UTR",comp),`:=`(ymin=ymin-3/10*model.height.diff,ymax=ymax+3/10*model.height.diff)
  ][grep("CDS",comp),`:=`(ymin=ymin-4/10*model.height.diff,ymax=ymax+4/10*model.height.diff)
  ][width!=0
  ][,bottom:=min(ymin)*0.975][]
  return(Feature.model)
}

add.mRNA.model <- function(feature.width,height,label=c("","5'UTR","CDS","3'UTR","")){
  componentStructure <- MakeModelStructure(feature.width = feature.width, height = height, label = label)
  list(geom_segment(data = componentStructure[grep("upstream|downstream",comp),],aes(x=xmin,xend=xmax,y=ymin,yend=ymax),color='black'),
       geom_rect(data = componentStructure[grep("UTR",comp),],aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill='grey70'),
       geom_rect(data = componentStructure[grep("CDS",comp),],aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill='grey50'),
       geom_text(data = componentStructure,aes(x=mid,y=label_pos,label=label)),
       geom_segment(data = head(componentStructure,-1),aes(x=xmax,xend=xmax,y=0,yend=bottom),linetype='dashed',color="grey")
  )
}

### Other function ####
# Subset representative models from GTF 
makeRepresentativeModelsGTF <- function(input,output,RepresentativeForms){
  Model     <-import.gff(con = input)
  Model.rep <-Model[mcols(Model)$transcript_id %in% RepresentativeForms]
  export.gff(Model.rep,con = output,format="gtf")
}

# Alignments To 5P ends
bedAlignemntTo5Pend <- function(file,MAPQ=NULL){
  # input file was converted from bam format via bedtools bamtobed
  # Unique hit is assigned 255 MAPQ if STAR aligner is used to map.
  for (file.i in unlist(file)){
    message("input:",file.i)
    NewFile <- sub(".bed.gz$",".5P.bed.gz",file.i)
    if (!file.exists(NewFile) ){
      Ailgnments <- fread(file.i,sep="\t",header = FALSE)
      if( !is.null(MAPQ)) {
        Ailgnments<-Ailgnments[V5>=MAPQ] # filter by MAPQ
      }
      Ailgnments[V6=="+",V3:=V2+1]
      Ailgnments[V6=="-",V2:=V3-1]
      fwrite(Ailgnments,file=NewFile,sep="\t",col.names=F)
      message("output:",NewFile) 
    }else{
      message(paste("Skip!","output file",NewFile, "is already existed.",sep=" ")) 
    }
  }
}
###### Main ########
setwd(dir = "R:/Qsync/XRNs/")
# TAIR10_GFF3_genes_transposons.gtf is made by gffread (http://ccb.jhu.edu/software/stringtie/gff.shtml#gffread)
# gffread -o TAIR10_GFF3_genes_transposons.gtf -T TAIR10_GFF3_genes_transposons.gff
GTF.file          = "./dataset/TAIR10_GFF3_genes_transposons.gtf"
GTF.rep.file      = "./dataset/TAIR10_representative_gene_models.gtf"
Rep_model_file    = "./dataset/TAIR10_representative_gene_models.txt"
txdb_file         = "./dataset/TAIR10_representative_gene_models.sqlite"
txGuitarTxdb.file = sub(".sqlite",".guitarTxdb",txdb_file)
# Subset representative models from Chr1 to Chr5
if (!file.exists(txdb_file)){
  Representative.forms <- fread(file = Rep_model_file,
                                sep="\t",
                                header = FALSE,
                                col.names = "ID",
                                skip = "AT")[grep("AT\\d+",ID)]
  makeRepresentativeModelsGTF(input = GTF.file,output = GTF.rep.file,RepresentativeForms = Representative.forms$ID)
  saveDb(makeTxDbFromGFF(file = GTF.rep.file,format = "gtf"),file = txdb_file)
}

# Build txGuitarTxdb
if( !file.exists(txGuitarTxdb.file) ){
  txdb <- loadDb(txdb_file)
  txGuitarTxdb <- makeGuitarTxdb(txdb,
                                 txpromoterLength = 100,
                                 txtailLength = 100,
                                 txMrnaComponentProp = c(1,1,1,1,1))
  saveRDS(txGuitarTxdb,txGuitarTxdb.file)
  rm(txGuitarTxdb)
}
# Alignments (bgzipped BED files) are converted into 5Pends and filtered by MAPQ
PARE.AlignemntFiles    <- list.files(pattern = "_MaxP.bed.gz", path = "./dataset/Bed", full.names = T)
bedAlignemntTo5Pend(PARE.AlignemntFiles,MAPQ = 255)
bed.Files              <- list.files(pattern = "_MaxP.5P.bed.gz", path = "./dataset/Bed", full.names = T)
names(bed.Files)       <- sub("_MaxP.5P.bed.gz","",basename(bed.Files))
bed.Files              <- as.list(bed.Files)

# Calculate the density of the bed files separately to prevent large ram usage
density.dir="PARE-density"
dir.create(density.dir, showWarnings = FALSE)
txMrnaComponentProp = c(1,2,4,2,1)
# record Component ratio to density output dir
fwrite(file = sprintf("./%s/ComponentProp.txt",density.dir),
       data.frame(txMrnaComponentProp=txMrnaComponentProp))

for ( ID in names(bed.Files) ){
    count <- GuitarPlotFast(txGuitarTxdb = txGuitarTxdb.file,
                            stBedFiles = bed.Files[ID],
                            txpromoterLength = 100,
                            txtailLength = 100,
                            headOrtail = TRUE,
                            enableCI = FALSE,
                            stSampleNum = 1,
                            mapFilterTranscript = TRUE,
                            pltTxType = c("mrna"),
                            txMrnaComponentProp = txMrnaComponentProp,
                            stGroupName = ID)
  fwrite(count$mrna$data,file=sprintf("./%s/%s.density",density.dir,ID),sep="\t")
}

# plot #
# Fetch and merge *.density files to make plot.
# density.dir = "PARE-density"
txMrnaComponentProp <- fread(file = sprintf("./%s/ComponentProp.txt",density.dir))$txMrnaComponentProp
DensityFiles  <- list.files(path = density.dir,pattern = "density",full.names = TRUE)
Density.dt    <- rbindlist(lapply(DensityFiles,fread,sep="\t",header=TRUE)
                           )[,`:=`(ID = factor(group,
                                               labels = c("WT","xrn2-1","xrn3-3","xrn4-6","fry1-6"),
                                               levels = c("WT","xrn2","xrn3","xrn4","fry1")),group=NULL)][]

ID.colors=setNames(c("#FF4040","#1874CD","#EEEE00","#00C5CD","#551A8B"),levels(Density.dt$ID))

Figure1E<-Density.plot(dt = rbind(Density.dt[ID %in% c("WT","xrn4-6","xrn2-1")][,group:="A"],
                                 Density.dt[ID %in% c("WT","xrn3-3","fry1-6")][,group:="B"]),
                      colorset = ID.colors,ncol = 2)+
  add.mRNA.model(feature.width = txMrnaComponentProp,
                 height = round(max(Density.dt$density)),
                 label=c("","5'UTR","CDS","3'UTR",""))+
  theme(legend.key.size = unit(1.5,"line"))
