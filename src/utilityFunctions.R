#   ##### Utility Functions ######
#HC: function to rename/reorder the columns for rbind of ssGSEAScore

createRNAseqProject <- function(project_name, meta_file, factorsToExclude) {
  prj <- ProjectSetUp$new(
    date                    = unlist(strsplit(x = as.character(Sys.time()), "\\s+"))[[1]],
    time                    = unlist(strsplit(x = as.character(Sys.time()), "\\s+"))[[2]],
    projectName             = project_name,
    annotationRDS           = "input/AnnotationRDS/annotation_ENSEMBL_gene.RDS",
    pcRDS                   = "input/AnnotationRDS/pc.other.HGNCTableFlat.rds",
    emRDS                   = "input/AnnotationRDS/EMGenes.RDS",
    tfRDS                   = "input/AnnotationRDS/TFs_no_epimachines.RDS",
    csRDS                   = "input/AnnotationRDS/CellSurface.RDS",
    cgaRDS                  = "input/AnnotationRDS/cancerGermlineAntigens.rds",
    ewsr1Fli1RDS            = "input/AnnotationRDS/EWSR1_FL1_DownstreamTargets.RDS",
    pax3Foxo1RDS             ="input/AnnotationRDS/PAX3_FOXO1_DownstreamTargets.RDS",
    
    BrainExpRDS             = "input/AnnotationRDS/VitalExpression/expressionTMM.RPKM.Brain.v2.RDS",
    HeartExpRDS             = "input/AnnotationRDS/VitalExpression/expressionTMM.RPKM.Heart.v2.RDS", 
    KidneyExpRDS            = "input/AnnotationRDS/VitalExpression/expressionTMM.RPKM.Kidney.v2.RDS", 
    LiverExpRDS             = "input/AnnotationRDS/VitalExpression/expressionTMM.RPKM.Liver.v2.RDS", 
    LungExpRDS              = "input/AnnotationRDS/VitalExpression/expressionTMM.RPKM.Lung.v2.RDS", 
    
    outputPrefix            = project_name,
    filterGenes             = TRUE,
    filterGeneMethod        = "bySum",
    factorName              = "DIAGNOSIS.Substatus.Tumor.Normal.Tissue",
    metadataFileRefCol      = "Sample.Biowulf.ID.GeneExp",
    metaDataFileName        = meta_file,
    outputdirRDSDir         = "output/GeneRDSOutput",
    outputdirTXTDir         = "output/GeneTXTOutput",
    gseaDir                 = "input/GSEA",
    plotsDir                = "output/Figures",
    plotsDataDir            = "output/FigureData",
    DiffGeneExpAnaDir       = "output/DiffExpResults",
    DiffGeneExpRDS          = "output/DiffGeneExpRDSOutput",
    factorsToExclude        = factorsToExclude
  )
  return(prj)
}

reOrderForssGSEA <- function(df, meta, ssGSEADF, dfCol, refCol) {
  ncol <- ncol(df)
  df$rn <- rownames(df)
  df$rn <- gsub("X0", "0", df$rn)
  df <- df %>% dplyr::inner_join(meta, by=c("rn"=dfCol))
  df <- df[,c(refCol, colnames(df)[1:ncol])]
  df <- df[match(colnames(ssGSEADF), df[,1]),]
  rownames(df) <- df[,1]
  rn <- df[,1]
  df <- df[,2:ncol(df)]
  df <- t(df)
  colnames(df) <- rn
  reOrderForssGSEA <- df
}

median.quartile <- function(x){
  out <- quantile(x, probs = c(0.5))
  names(out) <- c("y")
  return(out) 
}


## Flatten rcorr results in to a table
flattenCorrMatrix <- function(cormat, pmat, group_name="") {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut],
    group = group_name
  )
}

## Do Correlation for each diagnosis groups
rcorr_groups = function(annotDF = NULL, dataDF = NULL, GroupID =NULL, samplesID= NULL) {
  unique_groups <- unique(as.character(annotDF[,GroupID]))
  mergedDF <- do.call(rbind, lapply(unique_groups, function(x) {
    each_group_samples <- annotDF %>% dplyr::filter_(.dots =paste0("grepl(\"",x,"\",",GroupID,")") )
    each_group_samplesDF <- dataDF %>% dplyr::select(one_of(each_group_samples[,samplesID]))
    print(paste("group name", x,"  ", dim(each_group_samplesDF)[1], " ",dim(each_group_samplesDF)[2] ))
    each_group_samplesDF.corr <- rcorr(t(each_group_samplesDF), type = "spearman");
    return(flattenCorrMatrix(each_group_samplesDF.corr$r, each_group_samplesDF.corr$P, x))
  } ))
  return(mergedDF)
}

## coorelation plots
correlationPlots <- function(varName="", constName="", df=NA, customColorDF=NA, xlab="Log Total Clones"){
  
  print(paste(varName))
  df <- df %>% dplyr::filter_(.dots = paste0(varName, ">=", 5))
  customColorsVector <- setNames( as.character(customColorDF$Color), as.character(customColorDF$Diagnosis))
  corrTest <- cor.test(df[,constName], df[,varName], method = "spearman")
  if  ( corrTest$p.value < 2.2e-16 ) { corrTest$p.value = 2.2e-16 }
  plot <- ggplot(df, aes_string(x=constName, y=varName)) + 
    geom_smooth(method=lm,  fill="grey") +
    geom_point(aes(colour = factor(Diagnosis)), show.legend = T, size=2, shape=16) + 
    scale_colour_manual(values=customColorsVector) +
    theme_bw() +
    theme(axis.text=element_text(size=13)
          ,axis.title=element_text(size=13,face="bold")) +
    xlab(xlab) +
    ylab(paste("Erv expression", sep=" "))+
    # theme(legend.position = "none") +
    ggtitle(paste("Corr.Coeff = ", signif(corrTest$estimate[[1]],5), "\np-value = ", signif(corrTest$p.value,5), "", varName,sep=""))
  
  return(list(plot))
} 

ExonBedToExonGTF <- function(bedFile){
  bedFileExonDups <- bedFile %>% filter(V4 %in% bedFile$V4[which(duplicated(as.character(bedFile$V4)))])
  
  bedFile <- data.table(bedFile)
  bedFile[ , ExonIndex := 1:.N, by = V4 ]
  splitExonName <- strsplit(bedFile$V4, split="___")
  gene_id       <- sapply(splitExonName, function(x) { return(x[1]) } )
  exon_number   <- sapply(splitExonName, function(x) { return( paste(sort(unique(unlist(strsplit(x[2], ",")) )), collapse=",") ) } )
  transcript_id <- sapply(bedFile$V5, function(x) { return( paste(unique(unlist(strsplit(x, ",")) ), collapse=",") ) } )
  exon_id       <- paste(gene_id,"___",exon_number, ".", bedFile$ExonIndex, sep="")
  gene_name <- gene_id
  
  bedFileFinal <- data.frame(matrix(ncol = 9, nrow = nrow(bedFile)))
  colnames(bedFileFinal) <- c("Chr", "Anotation", "Feature", "Start", "End", "Strand1", "Strand2", "Strand3", "Info")
  bedFileFinal$Chr <- bedFile$V1
  bedFileFinal$Anotation <- rep("SeqCap_EZ_Panel", nrow(bedFile))
  bedFileFinal$Feature <- rep("exon", nrow(bedFile))
  bedFileFinal$Start <- bedFile$V2
  bedFileFinal$End <- bedFile$V3
  bedFileFinal$Strand1 <- rep(".", nrow(bedFile))
  bedFileFinal$Strand2 <- rep(".", nrow(bedFile))
  bedFileFinal$Strand3 <- rep(".", nrow(bedFile))
  bedFileFinal$Info <-   paste(
    
    paste("gene_id ","\"",gene_id,"\"",sep=""),
    paste("transcript_id ","\"",transcript_id,"\"",sep=""),
    paste("exon_number ","\"",exon_number,"\"",sep=""),
    paste("exon_id ","\"",exon_id,"\"",sep=""),
    paste("gene_name ","\"",gene_name,"\"",sep=""),
    sep=";"
  )
  write.table(bedFileFinal, "SeqCap_EZ_Panel_Neuro.target.hg19.merged.gtf", sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

mapableExonBases <- function(df){
  firstSeq = seq(1:nrow(df))
  GeneLength <- length ( unique( unlist ( lapply ( firstSeq, 
                                                   function(x, df) 
                                                   {
                                                     return(c(df[x,"Start"]:df[x,"End"]))
                                                   }, df=df 
  ) 
  )
  ) 
  )
  return(GeneLength)
}

makeAnnotationDFs <- function(AnnotationName="", Reffile=""){
  
  annotation_UCSC.1 <- read.csv(paste("C:/Users/sindiris/R Scribble/GTF/",Reffile, sep=""), sep="\t", header=FALSE, stringsAsFactors = FALSE)
  desc            <- lapply(annotation_UCSC.1$V9, function(x){ return(unlist(strsplit(as.character(x), split=";"))) } )
  annotation_UCSC.2 <- cbind(annotation_UCSC.1, t(as.data.frame(matrix(unlist(desc), nrow=length(unlist(desc[1]))))))
  annotation_UCSC.3 <- annotation_UCSC.2[order(annotation_UCSC.2[,1], annotation_UCSC.2[,10], annotation_UCSC.2[,12]), ]
  annotation_UCSC.4 <- subset(annotation_UCSC.3, annotation_UCSC.3[,3]=="exon")
  
  annotation_UCSC.Filt.5 <- annotation_UCSC.4[,c(1,4,5,7,10:14)]
  annotation_UCSC.Filt.6 <- apply(annotation_UCSC.Filt.5, 2, function(x) { return( gsub( "\\s+", "", x, perl=TRUE ) )} )
  annotation_UCSC.Filt.7 <- as.data.frame( apply(annotation_UCSC.Filt.6,2, function(x) { return( gsub("gene_id|transcript_id|exon_number|exon_id|gene_name","",x) ) } ))
  colnames(annotation_UCSC.Filt.7) <- c("Chr", "Start", "End", "Strand", "GeneID", "TranscriptID", "ExonNumb", "ExonID", "GeneName")
  annotation_UCSC.Filt.7$Start <- as.numeric(as.character(annotation_UCSC.Filt.7$Start))
  annotation_UCSC.Filt.7$End   <- as.numeric(as.character(annotation_UCSC.Filt.7$End))
  
  annotation_UCSC.exon <- ungroup(annotation_UCSC.Filt.7 %>% dplyr::mutate(ExonLength=End-Start+1 ) %>%
                                    dplyr::group_by(TranscriptID) %>% 
                                    dplyr::mutate(TransLength=sum(ExonLength) ) %>%
                                    dplyr::select_(.dots = c("Chr","Start","End","Strand", "GeneID", "GeneName", "TranscriptID", "TransLength","ExonID", "ExonLength"))) %>%
    dplyr::rename(Length=ExonLength)
  saveRDS(annotation_UCSC.exon, paste("annotation_",AnnotationName,"_exon.RDS", sep=""))
  
  annotation_UCSC.transcript <- annotation_UCSC.exon %>% dplyr::group_by(TranscriptID) %>% 
    dplyr::mutate(StartMin= min(Start), EndMax =max(End), Length=sum(Length) ) %>% 
    dplyr::select_(.dots = c("Chr","StartMin","EndMax","Strand", "GeneID", "GeneName", "TranscriptID", "Length")) %>% 
    dplyr::rename(Start=StartMin, End=EndMax) %>% 
    distinct(Chr, Start, End, Strand, GeneID, GeneName, Length)
  saveRDS(annotation_UCSC.transcript, paste("annotation_",AnnotationName,"_transcript.RDS", sep=""))
  
  annotation_UCSC.gene <- annotation_UCSC.exon %>%  dplyr::group_by(GeneID) %>% 
    dplyr::mutate(StartMin= min(Start), EndMax = max(End), GeneLength=mapableExonBases(data.frame(Start=Start, End=End)) ) %>%
    dplyr::select_(.dots = c("Chr","StartMin","EndMax","Strand", "GeneID", "GeneName", "GeneLength")) %>% 
    dplyr::rename(Start=StartMin, End=EndMax, Length=GeneLength) %>%
    distinct(Chr, Start, End, Strand, GeneID, GeneName, Length) %>%
    na.omit
  dupsGenes <- which(duplicated(annotation_UCSC.gene$GeneID) ==TRUE)
  if(AnnotationName=="UCSC")  annotation_UCSC.gene <- annotation_UCSC.gene[-dupsGenes,]
  saveRDS(annotation_UCSC.gene, paste("annotation_",AnnotationName,"_gene.RDS", sep=""))
}  

mergeMultipleLists <- function(groups) { as.character( unlist(sapply(groups, function(x){ 
  fileVectors <- list.files(x)
  fileList <- paste(x,"/",fileVectors, sep="")
  return( fileList )     })) ) 
}

zscore_All <- function(x){
  medX <- median(x)
  sdX <- sd(x, na.rm = FALSE)
  y <- (x - medX) / sdX
  return(y)
}

median_ALL <- function(x){
  medX <- median(x)
  y <- (x - medX)
  return(y)
}

zscore_WRT.Normal <- function(x){
  medX <- median(x[1:42])
  sdX <- sd(x[1:42], na.rm = FALSE)
  y <- (x - medX) / sdX
  return(y)
}

ttestFunc <- function(x, groupAvect="" , groupBvect=""){
  t_stat  <- t.test( x[groupAvect], x[groupBvect], alternative = c("two.sided") )
  outcome_list <- c("MedExpGroupA"=median(x[groupAvect]),  "MedExpGroupB"=median(x[groupBvect]), "TStat"=t_stat$statistic[["t"]], 
                    t_stat$p.value)
  return(outcome_list)
}

createBroadGCTFile <- function(x){
  RPKM_Data_Filt <- x
  metaDes                   <- matrix("", nrow = nrow(RPKM_Data_Filt), ncol = 2)  ; colnames(metaDes) <- c("Genes", "Description"); metaDes[,1] <- rownames(RPKM_Data_Filt)
  RPKM_Data_Filt_Meta       <- cbind(metaDes, RPKM_Data_Filt) 
  RPKM_Data_Filt_Meta       <- data.frame(lapply(RPKM_Data_Filt_Meta, as.character), stringsAsFactors=FALSE)
  metaSS                    <- matrix("", nrow = 3, ncol = ncol(RPKM_Data_Filt_Meta)); 
  metaSS[1,1]               <- "#1.2"; 
  metaSS[2,c(2,3)]          <- c(nrow(RPKM_Data_Filt_Meta), ncol(RPKM_Data_Filt_Meta)-2) ;
  metaSS[3,]                <- colnames(RPKM_Data_Filt_Meta)
  colnames( metaSS )        <- colnames(RPKM_Data_Filt_Meta)
  RPKM.Data.Filt.Meta.Broad <- rbind(as.data.frame(metaSS),RPKM_Data_Filt_Meta)
  return(RPKM.Data.Filt.Meta.Broad)
}

orderRowsOfDf <- function(targetOrder="", queryOrder="", df=""){
  matchRowOrder <- match( targetOrder, queryOrder )
  dfOrdered <- df[matchRowOrder,]
  return(dfOrdered)
}

orderColsOfDf <- function(targetOrder="", df=""){
  matchRowOrder <- match( targetOrder, queryOrder )
  dfOrdered <- df[,matchRowOrder]
  return(dfOrdered)
  
}

consolidateDF <- function(df, featureName="", funcName=""){
  dfConsolidated <- as.data.frame( df %>% dplyr::group_by_(.dots=c(featureName)) %>% dplyr::select(-c(1,2,3,5)) %>% 
                                     summarise_all(funs(!!funcName)) ) %>% data.frame() %>% 
    tibble::column_to_rownames( var=featureName )
  #%>% tibble::rownames_to_column()
  #%>% dplyr::select( -matches(featureName) )
  return(dfConsolidated)
}

GSVAFunction <- function(expressionFile, geneSet=geneSet, method="gsva", ssgsea.norm=FALSE, rnaseq=FALSE, annotationType=NULL){
  
  cytolyticDF       <- expressionFile[c("GZMA","GZMB","GZMH","GZMK", "GZMM", "PRF1"),] ; 
  
  Cytolytic <- apply(cytolyticDF,2,cytolyticScoreFun)
  gsva_es          <- gsva(as.matrix(expressionFile), geneSet, mx.diff=TRUE, method=method, ssgsea.norm=ssgsea.norm, rnaseq=rnaseq)
  
  if ( method=="gsva") { 
    gsva_es           <- rbind(gsva_es$es.obs,Cytolytic)
  } else { 
    gsva_es           <- rbind(gsva_es, Cytolytic) 
  }
  
  return(gsva_es)
}

GeneratePreRankGSEAcmd <- function(x, projectName=NULL, sample=NULL){
  label <- gsub(".gmt","",paste(sample,".",x,sep=""))
  cmd <- paste("java -cp /data/sindiris/Processing/GSEA/gsea2-2.2.2.jar -Xmx1024m xtools.gsea.GseaPreranked",
               " -gmx /data/sindiris/Processing/GSEA/",x,
               " -rnk /data/sindiris/Processing/GSEA/",projectName,"/rnk/",sample,
               " -chip /data/sindiris/Processing/GSEA/GENE_SYMBOL.chip ",
               " -collapse false",
               " -mode Max_probe",
               " -norm meandiv",
               " -nperm 1000",
               " -scoring_scheme weighted",
               " -rpt_label ",label,
               " -include_only_symbols  true",
               " -make_sets true",
               " -plot_top_x 20",
               " -rnd_seed timestamp",
               " -set_max 500",
               " -set_min 15",
               " -out /data/sindiris/Processing/GSEA/",projectName,
               " -gui false",
               sep="")
  return(cmd)
}

DESeq.getOffset <- function(y) {
  if (is.na(sizeFactors(y))) { stop("Call estimateSizeFactors first") }
  log(sizeFactors(y)) - mean(log(sizeFactors(y))) + mean(log(colSums(counts(y))))
}

hgncAnotate<- function(genelist="", columnName="gene_family", all=FALSE, outColNames="gene_family"){
  hgnc <- read.csv("T:/Sivasish_Sindiri/R_workspace/Cancer.Testes.Antigen/HGNC_GeneAnnotation.txt", sep = "\t", header = T,
                   stringsAsFactors = FALSE)
  if(!all){
    testResultDF<- do.call(rbind, lapply(genelist, function(x){ test <- hgnc %>% filter(symbol==x) %>% dplyr::select(one_of(columnName)) 
    if(nrow(test)==0){test<-"No desription found"}
    return(test)}))
  } else{
    testResultDF<- do.call(rbind, lapply(genelist, function(x){ test <- hgnc %>% filter(symbol==x) 
    if(nrow(test)==0){test<-"No desription found"}
    return(test)}))
  }
  
  colnames(testResultDF) <- outColNames
  return(testResultDF)
}

specify_decimal <- function(x, k) format(round(x, k), nsmall=k)

rangeZeroToOne <- function(X, exponential=FALSE, multFact=1){
  Y <- (X - min(X))/diff(range(X))
  if(exponential) {
    Y <- exp(Y*multFact)
  }
  return(Y)
}

correlation <- function(x,y,pval=FALSE){
  print(paste(length(as.numeric(unlist(x)))))
  print(paste(length(as.numeric(unlist(y)))))
  corrTest <- cor.test(as.numeric(unlist(x)), as.numeric(unlist(y)), method = "pearson",exact = FALSE)
  if(pval==TRUE){
    return(corrTest$p.value)
  } else {
    return(corrTest$estimate[["cor"]])
  }
  #return(cor(x, y))
}

cytolyticScoreFun <- function(x){
  x <- sum(log2(x+1))/6
  return(x)
}

factorizeColumn <- function( toFactor, asFactor){
  factorColumn <- factor(toFactor, levels=asFactor, ordered = TRUE )
  return(factorColumn)
}

getGeneNamesfromID <- function(annotation=NULL, GeneIdVector=NULL) {
  return( annotation %>% data.frame() %>% dplyr::filter(GeneID %in% GeneIdVector) %>% arrange(GeneID) %>% dplyr::select(matches("^GeneName$")) )
}

getGeneIDFromNames <- function(annotation=NULL, GeneNameVector=NULL) {
  return( annotation %>% data.frame() %>% dplyr::filter(GeneName %in% GeneNameVector) %>% arrange(GeneID) %>% dplyr::select(matches("^GeneID$")) )
}

getMetaDataStats <- function(metaDataDF=NULL, factorName=NULL, libType=NULL, CountCol=NULL, Sum=NULL, finalColName=NULL) {
  
  statDF <-  metaDataDF %>% group_by_(.dots= c(factorName, "Color", libType) ) %>% count_(var=as.name(CountCol)) %>% 
    dplyr::summarise(Count=n()) %>% 
    dplyr::group_by_(.dots= c(factorName)) %>% 
    dplyr::mutate( !!Sum := sum(Count)) %>% 
    spread_(libType, "Count") %>% 
    mutate_( .dots = setNames( list( interp(~paste(factorName ,"(", Sum , ")"), factorName=as.name(factorName), Sum=as.name(Sum) ) ), finalColName) ) %>% 
    data.frame() %>% mutate_at(c(3:6), funs(replace(., is.na(.), 0)))
  
  return(statDF)
}

getOncoprintOrder <- function(mat){
  # convert mat to mat_list
  get_type = function(mat) strsplit(mat, ";")[[1]]
  if(inherits(mat, "data.frame")) {
    mat = as.matrix(mat)
  }
  if(inherits(mat, "matrix")) {
    all_type = unique(unlist(lapply(mat, get_type)))
    all_type = all_type[!is.na(all_type)]
    all_type = all_type[grepl("\\S", all_type)]
    
    mat_list = lapply(all_type, function(type) {
      m = sapply(mat, function(x) type %in% get_type(x))
      dim(m) = dim(mat)
      dimnames(m) = dimnames(mat)
      m
    })
  }
  names(mat_list) = all_type
  
  arr = array(FALSE, dim = c(dim(mat_list[[1]]), length(all_type)), dimnames = c(dimnames(mat_list[[1]]), list(all_type)))
  for(i in seq_along(all_type)) {
    arr[, , i] = mat_list[[i]]
  }
  count_matrix = apply(arr, c(1, 2), sum)
  oncoprint_row_order = function() {
    order(rowSums(count_matrix), decreasing = TRUE)
  }
  
  oncoprint_column_order = function() {
    scoreCol = function(x) {
      score = 0
      for(i in 1:length(x)) {
        if(x[i]) {
          score = score + 2^(length(x)-i*1/x[i])
        }
      }
      return(score)
    }
    scores = apply(count_matrix[oncoprint_row_order(), ,drop = FALSE], 2, scoreCol)
    order(scores, decreasing=TRUE)
  }
  
  return(list(rowOder =oncoprint_row_order(), colOder=oncoprint_column_order() ))
}

pairedTesting <- function(df, ColsNamesToIgnore=NA, ColsNumbersToIgnore=NA, funcName=NA, corrMethod="pearson") {
  
  if(length(ColsNamesToIgnore) != 0 ) df %<>% dplyr::select_(.dots = setdiff(colnames(df),ColsNamesToIgnore))
  if(length(ColsNamesToIgnore) != 0 ) df %<>% dplyr::select_(.dots = setdiff(seq(1:ncol(df)),ColsNumbersToIgnore))
  
  ### Prepare combinations of test
  combos <- combn(ncol(df),2)
  
  ### Perform t.test & correraltion
  if(funcName=="t.test"){
    result <- do.call(rbind, apply(combos, 2, function(x) {
      Ttest <- t.test(df[, x[1]], df[, x[2]])
      CorrTest <- corr.test(df[, x[1], drop=FALSE], df[, x[2], drop=FALSE], method = corrMethod)
      out <- data.frame("var1" = colnames(df)[x[1]]
                        , "var2" = colnames(df[x[2]])
                        , "t.value" = sprintf("%.3f", Ttest$statistic)
                        , "ttest.p.value" = sprintf("%.3f", Ttest$p.value)
                        , "corr.Coefficient" = sprintf("%.3f", CorrTest$r[[1]]) 
                        , "corr.p.value" = sprintf("%.3f", CorrTest$p[[1]])
      )
      return(out)
    }))
    ttest.pvalue <- as.numeric(as.character(result$ttest.p.value)); corr.pvalue <- as.numeric(as.character(result$corr.p.value))
    result$holm.adjusted.ttest.p.val <- p.adjust(ttest.pvalue, method = "holm", n = length(ttest.pvalue))
    result$holm.adjusted.corr.p.val <- p.adjust(corr.pvalue, method = "holm", n = length(corr.pvalue))
  }
  return(result)
}

