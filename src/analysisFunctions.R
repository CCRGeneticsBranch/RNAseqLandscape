#   ##### Analysis Processing Functions ######

getCountObjRDS <- function(x){
  print(paste(x))
  featureCountRDS <- readRDS(x)
  return(featureCountRDS$count)
}

getCountObjTXT <- function(fileName, colNumb=1, rowNames=1){
  print(paste(fileName))
  featureCountTxt <- read.table(fileName, sep="\t", row.names = rowNames, header = 1);
  return(featureCountTxt[,colNumb, drop=FALSE])
}

makeFinalObj <- function(countObj=countObj, genesObj=genesObj, Geneslist="" ){
  targets <- gsub(" ", "", as.character(rownames(countObj)))
  annotation_Exon <- genesObj[match(targets, as.character(genesObj[,4])),]
  
  if(length(Geneslist) > 0 )
  {
    QueryGeneObj <- annotation_Exon[which(annotation_Exon$GeneName %in% Geneslist),]
  }
  else{
    QueryGeneObj <- annotation_Exon
  }
  
  Query_countObj <- countObj[rownames(QueryGeneObj), ]
  finalObj <- cbind(QueryGeneObj, Query_countObj)
  return(finalObj)
}

mergeCountsByObject <- function(x, type=""){
  print(paste(x))
  file_Dir_Gene = x
  fileName <- basename(file_Dir_Gene)
  GeneFiles             <- list.files(file_Dir_Gene) ; GeneFilesFormat <- gsub("",GeneFiles,pattern = paste("_counts.",type,".fc.RDS",sep=""))
  GeneFilesList         <- paste(file_Dir_Gene, "/", GeneFiles,sep="") ; length(GeneFilesList)
  
  fileList          <- GeneFilesList
  countObj          <- do.call(cbind,lapply(fileList,getCountObjRDS))
  
  saveRDS(countObj, paste(outputdirRDS, "/RawCount/", paste(fileName, type,"_RawCount.rds",sep=""), sep= ""))
  write.table(countObj, paste(outputdirTXT, "/RawCount/", paste(fileName, type,"_RawCount.txt",sep=""), sep= ""), sep="\t",row.names = TRUE, quote = FALSE)
}

getCountObjRSEM <- function(x, reqColName="", mergeColName=""){
  print(paste(x)) 
  
  sampleName <- gsub(".genes.results", "", basename(x))
  
  featureCountRDS <- fread(x, sep="\t", header = TRUE) %>% rename_(.dots=setNames(reqColName, sampleName))
  
  return(featureCountRDS[, mget(c("gene_id",sampleName))])
}

mergeCountsByObjectRSEM <- function(x, type="", reqColName="expected_count", mergeColName="gene_id"){
  print(paste(x))
  file_Dir_Gene = x
  fileName <- basename(file_Dir_Gene)
  GeneFiles         <- list.files(file_Dir_Gene) ; GeneFilesFormat <- gsub("",GeneFiles,pattern = paste(".",type,".results",sep=""))
  GeneFilesList     <- paste(file_Dir_Gene, "/", GeneFiles,sep="") ; length(GeneFilesList)
  
  fileList          <- GeneFilesList
  countObj          <- purrr::reduce(lapply(fileList,getCountObjRSEM, reqColName=reqColName, mergeColName=mergeColName), full_join, by=mergeColName)
  
  saveRDS(countObj, paste(outputdirRDS, "/RawCount/", paste(fileName, type,"_RawCount.rds",sep=""), sep= ""))
  write.table(countObj, paste(outputdirTXT, "/RawCount/", paste(fileName, type,"_RawCount.txt",sep=""), sep= ""), sep="\t",row.names = TRUE, quote = FALSE)
}

mergeCountsByObjectCCLE <- function(x, type="", extension="", colInterest=1){
  
  print(paste(x))
  file_Dir_Gene = x
  fileName <- basename(file_Dir_Gene)
  GeneFiles             <- list.files(file_Dir_Gene) ; GeneFilesFormat <- gsub("",GeneFiles,pattern = extension)
  GeneFilesList         <- paste(file_Dir_Gene, "/", GeneFiles,sep="") ; length(GeneFilesList)
  
  fileList          <- GeneFilesList
  countObj          <- do.call(cbind,lapply(fileList, getCountObjTXT, colNumb=colInterest, rowNames=1))
  colnames(countObj)<- GeneFilesFormat
  
  saveRDS(countObj, paste(outputdirRDS, "/RawCount/", paste(fileName, type,"_RawCount.rds",sep=""), sep= ""))
  write.table(countObj, paste(outputdirTXT, "/RawCount/", paste(fileName, type,"_RawCount.txt",sep=""), sep= ""), sep="\t",row.names = TRUE, quote = FALSE)
}

mergeDiffTestResults <- function(x, type="", saveDirPath="", extension="", colInterest=1){
  
  print(paste(x))
  file_Dir_Gene = x
  fileName <- basename(file_Dir_Gene) 
  dir.create(file.path(paste(saveDirPath,fileName,sep="/")))
  GeneFiles             <- list.files(file_Dir_Gene); GeneFiles <- GeneFiles[grep("*.txt", GeneFiles)]
  GeneFilesList         <- paste(file_Dir_Gene, "/", GeneFiles,sep="") ; length(GeneFilesList)
  
  countObj          <- do.call(cbind,lapply(GeneFilesList, getCountObjTXT, colNumb=colInterest, rowNames=1))
  
  write.table(countObj, paste(saveDirPath, paste(fileName, "/", type,"_MergedDiffExpResult.txt",sep=""), sep= "/"), sep="\t",row.names = TRUE, quote = FALSE)
}

tpm_default <- function(x){
  
  #Initialize Variables
  counts <- x$counts
  featureLength <- x$genes[,"Length"]
  meanFragmentLength <- rep(320, ncol(counts))
  
  # Ensure valid arguments.
  stopifnot(length(featureLength) == nrow(counts))
  stopifnot(length(meanFragmentLength) == ncol(counts))
  
  # Compute effective lengths of features in each library.
  effLen <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    featureLength - meanFragmentLength[i] + 1
  }))
  
  # Exclude genes with length less than the mean fragment length.
  idx <- apply(effLen, 1, function(x) min(x) > 1)
  counts <- counts[idx,]
  effLen <- effLen[idx,]
  featureLength <- featureLength[idx]
  
  # Process one column at a time.
  tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    rate = log(counts[,i]) - log(effLen[,i])
    denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))
  
  # Copy the row and column names from the original matrix.
  colnames(tpm) <- colnames(counts)
  rownames(tpm) <- rownames(counts)
  return(tpm)
}

fpkmToTpm <- function(fpkm){
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

keepProteinCoding <- function(geneMatrix=NULL, annotation=NULL, featureType="GeneID") {
  ###Keep only Protein Coding Genes
  PC <- read.table("C:/Users/sindiris/R Scribble/Annotation RDS/HGNC-protein-coding-List.txt", header = T, sep="\t")
  annotationPC <- annotation %>% filter(GeneName %in% PC$Genes)
  foundGenes <- which(row.names(geneMatrix) %in% annotationPC[,featureType] )
  geneMatrixNew <- geneMatrix %>% data.frame() %>% .[foundGenes,]
  annotationPC <- annotationPC %>% mutate(GeneID = factorizeColumn(annotationPC$GeneID, row.names(geneMatrixNew)))
  return(list(geneMatrixNew, annotationPC) )
}

normalize <- function(countObj, annotName = "", method="", fileName="", annotationRDS="", merge=T, saveFiles=F, condition="", proteinCodingOnly=TRUE){
  
  
  GeneDF_Norm_CPM.PC   <- data.frame()
  GeneDF_Norm_rpkm.PC  <- data.frame()
  GeneDF_Norm_tpm.PC   <- data.frame()
  
  Annotation <- data.frame(readRDS(annotationRDS))
  countObj   <- data.frame(countObj)
  
  if( annotName == "Exon")  { featureName = "ExonID"  }
  if( annotName == "Transcript"){ featureName = "TranscriptID" }
  if( annotName == "Gene"){ featureName = "GeneID" }
  
  rownames(Annotation) <- Annotation[,featureName]
  AnnotationNew <- Annotation %>% mutate(!!featureName := factorizeColumn(Annotation[[featureName]], row.names(countObj))) %>% 
    dplyr::arrange_(.dots=c(featureName))
  
  if (proteinCodingOnly) {
    countObjFinal <-  keepProteinCoding(geneMatrix = countObj, annotation=AnnotationNew, featureType=featureName)[[1]]
    AnnotationFinal<- keepProteinCoding(geneMatrix = countObj, annotation=AnnotationNew, featureType=featureName)[[2]] 
    AnnotationFinal <-  AnnotationFinal %>% drop_na()
    genesObj <- AnnotationFinal[,c(featureName, "Length")]
  } else {
    countObjFinal <-  countObj
    AnnotationFinal <- AnnotationNew
    genesObj <- AnnotationFinal[,c(featureName, "Length")]
  }
  
  if( annotName == "Exon") { 
    selectColumns <- c("Chr","Start","End","GeneName","TranscriptID","ExonID") 
  } else if(annotName == "Transcript"){ 
    selectColumns <- c("Chr","Start","End","GeneName","TranscriptID") 
  } else{ selectColumns <- c("Chr","Start","End","GeneName","GeneID") }
  
  ########################################### Choose the Method 
  print(paste("Chooseing the Normalization Method"))
  if( method=="EdgeR") {
    
    #### EdgeR
    ##  Make EdgeR Object
    colnames(genesObj) <- c("GeneID", "Length")
    GeneDF_EdgeR       <- DGEList(counts=countObjFinal, genes=genesObj, group = condition)
    ## Estimate Normalising Factors
    print(paste("Calculating scaling factor using calcNormFactors "))
    GeneDF.Norm  <- calcNormFactors(GeneDF_EdgeR) ; 
    print(paste("Estimating TPM, CPM and RPKM values"))
    ## Regularized Log Transformation using CPM, FPKM & TPM values
    GeneDF.CPM   <- as.data.frame(cpm(GeneDF.Norm,  normalized.lib.sizes = TRUE,log = FALSE)) 
    GeneDF.rpkm  <- as.data.frame(rpkm(GeneDF.Norm, normalized.lib.sizes = TRUE, log = FALSE))
    GeneDF.tpm   <- apply(rpkm(GeneDF.Norm, normalized.lib.sizes = TRUE), 2 , fpkmToTpm)
    #GeneDF.ScaledTpm <- t(t(GeneDF.tpm) / GeneDF.Norm$samples$norm.factors)
  }
  else if( method == "DESeq") {
    #### DESeq
    ##Make DESeq Object
    GeneDF_DESeq      <- DESeqDataSetFromMatrix(countData = as.matrix(countObjFinal), colData = DataFrame(condition), design = ~ condition)
    mcols(GeneDF_DESeq)$basepairs <- genesObj$Length-175+1
    ## Estimate SizeFactors
    GeneDF.Norm <- estimateSizeFactors(GeneDF_DESeq)
    ## Regularized Log Transformation using CPM, FPKM & TPM values
    GeneDF.fpm   <- as.data.frame(fpm(object = GeneDF.Norm,  robust = TRUE))
    GeneDF.rpkm  <- as.data.frame(fpkm(object = GeneDF.Norm,  robust = TRUE));
    GeneDF.tpm   <- apply(GeneDF.rpkm, 2 , fpkmToTpm)
    GeneDF.ScaledTpm <- t(t(GeneDF.tpm) / GeneDF.Norm$samples$norm.factors)
  }
  
  ########################################### Prepare final files
  print(paste("Prepare final files"))
  
  GeneDF_Norm_CPM  <- cbind(data.frame(AnnotationFinal[,selectColumns]), GeneDF.CPM )
  GeneDF_Norm_rpkm <- cbind(data.frame(AnnotationFinal[,selectColumns]), GeneDF.rpkm )
  GeneDF_Norm_tpm  <- cbind(data.frame(AnnotationFinal[,selectColumns]), GeneDF.tpm )
  #GeneDF_Norm_ScaledTpm <- cbind(data.frame(AnnotationFinal[,selectColumns]), GeneDF.ScaledTpm )
  
  # if (proteinCodingOnly) {
  #     GeneDF_Norm_CPM   <- keepProteinCoding(geneMatrix = GeneDF_Norm_CPM, annotation=AnnotationNew, featureType=featureName)[[1]]
  #     GeneDF_Norm_rpkm  <- keepProteinCoding(geneMatrix = GeneDF_Norm_rpkm, annotation=AnnotationNew, featureType=featureName)[[1]]
  #     GeneDF_Norm_tpm   <- keepProteinCoding(geneMatrix = GeneDF_Norm_tpm, annotation=AnnotationNew, featureType=featureName)[[1]]
  #     AnnotationFinal   <- keepProteinCoding(geneMatrix = countObj, annotation=AnnotationNew, featureType=featureName)[[2]]
  # }
  # 
  ########################################### Choose approprite folder and write files 
  if(merge) {
    RawCount=FPKM=CPM=TPM="/MergedFiles/"
  }else { RawCount="/RawCount/"; FPKM="/FPKM/" ; CPM="/CPM/" ; TPM="/TPM"}
  
  if(saveFiles)
  {
    print(paste("Writing final files"))
    saveRDS(countObj, paste(workdir,outputdirRDS, RawCount,fileName,"_Count_",annotName,".rds", sep= ""))
    saveRDS(GeneDF_Norm_rpkm, paste(workdir,outputdirRDS, FPKM,fileName,"_Norm_rpkm_",annotName,".rds", sep= ""))
    saveRDS(GeneDF_Norm_CPM, paste(workdir,outputdirRDS, CPM,fileName,"_Norm_cpm_",annotName,".rds", sep= ""))
    saveRDS(GeneDF_Norm_tpm, paste(workdir,outputdirRDS, TPM,fileName,"_Norm_tpm_",annotName,".rds", sep= ""))
    #saveRDS(GeneDF_Norm_ScaledTpm, paste(workdir,outputdirRDS, TPM,fileName,"_Norm_ScaledTPM_",annotName,".rds", sep= ""))
    
    write.table(countObj, paste(workdir,outputdirTXT, RawCount,fileName,"_Count_",annotName,".txt", sep= ""), sep="\t",row.names = TRUE, quote = FALSE)
    write.table(GeneDF_Norm_rpkm, paste(workdir,outputdirTXT, FPKM,fileName,"_Norm_rpkm_",annotName,".txt", sep= ""), sep="\t",row.names = FALSE, quote = FALSE)
    write.table(GeneDF_Norm_CPM, paste(workdir,outputdirTXT, CPM,fileName,"_Norm_cpm_",annotName,".txt", sep= ""), sep="\t",row.names = FALSE, quote = FALSE)
    write.table(GeneDF_Norm_tpm, paste(workdir,outputdirTXT, TPM,fileName,"_Norm_tpm_",annotName,".txt", sep= ""), sep="\t",row.names = FALSE, quote = FALSE)
    #write.table(GeneDF_Norm_ScaledTpm, paste(workdir,outputdirTXT, TPM,fileName,"_Norm_ScaledTPM_",annotName,".txt", sep= ""), sep="\t",row.names = FALSE, quote = FALSE)
  }
  
  #return(list(GeneDF.Norm, GeneDF_Norm_CPM, GeneDF_Norm_rpkm, GeneDF_Norm_tpm, GeneDF_Norm_ScaledTpm, AnnotationFinal))
  return(list(GeneDF.Norm, GeneDF_Norm_CPM, GeneDF_Norm_rpkm, GeneDF_Norm_tpm, AnnotationFinal))
  
}

performDiffGeneExp <- function(group1="", group2="Normal", group1Name="", group2Name="", metadataMapper.diffExpAnalysis="", mergeObjectsNoDup.diffExpAnalysis="",
                               selectGenesDF=NULL, maplot=FALSE, outDir=getwd(), method="exact", replicates=TRUE, proteinCodingOnly=FALSE,
                               groupfactorName = ""){
  
  ##Check if Dir exsists else make new one
  dir.create(file.path(outDir))
  group1.Filt <- paste(paste("^",group1, sep="") %>% paste(.,"$", sep=""), collapse = "|")
  group2.Filt <- paste(paste("^",group2, sep="") %>% paste(.,"$", sep=""), collapse = "|")
  
  subsetSamples <- metadataMapper.diffExpAnalysis %>% filter_(.dots=paste0("grepl(", "'", group1.Filt,"|",group2.Filt, "'" ,",", groupfactorName, ") & 
                                                                           !grepl(", "'", paste(germlineTissues, collapse="|"), "'",",", groupfactorName , ")" ))
  subsetSamples <- subsetSamples %>% mutate( DIAGNOSIS.Alias =factor(subsetSamples[,groupfactorName], levels=c(group1,group2))) %>% 
    dplyr::arrange(DIAGNOSIS.Alias); dim(subsetSamples)
  View(subsetSamples)
  
  filterMatrix  <- as.data.frame(mergeObjectsNoDup.diffExpAnalysis) %>% dplyr::select(one_of(as.character(subsetSamples$SAMPLE_ID.Alias))); print(paste(dim(filterMatrix)))
  
  filter_criteria.group1 <- interp(~y %in% x, .values=list(y = as.name(groupfactorName), x = group1))
  filter_criteria.group2 <- interp(~y %in% x, .values=list(y = as.name(groupfactorName), x = group2))
  
  group1Samples <- subsetSamples %>% filter_(filter_criteria.group1)
  group2Samples <- subsetSamples %>% filter_(filter_criteria.group2)
  
  group1Count <- group1Samples %>% summarise(total =n()) %>% as.numeric()
  group2Count <- group2Samples %>% summarise(total =n()) %>% as.numeric()
  
  if (group1Name != "" | group2Name != "") { 
    modelGroup = c(rep(group1Name, group1Count ), rep( group2Name, group2Count ))
  } else {
    modelGroup = c(rep( paste(group1, collapse="-"), group1Count ), rep( paste(group2, collapse="-"), group2Count ))
  }
  
  NormalizedDfs <- normalize(countObj=filterMatrix, annotName="Gene", method="EdgeR", annotationRDS=annotationRDS ,
                             fileName=mergeOutputFileName, merge=T, saveFiles=F, condition=modelGroup, proteinCodingOnly = proteinCodingOnly)
  
  groups <- unique(modelGroup)
  print(paste(groups))
  GeneDF.Norm = NormalizedDfs[[1]]
  GeneDF.FPKM = cbind(NormalizedDfs[[3]][,c(1:5)] , log2( NormalizedDfs[[3]][,-c(1:5)] + 1 ))
  rownames(GeneDF.FPKM) <- paste(GeneDF.FPKM[,5], GeneDF.FPKM[,4], sep="_")
  ### Using EdgeR Perform Differential Gene Expression Analsys
  
  designDisp                    <- model.matrix( ~modelGroup )
  
  if(replicates)  GeneDF_Dispersion             <- estimateDisp(GeneDF.Norm, design = designDisp )
  if(!replicates) GeneDF_Dispersion             <- estimateGLMCommonDisp(GeneDF.Norm, method="deviance",robust="TRUE",subset=NULL )
  
  if (method=="exact"){
    print(paste("Performing Differential Gene Expression Analsys using ExactTest"))
    GeneDF_DifferentialExpression <- exactTest(GeneDF_Dispersion, pair = c(unique(modelGroup)))
    GeneDF_DiffExp                <- GeneDF_DifferentialExpression$table
  } else {
    print(paste("Performing Differential Gene Expression Analsys NB-GLM"))
    fit <- glmFit(GeneDF_Dispersion, design = designDisp ) ; GeneDF_DifferentialExpression <- glmLRT(fit, coef=2)
    GeneDF_DiffExp               <- GeneDF_DifferentialExpression$table
  }
  
  GeneDF_DiffExp["AvglogFPKM"]        <- apply(GeneDF.FPKM[ , -c(1:5)] , 1, mean)
  if(replicates) {
    GeneDF_DiffExp[group1Name]        <- apply(GeneDF.FPKM[,as.character(group1Samples$SAMPLE_ID.Alias), drop=FALSE], 1, mean)
    GeneDF_DiffExp[group2Name]        <- apply(GeneDF.FPKM[,as.character(group2Samples$SAMPLE_ID.Alias), drop=FALSE], 1, mean)
  } else {
    GeneDF_DiffExp[group1Name]        <- GeneDF.FPKM[,as.character(group1Samples$SAMPLE_ID.Alias), drop=FALSE]
    GeneDF_DiffExp[group2Name]        <- GeneDF.FPKM[,as.character(group2Samples$SAMPLE_ID.Alias), drop=FALSE]
  }
  GeneDF_DiffExp["FDR"]   <- p.adjust(GeneDF_DiffExp$PValue, method="BH")
  
  
  print(paste("Saving Differential Gene Expression Analsys Results"))
  colnames(GeneDF_DiffExp)[c(1:4,7)] <- paste(group2Name, group1Name, colnames(GeneDF_DiffExp[c(1:4,7)]), sep = ".")
  #DiffExpSelectGenesWr <- cbind(data.frame("Genes"= rownames(DiffExpSelectGenes)), DiffExpSelectGenes)
  DiffExpSelectGenesWr <- tibble::rownames_to_column(GeneDF_DiffExp, var="GeneID") %>% data.frame()
  rownames(DiffExpSelectGenesWr) <- DiffExpSelectGenesWr$GeneID
  
  DiffExpSelectGenesWr <- left_join( DiffExpSelectGenesWr, GeneDF.FPKM[ , c(1:5)], by="GeneID")
  
  print(paste("1. Dimension of DiffExpSelectGenesWr " , dim(DiffExpSelectGenesWr)))
  
  ### Highlight and save the selected genes
  if( !is.null( selectGenesDF) ) {
    DiffExpSelectGenesWrFinal <- DiffExpSelectGenesWr %>% dplyr::filter(GeneName %in% c(selectGenesDF$GeneName))
    #DiffExpSelectGenesWr <- dplyr::full_join(selectGenesDF, DiffExpSelectGenesWr, by="GeneID")
    #DiffExpSelectGenesWrFinal <- DiffExpSelectGenesWr[complete.cases(DiffExpSelectGenesWr),]
    selectGenes=selectGenesDF$GeneID
  } else {
    DiffExpSelectGenesWrFinal <- DiffExpSelectGenesWr
    selectGenes=c()
  }
  
  print(paste("2. Dimension of DiffExpSelectGenesWrFinal " , dim(DiffExpSelectGenesWrFinal)))
  
  #### Append the MetaData
  #DiffExpSelectGenesWrFinal <- left_join( DiffExpSelectGenesWrFinal, GeneDF.FPKM[ , c(1,2,3,5)], by="GeneID")
  #print(paste("3. Dimension of DiffExpSelectGenesWrFinal " , dim(DiffExpSelectGenesWrFinal)))
  
  
  write.table(DiffExpSelectGenesWrFinal, file = paste(outDir, "/", paste(unique(modelGroup), collapse = "-VS-"), "selected.txt", sep=""), sep="\t", quote = FALSE, 
              row.names = F)
  #write.table(DiffExpSelectGenesWr, file = paste(outDir, "/", paste(unique(modelGroup), collapse = "-VS-"), "All.txt", sep=""), sep="\t", quote = FALSE, 
  #            row.names = F)
  
  ##ForMA plot
  if(maplot){
    print(paste("making MA Plot"))
    plot <- maPlot(GeneDF_DiffExp=GeneDF_DiffExp, minFCUP = 1.5, minFCDW = -1.5, group1=unique(modelGroup)[1], group2=unique(modelGroup)[2], 
                   selectGenes=selectGenes)
    return(plot)
  }
}

StringBeanPlot <- function(colList, Scores, orderOfFactor, orderOfSignature, standardize=TRUE){
  
  #function
  seqfunc <- function(x, start){ return(seq(start, start+1, length.out = x))}
  
  Scores$Diagnosis <- factor(Scores$Diagnosis, levels = orderOfFactor, ordered = TRUE)
  Scores  <- Scores %>% group_by(Diagnosis) %>% arrange(Diagnosis)
  
  counts <- Scores %>% group_by(Diagnosis) %>% summarise(n= n())
  xvals <- c()
  for(count in counts$n){ 
    xvalsNew <- seqfunc(count, start = 1);
    xvals <- c(xvals, xvalsNew)
  }
  
  Scores <- Scores %>% mutate_if(is.factor, as.character)
  Scores[colList]<- lapply(Scores[,colList], as.numeric)
  
  if(standardize==TRUE){
    Scores[,colList] <- apply(Scores[,colList], 2, zscore_All )
  }
  
  tidyScores <- Scores %>% tidyr::gather(Signatures, Scores, colList); 
  tidyScores$Signatures <- factor(tidyScores$Signatures, levels = orderOfSignature, ordered = TRUE)
  tidyScores <- tidyScores %>% dplyr::group_by(Signatures, Diagnosis, Scores) %>% arrange(Signatures, Diagnosis, Scores)
  tidyScores[,"SNONorm"] <- rep(xvals, length(orderOfSignature))
  
  ##Make median Segment
  medianY <- (tidyScores %>% group_by(Diagnosis, Signatures) %>% summarise(medianY=median(Scores)) %>% arrange(Signatures))$medianY
  medianX <- (tidyScores %>% group_by(Diagnosis, Signatures) %>% summarise(medianX=median(SNONorm))  %>% arrange(Signatures))$medianX
  segmentDF <- data.frame( xstart = medianX-0.05, ystart=medianY, xend=medianX+0.05, yend=medianY)
  segmentDF <- cbind(segmentDF, expand.grid(Diagnosis=unique(tidyScores$Diagnosis),Signatures=unique(tidyScores$Signatures)))
  #segmentDF <- data.frame( xstart = Xstart, ystart=medianY, xend=Xend, yend=medianY)
  
  summaryStats <- tidyScores %>% group_by(Diagnosis) %>% summarise(maxV = max(Scores), minV =min(Scores) )
  
  plot <- ggplot()+
    geom_point(data=tidyScores, aes(SNONorm, Scores, color = factor(Diagnosis)),show.legend = F, size=1.1)+
    #geom_violin(data=tidyScores, aes(SNONorm, Scores, color = factor(Diagnosis)),show.legend = F)+
    facet_grid(Signatures~Diagnosis ,switch = "both")+
    #ylim( -1.5, 1.5) + 
    #ylim( min(summaryStats$minV)-0.05,max(summaryStats$maxV)+0.05) + 
    #scale_colour_manual(values=c("dodgerblue4", "darkgoldenrod3", "firebrick4", "forestgreen")) +
    theme_bw() +
    geom_hline(yintercept=0, size=0.1) + 
    theme(axis.text.x = element_blank()
          ,axis.ticks.x =element_blank()
          ,axis.text.y = element_text(size=10, face="bold")
          ,axis.title.y=element_blank()
          ,axis.title.x=element_blank()
          ,panel.grid.major.x=element_blank()
          ,strip.text.y=element_text(size=14,face="bold")
          ,strip.text.x=element_text(size=14,face="bold", angle=90)
          ,strip.background=element_blank()
          ,panel.spacing=unit(.08, "lines")
          ,panel.spacing.y=unit(1, "lines")
          
    )+
    geom_segment(data = segmentDF, aes(x = xstart, xend = xend, y = ystart, yend = yend), size=1.5
                 #, linetype=2
                 , colour="red"
                 ,inherit.aes=FALSE
    )
  #labs(y="Standardized Infiltration Score")
  
  
  return(plot)
} 

OneVariablePlotSort <- function(colList, Scores, orderOfFactor, orderOfSignature, standardize=FALSE, logit =FALSE,
                                plotType="StringBean",customColors=NA, yLab="Score", summaryHlines =FALSE, 
                                sizeOfDots = 1, legendDisplay=TRUE){
  
  if (unique(is.na(customColors))) { customColors = setNames( StatsFinal$Color, StatsFinal$Diagnosis) }
  #function
  seqfunc   <- function(x, start){ return(seq(start, start+1, length.out = x))}
  xaxisSeq  <- function(x) {
    counts                  <- x %>% dplyr::group_by(Diagnosis) %>% dplyr::summarise(n= n())
    xvals                   <- c()
    for(count in counts$n){ 
      xvalsNew <- seqfunc(count, start = 1);
      xvals <- c(xvals, xvalsNew)
    }
    return(xvals)
  }
  drawStringBeanPlot <- function(x, tidyScoresPre){
    
    print(paste(x))
    tidyScores             <- tidyScoresPre %>% filter(Signatures == x) %>%  dplyr::group_by(Signatures,Diagnosis) %>% 
      dplyr::mutate(Med=median(Scores)) %>% 
      arrange(Signatures,Diagnosis,Scores) %>% 
      arrange(desc(Med)) %>% 
      ungroup() %>% 
      mutate( Diagnosis = factorizeColumn(Diagnosis, unique(as.character(Diagnosis) ) ),
              Color =  factorizeColumn(Color, unique(as.character(Color) ) ) ) %>% arrange(Diagnosis)
    tidyScores[,"SNONorm"] <- xaxisSeq(tidyScores)
    
    ##Make median Segment
    medianY <- (tidyScores %>% dplyr::group_by(Diagnosis, Signatures) %>% dplyr::summarise(medianY=median(Scores))  %>% dplyr::arrange(Diagnosis,Signatures))$medianY
    medianX <- (tidyScores %>% dplyr::group_by(Diagnosis, Signatures) %>% dplyr::summarise(medianX=median(SNONorm)) %>% dplyr::arrange(Diagnosis,Signatures))$medianX
    segmentDF <- data.frame( xstart = medianX-0.05, ystart=medianY, xend=medianX+0.15, yend=medianY)
    segmentDF <- cbind(segmentDF, expand.grid(Diagnosis=unique(tidyScores$Diagnosis),Signatures=unique(tidyScores$Signatures)))
    
    summaryStats <- tidyScores %>% group_by(Diagnosis) %>% summarise(maxV = max(Scores), minV =min(Scores) )
    scoreSummary <- summary(tidyScores$Scores)
    
    plot <- ggplot() +
      geom_point(data=tidyScores, aes(SNONorm, Scores, colour = factor(Diagnosis) ),show.legend = F, size=sizeOfDots) +
      scale_colour_manual(values=customColors  ) +
      #geom_violin(data=tidyScores, aes(SNONorm, Scores, color = factor(Diagnosis)),show.legend = F)+
      facet_grid(Signatures~Diagnosis ,switch = "both") +
      #ylim( min(summaryStats$minV)-0.05,max(summaryStats$maxV)+0.05) + 
      #ylim(-3,2.5) +
      labs( title= x ) +
      ylab( yLab ) +
      theme_bw() +
      geom_hline(yintercept=0, size=0.1) +
      theme( title = element_text(size=13, face="bold")
             ,axis.title.x = element_blank()
             ,axis.title.y=element_text(size=10, face="bold")
             ,axis.text.x = element_blank()
             ,axis.text.y = element_text(size=10, face="bold")
             ,axis.ticks.x =element_blank()
             ,strip.text.y= element_blank()
             ,strip.text.x=element_text(size=10,face="bold", angle=90, vjust=1)
             ,strip.background=element_blank()
             ,panel.grid.major.x=element_blank()
             ,panel.grid.minor.x=element_blank()
             ,panel.border = element_rect(colour = "black", fill=NA, size=0.0000000002, linetype = 1)
             ,panel.spacing = unit(0, "cm")
             ,strip.switch.pad.grid = unit(0, "cm")
      ) +
      geom_segment(data = segmentDF, aes(x = xstart, xend = xend, y = ystart, yend = yend), size=1.5
                   #, linetype=2
                   , colour="red"
                   ,inherit.aes=FALSE
      ) + 
      scale_y_continuous( expand = c(0, -0.05),
                          limits = c(min(tidyScores$Scores)-0.25, max(tidyScores$Scores)+0.25)) +
      scale_x_continuous( expand = c(0.1, 0))
    
    
    if(summaryHlines) { plot <- plot + 
      geom_hline(yintercept = scoreSummary[[2]], linetype="dashed", colour="#8888ff", size=0.4) +
      geom_hline(yintercept = scoreSummary[[4]], linetype="dashed", colour="#8888ff", size=0.4) + 
      geom_hline(yintercept = scoreSummary[[5]], linetype="dashed", colour="#8888ff", size=0.4)
    }     
    colnames(segmentDF)[3] <- x ; segmentDF$Diagnosis <- factor(segmentDF$Diagnosis, levels = orderOfFactor, ordered = TRUE);
    segmentDF <- segmentDF %>% arrange(Diagnosis)
    return(list(plot, segmentDF[,c(1,3)]) )
  }
  drawDensityPlot <- function(x, tidyScoresPre=NA , orderOfFactor=NA, customColors=NA, yLab =yLab){
    
    print(paste(x))
    tidyScores             <- tidyScoresPre %>% filter(Signatures == x) %>%  dplyr::group_by(Signatures,Diagnosis) %>% 
      dplyr::mutate(Med=median(Scores)) %>% arrange(Signatures,Diagnosis,Scores) %>% 
      arrange(desc(Med)) %>% 
      ungroup() #%>% 
    #mutate( Diagnosis = factorizeColumn(Diagnosis)) %>% arrange(Diagnosis)
    #mutate( Diagnosis = factorizeColumn(Diagnosis, orderOfFactor )) %>% arrange(Diagnosis)
    tidyScores[,"SNONorm"] <- xaxisSeq(tidyScores)
    
    ##Make median Segment
    medianY <- (tidyScores %>% dplyr::group_by(Diagnosis, Signatures) %>% dplyr::summarise(medianY=median(Scores)) %>% dplyr::arrange(Diagnosis,Signatures))$medianY
    medianX <- (tidyScores %>% dplyr::group_by(Diagnosis, Signatures) %>% dplyr::summarise(medianX=median(SNONorm))  %>% dplyr::arrange(Diagnosis,Signatures))$medianX
    segmentDF <- data.frame( xstart = medianX-0.05, ystart=medianY, xend=medianX+0.15, yend=medianY)
    segmentDF <- cbind(segmentDF, expand.grid(Diagnosis=unique(tidyScores$Diagnosis),Signatures=unique(tidyScores$Signatures)))
    
    summaryStats <- tidyScores %>% group_by(Diagnosis) %>% summarise(maxV = max(Scores), minV =min(Scores) )
    
    plot <- ggplot(data=tidyScores, aes(x = Scores, y = Diagnosis, height = ..density..)) +
      # to avoid overlaps of mountains , rel_min_height = 0.005, scale=0.9
      geom_density_ridges2(aes(fill = Diagnosis)) +
      scale_fill_manual(values=customColors) +
      guides(fill = legendDisplay ) +
      geom_vline(data=tidyScores, mapping=aes(xintercept=0), linetype = "dashed", colour = "maroon", size=1 ) +
      scale_y_discrete(expand = c(0.01, 0), limits = unique(rev(tidyScores$Diagnosis))) +
      scale_x_continuous(expand = c(0.01, 0)) +
      theme_ridges() + 
      theme(legend.title = element_text(size=15, face="bold") ) +
      labs( title= x ) +
      ylab("") +
      xlab(yLab)
    
    colnames(segmentDF)[3] <- x ; segmentDF$Diagnosis <- factor(segmentDF$Diagnosis, levels = orderOfFactor, ordered = TRUE);
    segmentDF <- segmentDF %>% arrange(Diagnosis)
    return(list(plot, segmentDF[,c(1,3)]) )
  }
  
  if(logit == TRUE) {
    Scores[,colList] <- apply(Scores[,colList,drop=FALSE] + 1 , 2, log2 )
  }
  
  if(standardize==TRUE){
    Scores[,colList] <- apply(Scores[,colList,drop=FALSE], 2, zscore_All )
  }
  
  tidyScoresPre <- Scores %>% tidyr::gather(Signatures, Scores, colList);
  mergeDF       <- merge(tidyScoresPre, StatsFinal, by.x="Diagnosis", by.y="Diagnosis", all.x=TRUE)
  tidyScoresPre  <- mergeDF[,c(1:4)] ; tidyScoresPre$Diagnosis <- factor(tidyScoresPre$Diagnosis, levels = orderOfFactor, ordered = TRUE)
  
  if( plotType =="StringBean") {
    plotList <- lapply(orderOfSignature, drawStringBeanPlot, tidyScoresPre)
  } else {
    plotList <- lapply(orderOfSignature, drawDensityPlot, tidyScoresPre=tidyScoresPre, orderOfFactor=orderOfFactor, customColors=customColors,
                       yLab =yLab)
  }
  return(plotList)
}   

maPlot <- function(GeneDF_DiffExp, minFCUP = 1.5, minFCDW = -1.5, group1="Normal", group2="Tumor", selectGenes=NULL, hgnc=""){
  
  logCPM3rdQ <- as.numeric( summary(GeneDF_DiffExp$logCPM)["3rd Qu."] ) 
  y.cept <- data.frame(point=c(1.5,-1.5), event = c("logFC >= 1.5", "logFC <= -1.5") )
  x.cept <- data.frame(point=c(0), event = c("Minimum CPM >= 1") )
  
  GeneDF_DiffExpPval <- GeneDF_DiffExp[which(GeneDF_DiffExp$FDR <= 0.05), ];
  print(paste( dim(GeneDF_DiffExpPval)))
  GeneDF_DiffExpPval$Color <- c("grey");    GeneDF_DiffExpPval$Genes <- rownames(GeneDF_DiffExpPval);
  GeneDF_DiffExpPval[which(GeneDF_DiffExpPval$logFC >=minFCUP), "Color"] <- "salmon1"
  GeneDF_DiffExpPval[which(GeneDF_DiffExpPval$logFC <=minFCDW), "Color"] <- "seagreen1"
  GeneDF_DiffExpPval[which(GeneDF_DiffExpPval$logCPM <=0 ), "Color"]      <- "grey"
  # GeneDF_DiffExpPval <- mutate(GeneDF_DiffExpPval, Color=if_else(Genes %in% hgnc$symbol, "salmon1", "mediumvioletred")); 
  # rownames(GeneDF_DiffExpPval) <- GeneDF_DiffExpPval$Genes
  # GeneDF_DiffExpPval[which(GeneDF_DiffExpPval$logCPM <=0 & GeneDF_DiffExpPval$logFC <=minFCUP & GeneDF_DiffExpPval$logFC >=minFCDW), "Color"]      <- "grey"
  
  if(is.null(selectGenes) ){
    
    GeneDF_DiffExpPval[which(GeneDF_DiffExpPval$logCPM >= logCPM3rdQ & GeneDF_DiffExpPval$logFC > minFCUP), "Color"] <- "indianred4"
    GeneDF_DiffExpPval[which(GeneDF_DiffExpPval$logCPM >= logCPM3rdQ & GeneDF_DiffExpPval$logFC < minFCDW), "Color"] <- "seagreen4"  
  }
  
  P <- ggplot(GeneDF_DiffExpPval, aes(logCPM, logFC, label=Genes)) + geom_point(colour=GeneDF_DiffExpPval$Color, size=1) +
    geom_hline(yintercept = 0, colour= "peachpuff4") +
    geom_hline(data=y.cept, mapping = aes(yintercept = point), colour= "peachpuff4") +
    geom_text(data=y.cept,  mapping = aes(y=point, x=10, label=event), size=4, angle=0, hjust=-0.4, vjust=0, fontface="bold") +
    geom_vline(data=x.cept, mapping = aes(xintercept = point), colour= "orange", size =1) +
    geom_text(data=x.cept,  mapping = aes(x=point, y=10, label=event), size=4, angle=90, vjust=-0.4, hjust=0, fontface="bold") +
    theme_bw() +
    xlim(-5, max(GeneDF_DiffExpPval$logCPM)+1) +
    ylim(-10, 15) +
    labs(title = paste("Mean-Average plot (MA-plot) for", group2, "vs", group1) ) + ylab("log-fold-change") +xlab("Mean Normalised Count")
  
  if(is.null(selectGenes) ){
    
    P <- P + geom_text_repel(aes(label=ifelse(logCPM >= logCPM3rdQ & logFC > minFCUP , as.character(Genes),'')),size=2,fontface = 'bold', color = 'navy') +
      geom_text_repel(aes(label=ifelse(logCPM >= logCPM3rdQ & logFC < minFCDW , as.character(Genes),'')),size=2,fontface = 'bold', color = 'navy')
    
  } else {
    selectGenesData = GeneDF_DiffExpPval %>% filter(Genes %in% selectGenes)
    P <- P + geom_point(data=selectGenesData, aes(logCPM, logFC), colour="navyblue", size=2)  +
      geom_text(data=selectGenesData, aes(label=Genes), hjust=-0.2, vjust=0, size=2,fontface = 'bold', color = 'navy')
  } 
  return(P)
}

volPlot <- function(GeneDF_DiffExp, minFCUP = 1.5, minFCDW = -1.5, group1="Normal", group2="Tumor", selectGenes=NULL, hgnc=""){
  
  #logPAdjust3rdQ <- as.numeric( summary(-log10(GeneDF_DiffExp$P.adjust))["3rd Qu."] )
  x.cept <- data.frame(point=c(1.5,-1.5), event = c("logFC >= 1.5", "logFC <= -1.5") )
  #y.cept <- data.frame(point=c(0), event = c(paste("Minimum FDR <=", logPAdjust3rdQ, sep="") ))
  
  GeneDF_DiffExpPval <- GeneDF_DiffExp
  GeneDF_DiffExpPval$Color <- c("grey");    GeneDF_DiffExpPval$Genes <- rownames(GeneDF_DiffExpPval);
  GeneDF_DiffExpPval[which(GeneDF_DiffExpPval$logFC >=minFCUP), "Color"] <- "salmon1"
  GeneDF_DiffExpPval[which(GeneDF_DiffExpPval$logFC <=minFCDW), "Color"] <- "seagreen1"
  
  P <- ggplot(GeneDF_DiffExpPval, aes(logFC, -log10(GeneDF_DiffExp$PValue), label=Genes)) + geom_point(colour=GeneDF_DiffExpPval$Color, size=1) +
    geom_vline(data=x.cept, mapping = aes(xintercept = point), colour= "peachpuff4") +
    geom_text(data=x.cept,  mapping = aes(x=point, y=4, label=event), size=4, angle=90, vjust=-0.4, hjust=0, fontface="bold") +
    #geom_vline(data=x.cept, mapping = aes(xintercept = point), colour= "orange", size =1) +
    theme_bw() +
    labs(title = paste("FDR vs logFC", group2, "vs", group1) ) + ylab("-log10(FDR)") +xlab("log fold change")
  
  
  if(!is.null(selectGenes) ){
    P <- P + geom_text(aes(label=ifelse(-log10(PValue) >= 4, as.character(Genes),'')),size=2,fontface = 'bold', color = 'navy')
  } else {
    selectGenesData = GeneDF_DiffExpPval %>% filter(Genes %in% selectGenes)
    P <- P  +
      geom_point(data=selectGenesData, aes(FDR, logFC), colour="navyblue", size=2)   +
      geom_text(data=selectGenesData, aes(label=Genes), hjust=-0.2, vjust=0, size=2,fontface = 'bold', color = 'navy')
  } 
  return(P)
}

hcPlot <- function(geneMatrix=NULL, expectedGroups=NULL){
  
  #### Trying 
  # P <- RPKM_Data_Filt_t %>% dist(method = "euclidean") %>% hclust( method="ward.D") %>% as.dendrogram %>% dendro_data(type = "rectangle") %>%
  #         set("labels_col", value = myPalette, k=expectedGroups) %>% 
  #         set("branches_k_color", value = myPalette, k = expectedGroups) %>%
  #         set("labels_cex", c(.9,1.2)) %>% 
  #         set("leaves_pch", 19)
  # ggplot(segment(P)) + 
  #   geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  #   coord_flip() + 
  #   scale_y_reverse(expand = c(0.2, 0)) + coord_polar(theta="x")
  #return(P)
}

makeRNKFiles <- function(DF, projectName=""){
  colNames <- colnames(DF)
  sapply(colNames[1], function(x){ 
    
    y <- DF[, x, drop=FALSE] %>% rownames_to_column("Name") %>% arrange_(.dots = paste0("desc(",x,")") )
    write.table(y, paste("./GSEA/rnk/", projectName,"/",x,".rnk",sep=""), sep = "\t",row.names = FALSE, quote = FALSE )})
}

PreRankedGSEA <- function(x, projectName=NULL, sample=NULL){
  
  label <- gsub(".gmt","",paste(sample,".",x,sep=""))
  cmd <- paste("java -cp /data/sindiris/Processing/GSEA/gsea2-2.2.2.jar -Xmx1024m xtools.gsea.GseaPreranked",
               " -gmx /data/sindiris/Processing/GSEA/",x,
               " -rnk /data/sindiris/Processing/GSEA/",projectName,"/rnk/",sample,
               " -chip /data/sindiris/Processing/GSEA/GENE_SYMBOL.chip ",
               " -collapse false",
               " -mode Max_probe",
               " -norm None",
               " -nperm 1000",
               " -scoring_scheme weighted",
               " -rpt_label ",label,
               " -include_only_symbols  true",
               " -make_sets true",
               " -plot_top_x 20",
               " -rnd_seed 149",
               " -set_max 500",
               " -set_min 15",
               " -out /data/sindiris/Processing/GSEA/",projectName,
               " -gui false",
               sep="")
  return(cmd)
}

NESorPvalGSEAPrerank <- function(x, colNumb = 1){
  
  files <- list.files(x, pattern = '^gsea_report_for_na_(pos|neg)_[1234567890]+\\.xls$', full.names = TRUE)
  df <- do.call(rbind, lapply(files, getCountObjTXT, colNumb=colNumb, rowNames=1)); print(paste(dim(df)))
  colnames(df) <- gsub("\\.\\/GSEA\\/prerankedGSEAOutput\\/\\/|.rnk.*.\\S+\\/gsea_report_for_na_(pos|neg)_\\S+.xls", "", files)[1]
  return(df)
}

filterSpecificCloneTypes  <- function(cloneData, cloneType){
  cloneDataFilt           <- cloneData %>% dplyr::filter(grepl(cloneType,v) | grepl("NF",v))
  return(cloneDataFilt)
}

correlationPlots <- function(varName="", constName="", df=NA){
  
  print(paste(varName))
  # plot <- ggscatter(df, x = constName, y =varName, 
  #                   add = "reg.line", conf.int = TRUE, 
  #                   cor.coef = TRUE, cor.method = "spearman",
  #                   xlab = constName, ylab = varName)
  
  corrTest <- cor.test(df[,constName], df[,varName], method = "spearman")
  if  ( corrTest$p.value < 2.2e-16 ) { corrTest$p.value = 2.2e-16 }
  plot <- ggplot(df, aes_string(x=constName, y=varName)) + 
    geom_smooth(method=lm,  fill="grey") +
    geom_point(aes(colour = factor(DIAGNOSIS.Substatus.Tumor.Normal.Tissue)), show.legend = T, size=4) + 
    scale_colour_manual(values=setNames( StatsFinal$Color, StatsFinal$DIAGNOSIS.Substatus.Tumor.Normal.Tissue)  ) +
    theme_bw() +
    theme(axis.text=element_text(size=13)
          ,axis.title=element_text(size=13,face="bold")) +
    xlab("Log Total Clones") +
    ylab(paste("Standardised Enrichment Score", sep=" "))+
    ggtitle(paste("Corr.Coeff = ", signif(corrTest$estimate[[1]],5), "\np-value = ", signif(corrTest$p.value,5), "                                            ", varName,sep=""))
  
  return(list(plot))
}   

makeEntropyInput <- function(filename, inputDir="", outputDir="", cloneType="") {
  outfileName <- paste0(outputDir, gsub("convert.|.clones.txt","",filename ), ".Entropy.txt")
  exomeData <- read.csv( paste(inputDir, filename, sep=""), sep="\t", header = TRUE )
  
  if(cloneType != "") {
    exomeData               <- filterSpecificCloneTypes(cloneData = exomeData, cloneType = cloneType)
  }
  
  if(nrow(exomeData)>0){
    exomeDataEntropy <- data.frame(VJcombo=paste(exomeData$v,exomeData$j,sep="."), Counts=exomeData$count, Vcassette=exomeData$v, 
                                   Jcassette=exomeData$j, aaCDR3_filtered = exomeData$cdr3aa, ntCDR3= exomeData$cdr3nt)
  } else {
    exomeDataEntropy <- emptyDFEntropy
  }
  write.table(exomeDataEntropy, outfileName, sep = "\t", row.names = FALSE, quote = FALSE)
}

immunoseqv2 <- function(x) {
  filename = x
  #outfileName <- paste0("input/MiXCR/immunoseqv2EntropyNitin/", gsub("convert.|.clones.txt","",filename ), ".Entropy.txt")
  exomeData <- read.csv( paste("input/MiXCR/immunoseqv2/", filename, sep=""), sep="\t", header = TRUE, stringsAsFactors = FALSE )
  exomeDataFilt <- exomeData[,1:47] %>% dplyr::filter(sequenceStatus %in% c("In"))
  colnames(exomeDataFilt)[3] <- "count..templates.reads."
  
  print(paste("input/MiXCR/immunoseqv2/", filename, sep=""))
  print(dim(exomeDataFilt))
  
  if(nrow(exomeDataFilt)>0){
    
    emptyVGeneName <- which(exomeDataFilt$vGeneName == "")
    emptyJGeneName <- which(exomeDataFilt$jGeneName == "")
    exomeDataFilt[emptyVGeneName, c("vGeneName")] <- sapply(exomeDataFilt[emptyVGeneName,c("vGeneNameTies")], function(x){  unlist(strsplit(as.character(x), ","))[1] })
    exomeDataFilt[emptyJGeneName, c("jGeneName")] <- sapply(exomeDataFilt[emptyJGeneName,c("jGeneNameTies")], function(x){  unlist(strsplit(as.character(x), ","))[1] })
    
    
    
    exomeDataEntropy <- data.frame(VJcombo=paste(exomeDataFilt$vGeneName,exomeDataFilt$jGeneName,sep="."), 
                                   Counts=exomeDataFilt$count..templates.reads., 
                                   Vcassette=exomeDataFilt$vGeneName, 
                                   Jcassette=exomeDataFilt$jGeneName, 
                                   aaCDR3_filtered = exomeDataFilt$aminoAcid, 
                                   ntCDR3= exomeDataFilt$nucleotide)
  } else {
    exomeDataEntropy <- emptyDFEntropy
  }
  #write.table(exomeDataEntropy, outfileName, sep = "\t", row.names = FALSE, quote = FALSE)
  
}

extractEntropyResults_ShenLabPaper <- function(x,value){
  print(x)
  exomeData <- read.csv( paste(folderName, x, sep=""), sep="\t", header = TRUE )
  print(dim(exomeData))
  if(nrow(exomeData) == 0){
    exomeData <- emptyDFEntropyResults
  }
  return(exomeData[value,])
}
