library(R6)
library(assertthat)
library(microbenchmark)

# ### Project SetUp Class ####
ProjectSetUp <- R6Class(
  classname = "ProjectSetUP",
  portable = TRUE,
  private  = list(
    
    ## Check validity of directories
    checkDirExists = function() {
      
      projectDirPath = self$workDir
      ## Check & create if the project directory exists ?? 
      if(!dir.exists(file.path(self$workDir))) {
        warning(paste(self$workDir, "No Such directory exsist !!"))
        paste0("Creating working dir for this project ", projectDirPath)
        dir.create(self$workDir, showWarnings = F)
      }
      
      ## Create all the required dirs for the project.
      dir.create(projectDirPath, showWarnings = F)
      self$fileDirs <- paste(projectDirPath, c(self$plotsDir, self$plotsDataDir, self$DiffGeneExpAnaDir, self$DiffGeneExpRDS), sep="/")
      lapply(self$fileDirs, function(x){ if(!dir.exists(x)) dir.create(x, showWarnings = F) })
    },
    ## Generate filters to exclude given list
    generateFiltersToExclude = function(factorsToExclude=NA){
      
      lapply(factorsToExclude, function(x){
        subList=x
        sublistNames=names(x)
        paste0("!(",
               paste0(sapply(sublistNames, function(y){
                 
                 #print(subList[y])
                 filter <- paste0(y," %in% ","\"",unlist(unname(subList[y])), "\"" )
                 return(filter)
               }), collapse = " & " ),
               ")"
        )
      })
    },
    ## Read metadata file
    readMetaData = function() {
      self$metaDataDF <- read.csv(paste0(self$workDir,"/",self$metaDataFileName), sep="\t", header = T);
      print(paste0("Dimension of metadata is ", paste(dim(self$metaDataDF))))
      if(!is.na(self$factorsToExclude)){
        filters <- private$generateFiltersToExclude(factorsToExclude=self$factorsToExclude)
        print(paste(filters))
        test <- sapply(filters, function(x){
          self$metaDataDF <- dplyr::filter_(self$metaDataDF, .dots=list(x))
        })
      }
      print(paste0("Dimension of metadata after applying parameter 'factorsToExclude' is ", paste(dim(self$metaDataDF)[1])))
      print("Make R valid names in the metadata file and storing it as validMetaDataDF")
      tempDF <- as.data.frame(lapply(self$metaDataDF[,c("Patient.ID","Sample.Data.ID","Sample.Biowulf.ID.GeneExp","Sample.ID.Alias")], make.names)) %>% data.frame()
      self$validMetaDataDF <- cbind(tempDF, self$metaDataDF[,which(!names(self$metaDataDF) %in% c("Patient.ID","Sample.Data.ID","Sample.Biowulf.ID.GeneExp","Sample.ID.Alias"))])
    },
    ## Read reference annotation file ENSEMBL or UCSC
    readAnnotation = function() {
      self$annotationDF <- readRDS(self$annotationRDS) %>% as.data.frame()
    },
    ## Read ProteinCoding annotation file
    readProteinCoding = function(){
      self$pcDF <- as.data.frame( readRDS(self$pcRDS) )
      colnames(self$pcDF)[c(1,2)] <- c("GeneName","GeneID")
    },
    ## Read TranscriptionFactor annotation file
    readTranscriptionFactor = function(){
      self$tfDF <- readRDS(self$tfRDS)
    },
    ## Read ExhaustionMarker annotation file
    readExhaustionMarker = function(){
      self$emDF <- readRDS(self$emRDS)
    },
    ## Read CellSurface annotation file
    readCellSurface = function(){
      self$csDF <- readRDS(self$csRDS)  %>% dplyr::filter(NewCount >= 5)
    },
    ## Read CancerGermlineAntigens annotation file
    readCancerGermlineAntigens = function(){
      self$cgaDF <- readRDS(self$cgaRDS) 
    },
    ## Read pax3Foxo1Targets annotation file
    readpax3Foxo1Targets = function(){
      self$pax3Foxo1DF <- readRDS(self$pax3Foxo1RDS) 
    },
    ## Read ewsr1Fli1Targets annotation file
    readewsr1Fli1Targets = function(){
      self$ewsr1Fli1DF <- readRDS(self$ewsr1Fli1RDS) 
    },
    ## Read BrainExp annotation file
    readBrainExp = function(){
      self$BrainExpDF <- readRDS(self$BrainExpRDS) 
    },
    ## Read HeartExp annotation file
    readHeartExp = function(){
      self$HeartExpDF <- readRDS(self$HeartExpRDS) 
    },
    ## Read KidneyExp annotation file
    readKidneyExp = function(){
      self$KidneyExpDF <- readRDS(self$KidneyExpRDS) 
    },
    ## Read LiverExp annotation file
    readLiverExp = function(){
      self$LiverExpDF <- readRDS(self$LiverExpRDS) 
    },
    ## Read LungEx annotation file
    readLungExpRDS = function(){
      self$LungExpDF <- readRDS(self$LungExpRDS) 
    },
    ## Get Color map
    getFactorColorMap = function(){
      customColorsDF <- self$metaDataDF %>% dplyr::filter(! LIBRARY_TYPE %in% c("CellLine", "Normal"))
      customColorsDF <-  customColorsDF[,c(self$factorName, "Color")] %>% data.frame() %>% dplyr::distinct()
      customColorsDF[,self$factorName]  <- gsub(".Tumor", "", customColorsDF[,self$factorName])
      colnames(customColorsDF)[1] <- "Diagnosis";
      self$customColorsDF <- customColorsDF
    },
    ## Get Color map
    getFactorColorMapAll = function(){
      customColorsDFAll <- self$metaDataDF %>% dplyr::filter(! LIBRARY_TYPE %in% c("Normal"))
      customColorsDFAll <-  customColorsDFAll[,c(self$factorName, "Color.Substatus")] %>% data.frame() %>% dplyr::distinct()
      # customColorsDFAll[,self$factorName]  <- gsub(".Tumor", "", customColorsDFAll[,self$factorName])
      colnames(customColorsDFAll)[1] <- "Diagnosis";
      self$customColorsDFAll <- customColorsDFAll
    }
    ## Make a Tree Map
  ),
  public = list(
    date                    = NULL,
    time                    = NULL,
    projectName             = NULL,
    annotationRDS           = NULL,
    outputPrefix            = NULL,
    filterGenes             = NULL,
    filterGeneMethod        = NULL,
    factorName              = NULL,
    metadataFileRefCol      = NULL,
    metaDataFileName        = NULL,
    workDir                 = NULL,
    outputdirRDSDir         = NULL,
    outputdirTXTDir         = NULL,
    gseaDir                 = NULL,
    plotsDir                = NULL,
    plotsDataDir            = NULL,
    DiffGeneExpAnaDir       = NULL,
    DiffGeneExpRDS          = NULL,
    metaDataDF              = NULL,
    annotationDF            = NULL,
    pcDF                    = NULL,
    pcRDS                   = NULL,
    tfDF                    = NULL,
    tfRDS                   = NULL,
    csDF                    = NULL, 
    csRDS                   = NULL,
    emDF                    = NULL, 
    emRDS                   = NULL,
    cgaDF                   = NULL,
    cgaRDS                  = NULL,
    pax3Foxo1DF             = NULL,
    pax3Foxo1RDS            = NULL,
    ewsr1Fli1DF             = NULL,
    ewsr1Fli1RDS            = NULL,
  
    BrainExpDF              = NULL,               
    BrainExpRDS             = NULL,
    HeartExpDF              = NULL,
    HeartExpRDS             = NULL,
    KidneyExpDF             = NULL,
    KidneyExpRDS            = NULL,
    LiverExpDF              = NULL,
    LiverExpRDS             = NULL,
    LungExpDF               = NULL,
    LungExpRDS              = NULL,
    
    factorsToExclude        = NULL,
    fileDirs                = NULL, 
    validMetaDataDF         = NULL,
    customColorsDF          = NULL,
    customColorsDFAll       = NULL,
    initialize              = function(date = NA, time =NA, projectName = NA, annotationRDS = NA, outputPrefix = NA,
                                       filterGenes = NA, filterGeneMethod = NA, factorName = NA, metaDataFileName = NA, 
                                       workDir = NA, outputdirRDSDir = NA, outputdirTXTDir = NA,gseaDir = NA, plotsDir = NA, 
                                       plotsDataDir = NA, DiffGeneExpAnaDir = NA, DiffGeneExpRDS = NA, pcRDS = NA,tfRDS=NA, emRDS=NA,
                                       csRDS =NA, cgaRDS=NA, pax3Foxo1RDS=NA, ewsr1Fli1RDS=NA, factorsToExclude=NA,
                                       BrainExpRDS=NA, HeartExpRDS=NA, KidneyExpRDS=NA, LiverExpRDS=NA, LungExpRDS=NA, metadataFileRefCol=NA) {
      
      self$date <- date
      self$time <- time
      self$projectName <- projectName
      self$annotationRDS <- annotationRDS
      self$outputPrefix <- outputPrefix
      self$filterGenes <- filterGenes
      self$filterGeneMethod <- filterGeneMethod
      self$factorName <- factorName
      self$metadataFileRefCol <- metadataFileRefCol
      self$metaDataFileName <- metaDataFileName
      self$workDir <- "./"
      self$outputdirRDSDir <- outputdirRDSDir
      self$outputdirTXTDir <- outputdirTXTDir
      self$gseaDir <- gseaDir
      self$plotsDir <- plotsDir
      self$plotsDataDir <- plotsDataDir
      self$DiffGeneExpAnaDir <- DiffGeneExpAnaDir
      self$DiffGeneExpRDS <- DiffGeneExpRDS
      
      self$pcRDS <- pcRDS
      self$tfRDS <- tfRDS
      self$emRDS <- emRDS
      self$csRDS <- csRDS
      self$cgaRDS <- cgaRDS
      
      self$pax3Foxo1RDS <- pax3Foxo1RDS
      self$ewsr1Fli1RDS <- ewsr1Fli1RDS
      
      self$BrainExpRDS <- BrainExpRDS
      self$HeartExpRDS <- HeartExpRDS
      self$KidneyExpRDS <- KidneyExpRDS
      self$LiverExpRDS <- LiverExpRDS
      self$LungExpRDS <- LungExpRDS
      
      
      self$factorsToExclude <- factorsToExclude
      private$checkDirExists()
      private$readMetaData()
      private$readAnnotation()
      private$getFactorColorMap()
      private$getFactorColorMapAll()
      if (!is.na(pcRDS)){ private$readProteinCoding() }
      if (!is.na(tfRDS)){ private$readTranscriptionFactor() }
      if (!is.na(emRDS)){ private$readExhaustionMarker() }
      if (!is.na(csRDS)){ private$readCellSurface() }
      if (!is.na(cgaRDS)){ private$readCancerGermlineAntigens() }
      
      if (!is.na(pax3Foxo1RDS)){ private$readpax3Foxo1Targets() }
      if (!is.na(ewsr1Fli1RDS)){ private$readewsr1Fli1Targets() }
      
      if (!is.na(BrainExpRDS)){ private$readBrainExp() }
      if (!is.na(HeartExpRDS)){ private$readHeartExp() }
      if (!is.na(KidneyExpRDS)){ private$readKidneyExp() }
      if (!is.na(LiverExpRDS)){ private$readLiverExp() }
      if (!is.na(LungExpRDS)){ private$readLungExpRDS() }
      
      print(paste0("Changing to working directory ", self$workDir))
      
    }
  ) 
)

# ### Utilities Class ####

CoreUtilities <- R6Class(
  classname = "CoreUtilities",
  portable  = TRUE,
  ## Read CSV or TXT files
  private   = list(
    readTXTFiles  = function(x, fileSuffix=NA, colNameSelect=NA, primaryID=NA ){
      if(!is.na(fileSuffix)) 
      {
        sampleName <- gsub(fileSuffix, "", basename(x)) 
      } else { 
        sampleName <- basename(x) 
      }
      ## print(paste0("I am in readTXTFiles ", sampleName, "  ", colNameSelect))
      print(x)
      rdsObj <- fread(x, sep="\t", header = TRUE) %>% rename_(.dots=setNames(colNameSelect, sampleName))
      ## Operation to include TCGA datasets.
      rdsObj[,c(primaryID) := sapply(.SD,gsub,pattern = "\\.[0-9]*",replacement = ""), .SDcols=c(primaryID)]
      return(rdsObj[, c(primaryID,sampleName), with=FALSE])
      
    },
    ## Merge CSV or TXT files
    mergeTXTFiles = function( x, fileSuffix=NA, colNameSelect=NA, primaryID=NA ){
      
      print(paste(" Total files in the input folder ", length(x), "Total libraries in the metadata ", length( self$allFileList)))
      #Df_results <- data.frame( basename(x), result = grepl(paste(self$allFileList, collapse = "|"),basename(x)))
      #selectedFileList <- x[which(Df_results$result == TRUE)]
      selectedFileList <- x[which(basename(x) %in% self$allFileList)]
      ## print(paste0('Which files in input folder are IN the meta data ',length(selectedFileList), paste(basename(selectedFileList), collapse=";")))
      notselectedFileListFolder      <- x[which(!basename(x)  %in%  self$allFileList )];
      ## print(paste0('Which files in the input folder are NOT IN the metadata ',length(notselectedFileListFolder), paste(basename(notselectedFileListFolder), collapse=";")))
      notselectedFileListMeta        <- self$allFileList[which(!self$allFileList  %in% basename(x))];
      print(paste0('Which files in the metadata are NOT IN the folder ',length(notselectedFileListMeta), paste(notselectedFileListMeta, collapse=";")))
      
      if( length(notselectedFileListMeta) >= 1) {                
        cat( paste(" Following files are not present in input folder.Please record them in metadata file or remove from the input file folder"
                   ,paste(basename(notselectedFileListMeta), collapse = "\n")))
        # cat("\n\n")
        # cat( paste(" Following files are not present in metadata File .Please record them in metadata file or remove from the input file folder"
        #       ,paste(basename(notselectedFileListFolder), collapse = "\n")))
        stop("please check the above error")
      }
      
      print(paste0("Selecting ", length(selectedFileList), " files out of ", length(x), " from the given folder"))
      View(data.frame("Selected"=basename(selectedFileList)))
      
      dataMatrixLists     <- lapply(selectedFileList, private$readTXTFiles, fileSuffix=fileSuffix,colNameSelect=colNameSelect, primaryID=primaryID)
      dataMatrix          <- purrr::reduce(dataMatrixLists, full_join, by=primaryID)
      
      return(dataMatrix)
    }
  ),
  public    = list(
    workDir                = NULL,
    projectName            = NULL,
    allFileList            = NULL,
    project                = NULL,
    initialize             = function(ProjectSetUpObject = NA ){
      
      assert_that("ProjectSetUP" %in% class(ProjectSetUpObject), 
                  msg="Please setup Project using ProjectSetUp class !!\nProjectSetUpObject cannot be NA !!")
      self$workDir     <- ProjectSetUpObject$workDir
      self$projectName <- ProjectSetUpObject$projectName
      
    },
    ## Get merged matrix
    getMergedMatrix = function(dir = NA, fileFormat = NA, colNameSelect = NA, colIndexSelect = NA, isRowNames = FALSE, rowNamesColInFile = NA,
                               fileSuffix=NA, primaryID=NA, metadata=NA, metadataFileRefCol=NA){
      
      if( !is.na(colNameSelect) & !is.na(colIndexSelect)) stopifnot("Use either colNameSelect or colIndexSelect but not both")
      if( is.na(colNameSelect) & is.na(colIndexSelect)) stopifnot("Use either colNameSelect or colIndexSelect but not both")
      assert_that( !is.na(fileFormat)  , msg ="Please provide file format to parse. Following are supported, \"rds\", \"csv or text\", \"feather\" ")
      #
      if( is.na(colNameSelect) & is.na(colIndexSelect)) {
        print("First column of every file will be selected.Or use either \"colNameSelect\" or \"colIndexSelect\"")
        colIndexSelect = 1
      }
      if ( is.na(colNameSelect) ) { 
        colNameSelect  =  colIndexSelect
      } 
      if(isRowNames) assert_that(!is.na(rowNamesColInFile), msg= "Please provide row names index , \"rowNamesColInFile\" can't be NA ")
      
      self$allFileList = paste0(as.character(metadata[,metadataFileRefCol]),fileSuffix)
      
      dirs                <- list.dirs(paste0(self$workDir,"/",dir))[-1]
      folderNames         <- basename(dirs)
      
      switch(fileFormat,
             
            "rds"     = {
               mergedDataList        <- lapply( lapply(dirs, list.files, full.names=T),
                                              private$mergeRDSFiles, fileFormat = fileFormat, colInterest = colIndexSelect, fileSuffix = fileSuffix,
                                                       colNameSelect = colNameSelect )
               rowNames              <- rownames(mergedDataList[[1]])
               mergedData            <- dplyr::bind_cols( mergedDataList)
               rownames(mergedData)  <- rowNames
               return(mergedData)
                    },
             "txt"    = {
               mergedDataList        <- lapply( lapply(dirs, list.files, full.names=T),
                                                private$mergeTXTFiles, fileSuffix = fileSuffix, colNameSelect = colNameSelect, primaryID=primaryID )
               mergedData            <-  purrr::reduce( mergedDataList,  full_join, by=primaryID )
               mergedData            <-  tibble::column_to_rownames(mergedData, var=primaryID)
               return(mergedData)
             }
      )
    },
    ## Consolidate the matrix
    consolidateDF = function(df = NA, featureName=NA, funcName=NA, colsToExclude = NA){
      assert_that( funcName %in% c("sum", "max", "mean"), msg= "For now duplicate featureName can be consolidated only by sum, max or mean" )
      ### Get character or factor columns
      listColsToExclude = names(which(sapply(df, class) == "character" | sapply(df, class) == "factor"))
      assert_that(length(listColsToExclude) != 0, msg="Character and Factor columns found. Please make sure the object has all numeric columns !!")
      if(!is.na(colsToExclude)) assert_that(colsToExclude %in% colnames(df), msg = "colsToExclude not found in the data frame !!")
      
      assert_that(!is.na(featureName), msg="featureName cannot be NA")
      
      if(funcName == "sum")  sumFunc <- function(x)  sum(x)
      if(funcName == "max")  sumFunc <- function(x)  max(x)
      if(funcName == "mean") sumFunc <- function(x)  mean(x)
      
      if(is.na(colsToExclude)) df[, lapply(.SD, sumFunc), by=featureName]
      else df[, lapply(.SD, sumFunc), by=featureName, .SDcols = -colsToExclude]
    },
    ## matchAndChangeColNames columns
    matchAndChangeColNames = function(){
      metaDT <-  setDT(rnaseqProject$metaDataDF)
      setDT(metaDT)
    },
    ## Factorise columns
    factorizeColumn = function( toFactor, asFactor){
      factorColumn <- factor(toFactor, levels=asFactor, ordered = TRUE )
      return(factorColumn)
    },
    ## Keep Protein Coding
    keepProteinCoding = function(geneMatrix=NULL, annotation=NULL, featureType="GeneID") {
      ###Keep only Protein Coding Genes
      # PC <- read.table("C:/Users/sindiris/R Scribble/Annotation RDS/HGNC-protein-coding-List.txt", header = T, sep="\t")
      # annotationPC <- annotation %>% filter(GeneName %in% PC$Genes)
      # foundGenes <- which(row.names(geneMatrix) %in% annotationPC[,featureType] )
      # geneMatrixNew <- geneMatrix %>% data.frame() %>% .[foundGenes,]
      # annotationPC <- annotationPC %>% mutate(GeneID = self$factorizeColumn(annotationPC$GeneID, row.names(geneMatrixNew)))
      # return(list(geneMatrixNew, annotationPC) )
      
      foundGenes <- which(row.names(geneMatrix) %in% rnaseqProject$pcDF[,featureType] )
      geneMatrixNew <- geneMatrix %>% data.frame() %>% .[foundGenes,]
      annotationPCNew <- annotation %>% filter(GeneID %in% rnaseqProject$pcDF[,featureType] )
      
      return(list(geneMatrixNew, annotationPCNew) )
      
    },
    ## convert FPKM to TPM
    fpkmToTpm = function(fpkm){
      exp( (log(fpkm) - log(sum(fpkm))) + log(1e6))
    },
    ## subsetMetaData
    subsetMetaData = function(colnamesDF=NA){
      df <- dplyr::left_join(colnamesDF,rnaseqProject$validMetaDataDF, by="Sample.Biowulf.ID.GeneExp") %>% dplyr::filter(complete.cases(.))
      if( ncol(df) == ncol(rnaseqProject$validMetaDataDF) | nrow(df) == nrow(rnaseqProject$validMetaDataDF ) ) {
        rnaseqProject$validMetaDataDF <- df
      } else {
        stop("There are possible multiple issues in the provided metadata input file. For example, syntatically incorrect entries/names in metadata")
      }
    },
    ## Annotate a gene expression df with "GeneID" as primary key
    featureNameAnot = function(annotationDF=NA, querryDF=NA, identifier=NA){
      annotDF <- dplyr::left_join( annotationDF, querryDF, by=identifier)
      print(paste("Annotating expression DF"))
      print(dim(annotDF))
      return(annotDF)
    },
    ## perform Hirarchial clustering
    performClustering = function( df = NA) {
      RPKM_Data_Filt_t=t(df)
      hc<-hclust(dist(RPKM_Data_Filt_t,"euclidean"),"ward.D")
      pdf(paste0(rnaseqProject$workDir,"/",rnaseqProject$plotsDir,"/",rnaseqProject$projectName,"_hc.pdf"),height=15,width=15)
      plot(as.phylo(hc),type = "fan", label.offset=1, cex=1)
      dev.off()
    },
    ## Standardise a matrix
    zscore_All = function( x = NA) {
      medX <- median(x)
      sdX <- sd(x, na.rm = FALSE)
      y <- (x - medX) / sdX
      return(y)
    },
    ## Create Broad ssGSEA input file
    createBroadGCTFile = function(x = NA ){
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
    },
    ## Parse Broad ssGSEA input file
    parseBroadGTCOutFile = function(fileName = NA){
      expressionTMM.RPKM.ssGSEA.output.pre <- read.csv(fileName,sep="\t",header = FALSE, stringsAsFactors = FALSE )[ -c(1:2), ]
      colnames( expressionTMM.RPKM.ssGSEA.output.pre ) <- as.character(unname(unlist(expressionTMM.RPKM.ssGSEA.output.pre[1,])))
      expressionTMM.RPKM.ssGSEA.output.pre <- expressionTMM.RPKM.ssGSEA.output.pre[-c(1),-c(2)] %>% tibble::remove_rownames() %>% 
        tibble::column_to_rownames(var = "Name")
      expressionTMM.RPKM.ssGSEA.output <- data.frame(lapply(expressionTMM.RPKM.ssGSEA.output.pre,as.numeric))
      rownames( expressionTMM.RPKM.ssGSEA.output ) <- rownames(expressionTMM.RPKM.ssGSEA.output.pre)
      return(expressionTMM.RPKM.ssGSEA.output)
    },
    ## Geometric mean
    calculateGeoMean = function(x){
      ## gemoMean
      y <- exp(mean(log(x[is.finite(log(x+0.01))]),na.rm=T))
      ## arithmatic mean
      ## y <- sum(log2(x+1))/6
      return(y)
    },
    ## calculate cytolytic scores
    cytolyticScore = function(expDF = NA ) {
      cytolyticDF       <- expDF[c("GZMA","GZMB","GZMH","GZMK", "GZMM", "PRF1"),] ; 
      CytolyticScores   <- apply(cytolyticDF,2, self$calculateGeoMean) %>% data.frame() %>% t()
      rownames(CytolyticScores) <- "CytolyticScore"
      return(CytolyticScores)
    },
    ## make one variable plots
    OneVariablePlotSort = function(colList=NA, Scores=NA, orderOfFactor=NA, orderOfSignature=NA, standardize=FALSE, logit =FALSE,
                                   plotType="StringBean",customColorDF=NA, yLab="Score", summaryHlines =FALSE, 
                                   sizeOfDots = 1, legendDisplay=TRUE, logBase=2){
      
      #if (unique(is.na(customColors))) { customColors = setNames( StatsFinal$Color, StatsFinal$Diagnosis) }
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
      drawStringBeanPlot <- function(x, tidyScoresPre, customColors=NA, yLab =yLab){
        
        print(paste(x))
        tidyScores             <- tidyScoresPre %>% filter(orderOfSignature == x) %>%  dplyr::group_by(orderOfSignature,Diagnosis) %>% 
          dplyr::mutate(Med=median(Scores)) %>% arrange(orderOfSignature,Diagnosis,Scores) %>% 
          arrange(desc(Med)) %>% 
          ungroup() %>%
          mutate( Diagnosis = self$factorizeColumn(Diagnosis, unique(Diagnosis) ) ) %>% arrange(Diagnosis)
        #mutate( Diagnosis = factorizeColumn(Diagnosis, orderOfFactor )) %>% arrange(Diagnosis)
        tidyScores[,"SNONorm"] <- xaxisSeq(tidyScores)      
        
        ##Make median Segment
        medianY <- (tidyScores %>% dplyr::group_by(Diagnosis, orderOfSignature) %>% dplyr::summarise(medianY=median(Scores)) %>% dplyr::arrange(Diagnosis,orderOfSignature))$medianY
        medianX <- (tidyScores %>% dplyr::group_by(Diagnosis, orderOfSignature) %>% dplyr::summarise(medianX=median(SNONorm))  %>% dplyr::arrange(Diagnosis,orderOfSignature))$medianX
        segmentDF <- data.frame( xstart = medianX-0.05, ystart=medianY, xend=medianX+0.15, yend=medianY)
        segmentDF <- cbind(segmentDF, expand.grid(Diagnosis=unique(tidyScores$Diagnosis),orderOfSignature=unique(tidyScores$orderOfSignature)))
        
        summaryStats <- tidyScores %>% group_by(Diagnosis) %>% summarise(maxV = max(Scores), minV =min(Scores) )
        #View(summaryStats)
        plot <- ggplot() +
          geom_point(data=tidyScores, aes(SNONorm, Scores, colour = factor(Diagnosis) ),show.legend = F, size=sizeOfDots) +
          scale_colour_manual(values=customColors  ) +
          #geom_violin(data=tidyScores, aes(SNONorm, Scores, color = factor(Diagnosis)),show.legend = F)+
          facet_grid(orderOfSignature~Diagnosis ,switch = "both") +
          #ylim( min(summaryStats$minV)-0.05,max(summaryStats$maxV)+0.05) + 
          #ylim(-3,2.5) +
          labs( title= x ) +
          ylab( yLab ) +
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
                 ,panel.grid.major.y=element_blank()
                 ,panel.grid.minor.x=element_blank()
                 ,panel.border = element_rect(colour = "lightgrey", fill=NA, size=0.0000000002, linetype = 1)
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
        tidyScores             <- tidyScoresPre %>% filter(orderOfSignature == x) %>%  dplyr::group_by(orderOfSignature,Diagnosis) %>% 
          dplyr::mutate(Med=median(Scores)) %>% arrange(orderOfSignature,Diagnosis,Scores) %>% 
          arrange(desc(Med)) %>% 
          ungroup() %>%
          mutate( Diagnosis = self$factorizeColumn(Diagnosis, unique(Diagnosis) ) ) %>% arrange(Diagnosis)
        #mutate( Diagnosis = factorizeColumn(Diagnosis, orderOfFactor )) %>% arrange(Diagnosis)
        tidyScores[,"SNONorm"] <- xaxisSeq(tidyScores)
        
        ##Make median Segment
        medianY <- (tidyScores %>% dplyr::group_by(Diagnosis, orderOfSignature) %>% dplyr::summarise(medianY=median(Scores)) %>% dplyr::arrange(Diagnosis,orderOfSignature))$medianY
        medianX <- (tidyScores %>% dplyr::group_by(Diagnosis, orderOfSignature) %>% dplyr::summarise(medianX=median(SNONorm))  %>% dplyr::arrange(Diagnosis,orderOfSignature))$medianX
        segmentDF <- data.frame( xstart = medianX-0.05, ystart=medianY, xend=medianX+0.15, yend=medianY)
        segmentDF <- cbind(segmentDF, expand.grid(Diagnosis=unique(tidyScores$Diagnosis),orderOfSignature=unique(tidyScores$orderOfSignature)))
        
        summaryStats <- tidyScores %>% group_by(Diagnosis) %>% summarise(maxV = max(Scores), minV =min(Scores) )
        plot <- ggplot(data=tidyScores, aes(x = Scores, y = Diagnosis, height = ..density..)) +
          # to avoid overlaps of mountains , rel_min_height = 0.005, scale=0.9
          geom_density_ridges2(aes(fill = Diagnosis)) +
          scale_fill_manual(values=customColors, guide=FALSE) +
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
        Scores[,colList] <- apply(Scores[,colList,drop=FALSE] + 1 , 2, function(x) {return(logb(x,logBase))} )
      }
      
      if(standardize==TRUE){
        Scores[,colList] <- apply(Scores[,colList,drop=FALSE], 2, self$zscore_All )
      }
      
      tidyScoresPre <- Scores %>% tidyr::gather(orderOfSignature, Scores, colList);
      mergeDF       <-  merge(tidyScoresPre, customColorDF, by.x="Diagnosis", by.y="Diagnosis", all.x=TRUE)
      tidyScoresPre  <- mergeDF[,c(1:4)] ; # tidyScoresPre$Diagnosis <- factor(tidyScoresPre$Diagnosis, levels = orderOfFactor, ordered = TRUE)
      
      ## Set customColors
      customColorsVector <- setNames( as.character(customColorDF$Color), as.character(customColorDF$Diagnosis))
      
      if( plotType =="StringBean") {
        plotLists <- lapply(orderOfSignature, drawStringBeanPlot, tidyScoresPre, customColors=customColorsVector,yLab =yLab)
      } else {
        plotLists <- lapply(orderOfSignature, drawDensityPlot, tidyScoresPre=tidyScoresPre, orderOfFactor=orderOfFactor, customColors=customColorsVector,
                            yLab =yLab)
      }
      return(plotLists)
    },
    OneVariablePlotSortMean = function(colList=NA, Scores=NA, orderOfFactor=NA, orderOfSignature=NA, standardize=FALSE, logit =FALSE,
                                   plotType="StringBean",customColorDF=NA, yLab="Score", summaryHlines =FALSE, 
                                   sizeOfDots = 1, legendDisplay=TRUE, logBase=2){
      
      #if (unique(is.na(customColors))) { customColors = setNames( StatsFinal$Color, StatsFinal$Diagnosis) }
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
      drawStringBeanPlot <- function(x, tidyScoresPre, customColors=NA, yLab =yLab){
        
        print(paste(x))
        tidyScores             <- tidyScoresPre %>% filter(orderOfSignature == x) %>%  dplyr::group_by(orderOfSignature,Diagnosis) %>% 
          dplyr::mutate(Med=mean(Scores)) %>% arrange(orderOfSignature,Diagnosis,Scores) %>% 
          arrange(desc(Med)) %>% 
          ungroup() %>%
          mutate( Diagnosis = self$factorizeColumn(Diagnosis, unique(Diagnosis) ) ) %>% arrange(Diagnosis)
        #mutate( Diagnosis = factorizeColumn(Diagnosis, orderOfFactor )) %>% arrange(Diagnosis)
        tidyScores[,"SNONorm"] <- xaxisSeq(tidyScores)      
        
        ##Make mean Segment
        meanY <- (tidyScores %>% dplyr::group_by(Diagnosis, orderOfSignature) %>% dplyr::summarise(meanY=mean(Scores)) %>% dplyr::arrange(Diagnosis,orderOfSignature))$meanY
        meanX <- (tidyScores %>% dplyr::group_by(Diagnosis, orderOfSignature) %>% dplyr::summarise(meanX=mean(SNONorm))  %>% dplyr::arrange(Diagnosis,orderOfSignature))$meanX
        segmentDF <- data.frame( xstart = meanX-0.05, ystart=meanY, xend=meanX+0.15, yend=meanY)
        segmentDF <- cbind(segmentDF, expand.grid(Diagnosis=unique(tidyScores$Diagnosis),orderOfSignature=unique(tidyScores$orderOfSignature)))
        
        summaryStats <- tidyScores %>% group_by(Diagnosis) %>% summarise(maxV = max(Scores), minV =min(Scores) )
        scoreSummary <- summary(tidyScores$Scores)
        
        plot <- ggplot() +
          geom_point(data=tidyScores, aes(SNONorm, Scores, colour = factor(Diagnosis) ),show.legend = F, size=sizeOfDots) +
          scale_colour_manual(values=customColors  ) +
          #geom_violin(data=tidyScores, aes(SNONorm, Scores, color = factor(Diagnosis)),show.legend = F)+
          facet_grid(orderOfSignature~Diagnosis ,switch = "both") +
          #ylim( min(summaryStats$minV)-0.05,max(summaryStats$maxV)+0.05) + 
          #ylim(-3,2.5) +
          labs( title= x ) +
          ylab( yLab ) +
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
                 ,panel.grid.major.y=element_blank()
                 ,panel.grid.minor.x=element_blank()
                 ,panel.border = element_rect(colour = "lightgrey", fill=NA, size=0.0000000002, linetype = 1)
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
        tidyScores             <- tidyScoresPre %>% filter(orderOfSignature == x) %>%  dplyr::group_by(orderOfSignature,Diagnosis) %>% 
          dplyr::mutate(Med=mean(Scores)) %>% arrange(orderOfSignature,Diagnosis,Scores) %>% 
          arrange(desc(Med)) %>% 
          ungroup() %>%
          mutate( Diagnosis = self$factorizeColumn(Diagnosis, unique(Diagnosis) ) ) %>% arrange(Diagnosis)
        #mutate( Diagnosis = factorizeColumn(Diagnosis, orderOfFactor )) %>% arrange(Diagnosis)
        tidyScores[,"SNONorm"] <- xaxisSeq(tidyScores)
        
        ##Make mean Segment
        meanY <- (tidyScores %>% dplyr::group_by(Diagnosis, orderOfSignature) %>% dplyr::summarise(meanY=mean(Scores)) %>% dplyr::arrange(Diagnosis,orderOfSignature))$meanY
        meanX <- (tidyScores %>% dplyr::group_by(Diagnosis, orderOfSignature) %>% dplyr::summarise(meanX=mean(SNONorm))  %>% dplyr::arrange(Diagnosis,orderOfSignature))$meanX
        segmentDF <- data.frame( xstart = meanX-0.05, ystart=meanY, xend=meanX+0.15, yend=meanY)
        segmentDF <- cbind(segmentDF, expand.grid(Diagnosis=unique(tidyScores$Diagnosis),orderOfSignature=unique(tidyScores$orderOfSignature)))
        
        summaryStats <- tidyScores %>% group_by(Diagnosis) %>% summarise(maxV = max(Scores), minV =min(Scores) )
        
        plot <- ggplot(data=tidyScores, aes(x = Scores, y = Diagnosis, height = ..density..)) +
          # to avoid overlaps of mountains , rel_min_height = 0.005, scale=0.9
          geom_density_ridges2(aes(fill = Diagnosis)) +
          scale_fill_manual(values=customColors, guide=FALSE) +
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
        Scores[,colList] <- apply(Scores[,colList,drop=FALSE] + 1 , 2, function(x) {return(logb(x,logBase))} )
      }
      
      if(standardize==TRUE){
        Scores[,colList] <- apply(Scores[,colList,drop=FALSE], 2, self$zscore_All )
      }
      
      tidyScoresPre <- Scores %>% tidyr::gather(orderOfSignature, Scores, colList);
      mergeDF       <-  merge(tidyScoresPre, customColorDF, by.x="Diagnosis", by.y="Diagnosis", all.x=TRUE)
      tidyScoresPre  <- mergeDF[,c(1:4)] ; # tidyScoresPre$Diagnosis <- factor(tidyScoresPre$Diagnosis, levels = orderOfFactor, ordered = TRUE)
      
      ## Set customColors
      customColorsVector <- setNames( as.character(customColorDF$Color), as.character(customColorDF$Diagnosis))
      
      if( plotType =="StringBean") {
        plotLists <- lapply(orderOfSignature, drawStringBeanPlot, tidyScoresPre, customColors=customColorsVector,yLab =yLab)
      } else {
        plotLists <- lapply(orderOfSignature, drawDensityPlot, tidyScoresPre=tidyScoresPre, orderOfFactor=orderOfFactor, customColors=customColorsVector,
                            yLab =yLab)
      }
      return(plotLists)
    },
    ## memo sort
    memoSort = function(M = NA) {
      geneOrder <- sort(rowSums(M), decreasing=TRUE, index.return=TRUE)$ix;
      scoreCol <- function(x) {
        score <- 0;
        for(i in 1:length(x)) {
          if(x[i]) {
            score <- score + 2^(length(x)-i);
          }
        }
        return(score);
      }
      scores <- apply(M[geneOrder, ], 2, scoreCol);
      sampleOrder <- sort(scores, decreasing=TRUE, index.return=TRUE)$ix;
      return(M[geneOrder, sampleOrder]);
    },
    ## TCR analysis functions
    ## Filter specific clones
    filterSpecificCloneTypes  = function(cloneData, cloneType){
      cloneDataFilt           <- cloneData %>% dplyr::filter(grepl(cloneType,v) | grepl("NF",v))
      return(cloneDataFilt)
    },
    ## make correlation plots for immune scores and clone count
    correlationPlots = function(varName="", constName="", df=NA){
      
      print(paste(varName))
      # plot <- ggscatter(df, x = constName, y =varName, 
      #                   add = "reg.line", conf.int = TRUE, 
      #                   cor.coef = TRUE, cor.method = "spearman",
      #                   xlab = constName, ylab = varName)
      
      corrTest <- cor.test(df[,constName], df[,varName], method = "spearman")
      ## if  ( corrTest$p.value < 2.2e-16 ) { corrTest$p.value = 2.2e-16 }
      plot <- ggplot(df, aes_string(x=constName, y=varName)) + 
        geom_smooth(method=lm,  fill="grey") +
        geom_point(aes(colour = factor(Diagnosis)), show.legend = T, size=3, shape=16) + 
        scale_colour_manual(values=setNames( StatsFinal$Color, StatsFinal$Diagnosis)  ) +
        theme_bw() +
        theme(axis.text=element_text(size=13)
              ,axis.title=element_text(size=13,face="bold")) +
        xlab("Log Total Clones") +
        ylab(paste("Standardised Enrichment Score", sep=" "))+
        ggtitle(paste("Corr.Coeff = ", signif(corrTest$estimate[[1]],5), "\np-value = ", signif(corrTest$p.value,5), "                                            ", varName,sep=""))
      
      return(list(plot))
    }   ,
    ## Prepare Input for entropy and clonality
    makeEntropyInput = function(filename, inputDir="", outputDir="") {
      outfileName <- paste0(outputDir, gsub("convert.|.clones.txt","",filename ), ".Entropy.txt")
      exomeData <- read.csv( paste(inputDir, filename, sep=""), sep="\t", header = TRUE )
      if(nrow(exomeData)>0){
        exomeDataEntropy <- data.frame(VJcombo=paste(exomeData$v,exomeData$j,sep="."), Counts=exomeData$count, Vcassette=exomeData$v, 
                                       Jcassette=exomeData$j, aaCDR3_filtered = exomeData$cdr3aa, ntCDR3= exomeData$cdr3nt)
      } else {
        exomeDataEntropy <- emptyDFEntropy
      }
      write.table(exomeDataEntropy, outfileName, sep = "\t", row.names = FALSE, quote = FALSE)
    },
    ## Convert immunoseq data to compatible file format
    immunoseqv2 = function(x) {
      filename = x
      outfileName <- paste0("./immunoseqv2EntropyNitin/", gsub("convert.|.clones.txt","",filename ), ".Entropy.txt")
      exomeData <- read.csv( paste("./immunoseqv2/", filename, sep=""), sep="\t", header = TRUE, stringsAsFactors = FALSE )
      exomeDataFilt <- exomeData[,1:47] %>% dplyr::filter(sequenceStatus %in% c("In"))
      colnames(exomeDataFilt)[3] <- "count..templates.reads."
      
      print(paste("./immunoseqv2/", filename, sep=""))
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
      write.table(exomeDataEntropy, outfileName, sep = "\t", row.names = FALSE, quote = FALSE)
      
    },
    ## Makepreranke files for GSEA/ssGSEA
    getCountObjTXT = function(fileName, colNumb=1, rowNames=1){
      print(paste(fileName))
      featureCountTxt <- read.table(fileName, sep="\t", row.names = rowNames, header = 1);
      return(featureCountTxt[,colNumb, drop=FALSE])
    },
    NESorPvalGSEAPrerank = function(x, colNumb = 1){
      
      files <- list.files(x, pattern = '^gsea_report_for_na_(pos|neg)_[1234567890]+\\.xls$', full.names = TRUE)
      df <- do.call(rbind, lapply(files, self$getCountObjTXT, colNumb=colNumb, rowNames=1)); print(paste(dim(df)))
      colnames(df) <- gsub("\\.\\/GSEA\\/prerankedGSEAOutput\\/\\/|.rnk.*.\\S+\\/gsea_report_for_na_(pos|neg)_\\S+.xls", "", files)[1]
      return(df)
    },
    makeRNKFiles = function(DF, projectName=""){
      colNames <- colnames(DF)
      sapply(colNames[1], function(x){ 
        
        y <- DF[, x, drop=FALSE] %>% rownames_to_column("Name") %>% arrange_(.dots = paste0("desc(",x,")") )
        write.table(y, paste("./GSEA/rnk/", projectName,"/",x,".rnk",sep=""), sep = "\t",row.names = FALSE, quote = FALSE )})
    },
    PreRankedGSEA = function(x, projectName=NULL, sample=NULL){
      
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
    },
    flattenCorrMatrix = function(cormat, pmat) {
      ut <- upper.tri(cormat)
      data.frame(
        row = rownames(cormat)[row(cormat)[ut]],
        column = rownames(cormat)[col(cormat)[ut]],
        cor  =(cormat)[ut],
        p = pmat[ut]
      )
    },
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
    },
    validfMatrix = function(df=NA){
      tempDF <- as.data.frame(lapply(df, make.names)) %>% data.frame()
      return(tempDF)
    }
  )
)

# ### GeneExpression normalization using multiple tools ####
# "keepProteinCoding(), fpkmToTpm()" function in another class

GeneExpNormalization <- R6Class(
  classname = "GeneExpNormalization",
  inherit   = CoreUtilities, 
  portable  = TRUE,
  private   = list(
    GeneDF               = NULL,
    GeneDFNorm           = NULL,
    miniAnnotationDF     = NULL,
    countObj             = NULL,
    featureType          = NULL,
    packageRNAseq        = NULL,
    annotationDF         = NULL,
    design               = NULL,
    proteinCodingOnly    = NULL,
    corUtilsFuncs        = NULL,
    saveFiles            = NULL,
    setUp                = function(){
      
      if( private$featureType == "Exon")      { private$featureType = "ExonID"       }
      if( private$featureType == "Transcript"){ private$featureType = "TranscriptID" }
      if( private$featureType == "Gene")      { private$featureType = "GeneID"       }
      
      private$annotationDF %<>% mutate(!!private$featureType := super$factorizeColumn(private$annotationDF[[private$featureType]], 
                                                                                      row.names(private$countObj))) %>% 
        dplyr::arrange_(.dots=c(private$featureType))
      
      if(private$proteinCodingOnly) {
        keepProteinCodingOut     <- super$keepProteinCoding(geneMatrix = private$countObj, 
                                                            annotation=private$annotationDF, featureType=private$featureType)
        private$countObj         <- keepProteinCodingOut[[1]]
        private$annotationDF     <- keepProteinCodingOut[[2]]
        print(paste0(" dim of gene exp ", dim(private$countObj) ))
        print(paste0(" dim of annotation ", dim(private$annotationDF) ))
      }
      
      private$miniAnnotationDF <- private$annotationDF[,c(private$featureType, "Length")]
      colnames(private$miniAnnotationDF) <- c("GeneID", "Length")
      
      switch(private$packageRNAseq,
             
             "edgeR"= {
               
               ##  Make EdgeR Object
               private$GeneDF                  <- DGEList(counts=private$countObj, genes=private$miniAnnotationDF, 
                                                          group = as.factor(as.character(private$design)) )
               ## Estimate Normalising Factors
               private$GeneDFNorm              <- calcNormFactors(private$GeneDF)
               
             },
             "deseq2"= {
               
               ##  Make deseq2 Object
               condition                       <- private$design
               private$GeneDF                  <- DESeqDataSetFromMatrix(countData = private$countObj, colData = DataFrame(condition), design = ~condition)
               mcols(private$GeneDF)$basepairs <- private$miniAnnotationDF$Length-174
               ## Estimate SizeFactors
               private$GeneDFNorm              <- estimateSizeFactors(private$GeneDF)
             }
      )
      
    }
  ),
  public    = list(
    
    initialize           = function(countObj = NA, featureType = NA, packageRNAseq = NA, annotationDF = NA, design = NA,
                                    proteinCodingOnly = TRUE, corUtilsFuncs = NA, saveFiles=FALSE ) 
    {
      
      assert_that( class(countObj)[1] == "matrix", msg = "Please provide raw count in matrix format")
      assert_that( nrow(countObj) == nrow(annotationDF), msg = "Mismatch  of gene/transcript entries in between countObj vs annotationDF " )
      assert_that( is.logical(proteinCodingOnly) , msg = "proteinCodingOnly parameter should be boolean")
      
      private$countObj <- countObj
      private$featureType <- featureType
      private$packageRNAseq <- packageRNAseq
      private$annotationDF  <- annotationDF
      private$design        <- design
      private$proteinCodingOnly <- proteinCodingOnly
      private$corUtilsFuncs  <- corUtilsFuncs
      private$saveFiles      <- saveFiles
      private$setUp()
    },
    edgeRMethod      = function(x, logtransform=FALSE, zscore=FALSE) {
      assert_that(private$packageRNAseq == "edgeR", msg = paste0("GeneExpNormalization Object was created for ",private$packageRNAseq,
                                                                 ". Please create a new GeneExpNormalization Object for edgeR"))
      
      assert_that(x %in% c("CPM", "TMM-RPKM", "TPM", "NormFactorDF", "RawCounts"), msg = "This function can only generate \"CPM\", \"TMM-RPKM\", \"TPM\" values ")
      
      if(x == "NormFactorDF") return(private$GeneDFNorm)
      if(x == "RawCounts") { 
        rawCounts <-  private$GeneDFNorm$counts %>% data.frame() %>% tibble::rownames_to_column(var="GeneID")
        return( private$corUtilsFuncs$featureNameAnot(querryDF=rawCounts, identifier="GeneID", annotationDF=private$annotationDF)  ) 
      }
      if(x == "CPM" )         { 
        cpmDF <- as.data.frame(cpm(private$GeneDFNorm,  normalized.lib.sizes = TRUE,log = FALSE)) %>% tibble::rownames_to_column(var="GeneID")
        return( private$corUtilsFuncs$featureNameAnot(querryDF=cpmDF, identifier="GeneID", annotationDF=private$annotationDF) )  
      }
      if(x == "TMM-RPKM" )    { 
        rpkmDF =as.data.frame(rpkm(private$GeneDFNorm, normalized.lib.sizes = TRUE, log = FALSE))
        if (isTRUE(logtransform) & isTRUE(zscore)) {
          
          print("Log transforming and standardising")
          rpkmDF <- log2(rpkmDF+1) 
          rpkmDF <- apply(rpkmDF, 1, function(x){
            medX <- median(x)
            sdX <- sd(x, na.rm = FALSE)
            y <- (x - medX) / sdX
            return(y)
          }) %>% t() %>% data.frame()
          
        } else if ( isTRUE(logtransform) &  !isTRUE(zscore)  ){
          
          print("only Log transforming ")
          rpkmDF <- log2(rpkmDF+1) 
          
        }
        rpkmDF <- rpkmDF %>% tibble::rownames_to_column(var="GeneID")
        
        return( private$corUtilsFuncs$featureNameAnot(querryDF=rpkmDF, identifier="GeneID", annotationDF=private$annotationDF) ) 
      }
      if(x == "TPM" )         { 
        tpmDF <- apply(rpkm(private$GeneDFNorm, normalized.lib.sizes = TRUE), 2 , super$fpkmToTpm) %>% tibble::rownames_to_column(var="GeneID")
        return( private$corUtilsFuncs$featureNameAnot(querryDF=tpmDF, identifier="GeneID", annotationDF=private$annotationDF) )  
      }
    },
    deseq2           = function(x) {
      
      assert_that(private$packageRNAseq == "deseq2", msg = paste0("GeneExpNormalization Object was created for ",private$packageRNAseq,
                                                                  ". Please create a new GeneExpNormalization Object for deseq2"))
      
      assert_that(x %in% c("FPM", "TMM-RPKM", "TPM", "RLOG", "VST"), 
                  msg = "This function can only generate \"CPM\", \"TMM-RPKM\", \"TPM\", \"RLOG\", \"VST\" values ")
      
      
      if(x == "NormFactorDF") return(private$GeneDFNorm)
      if(x == "FPM" )         return( as.data.frame(fpm(object = private$GeneDFNorm,  robust = TRUE))    )
      if(x == "TMM-RPKM" )   return( as.data.frame(fpkm(object = private$GeneDFNorm,  robust = TRUE))   )
      if(x == "TPM" )         return( apply(as.data.frame(fpkm(object = private$GeneDFNorm,  robust = TRUE)), 2 , fpkmToTpm))
      if(x == "RLOG" )        return( rlog(private$GeneDFNorm, blind = TRUE)              )
      if(x == "VST" )         return( vst(private$GeneDFNorm, blind = TRUE))
    }
  )
)

# ### Differntial Gene Expression Analysis ####

DifferentialGeneExp <- R6Class(
  classname = "DifferentialGeneExp",
  portable  = TRUE,
  private   = list
  (
  
  countObj            = NULL,
  metadataDF          = NULL,
  group1Flatten       = NULL,
  group2Flatten       = NULL,
  GeneDF_DiffExp      = NULL,
  pairGroup1Name      = NULL,
  group1Samples       = NULL,
  pairGroup2Name      = NULL,
  group2Samples       = NULL,
  fileDirs            = NULL,
  corUtilsFuncs        = NULL,
  setUP               = function(x) {},
  changeGroupName     = function(vectorOfNames, toChangeName, inVectorName){
    return(gsub( paste0(paste0("^",vectorOfNames) %>% paste0(.,"$"), collapse = "|"),toChangeName,inVectorName))
  },
  flattenGroup        = function(x){
    unlist(
      lapply(x, function(y){  
        if(y$each){
          return( as.list(sapply(unlist(unname(y[1])), function(y){ return(y) }) ) )
        } else {
          return(y[1])
        }
      }), recursive = FALSE)
  },
  makePairs           = function(x, y){
    pairs = unlist(lapply(x, function(X) {
      lapply(y, function(Y) {
        c(X, Y)
      })
    }), recursive=FALSE)
    return(pairs)
  },
  DiffGeneExp         = function( group1Count = NA, group2Count = NA, pairGroup1Name = NA,pairGroup2Name = NA, DGEobj = NA  ){
    
    ## Use the appropriate method depending upon the size of the replicates
    if( group1Count > 1 & group2Count > 1 )   {
      
      print("Groups have replicates")
      ## Generate model & design for differential gene expression (edgeR only for now
      
      #Using exactTest()
      modelGroup   <- factor( c( rep(pairGroup1Name, group1Count),rep(pairGroup2Name, group2Count)) )
      modelDesign  <- model.matrix( ~modelGroup )
      ## Estimate Dispersion using estimateDisp() 
      GeneDF_Dispersion  <- estimateDisp(DGEobj, design = modelDesign )
      print("Predicting differntially expressed genes using exactTest()")
      print(paste0(pairGroup1Name,"  ",group1Count, "  ", pairGroup2Name,"  ",group2Count))
      private$GeneDF_DiffExp <- exactTest(GeneDF_Dispersion, pair = c(unique(modelGroup)))$table
      
      # # Using Quasilikelihood ratio test() **************** ISSUE groups getting sorted ************
      # modelGroup   <- factor( c( rep(pairGroup1Name, group1Count), rep(pairGroup2Name, group2Count)) )
      # modelDesign  <- model.matrix( ~modelGroup )
      # print("Predicting differntially expressed genes using Quasi Likelihood ratio test glmQLFit() & glmQLFTest()")
      # GeneDF_Dispersion  <- estimateDisp(DGEobj, design = modelDesign )
      # fit                        <- glmQLFit(GeneDF_Dispersion, design = modelDesign )
      # private$GeneDF_DiffExp     <- glmQLFTest(fit, coef=2)$table
      
      #print("Predicting differntially expressed genes using LimmaVoom")
      # GeneDF_Dispersion_v   <- voom(DGEobj, modelDesign)
      # fit_v                 <- lmFit(GeneDF_Dispersion_v, modelDesign)
      # fit_v                 <- eBayes(fit_v)
      # modelGroup_v   <- factor( c(rep(pairGroup2Name, group2Count), rep(pairGroup1Name, group1Count)) )
      # modelDesign_v  <- model.matrix( ~modelGroup_v )
      # GeneDF_Dispersion_vwts  <- voomWithQualityWeights(DGEobj, design=modelDesign_v, normalization="none", plot=FALSE)
      # fit_vwts                <- eBayes(lmFit(GeneDF_Dispersion_vwts,design=modelDesign_v))
      # GeneDF_DiffExp_V        <- topTable(fit_vwts,coef=2,number=length(fit_vwts$genes[,1]),sort.by="none")
      # print(head(GeneDF_DiffExp_V))
    }
    else  {
      ## Generate model & design for differential gene expression (edgeR only for now)
      modelGroup                 <- factor(c(rep(pairGroup1Name, group1Count), rep( pairGroup2Name, group2Count   )))
      modelDesign                <- model.matrix( ~modelGroup )
      ## Estimate Dispersion using estimateDisp() 
      print("Predicting differntially expressed genes using  estimateGLMCommonDisp() and exact test ")
      GeneDF_Dispersion          <- estimateGLMCommonDisp(DGEobj, method="deviance",robust="TRUE",subset=NULL )
      private$GeneDF_DiffExp     <- exactTest(GeneDF_Dispersion, pair = c(unique(modelGroup)))$table
    }
    
    private$GeneDF_DiffExp["FDR"]   <- p.adjust(private$GeneDF_DiffExp$PValue, method="BH")
    private$GeneDF_DiffExp          <- tibble::add_column(private$GeneDF_DiffExp, GeneID = rownames(private$GeneDF_DiffExp), .before = 1)
    
    return(private$GeneDF_DiffExp)
  },
  MeanGroupExpression = function(group1Count = NA, group2Count = NA, geneExpMatrix = NA){
    
    ## GeneExpression Mean
    if(group1Count > 1) {
      private$GeneDF_DiffExp[private$pairGroup1Name]        <- apply(geneExpMatrix[,as.character(private$group1Samples[,1]), drop=FALSE], 1, mean)
      #group1MeanExp                                        <- apply(geneExpMatrix[,as.character(private$group1Samples[,1]), drop=FALSE], 1, mean)
      #print(paste(head(group1MeanExp)))
    } else {
      private$GeneDF_DiffExp[private$pairGroup1Name]        <- geneExpMatrix[,as.character(private$group1Samples[,1]), drop=FALSE]
    }
    
    if(group2Count > 1) {
      private$GeneDF_DiffExp[private$pairGroup2Name]        <- apply(geneExpMatrix[,as.character(private$group2Samples[,1]), drop=FALSE], 1, mean)
      #group2MeanExp                                        <- apply(geneExpMatrix[,as.character(private$group2Samples[,1]), drop=FALSE], 1, mean)
      #print(paste(head(group2MeanExp)))
      
    } else {
      private$GeneDF_DiffExp[private$pairGroup2Name]        <- geneExpMatrix[,as.character(private$group2Samples[,1]), drop=FALSE]
    }
    
    private$GeneDF_DiffExp["FDR"]   <- p.adjust(private$GeneDF_DiffExp$PValue, method="BH")
    private$GeneDF_DiffExp["myLogFC"]   <- log2(private$GeneDF_DiffExp[private$pairGroup2Name]+1) - log2(private$GeneDF_DiffExp[private$pairGroup1Name]+1)
    return(private$GeneDF_DiffExp)
  },
  executeDiffGeneExp  = function( pairedList = NA){
    
    ## Get the group names
    pairGroup1List <- private$group1Flatten[  pairedList[1] ]
    private$pairGroup1Name <- names(pairGroup1List)
    pairGroup1 <- unname(unlist(pairGroup1List))
    
    pairGroup2List <- private$group2Flatten[  pairedList[2] ]
    private$pairGroup2Name <- names(pairGroup2List)
    pairGroup2 <- unname(unlist(pairGroup2List))
    
    ## Make Custom Filters
    group1.group2.Filt   <- interp(~y %in% x, .values=list(y = as.name(self$groupColumnName), x = c(pairGroup1, pairGroup2) ))
    
    factorizeGroupColumn <- interp( ~factor(groupColumnName, levels = val, ordered = TRUE), 
                                    groupColumnName=as.name(self$groupColumnName), val=c(pairGroup1, pairGroup2) )
    
    ## Filter the metadata table for group1 and group2
    subSetSamples <- private$metadataDF %>% 
      filter_(group1.group2.Filt) %>% 
      mutate_(.dots = setNames( list(factorizeGroupColumn) , self$groupColumnName )) %>%
      dplyr::arrange_(.dots = self$groupColumnName)
    
    private$group1Samples <- subSetSamples %>% filter_(interp(~y %in% x, .values=list(y = as.name(self$groupColumnName), x = pairGroup1))) %>%
      select_(.dots=list(self$samplesColumnName))
    private$group2Samples <- subSetSamples %>% filter_(interp(~y %in% x, .values=list(y = as.name(self$groupColumnName), x = pairGroup2))) %>%
      select_(.dots=list(self$samplesColumnName))
    
    ## Subset Count matrix
    pairCountMatrix <- private$countObj %>% data.frame() %>% dplyr::select_(.dots=as.character(subSetSamples[,self$samplesColumnName]))
    
    ## Getting count of each group
    groupCount  <- table(subSetSamples[,self$groupColumnName])
    group1Count <- sum(groupCount[pairGroup1])
    group2Count <- sum(groupCount[pairGroup2])
    
    ## Generate model & design for differential gene expression (edgeR only for now)
    modelGroup   <- factor(c(rep(private$pairGroup1Name, group1Count), rep( private$pairGroup2Name, group2Count   )))
    modelDesign  <- model.matrix( ~modelGroup )
    
    ## Generate Gene expression in the desired Unit.
    expressionObj <- GeneExpNormalization$new(
      countObj          = as.matrix(pairCountMatrix), 
      featureType       = self$featureType, 
      packageRNAseq     = self$packageRNAseq, 
      annotationDF      = rnaseqProject$annotationDF, 
      design            = modelDesign[,2], 
      proteinCodingOnly = FALSE,
      corUtilsFuncs     = private$corUtilsFuncs
    )
    geneExpression = expressionObj$edgeRMethod(self$expressionUnit)
    
    ## perform differential gene expression.
    private$GeneDF_DiffExp <- private$DiffGeneExp(group1Count = group1Count, group2Count = group2Count, DGEobj = expressionObj$edgeRMethod("NormFactorDF"), 
                                                  pairGroup1Name = private$pairGroup1Name, pairGroup2Name = private$pairGroup2Name)
    
    ## calculat mean expression
    private$GeneDF_DiffExp                <- private$MeanGroupExpression(group1Count = group1Count, group2Count = group2Count, geneExpMatrix = geneExpression)
    private$GeneDF_DiffExp["AvglogFPKM"]  <- apply(geneExpression[,-c(1:7)] , 1, mean)
    
    ## Annotate Genes
    private$GeneDF_DiffExp                <- corUtilsFuncs$featureNameAnot(querryDF=private$GeneDF_DiffExp, identifier="GeneID", annotationDF = rnaseqProject$annotationDF)
    
    ## Filter Genes
    folderName <- paste0(private$pairGroup1Name, "_", private$pairGroup2Name)
    private$filterGenes(filterName="all", folderName = folderName )
    if( self$subsetGenes ){
      private$filterGenes(filterName="ProteinCoding"         , folderName = folderName )
      private$filterGenes(filterName="CellSurface"           , folderName = folderName )
      private$filterGenes(filterName="TranscriptionFactor"   , folderName = folderName )
      private$filterGenes(filterName="CancerGermlineAntigen" , folderName = folderName )
      private$filterGenes(filterName="ExhaustionMarkers" , folderName = folderName )
    }
    return(private$GeneDF_DiffExp)
  },
  appendAnnotation = function(df=NA, annottype=NA){
    GeneDF_DiffExp <- df %>% mutate(
      ProteinCoding = ifelse(GeneName %in% rnaseqProject$pcDF[,"GeneName"], "Y", "N"),
      CellSurface   = ifelse(GeneName %in% rnaseqProject$csDF[,"GeneName"], "Y", "N"),
      ExhaustionMarkers   = ifelse(GeneName %in% rnaseqProject$emDF[,"GeneName"], "Y", "N"),
      TranscriptionFactor     = ifelse(GeneName %in% rnaseqProject$tfDF[,"GeneName"], "Y", "N"),
      CancerGermlineAntigen   = ifelse(GeneName %in% rnaseqProject$cgaDF[,"GeneName"], "Y", "N"),
      PAX3FOXO1     = ifelse(GeneName %in% rnaseqProject$pax3Foxo1DF[,"GeneName"], "Y", "N"),
      EWSR1FL1      = ifelse(GeneName %in% rnaseqProject$ewsr1Fli1DF[,"GeneName"], "Y", "N")
    )
    GeneDF_DiffExp <- Reduce(function(x,y) merge(x,y,by=c("GeneID", "GeneName"),all=TRUE) ,list(GeneDF_DiffExp,
                                                                                                rnaseqProject$BrainExpDF[  which(GeneDF_DiffExp$GeneName %in% rnaseqProject$BrainExpDF$GeneName ), c("GeneID", "GeneName", "Brain.MeanExp") ],
                                                                                                rnaseqProject$HeartExpDF[  which(GeneDF_DiffExp$GeneName %in% rnaseqProject$HeartExpDF$GeneName ), c("GeneID", "GeneName", "Heart.MeanExp") ],
                                                                                                rnaseqProject$KidneyExpDF[ which(GeneDF_DiffExp$GeneName %in% rnaseqProject$KidneyExpDF$GeneName), c("GeneID", "GeneName", "Kidney.MeanExp")],
                                                                                                rnaseqProject$LiverExpDF[  which(GeneDF_DiffExp$GeneName %in% rnaseqProject$LiverExpDF$GeneName ), c("GeneID", "GeneName", "Liver.MeanExp") ], 
                                                                                                rnaseqProject$LungExpDF[   which(GeneDF_DiffExp$GeneName %in% rnaseqProject$LungExpDF$GeneName  ), c("GeneID", "GeneName", "Lung.MeanExp")  ]
    ))
    
    if(!is.na(annottype) & !annottype %in% c("all")) {
      GeneDF_DiffExp <- GeneDF_DiffExp %>% dplyr::filter_(.dots=paste0(annottype, " == \"Y\""))
    } 
    
    return(GeneDF_DiffExp)
  },
  ## Annotate a gene expression df with "GeneID" as primary key
  # featureNameAnot = function(querryDF=NA, identifier=NA){
  #   annotDF <- dplyr::left_join(rnaseqProject$annotationDF, querryDF, by=identifier)
  #   print(paste("Annotating expression DF"))
  #   print(dim(annotDF))
  #   return(annotDF)
  # },
  filterGenes = function(filterName= NA , folderName = NA){
    
    rdsDir  <- paste0( private$fileDirs[4],"/", folderName, "/" )
    txtDir  <- paste0( private$fileDirs[3],"/", folderName, "/" )
    
    if(!dir.exists(rdsDir)) { dir.create(rdsDir) }
    if(!dir.exists(txtDir)) { dir.create(txtDir) }
    
    rdsfile <- paste0(rdsDir, paste0(folderName,"_",filterName,".rds") )
    txtFile <- paste0(txtDir, paste0(folderName,"_",filterName,".txt") )
    
    switch(filterName,
           
           "ProteinCoding"= {
             GeneDF_DiffExp <- dplyr::left_join(private$GeneDF_DiffExp, rnaseqProject$pcDF, by=c("GeneID","GeneName"))
             GeneDF_DiffExp <- private$appendAnnotation(df=GeneDF_DiffExp, annottype ="ProteinCoding" )
           },
           "CellSurface"= {
             GeneDF_DiffExp <- dplyr::left_join(private$GeneDF_DiffExp, rnaseqProject$pcDF,by=c("GeneID","GeneName"))
             GeneDF_DiffExp <- private$appendAnnotation(df=GeneDF_DiffExp, annottype ="CellSurface" )
           },
           "TranscriptionFactor"= {
             GeneDF_DiffExp <- dplyr::left_join(private$GeneDF_DiffExp, rnaseqProject$pcDF,by=c("GeneID","GeneName"))
             GeneDF_DiffExp <- private$appendAnnotation(df=GeneDF_DiffExp, annottype ="TranscriptionFactor" )
           },
           "CancerGermlineAntigen"= {
             GeneDF_DiffExp <- dplyr::left_join(private$GeneDF_DiffExp, rnaseqProject$pcDF, by=c("GeneID","GeneName"))
             GeneDF_DiffExp <- private$appendAnnotation(df=GeneDF_DiffExp, annottype ="CancerGermlineAntigen" )
           },
           "ExhaustionMarkers"= {
             GeneDF_DiffExp <- dplyr::left_join(private$GeneDF_DiffExp, rnaseqProject$pcDF, by=c("GeneID","GeneName"))
             GeneDF_DiffExp <- private$appendAnnotation(df=GeneDF_DiffExp, annottype ="ExhaustionMarkers" )
           },
           "all" = {
             GeneDF_DiffExp <- private$GeneDF_DiffExp
             GeneDF_DiffExp <- dplyr::left_join(GeneDF_DiffExp, rnaseqProject$pcDF, by=c("GeneID","GeneName"))
             GeneDF_DiffExp <- private$appendAnnotation(df=GeneDF_DiffExp, annottype ="all" )
           }
    )
    saveRDS(GeneDF_DiffExp, rdsfile)
    write.table(x = GeneDF_DiffExp, file = txtFile, sep="\t", row.names = FALSE, quote=FALSE)
  }
  ),
  public    = list
  (
  
  group1            = NULL,
  group2            = NULL,
  packageRNAseq     = NULL,
  factorsExclude    = NULL,
  groupColumnName   = NULL,
  samplesColumnName = NULL,
  pairedList        = NULL,
  diffGeneExpList   = NULL,
  expressionUnit    = NULL,
  featureType       = NULL,
  writeFiles        = NULL,
  subsetGenes       = NULL,
  initialize        = function(countObj = NA, metadataDF = NA, packageRNAseq = NA, expressionUnit = NA, featureType = NA,
                               group1   = NA, group2   = NA, groupColumnName = NA, samplesColumnName = NA,
                               writeFiles = FALSE, fileDirs = NA, corUtilsFuncs = NA , subsetGenes= TRUE ){
    
    private$countObj           <- countObj
    private$metadataDF         <- metadataDF
    self$expressionUnit        <- expressionUnit
    self$packageRNAseq         <- packageRNAseq
    self$featureType           <- featureType
    self$group1                <- group1
    self$group2                <- group2
    self$groupColumnName       <- groupColumnName
    self$samplesColumnName     <- samplesColumnName
    self$writeFiles            <- writeFiles
    self$subsetGenes           <- subsetGenes
    private$fileDirs           <- fileDirs
    private$corUtilsFuncs      <- corUtilsFuncs
    
    ## Get the group names
    private$group1Flatten = private$flattenGroup(self$group1)
    private$group2Flatten = private$flattenGroup(self$group2)
    
    ## Pair elements from two groups
    self$pairedList    = private$makePairs(names(private$group1Flatten), names(private$group2Flatten))
    
    ## printing
    print("printing Group1")
    print(private$group1Flatten)
    print("printing Group2")
    print(private$group2Flatten)
    
  },
  performDiffGeneExp   = function(){
    self$diffGeneExpList <- lapply(self$pairedList, private$executeDiffGeneExp)
    return(self$diffGeneExpList)
  }
  )
)
