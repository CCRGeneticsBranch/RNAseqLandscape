rm(list=ls())

## Source all classes and packages ####
project_name <- basename(getwd())
source("src/utilityPackages.R")
source("src/statisticalPackages.R")
source("src/utilityFunctions.R")
source("src/analysisFunctions.R")
source("src/class.R")

dir.create("output", showWarnings = F)
## Instantiate a new Object of type ProjectSetUp ####

meta_file <- "input/MetadataMapper.v10.txt"
meta_tcga_file <- "input/MetadataMapper.v10.tcga.txt"
meta_tcr_file <- "input/MetadataMapper.v10.TCR.txt"

rnaseqProject <- createRNAseqProject(project_name, meta_file, list("CellLine"=list("LIBRARY_TYPE"="CellLine")))

## Add utility functions to the project ####
corUtilsFuncs <- CoreUtilities$new(  ProjectSetUpObject = rnaseqProject )

# Make a Tree Map ####
StatsFinal <-  rnaseqProject$metaDataDF %>% group_by_(.dots= c("DIAGNOSIS.Alias","DIAGNOSIS.Alias.TreeMap",rnaseqProject$factorName, "Color", "LIBRARY_TYPE.TreeMap") ) %>% 
  count_(var=as.name("Sample.ID")) %>% dplyr::summarise(Count=n()) %>% 
  dplyr::group_by_(.dots= c(rnaseqProject$factorName)) %>%  
  dplyr::mutate( SampleSum := sum(Count)) %>% 
  spread_("LIBRARY_TYPE.TreeMap", "Count") %>% 
  mutate_( .dots = setNames( list( interp(~paste(rnaseqProject$factorName ,"(", Sum , ")"), 
                                          factorName=as.name(rnaseqProject$factorName), Sum=as.name("SampleSum") ) ), "LegendSampleSum") ) %>% 
  data.frame() %>% distinct(DIAGNOSIS.Alias.TreeMap,SampleSum, .keep_all = TRUE)

StatsFinal <- StatsFinal %>% dplyr::filter(DIAGNOSIS.Alias != "NS")
StatsFinal[,"LegendSampleSum"] <- paste(StatsFinal[,"DIAGNOSIS.Alias.TreeMap"],"( ",StatsFinal[,"SampleSum"], " )",sep="")
pdf(file=paste0(rnaseqProject$plotsDir, "/Diagnosis Tree Map All.pdf"), height=8, width= 10)
treemap(dtf=data.frame(StatsFinal), index=c("DIAGNOSIS.Alias", "DIAGNOSIS.Alias.TreeMap"),
        vSize="SampleSum",
        type="categorical",
        vColor="LegendSampleSum",
        palette = as.character(StatsFinal$Color),
        fontcolor.labels=c("black"),
        bg.labels=c("#CCCCCCDC"),
        algorithm = "squarified",
        inflate.labels=F,
        fontsize.labels = 10,
        fontsize.legend = 10,
        border.lwds=0.9,
        title = "Samples Map",
        title.legend = "Histology"
)
dev.off()

## Generate expression matrix ####
rm(mergeObjectsNoDup)
mergeObjectsNoDup_data <- readRDS("input/GeneRDS/RawCount/All.samples.Tumor.Normal.RDS")
## Tumor Normal and Cellline
## mergeObjectsNoDup_data <- readRDS("RNASeq.RSEM/GeneRDSOutput/RawCount/All.samples.Tumor.Normal.Celline.RDS")

design <- rnaseqProject$metaDataDF


# ## Rearrange the design matrix and data matrix
design$LIBRARY_TYPE <- factor(design$LIBRARY_TYPE, levels =c("Tumor", "Normal"), ordered = TRUE)
design$Target.Status.Life <- factor(design$Target.Status.Life, levels = c("Alive", "Dead", ""), ordered = TRUE)
design$Target.Status.Risk <- factor(design$Target.Status.Risk, levels = c("Low.Risk", "Intermediate.Risk", "High.Risk"), ordered = TRUE)
design <- design %>% arrange( LIBRARY_TYPE, DIAGNOSIS.Substatus.Tumor.Normal.Tissue,Target.Status.Life,Target.Status.Risk )
mergeObjectsNoDup <- mergeObjectsNoDup_data %>% dplyr::select(one_of(c("rn",as.character(design[,c(rnaseqProject$metadataFileRefCol)])))); 
#mergeObjectsNoDup <- mergeObjectsNoDup_data %>% dplyr::select(one_of(as.character(design[,c(rnaseqProject$metadataFileRefCol)]))); 

## Check if designmatrix and count matrix have same order of columns
#HC: ignore the first "rn" column
#View(data.frame(count_names=colnames(mergeObjectsNoDup)[2:ncol(mergeObjectsNoDup)], design_names=design[,rnaseqProject$metadataFileRefCol]))

#HC: check if colnames are the same
if (all(colnames(mergeObjectsNoDup)[2:ncol(mergeObjectsNoDup)] == design[,rnaseqProject$metadataFileRefCol])) {
  print("Column matched!")
} else {
  print("Column not matched!")
}
## Evaluate presence of duplicate features (genes) and consolidate them ####
#mergeObjectsNoDup.pre <- setDT(mergeObjectsNoDup, keep.rownames = TRUE) %>% dplyr::rename(GeneID = rn) 
mergeObjectsNoDup.pre <- mergeObjectsNoDup %>% dplyr::rename(GeneID = rn) 
design$LIBRARY_TYPE <- factor(design$LIBRARY_TYPE, levels =c("Tumor", "CellLine", "Normal"), ordered = TRUE)
mergeObjectsNoDup.pre <- dplyr::left_join(rnaseqProject$annotationDF[,c("GeneID", "GeneName")], mergeObjectsNoDup.pre, by="GeneID") %>% 
  data.table()
mergeObjectsConso     <- corUtilsFuncs$consolidateDF(mergeObjectsNoDup.pre[,-c("GeneID")], funcName = "max", featureName = "GeneName")
mergeObjectsConso     <- dplyr::full_join(mergeObjectsConso, rnaseqProject$annotationDF[,c("GeneID", "GeneName")], by="GeneName") %>%  
  data.table()
mergeObjectsConso     <- subset(mergeObjectsConso,!duplicated(mergeObjectsConso$GeneName))
mergeObjectsConso     <- mergeObjectsConso[complete.cases(mergeObjectsConso), ]; dim(mergeObjectsConso)
mergeObjectsConso     <- mergeObjectsConso[,-c("GeneName")]         %>% 
  data.frame()                               %>% 
  tibble::column_to_rownames(var = "GeneID") %>% 
  as.matrix() ; dim(mergeObjectsConso)
## matching above data frame with the annotationDF
rnaseqProject$annotationDF <- rnaseqProject$annotationDF %>% dplyr::filter(GeneID %in% rownames(mergeObjectsConso)); dim(rnaseqProject$annotationDF)

## Subset metaDataDF by the number of samples in the folder ####
colnamesDF           <- data.frame( "Sample.Biowulf.ID.GeneExp"= colnames(mergeObjectsConso))
corUtilsFuncs$subsetMetaData(colnamesDF=colnamesDF)

## Instantiate a new Object of type GeneExpNormalization ####
expressionObj        <- GeneExpNormalization$new(
  
  countObj          = as.matrix(mergeObjectsConso), 
  featureType       = "Gene", 
  packageRNAseq     = "edgeR", 
  annotationDF      = rnaseqProject$annotationDF, 
  design            = design[,rnaseqProject$factorName], 
  #design           = newMetaDataDF[,rnaseqProject$factorName],
  proteinCodingOnly = FALSE,
  corUtilsFuncs     = corUtilsFuncs
)

## Get expression in desired units ####
### RawCounts
#expressionTMM.Counts          = expressionObj$edgeRMethod("RawCounts")
## Normalised counts
#expressionTMM.NormDF         = expressionObj$edgeRMethod("NormFactorDF")

### RPKM
expressionTMM.RPKM            = expressionObj$edgeRMethod("TMM-RPKM", logtransform = TRUE, zscore = FALSE)
designMatrix                  <- corUtilsFuncs$validfMatrix(df = design)
### Zscore ###
expressionTMM.RPKM.zscore <- expressionObj$edgeRMethod("TMM-RPKM", logtransform = TRUE, zscore = TRUE)
## Replace NA with 0
expressionTMM.RPKM.zscore[is.na(expressionTMM.RPKM.zscore)] <- 0

## Arrange data by histology and Library type
expressionTMM.RPKM.arr <- expressionTMM.RPKM %>% dplyr::select(one_of("Chr","Start","End","Strand","GeneID","GeneName","Length",
                                                                      as.character(gsub("-",".",designMatrix[,rnaseqProject$metadataFileRefCol]))))

## Add additional annotations (sample Id alias) ####
AliasNames_df                 <- dplyr::left_join( data.frame("Sample.Biowulf.ID.GeneExp"=colnames(expressionTMM.RPKM.arr)), 
                                                   designMatrix[,c(rnaseqProject$metadataFileRefCol,rnaseqProject$factorName,"Sample.ID.Alias", 
                                                                   "Sample.Data.ID", "DIAGNOSIS.Alias","Annotation_Target_Khanlab_jun", 
                                                                   "SampleID..In.Paper")] )
AliasColnames                 <- c(as.character(AliasNames_df[c(1:7),1]), as.character(AliasNames_df[-c(1:7),7])); 

## Check if designmatrix and count matrix have same order of columns
#View(data.frame(count_names=colnames(expressionTMM.RPKM.arr)[-c(1:7)], 
#                design_names=designMatrix[,rnaseqProject$metadataFileRefCol],
#                AliasColnames = as.character(AliasNames_df[-c(1:7),7])))

## Perform Sanity Check for the above operations #####
stopifnot( length(colnames(expressionTMM.RPKM.arr)) == length(AliasColnames) )
colnames(expressionTMM.RPKM.arr)  <- AliasColnames

# ## Arrange data by histology and Library type (Zscore)
#HC: uncomment this so that we can continue
 expressionTMM.RPKM.arr.zscore <- expressionTMM.RPKM.zscore  %>% dplyr::select(one_of("Chr","Start","End","Strand","GeneID","GeneName","Length",
                                                                                     as.character(gsub("-",".",designMatrix[,rnaseqProject$metadataFileRefCol]))))

## Add additional annotations (sample Id alias) ####
AliasNames_df                 <- dplyr::left_join( data.frame("Sample.Biowulf.ID.GeneExp"=colnames(expressionTMM.RPKM.arr.zscore)), 
                                                   designMatrix[,c(rnaseqProject$metadataFileRefCol,rnaseqProject$factorName,"Sample.ID.Alias", 
                                                                   "Sample.Data.ID", "DIAGNOSIS.Alias","Annotation_Target_Khanlab_jun")] )
AliasColnames                 <- c(as.character(AliasNames_df[c(1:7),1]), as.character(AliasNames_df[-c(1:7),6]));

## Perform Sanity Check for the above operations #####
stopifnot( length(colnames(expressionTMM.RPKM.arr.zscore)) == length(AliasColnames) )
colnames(expressionTMM.RPKM.arr.zscore)  <- AliasColnames

### Performing ssGSEA output analysis. ( Plotting the scores across histology ) ##########

### Prepare input for ssGSEA broad gene pattern
expressionTMM.RPKM.GSEA.Input <- expressionTMM.RPKM.arr[, -c(1:7)]; 
rownames(expressionTMM.RPKM.GSEA.Input) <- expressionTMM.RPKM.arr$GeneName
expressionTMM.RPKM.GSEA.print = corUtilsFuncs$createBroadGCTFile(expressionTMM.RPKM.GSEA.Input)

#******************* HC: we should run all the way to here *************************
#*
## Only For TCGA+Khanlab dataSet 
#HC: score with TCGA

rnaseqTCGAProject <- createRNAseqProject(project_name, meta_tcga_file, list("CellLine"=list("LIBRARY_TYPE"="CellLine")))

khanlab.TCGA.geneExp <- readRDS("input/GeneRDS/RPKM_Data_Filt_Consolidated.GeneNames.all.TCGA.Khanlab.pc.log22019-03-19.rds")
expressionTMM.TCGA.RPKM.GSEA.Input <- khanlab.TCGA.geneExp[, -c(1:7)]; rownames(expressionTMM.TCGA.RPKM.GSEA.Input) <- khanlab.TCGA.geneExp[,6]

## Read the ssGSEA output
ssGSEAScores            <- corUtilsFuncs$parseBroadGTCOutFile("input/GSEA/RPKM_Data_Filt_Consolidated.GeneNames.all.Khanlab.pc.log2.2019-06-14.PROJ.gct")
ssGSEAScoresTCGA            <- corUtilsFuncs$parseBroadGTCOutFile("input/GSEA/RPKM.TMM.RPKM.GSEA.Input.All.TCGA.Khanlab.2019-03-19.PROJ.gct")

## Add custom expression like cytolytic scre and HLA gene expression to the ssGSEA Outpuut file.
cytolyticScore          <- corUtilsFuncs$cytolyticScore(expressionTMM.TCGA.RPKM.GSEA.Input)
HLA_cytolyticScore      <- rbind(expressionTMM.TCGA.RPKM.GSEA.Input[c("HLA-A", "HLA-B", "HLA-C"),], cytolyticScore)
#View(data.frame(colnames(HLA_cytolyticScore), colnames(ssGSEAScoresTCGA)))
ssGSEAScores.HLA.Cyto   <- rbind(ssGSEAScoresTCGA,HLA_cytolyticScore)

#HC: column order is different. need to make them consistent

## Plot the one variable plot

## Filter specified Diagnosis
factorsToExclude              = paste(c("NS.", "YST", "Teratoma"), collapse = "|")
selected.metadata              <- rnaseqTCGAProject$validMetaDataDF  %>% 
  filter_(  .dots = paste0("!grepl(", "'", factorsToExclude , "'" ,",", rnaseqTCGAProject$factorName, ")")) %>% 
  dplyr::select_( .dots=c(rnaseqTCGAProject$metadataFileRefCol, rnaseqTCGAProject$factorName ) )

ssGSEAScores.HLA.Cyto.Selected <- ssGSEAScores.HLA.Cyto %>% dplyr::select(one_of(as.character(selected.metadata[, rnaseqTCGAProject$metadataFileRefCol])))
dim(ssGSEAScores.HLA.Cyto.Selected)

## sanity check Checking metadata vs data ##
stopifnot( ncol(ssGSEAScores.HLA.Cyto.Selected) == length(as.character(selected.metadata$Sample.Biowulf.ID.GeneExp)) )

## Preparing the expression matrix for string plot, by appending metadata
Scores <- cbind(t(ssGSEAScores.HLA.Cyto.Selected), selected.metadata[,rnaseqProject$factorName, drop=FALSE]) %>% 
  dplyr::rename_(.dots = setNames( list(rnaseqProject$factorName), list("Diagnosis") )) #%>%
#dplyr::mutate(Diagnosis = factor(Diagnosis, ordered = TRUE, levels = orderOfFactor))

### Setting up variables for  string plot
## Set the order of Diagnosis to appear
orderOfFactor    <- unique(Scores$Diagnosis)
## Set the order of signature to appear
orderOfSignature <- colnames(Scores)[-ncol(Scores)]
## Total list of signatures
colList          <- c(1:(ncol(Scores)-1))
## Generate custom colors
customColorDF    <- rnaseqTCGAProject$customColorsDFAll

### Filter the score matrix for diffent categories in the plot

### For all nothing to be changed
### For only khanlab filter the tidy score matrix & color matrix

#Scores <- Scores %>% filter(!grepl("TCGA.",Diagnosis)) %>% filter(!grepl(".CellLine",Diagnosis))
#customColorDF <- customColorDF %>% filter(!grepl("TCGA.",Diagnosis)) %>% filter(!grepl(".CellLine",Diagnosis))

## Plot the onevariable plot
plotLists        <- corUtilsFuncs$OneVariablePlotSort(colList, Scores=Scores, orderOfFactor, orderOfSignature, standardize =FALSE, logit =FALSE, plotType = "density",
                                                      yLab = "Standardised enrichment score", legendDisplay = FALSE, customColorDF = customColorDF )
## Save the plots
EnrischmentScorePlots <- lapply(plotLists, function(l) l[[1]])
SBName                <- paste(rnaseqTCGAProject$workDir, rnaseqTCGAProject$plotsDir,"TMM-RPKM.ssGSEA.enrichmentScores.all.pc.Khanlab.TCGA.log.pdf",sep="/")
ggsave(SBName, marrangeGrob(EnrischmentScorePlots,ncol=2,nrow=1 ), width = 20, height = 15 )

### For Everything except Cellline

Scores <- Scores %>%  filter(!grepl(".CellLine",Diagnosis))
customColorDF <- customColorDF  %>% filter(!grepl(".CellLine",Diagnosis))

## Plot the onevariable plot
plotLists        <- corUtilsFuncs$OneVariablePlotSort(colList, Scores=Scores, orderOfFactor, orderOfSignature, standardize =TRUE, logit =FALSE, plotType = "density",
                                                      yLab = "Standardised enrichment score", legendDisplay = FALSE, customColorDF = customColorDF )
## Save the plots
EnrischmentScorePlots <- lapply(plotLists, function(l) l[[1]])
SBName                <- paste(rnaseqTCGAProject$workDir, rnaseqTCGAProject$plotsDir,"Figure1B&S1.pdf",sep="/")
ggsave(SBName, marrangeGrob(EnrischmentScorePlots,ncol=2,nrow=1 ), width = 20, height = 15 )

## Plot to do percent samples enriched across cancer types

## Using  scores
# dropSignatures    <- c("Macrophages_M0","Macrophages_M1", "Macrophages_M2","Dendritic_cells_activated")
factorsToExclude  <- paste(c("NS", "YST", "Teratoma"), collapse = "|")

## Read and parse the ssGSEA Output from Broad GenePattern
Scores            <- corUtilsFuncs$parseBroadGTCOutFile("input/GSEA/RPKM_Data_Filt_Consolidated.GeneNames.all.pc.log2.2019-01-31.PROJ.KeggSig.gct")
## Standardizing the raw score to amplify the difference.
ScoresZscore      <- apply(Scores[c(1:24,43),],1, corUtilsFuncs$zscore_All)                  
#ScoresZscore      %<>% data.frame() %<>% tibble::rownames_to_column(var=rnaseqProject$metadataFileRefCol) # %<>% dplyr::select(-one_of(dropSignatures ))
ScoresZscore <- ScoresZscore %>% data.frame() %>% tibble::rownames_to_column(var=rnaseqProject$metadataFileRefCol) # %<>% dplyr::select(-one_of(dropSignatures ))

## Preparing the data and prepare for heatmap
rnaseqProject$metaDataDF <- as.data.frame( apply(rnaseqProject$metaDataDF, 2, as.character), stringsAsFactors = FALSE )
ScoresForGather   <- tidyr::gather(ScoresZscore, key="GeneSet", value="Score", -!!rnaseqProject$metadataFileRefCol )
ScoresForGather   <- dplyr::left_join(ScoresForGather, 
                                      rnaseqProject$metaDataDF[,c(rnaseqProject$metadataFileRefCol, rnaseqProject$factorName)], 
                                      by=rnaseqProject$metadataFileRefCol) %>% 
  dplyr::filter_(  .dots = paste0("!grepl(", "'", factorsToExclude , "'" ,",", rnaseqProject$factorName, ")")) %>% 
  dplyr::rename_(.dots = setNames(list(rnaseqProject$factorName),c("Diagnosis")) )

## Order of genesets
genesets <- c("ImmuneSignature",
              "StromalSignature",
              "Kegg_Antigen_processing_and_presentation",
              "T.cells_CD8",
              "T.cells_CD4_naive",
              "T.cells_CD4_memory_resting",
              "T.cells_CD4_memory_activated",
              "T.cells_follicular_helper",
              "T.cells_regulatory",
              "T.cells_gamma_delta",
              "NK.cells_activated",
              "NK.cells_resting",
              "B.cells_naive",
              "B.cells_memory",
              "Plasma_cells",
              "Monocytes",
              "Dendritic_cells_resting",
              "Macrophages_M0",
              "Macrophages_M1",
              "Macrophages_M2",
              "Dendritic_cells_activated",
              "Neutrophils",
              "Eosinophils",
              "Mast_cells_resting",
              "Mast_cells_activated"
)

Diagnosis <- c("WT", "SS", "CCSK", "EWS",  "RMS.FN", "RMS.FP", "NB.MYCN.A","NB.Unknown", "DSRCT", "NB.MYCN.NA",  "OS", "UDS", 
               "ML", "HBL", "ASPS")
ScoresForGather$GeneSet <- factor(ScoresForGather$GeneSet, levels = genesets, ordered = TRUE)
ScoresForGather$Diagnosis <- factor(ScoresForGather$Diagnosis, levels = Diagnosis, ordered = TRUE)

ScoresForGatherPercent         <- ScoresForGather                                             %>% 
  dplyr::group_by(Diagnosis, GeneSet)                         %>% 
  dplyr::mutate(TotalCount = n(), Enriched = sum(Score > 0 )) %>% 
  dplyr::mutate(SamplePercent = (Enriched/TotalCount)*100 )   
ScoresForGatherUnique          <- ScoresForGatherPercent[,c(2,4,7)]                           %>% 
  ungroup()                                                   %>% 
  distinct()
ScoresForSpread                <- tidyr::spread( ScoresForGatherUnique, Diagnosis, SamplePercent ) %>% t() 
colnames(ScoresForSpread)      <- ScoresForSpread[1,]; 
ScoresForSpreadHeat            <- ScoresForSpread[-1,]
ScoresForSpreadHeat            <- t(apply(ScoresForSpreadHeat, 1, as.numeric))
colnames(ScoresForSpreadHeat)  <- colnames(ScoresForSpread)
#write.table(ScoresForSpread, "RNASeq/PlotData/ScoresperDiagSigEnrichZscore.txt", sep="\t", quote = F, col.names = T, row.names = F)

## Open the PDF File
pdf( paste(rnaseqProject$workDir, rnaseqProject$plotsDir,"Figure1A.pdf", sep="/"), height=10, width = 20)

### Using two different packages to plot
## Using SuperheatMap package
# superheat(ScoresForSpreadHeat,
#           bottom.label.text.angle=90,
#           title.size = 6,
#           heat.pal = c( "#4FFC07","#273746", "#F92908"),
#           pretty.order.rows = T,
#           pretty.order.cols = T
# )

## Using fheatmap (presently using)
breaks <- seq(min(ScoresForSpreadHeat),max(ScoresForSpreadHeat), by=0.1)
#Yellow blue grey
#matrix_color_vector <- colorpanel(n=length(breaks)-1,low="#F4D03F",mid="#273746",high="#5DADE2")
# #green black red
# matrix_color_vector <- colorpanel(n=length(breaks)-1,low="#4FFC07",mid="#0B0B0B",high="#F92908")
#black red
matrix_color_vector <- colorpanel(n=length(breaks)-1,low="#4FFC07",mid="#273746",high="#F92908")

fheatmap(t(ScoresForSpreadHeat), display_tree_col = F,cluster_rows = F, cluster_cols = F, mat_color = matrix_color_vector,
         row_fontsize = 5, col_fontsize = 5, cell_border = F, cell_border_col = "#A6ACAF",seed = 10,
         title = "Percent samples enriched across cancer types")

dev.off()

### Perform correlation analysis  for DDR genes ####
## Read the ssGSEA output
ssGSEAScores            <- corUtilsFuncs$parseBroadGTCOutFile("input/GSEA/RPKM_Data_Filt_Consolidated.GeneNames.all.pc.log2.2019-01-31.PROJ.KeggSig.gct")
## Read the DDR genes file
DDRGenes <- read.table("input/AnnotationRDS/276_DDR_Genes.v2.txt", sep="\t", header = T, stringsAsFactors = FALSE)
## sanity check if any DDR genes are missing form the expression matrix
excludeGenesIndx <- which(! DDRGenes$DDRGenes %in% rownames(expressionTMM.RPKM.GSEA.Input))
DDRGenesPresent <- DDRGenes$DDRGenes[-excludeGenesIndx]; length(DDRGenesPresent)
## Get gene expression the present genes
DDRGenesGeneExp     <- expressionTMM.RPKM.GSEA.Input[DDRGenesPresent,]; dim(DDRGenesGeneExp)
## Sanity check for NA or Inf
indx <- apply(DDRGenesGeneExp, 1, function(x) any(is.na(x) | is.infinite(x)))
rownames(DDRGenesGeneExp)[indx];
DDRGenesGeneExpFinal <-  DDRGenesGeneExp[complete.cases(DDRGenesGeneExp), ]; dim(DDRGenesGeneExpFinal)
## bind geneexpression matrix to immunescore matrix and perform correlation
#HC: rename and reorder 
DDRGenesGeneExpFinal_tmp <- reOrderForssGSEA(as.data.frame(t(DDRGenesGeneExpFinal)), rnaseqProject$validMetaDataDF, ssGSEAScores, "SampleID..In.Paper", rnaseqProject$metadataFileRefCol)

ssGSEAScores.DDRGenes  <- rbind(ssGSEAScores,DDRGenesGeneExpFinal_tmp)
ssGSEAScores.DDRGenes.corr <- rcorr(t(ssGSEAScores.DDRGenes),  type = "spearman");
## Plot the heatmap
col<- colorRampPalette(c("blue", "white", "red"))(20)
pdf("output/Figures/ImmuneScoreVSDDRgeneExp.pdf", height = 55, width = 55)
heatmap(x = ssGSEAScores.DDRGenes.corr$r, col = col, symm = TRUE, keep.dendro = FALSE)
dev.off()
## slicing out the interesting part
ssGSEAScores.DDRGenes  <- rbind(ssGSEAScores[c("T.regulatory_PMID_30127393_neg", "T.regulatory_PMID_30127393_pos"),],DDRGenesGeneExpFinal_tmp)
ssGSEAScores.DDRGenes.corr <- rcorr(t(ssGSEAScores.DDRGenes),  type = "spearman"); 
pdf("output/Figures/ImmuneScoreVSDDRgeneExpSelected.pdf", height = 55, width = 55)
heatmap(x = ssGSEAScores.DDRGenes.corr$r, col = col, symm = TRUE, keep.dendro = FALSE)
dev.off()
## Writing data to files
write.table(ssGSEAScores.DDRGenes, "output/FigureData/ssGSEAScores.DDRGenes.txt", quote = FALSE, sep = "\t")
write.table(ssGSEAScores.DDRGenes.corr$p, "output/FigureData/ImmuneSoreVSDDRgexp.p.txt", quote = FALSE)

### Perform correlation analysis  for CTA genes ####
## Read the ssGSEA output
ssGSEAScores            <- corUtilsFuncs$parseBroadGTCOutFile("input/GSEA/RPKM_Data_Filt_Consolidated.GeneNames.all.Khanlab.pc.log2.2019-06-14.PROJ.kegg.immune.rm.dups.gct")
## Read the DDR genes file
rm(CGAGenes)
CGAGenes <- as.character(rnaseqProject$csDF$GeneName)
## sanity check if any DDR genes are missing form the expression matrix
excludeGenesIndx <- which(! CGAGenes %in% rownames(expressionTMM.RPKM.GSEA.Input)); length(excludeGenesIndx)
DDRGenesabsent <- CGAGenes[excludeGenesIndx]; as.character(DDRGenesabsent)
DDRGenesPresent <- CGAGenes[-excludeGenesIndx]; length(DDRGenesPresent)
## Get gene expression the present genes
DDRGenesGeneExp     <- expressionTMM.RPKM.GSEA.Input[DDRGenesPresent,]; dim(DDRGenesGeneExp)
## Sanity check for NA or Inf
indx <- apply(DDRGenesGeneExp, 1, function(x) any(is.na(x) | is.infinite(x)))
rownames(DDRGenesGeneExp)[indx];
DDRGenesGeneExpFinal <-  DDRGenesGeneExp[complete.cases(DDRGenesGeneExp), ]; dim(DDRGenesGeneExpFinal)
## bind geneexpression matrix to immunescore matrix and perform correlation
cytolyticScore          <- corUtilsFuncs$cytolyticScore(expressionTMM.RPKM.GSEA.Input)

cytolyticScore_tmp <- reOrderForssGSEA(as.data.frame(t(cytolyticScore)), rnaseqProject$validMetaDataDF, ssGSEAScores, "SampleID..In.Paper", rnaseqProject$metadataFileRefCol)
ssGSEAScores.CGAGenes  <- rbind(ssGSEAScores[4,],cytolyticScore_tmp,DDRGenesGeneExpFinal_tmp) %>% data.frame()

## Remove Normal samples
Tumor_samples_annot <- AliasNames_df %>% dplyr::filter(!grepl('Teratoma|NS',DIAGNOSIS.Substatus.Tumor.Normal.Tissue))
Tumor_samples_annot <- Tumor_samples_annot[complete.cases(Tumor_samples_annot),]
ssGSEAScores.CGAGenes.Tumor <- ssGSEAScores.CGAGenes %>% dplyr::select(one_of(as.character(Tumor_samples_annot$Sample.Biowulf.ID.GeneExp)))

Tumor_samples_annot$DIAGNOSIS.Substatus.Tumor.Normal.Tissue <- gsub('NS.*','NS',Tumor_samples_annot$DIAGNOSIS.Substatus.Tumor.Normal.Tissue)
## Perform rcorr for each tumor groups 
correlationDF <- rcorr_groups(Tumor_samples_annot, ssGSEAScores.CGAGenes.Tumor, rnaseqProject$factorName, "Sample.Biowulf.ID.GeneExp")
correlationDF_cytolytic <- correlationDF %>% filter(grepl('CytolyticScore|T.cells_CD8', row))
write.table(correlationDF_cytolytic, "output/FigureData/CD8.cytolytic.cs.each.Tumor.corr.txt", quote = FALSE, sep="\t")


### Perform annotation for the above
table_diff_Exp_annot <- read.table("input/CD8/CD8.annotation/TranscriptionFactor.Summarised.traditionalRank.Dexp.txt",
                                   sep="\t", header = T, row.names = 1)
colnames(table_diff_Exp_annot) <- gsub("NormalsStatus", '', colnames(table_diff_Exp_annot))

cor.data <- read.table("input/CD8/CD8.annotation/CD8.cytolytic.tf.each.Tumor.corr.txt",
                       sep="\t", header = T, row.names = 1)

diff_exp_annotate <- function(x, lookupDF=NA){
  value = lookupDF[x["column"],x["group"]]
  if(is.null(value)){
    value= NA
  }
  x$DiffExp <- value
  return(data.frame(x))
}

annotatedList <- apply(cor.data, 1, diff_exp_annotate, lookupDF=table_diff_Exp_annot)
annotatedDF <- do.call(rbind.data.frame,annotatedList)
write.table(annotatedDF, 
            "output/FigureData/CD8.cytolytic.tf.each.Tumor.corr.diff.exp.txt", quote = FALSE, sep="\t")


#HC: Fig2
#### Perform correlation analysis between Exhaustion markers and immune signatures ####
## Read the ssGSEA output
ssGSEAScores            <- corUtilsFuncs$parseBroadGTCOutFile("input/GSEA/RPKM_Data_Filt_Consolidated.GeneNames.all.Khanlab.pc.log2.2019-06-14.PROJ.kegg.immune.rm.dups.gct")
## Read the DDR genes file
EMGenesPre <- read.table("input/AnnotationRDS/Exhaustion_markers_genes.txt", sep="\t", header = T, stringsAsFactors = FALSE)
EMGenes <- data.frame(Exhaustion_Marker_Genes = c(EMGenesPre$Exhaustion_Marker_Genes, 'C10orf54') )
## sanity check if any EM genes are missing form the expression matrix
excludeGenesIndx <- which(!EMGenes$Exhaustion_Marker_Genes %in% rownames(expressionTMM.RPKM.GSEA.Input))
EMGenesPresent <- EMGenes[!excludeGenesIndx,]
if( length(EMGenesPresent) != 0) {
  print("Some Genes are not found !!")
  EMGenesPresent <- EMGenes$Exhaustion_Marker_Genes[-excludeGenesIndx]; length(EMGenesPresent)
} else {
  print("All Genes are found !!")
  EMGenesPresent <- EMGenes$Exhaustion_Marker_Genes; length(EMGenesPresent)
}
## Get gene expression the present genes
EMGenesGeneExp     <- expressionTMM.RPKM.GSEA.Input[EMGenesPresent,]; dim(EMGenesGeneExp)
## Sanity check for NA or Inf
index <- apply(EMGenesGeneExp, 1, function(x) any(is.na(x) | is.infinite(x)))
rownames(EMGenesGeneExp)[index];
EMGenesGeneExpFinal <-  EMGenesGeneExp[complete.cases(EMGenesGeneExp), ]; dim(EMGenesGeneExpFinal)
## bind geneexpression matrix to immunescore matrix and perform correlation

EMGenesGeneExpFinal <- reOrderForssGSEA(as.data.frame(t(EMGenesGeneExpFinal)), rnaseqProject$validMetaDataDF, ssGSEAScores, "SampleID..In.Paper", rnaseqProject$metadataFileRefCol)

ssGSEAScores.EMGenes  <- rbind(ssGSEAScores, EMGenesGeneExpFinal); dim(ssGSEAScores.EMGenes)

# ## For making heatmap
# ssGSEAScores.EMGenes  <- rbind(ssGSEAScores[c("T.cells_CD4_memory_activated", "T.cells_CD8"),], EMGenesGeneExpFinal); dim(ssGSEAScores.EMGenes)

## Get correlation between Gene expression and Immune signature irrespective of diagnosis
ssGSEAScores.EMGenes.T <- t(ssGSEAScores.EMGenes)
ssGSEAScores.EMGenes.corr <- rcorr(ssGSEAScores.EMGenes.T,  type = "spearman");

## Get correlation Immune gene signatures irrespective of diagnosis
ssGSEAScores.T <- t(ssGSEAScores)
ssGSEAScores.corr <- rcorr(ssGSEAScores.T,  type = "spearman");

## coorelation plot
ssGSEAScores.T.df <- data.frame(ssGSEAScores.T)
ssGSEAScores.T.df[,"diff"] <- abs(ssGSEAScores.T[,"Kegg_Antigen_processing_and_presentation"]-ssGSEAScores.T[,"ImmuneSignature"])

## Plot heatmap for the above
col<- colorRampPalette(c("blue", "white", "red"))(100)
pdf("output/Figures/ImmuneScoreVSDDRgeneExp.v2.pdf", height = 55, width = 55)
#heatmap(x = ssGSEAScores.EMGenes.corr$r, col = col, symm = TRUE, keep.dendro = FALSE)
pheatmap(ssGSEAScores.corr$r, col = col)
dev.off()

## Add annotation
ssGSEAScores.EMGenes.T <- t(ssGSEAScores.EMGenes)
ssGSEAScores.EMGenes.T <- data.frame(ssGSEAScores.EMGenes.T) %>%  tibble::rownames_to_column(var="Sample.Biowulf.ID.GeneExp")
ssGSEAScores.EMGenes.T.annot.all <- dplyr::left_join(ssGSEAScores.EMGenes.T, AliasNames_df[,c(1,2,3)], 
                                                     by="Sample.Biowulf.ID.GeneExp") %>%
  dplyr::rename(Diagnosis=DIAGNOSIS.Substatus.Tumor.Normal.Tissue)
ssGSEAScores.EMGenes.T.annot <- ssGSEAScores.EMGenes.T.annot.all %>% filter( !grepl("NS", Diagnosis) )  %>%
  tibble::column_to_rownames(var="Sample.Biowulf.ID.GeneExp")
## Make Correlation matrix
rownames <- rep( EMGenes$Exhaustion_Marker_Genes , length( unique(ssGSEAScores.EMGenes.T.annot$Diagnosis)) ); length(rownames)
gp = dplyr::group_by(ssGSEAScores.EMGenes.T.annot, Diagnosis)
gpTowrite <- gp ; rownames(gpTowrite) <- rownames(ssGSEAScores.EMGenes.T.annot)

CorrDF.Out.R <- dplyr::do(gp, data.frame(Cor=t(corr.test(.[,1:43], .[,44:66], method = "spearman")$r))) %>% data.frame() %>%  
  mutate_all( funs_( interp( ~replace(., is.na(.),0) ) ) ); dim(CorrDF.Out.R); head(CorrDF.Out.R)
CorrDF.Out.R <- tibble::add_column(CorrDF.Out.R, GeneNames= rownames, .after=1) 
CorrDF.Out.R$GeneNames <- factor(CorrDF.Out.R$GeneNames, levels = sort(unique(CorrDF.Out.R$GeneNames)))
write.table(CorrDF.Out.R, paste0("output/FigureData/CorrDF.Out.R.EM.spearman.v2.txt"), sep="\t", row.names = F, quote = FALSE )
CorrDF.Out.P <- dplyr::do(gp, data.frame(Cor=t(corr.test(.[,1:43], .[,44:66], method = "spearman")$p))) %>% data.frame()
CorrDF.Out.P <- tibble::add_column(CorrDF.Out.P, GeneNames= rownames, .after=1)
write.table(CorrDF.Out.P, paste("output/FigureData/CorrDF.Out.P.EM.spearnman.v2.txt", sep=""), sep="\t", row.names = F, quote = FALSE )      

## Post Analysis; Merge coefficient and p values
CorrDF.Out.R.TC8.TC4 <- CorrDF.Out.R %>% dplyr::select(matches("Diagnosis|GeneNames|CD8|CD4"))
CorrDF.Out.P.TC8.TC4 <- CorrDF.Out.P %>% dplyr::select(matches("Diagnosis|GeneNames|CD8|CD4"))      
CorrDF.Out.R.P.TC8.TC4 <- merge(CorrDF.Out.R.TC8.TC4, CorrDF.Out.P.TC8.TC4, by=c("Diagnosis", "GeneNames"))
CorrDF.Out.R.P.TC8.TC4$Diagnosis <- factor(as.character(CorrDF.Out.R.P.TC8.TC4$Diagnosis),ordered = TRUE,
                                           #levels = sort(unique( CorrDF.Out.R.P.TC8.TC4$Diagnosis )) )
                                           levels = c("NB.MYCN.NA","NB.MYCN.A","EWS","DSRCT","OS","RMS.FP","RMS.FN",
                                                      "Teratoma","SS","CCSK","NB.Unknown","ASPS","HBL","WT","ML","UDS","YST") )
## Filter correlation based on the filter
CorrDF.Out.R.P.CD8.CD4_MemoryAct <- CorrDF.Out.R.P.TC8.TC4 %>% 
  mutate(TC8.TC4.MemoryAct=ifelse( abs( Cor.T.cells_CD8.x >= 0.3 & Cor.T.cells_CD8.y <= 0.05 ) |
                                     abs( Cor.T.cells_CD4_memory_activated.x >= 0.3 & Cor.T.cells_CD4_memory_activated.y <= 0.05 ), 1, 0 ) ) %>%
  dplyr::select(Diagnosis, GeneNames, TC8.TC4.MemoryAct) 
CorrDF.Out.R.P.CD8.CD4_MemoryAct.Spread <- CorrDF.Out.R.P.CD8.CD4_MemoryAct %>% tidyr::spread(Diagnosis,TC8.TC4.MemoryAct)
CD8.CD4_MemoryAct.Spread <- CorrDF.Out.R.P.CD8.CD4_MemoryAct.Spread  %>% dplyr::mutate(Count=rowSums(.[2:ncol(CorrDF.Out.R.P.CD8.CD4_MemoryAct.Spread)]))
CD8.CD4_MemoryAct.Spread <- tibble::add_column( CD8.CD4_MemoryAct.Spread, Legend=paste(CD8.CD4_MemoryAct.Spread$GeneNames, 
                            "(", CD8.CD4_MemoryAct.Spread$Count, ")"), .after=1) %>% data.frame() %>%  dplyr::select(-contains("GeneNames"))

## Filter correlation based on the filter plot correlation values only.
CorrDF.Out.R.P.CD8.heatmap <- CorrDF.Out.R.P.TC8.TC4 %>% 
  mutate(TC8.TC4.MemoryAct=ifelse( abs( Cor.T.cells_CD8.y <= 0.05 ), Cor.T.cells_CD8.x, 0 ) ) %>%
  dplyr::select(Diagnosis, GeneNames, TC8.TC4.MemoryAct) 
CorrDF.Out.R.P.CD8.heatmap.Spread <- CorrDF.Out.R.P.CD8.heatmap %>% tidyr::spread(Diagnosis,TC8.TC4.MemoryAct)
CorrDF.Out.R.P.CD8.heatmap.Spread <- CorrDF.Out.R.P.CD8.heatmap.Spread  %>% dplyr::mutate(Count=rowSums(.[2:ncol(CorrDF.Out.R.P.CD8.heatmap.Spread)]))
CorrDF.Out.R.P.CD8.heatmap.Spread <- tibble::add_column( CorrDF.Out.R.P.CD8.heatmap.Spread, Legend=paste(CorrDF.Out.R.P.CD8.heatmap.Spread$GeneNames), .after=1) %>% 
  data.frame() %>% dplyr::select(-contains("GeneNames"))

## Plot the heatmap
CorrDF <- CorrDF.Out.R.P.CD8.heatmap.Spread
#CorrDF <- CorrDF.Out.R.P.CD8.CD4_MemoryAct
CorrDF <- CorrDF %>% 
  dplyr::arrange(Count) %>% 
  dplyr::filter( Count > 0 ) %>% 
  dplyr::select(-one_of("Count")) %>% 
  dplyr::select(-matches("ML|UDS|YST")) %>%
  tibble::column_to_rownames(var="Legend")
# CD8.CD4_MemoryAct.Spread %<>%   dplyr::arrange(desc(Legend)) %<>% tibble::column_to_rownames(var="Legend")
# CD8.CD4_MemoryAct.Spread %<>% column_to_rownames(var="Legend")

colfunc <- colorRampPalette(c("#f2f2f2", "firebrick3"))
pdf("output/Figures/Figure2B.pdf")
pheatmap(t(CorrDF), 
         #color =c("#e0e0d1", "#004080"), 
         #color =c("#f2f2f2", "#004080"), 
         color = colfunc(10),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         border_color = "grey", 
         #border_color = NA,
         treeheight_row = 0,
         fontsize = 12,
         legend = TRUE,
         cellheight = 15,
         cellwidth=15,
         main="CD8 correlation with EM")
dev.off()


### Perform Differential gene expression analysis ####
#HC start here

## Control groups ##
Brain          <- c("NS.cerebellum","NS.cerebrum")
Heart          <- c("NS.heart", "NS.heart")
Kidney         <- c("NS.kidney")
Liver          <- c("NS.liver")
Lung           <- c("NS.lung")
germline       <- c("NS.testis","NS.ovary")
vitalNormals   <- c("NS.heart","NS.kidney","NS.liver","NS.lung")
vital.Brain.Normals   <- c("NS.cerebellum","NS.cerebrum", "NS.heart","NS.kidney","NS.liver","NS.lung")
othersNormals  <- c("NS.adrenalgland","NS.bladder","NS.colon","NS.ileum","NS.ovary","NS.pancreas","NS.prostate", 
                    "NS.skeletalmuscle","NS.spleen", "NS.stomach","NS.testis", "NS.ureter", "NS.uterus")
Normals        <- c("NS.adrenalgland","NS.bladder","NS.cerebellum","NS.cerebrum","NS.colon","NS.heart",
                    "NS.ileum","NS.kidney","NS.liver","NS.lung","NS.ovary","NS.pancreas","NS.prostate", 
                    "NS.skeletalmuscle","NS.spleen", "NS.stomach","NS.testis", "NS.ureter", "NS.uterus")
NormalsNoGermLine <- c("NS.adrenalgland","NS.bladder","NS.cerebellum","NS.cerebrum","NS.colon","NS.heart",
                       "NS.ileum","NS.kidney","NS.liver","NS.lung","NS.pancreas","NS.prostate", 
                       "NS.skeletalmuscle","NS.spleen", "NS.stomach", "NS.ureter", "NS.uterus")

tumorSubStatus.polyA    <- c("RMS.FP" , "RMS.FN", "EWS" ,"ASPS", "DSRCT", "HBL", "ML", "NB.MYCN.NA","NB.MYCN.A", "NB.Unknown", "OS", 
                             "SS", "Teratoma" ,"UDS" ,"YST")
tumorSubStatus.ribozero <-  c("WT" ,"CCSK")
Tumors                  <-  c("ASPS","DSRCT", "EWS" ,"HBL", "ML", "NB" ,"OS", "RMS", "SS", "Teratoma" ,"UDS" ,"YST","WT", "CCSK")

#HC: no WT and CCSK here
### Intantiate a new Differential Gene Expression Object ####
dgeObj  <- DifferentialGeneExp$new(
  countObj          = expressionObj$edgeRMethod("NormFactorDF")$counts,
  group1            = list(list("Normals"=NormalsNoGermLine,each=FALSE)),
  group2            = list(list("Tumor"=tumorSubStatus.polyA, each=TRUE)),
  packageRNAseq     = "edgeR",
  groupColumnName   = rnaseqProject$factorName,
  metadataDF        = rnaseqProject$metaDataDF,
  samplesColumnName = "Sample.Biowulf.ID.GeneExp",
  expressionUnit    = "TMM-RPKM",
  featureType       = "Gene",
  writeFiles        = TRUE,
  fileDirs          = rnaseqProject$fileDirs,
  subsetGenes       = TRUE,
  corUtilsFuncs     = corUtilsFuncs 
)
 

DiffExpObj <- dgeObj$performDiffGeneExp()

### Filtering ####
# Step 0   Define fucnctions ####

getCountObjTXT <- function(fileName, colNumb=1, rowNames=1){
  print(paste(fileName))
  featureCountTxt <- read.csv(fileName, sep="\t", row.names = "GeneName", header = 1);
  return(featureCountTxt[,colNumb, drop=FALSE])
}

mergeDiffTestResults <- function(x, type="", saveDirPath="", extension="", colInterest=1, rowNamesCol =1,
                                 fileSuffix=".txt"){
  
  print(paste(x))
  file_Dir_Gene = x
  fileName <- basename(file_Dir_Gene) 
  dir.create(file.path(paste(saveDirPath,fileName,sep="/")))
  GeneFiles             <- list.files(file_Dir_Gene); GeneFiles <- GeneFiles[grep(fileSuffix, GeneFiles)]
  GeneFilesList         <- paste(file_Dir_Gene, "/", GeneFiles,sep="") ; length(GeneFilesList)
  
  countObj          <- do.call(cbind,lapply(GeneFilesList, getCountObjTXT, colNumb=colInterest, rowNames=rowNamesCol))
  countObj_print    <- countObj %>% tibble::rownames_to_column(var="GeneName")
  
  write.table(countObj_print, paste(saveDirPath, paste(fileName, "/", type, ".DiffExp.txt",sep=""), sep= "/"), sep="\t",
              row.names = FALSE, quote = FALSE)
}
# Step 1.  Set the filters and annotation ####

## Javed's Filter for all three categories
#group2FPKM.T = 1 ; group1FPKM.T = 1;  PValue.T = 0.05 ; logFoldDiff.T = 0 ; FDR_value.T = 0.05 ; vitalFPKM.T = 0
group2FPKM.T = 5 ; group1FPKM.T = 1;  PValue.T = 0.00001 ; logFoldDiff.T = 4 ; FDR_value.T = 0.05 ; vitalFPKM.T = 1
Zscored.logFC = 0.25; Zscore.group2 = 0.5; group2FPKM = 5; group1FPKM = 1;  PValue = 0.001; logFC = 4; FDR = 0.05

file.copy("input/DiffExpResults/Normals_CCSK", "output/DiffExpResults", recursive=TRUE)
file.copy("input/DiffExpResults/Normals_WT", "output/DiffExpResults", recursive=TRUE)

selectedGeneLists <- c("CancerGermlineAntigen", "CellSurface", "TranscriptionFactor")

for (selectedGeneList in selectedGeneLists) {
  group_order <- read.table("input/DE_heatmap_order.txt")
  selected <- read.table(paste0("input/", selectedGeneList, "_GeneList.txt"))
  
  MergedDiffExpResultDir <- paste0("output/MergedDiffExpResults/",selectedGeneList)
  
  # Step 2.  Perform Merging of differential expression file across groups ####
  dir.create("output/MergedDiffExpResults", showWarnings = F)
  dir.create(MergedDiffExpResultDir, showWarnings = F)
  ConditionGroup <- c(unique(sapply(dgeObj$pairedList, function(x){ return(paste(x[1],x[2],sep = "_"))  })), c("Normals_WT", "Normals_CCSK") )
  #ConditionGroup <- c(unique(sapply(dgeObj$pairedList, function(x){ return(paste(x[1],x[2],sep = "_"))  })))
  groups <- list.dirs(paste("output/DiffExpResults/", sep=""))[-1]; groups[1]
  output <- sapply(groups, mergeDiffTestResults, type="Gene", colInterest=c(7,9,10,11,12, 15:28), rowNamesCol = 2,
                   fileSuffix=paste0(selectedGeneList,".txt"), saveDirPath=MergedDiffExpResultDir)
  
  # Step 3.  Core Function and save files ####
  allTumorStats <- do.call(cbind, lapply(ConditionGroup, function(x){
    diff_file <- paste(MergedDiffExpResultDir,"/",x,"/Gene.DiffExp.txt",sep="")
    print(diff_file)
    tumorData <- read.csv( diff_file, sep="\t", header = T, stringsAsFactors = FALSE ) 
    
    ## Actual filtering
    groupsCompare <- unlist(strsplit(x, "_"))
    print(groupsCompare)
    filterDFByColNames <- c("logFC",           groupsCompare[2], groupsCompare[1])
    newColNames <- paste0("Zscored.",c("logFC",  groupsCompare[2], groupsCompare[1]))
    
    ##Zscoreing matrix
    tumorDataPvalue        <- tumorData ; print(dim(tumorData))
    tumorDataPvalue_Zscore <- apply(tumorDataPvalue[,c("logFC",              groupsCompare[2], groupsCompare[1])],2,corUtilsFuncs$zscore_All)
    colnames(tumorDataPvalue_Zscore) <- newColNames;
    tumorDataPvalue_Zscore <- cbind(tumorDataPvalue[,c("GeneName"),drop = FALSE], tumorDataPvalue_Zscore ) %>% data.frame()
    tumorAllData <- left_join(tumorDataPvalue_Zscore, tumorDataPvalue, by="GeneName") ; 
    
    print(paste("Dim of", diff_file))
    dim(tumorAllData)
    
    ## Zscore Ranking filter
    tumorAllData.zscore <- tumorAllData %>% 
      dplyr::filter_(.dots=paste0( 
        groupsCompare[2]," >= ", group2FPKM ,
        " &  Zscored.logFC   >= ", Zscored.logFC,
        " &   logFC >", logFC,
        " & ", paste0("Zscored.",groupsCompare[2]), " >= ", Zscore.group2 )) %>%
      dplyr::arrange_(.dots = paste0("desc(","Zscored.",groupsCompare[2], ")" ) )
    
    ## Complete filtered gene List with zscoring filter
    ## write.table(tumorAllData.zscore, paste(MergedDiffExpResultDir,"/",x,"/",x,".filteredgenes.ZscoringRank.txt",sep=""), sep="\t", row.names = FALSE, quote = FALSE)
    
    ## Complete filtered gene List (Binary) with zscoring filter
    selectGenes <- tumorAllData.zscore %>%  dplyr::select(GeneName)
    tumorData.zscoreRanking <- tumorDataPvalue
    tumorData.zscoreRanking["status"] <- 0
    print(head(tumorData.zscoreRanking))
    statusDF.zscoreRanking <- tumorData.zscoreRanking %>% mutate(status=ifelse(GeneName %in% selectGenes$GeneName, 1, 0)) %>% dplyr::select(GeneName, status) %>% plyr::rename(c('status'=paste(groupsCompare[2],groupsCompare[1],"Status", sep="")))
    statusDF.zscoreRanking <- statusDF.zscoreRanking[, !duplicated(colnames(statusDF.zscoreRanking))] 
    
    # write.table(statusDF.zscoreRanking, paste(MergedDiffExpResultDir,"/",selectedGeneList,".Summarised.zscoreRank.txt",sep=""), sep="\t", row.names = FALSE, quote = FALSE)
    
    
    ## Javed's Filtering
    print(colnames(tumorAllData))
    tumorAllData.filt <- tumorAllData %>% dplyr::filter_(.dots=paste0( groupsCompare[2], " >= ", group2FPKM.T ,
                                                                       " &  logFC >=", logFoldDiff.T,
                                                                       " &  PValue   <=", PValue.T ,
                                                                       " &  Brain.MeanExp  < ", vitalFPKM.T ,
                                                                       " &  Heart.MeanExp  < ", vitalFPKM.T   )) %>%
      dplyr::arrange_(.dots = paste0("desc(","Zscored.",groupsCompare[2], ")" ) )
    
    ## Complete filtered gene List with traditional filter
    write.table(tumorAllData.filt, paste(MergedDiffExpResultDir,"/",x,"/",x,".filteredgene.traditionalRank.txt",sep=""), sep="\t", row.names = FALSE, quote = FALSE)
    
    ## Complete filtered gene List (Binary) with traditional filter 
    selectGenes <- tumorAllData.filt %>%  dplyr::select(GeneName)
    tumorData.traditionalRanking <- tumorDataPvalue
    tumorData.traditionalRanking["status"] <- 0
    head(tumorData.traditionalRanking)
    statusDF.traditionalRanking <- tumorData.traditionalRanking %>% mutate(status=ifelse(GeneName %in% selectGenes$GeneName, 1, 0)) %>% dplyr::select(GeneName, status) %>% plyr::rename(c('status'=paste(groupsCompare[2],groupsCompare[1],"Status", sep="")))
    statusDF.traditionalRanking <- statusDF.traditionalRanking[, !duplicated(colnames(statusDF.traditionalRanking))] 
    
    write.table(statusDF.traditionalRanking, paste(MergedDiffExpResultDir,"/",selectedGeneList,".Summarised.traditionalRank.txt",sep=""), sep="\t", row.names = FALSE, quote = FALSE)
    
    
    return(list(statusDF.zscoreRanking,statusDF.traditionalRanking))
  }))
  
  allInOneFile <- do.call(rbind, lapply(ConditionGroup, function(x){
    tumorData <- read.csv( paste(MergedDiffExpResultDir,"/",x,"/Gene.DiffExp.txt",sep=""), sep="\t", header = T, stringsAsFactors = FALSE ) 
    colnames(tumorData)[6] <- "Tumor"
    tumorData$Group <- x
    return(tumorData)
  })); dim(allInOneFile)
  write.table(allInOneFile, paste(MergedDiffExpResultDir,"/",selectedGeneList,".allSamples.txt",  sep=""),
              sep="\t", row.names = FALSE, quote = FALSE)
  
  #tumorStatusDF <- allTumorStats[, !duplicated(colnames(allTumorStats))]  %>% mutate(RowSum= rowSums(.[-1]))
  
  # Step 4.  Merge multiple DFs, memo-sort each DF and plot ####
  #HC: we use traditional
  allTumorMergedStats <- lapply(1:nrow(allTumorStats), function(x){
    mergedDF <- do.call(cbind, allTumorStats[x,])
    mergedDF <- mergedDF[, !duplicated(colnames(mergedDF))]        %>% 
      tibble::column_to_rownames(var="GeneName")
    mergedDF.Memo        <- corUtilsFuncs$memoSort(M=mergedDF)
    mergedDF.Memo$RowSum <- apply(mergedDF.Memo, 1, function(x) sum(x!=0))
    mergedDF.Memo <- mergedDF.Memo %>% tibble::rownames_to_column(var="GeneName")
    return(mergedDF.Memo)
  })
  
  # Step 5.  Save the files ####
  # write.table(allTumorMergedStats[[1]], paste(MergedDiffExpResultDir,"/",selectedGeneList,".Summarised.ZscoreRank.Dexp.txt",  sep=""),
  #             sep="\t", row.names = FALSE, quote = FALSE)
  write.table(allTumorMergedStats[[2]], paste(MergedDiffExpResultDir,"/",selectedGeneList,".Summarised.traditionalRank.Dexp.txt",sep=""),
              sep="\t", row.names = FALSE, quote = FALSE)
  
  # Step 6.  Select rows for heatmap ####
  #allTumorStatsFinal <- read.table(paste(MergedDiffExpResultDir,"/",selectedGeneList,".Summarised.ZscoreRank.Dexp.txt",sep=""),sep="\t", header = TRUE)
  allTumorStatsFinal <- read.table(paste(MergedDiffExpResultDir,"/",selectedGeneList,".Summarised.traditionalRank.Dexp.txt",sep=""),sep="\t", header = TRUE)
  CTA.Filt <- allTumorStatsFinal %>% filter(RowSum >= 1) %>% 
    dplyr::arrange(-RowSum) %>% t() %>% data.frame()
  colnames(CTA.Filt) <- as.character(unlist(CTA.Filt[c("GeneName"),]))
  CTA.Filt.sorted <- CTA.Filt[-1,]
  CTA.Filt.sorted <- CTA.Filt.sorted %>% tibble::rownames_to_column("Diagnosis")
  CTA.Filt.sorted <- CTA.Filt.sorted %>% dplyr::arrange_(.dots = list(paste0("desc(",colnames(CTA.Filt.sorted)[2], ")")))
  rownames(CTA.Filt.sorted) <- CTA.Filt.sorted[,1]
  CTA.Filt.sorted <- CTA.Filt.sorted[-1,]
  CTA.Filt.sorted <- CTA.Filt.sorted[,-1]
  #filter(, GeneName %in% c("CD99", "FGFR4", "ALK", "GPC2", "MYCN", "MYOG", "MYOD1", "IGF2", "CTAG1B"))
  #View(CTA.Filt.sorted);dim(CTA.Filt.sorted)
  
  # Step 7.  Plot the heatmap ####
  rownames(CTA.Filt.sorted) <- gsub("NormalsStatus", "", rownames(CTA.Filt.sorted))
  indx <- sapply(CTA.Filt.sorted, is.character)
  CTA.Filt.sorted[indx] <- lapply(CTA.Filt.sorted[indx], function(x) as.numeric(as.character(x)))
  
  #HC: filter and re-order the matrix
  
  CTA.Filt.sorted <- CTA.Filt.sorted[match(group_order$V1, rownames(CTA.Filt.sorted)),]
  CTA.Filt.sorted <- as.data.frame(CTA.Filt.sorted) %>% dplyr::select(any_of(selected$V1))
  
  h = 6
  w = 12
  if (selectedGeneList == "CellSurface") {
    CTA.Filt.sorted <- t(CTA.Filt.sorted)
    h = 18
    w = 8
  }
  if (selectedGeneList == "TranscriptionFactor") {
    CTA.Filt.sorted <- t(CTA.Filt.sorted)
    h = 15
    w = 8
  }
  
  dim(CTA.Filt.sorted)
  
  pdf( paste(rnaseqProject$workDir, rnaseqProject$plotsDir, 
             paste0("Figure4_", selectedGeneList, ".pdf"),  sep="/"), height = h, width = w)
  pheatmap(CTA.Filt.sorted, 
           #color =c("#e0e0d1", "#004080"), 
           color =c("#f2f2f2", "#004080"), 
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           border_color = "grey", 
           #border_color = NA,
           angle_col = 45,
           treeheight_row = 0,
           fontsize = 10,
           cellwidth = 12,
           cellheight = 12,
           legend = FALSE )
  dev.off()
}
### Figure 4-ABC ends here

### Violin plot for exhaustion markers ####
toPlotDF <- readRDS("input/toplotDF.rds")
RPKM.Data.Exhaustion   <- expressionTMM.RPKM.arr %>% dplyr::filter(GeneName %in% c("C10orf54",as.character(rnaseqProject$emDF$GeneName))) %>% 
  dplyr::arrange(GeneName)
Exhaustion.Transpose           <- as.data.frame(t(RPKM.Data.Exhaustion[,-c(1:7)]))
colnames(Exhaustion.Transpose) <- RPKM.Data.Exhaustion$GeneName
#HC: use SampleID.in.paper and replace X0 with O
var_name <- "Sample.Biowulf.ID.GeneExp"
var_name <- "SampleID..In.Paper"
Exhaustion.Transpose <- Exhaustion.Transpose %>% tibble::rownames_to_column(var=var_name)
Exhaustion.Transpose[[var_name]] <- gsub("X0", "0", Exhaustion.Transpose[[var_name]])
Exhaustion.Transpose.diag <- dplyr::left_join(Exhaustion.Transpose, rnaseqProject$metaDataDF[,c(var_name, "DIAGNOSIS.Substatus.Tumor.Normal.Tissue")], 
                                              by=var_name) %>% dplyr::rename(Diagnosis=DIAGNOSIS.Substatus.Tumor.Normal.Tissue)
Rm.Normal.Exhaustion.Transpose <- Exhaustion.Transpose.diag %>% filter(!grepl("^NS.*", Exhaustion.Transpose.diag$Diagnosis))
finalExhaustionMatrix <-  melt(Rm.Normal.Exhaustion.Transpose[,-1], id.var = "Diagnosis")

finalExhaustionMatrix.tidy <- finalExhaustionMatrix %>% dplyr::group_by(Diagnosis, variable) %>% 
  dplyr::mutate(Med=median(value)) %>% arrange(Diagnosis, variable, value) %>% 
  arrange(desc(Med)) %>% 
  ungroup() %>% 
  mutate( Diagnosis.Marker = factor(paste(Diagnosis,variable,sep="."), levels= unique(paste(Diagnosis,variable,sep=".")),
                                    order = TRUE) ) #%>% 
#arrange(Diagnosis.Marker)

pdf("output/Figures/Figure2A&S3.pdf",height=25,width=20)
customColorsVector <- setNames(unique(as.character(toPlotDF$Color.Jun)),unique(as.character(toPlotDF$Diagnosis)) )
ggplot(finalExhaustionMatrix.tidy, aes(x=Diagnosis.Marker, y=log2(value+1) , fill=Diagnosis)) + 
  ##ggplot(data, aes(x=Group, y=log2(ENSG00000182752))) + 
  #geom_boxplot(varwidth = TRUE,notch = FALSE) + 
  geom_violin(scale = "width",trim = T) + 
  stat_summary(fun=median.quartile,geom='point') +
  #geom_jitter(width=0.1) +
  theme_bw() + 
  ylab( paste("log2(FPKM)") ) +
  #ylim(0,10) +
  xlab( "Diagnosis" ) +
  #geom_hline(yintercept=0, size=0.1) + 
  scale_fill_manual(values = customColorsVector) + 
  theme( title = element_text(size=13, face="bold")
         ,axis.title.x = element_text(size=13, face="bold")
         ,axis.title.y = element_text(size=13, face="bold")
         ,axis.text.x = element_text(size=10, face="bold", angle=90, vjust=1)
         ,axis.text.y = element_text(size=10, face="bold")
         ,axis.ticks.x =element_blank()
         ,strip.text.y= element_blank()
         ,strip.text.x=element_text(size=13,face="bold")
         ,strip.background=element_blank()
         ,panel.grid.major.x=element_blank()
         ,panel.grid.minor.x=element_blank()
         ,panel.border = element_rect(colour = "black", fill=NA, size=0.0000000002, linetype = 2)
         ,panel.spacing = unit(0, "cm")
         ,strip.switch.pad.grid = unit(0, "cm")
  ) + facet_wrap( ~ variable, scales="free", nrow = 7) +
  scale_x_discrete(labels=setNames(as.character(finalExhaustionMatrix.tidy$Diagnosis), finalExhaustionMatrix.tidy$Diagnosis.Marker))
dev.off()


#HC: Figure S7 needs cellline
rnaseqProjectAll <- createRNAseqProject(project_name, meta_file, list('None'=list("LIBRARY_TYPE"="")))

### Violin plot for HLA-A,HLA-B,HLA-C ####################################################

#RPKM.Data.Exhaustion   <- expressionTMM.RPKM %>% dplyr::filter(GeneName %in% c("HLA-A", "HLA-B", "HLA-C")) %>% dplyr::arrange(GeneName)

metaDataDFAll <- rnaseqProjectAll$metaDataDF

#HC: get expression data including cell lines

## Tumor Normal and Cellline
mergeObjectsNoDupAll_data <- readRDS("input/GeneRDS/RawCount/All.samples.Tumor.Normal.Celline.RDS")
mergeObjectsNoDupAll <- data.frame("GeneID"=rownames(mergeObjectsNoDupAll_data))
mergeObjectsNoDupAll <- cbind(mergeObjectsNoDupAll, mergeObjectsNoDupAll_data)
mergeObjectsNoDupAll <- dplyr::left_join(rnaseqProjectAll$annotationDF[,c("GeneID", "GeneName")], mergeObjectsNoDupAll, by="GeneID") %>% data.table()
mergeObjectsConsoAll <- corUtilsFuncs$consolidateDF(mergeObjectsNoDupAll[,-c("GeneID")], funcName = "max", featureName = "GeneName")
mergeObjectsConsoAll <- dplyr::full_join(mergeObjectsConsoAll, rnaseqProjectAll$annotationDF[,c("GeneID", "GeneName")], by="GeneName") %>% data.table()
mergeObjectsConsoAll <- subset(mergeObjectsConsoAll,!duplicated(mergeObjectsConsoAll$GeneName))
mergeObjectsConsoAll <- mergeObjectsConsoAll[complete.cases(mergeObjectsConsoAll), ]; dim(mergeObjectsConsoAll)
mergeObjectsConsoAll <- mergeObjectsConsoAll[,-c("GeneName")] %>% data.frame() %>% tibble::column_to_rownames(var = "GeneID") %>% 
  as.matrix() ; dim(mergeObjectsConsoAll)

rnaseqProjectAll$annotationDF <- rnaseqProjectAll$annotationDF %>% dplyr::filter(GeneID %in% rownames(mergeObjectsConsoAll));
design <- rnaseqProjectAll$metaDataDF
expressionObjAll        <- GeneExpNormalization$new(
  countObj          = as.matrix(mergeObjectsConsoAll), 
  featureType       = "Gene", 
  packageRNAseq     = "edgeR", 
  annotationDF      = rnaseqProjectAll$annotationDF, 
  design            = design[,rnaseqProjectAll$factorName], 
  proteinCodingOnly = FALSE,
  corUtilsFuncs     = corUtilsFuncs
)

expressionTMM.RPKM.All = expressionObjAll$edgeRMethod("TMM-RPKM", logtransform = TRUE, zscore = FALSE)

plotMultiGG <- function(ps, filename, nrow=3, ncol=3, height=12, width=18) {
  nplots <- nrow * ncol
  npages <- length(ps) %/% (nplots)
  if (length(ps) %% (nplots) > 0 )
    npages <- npages + 1
  pdf(filename, height=height, width=width)
  for (i in c(1:npages)) {
    sidx <- (nplots*(i-1)+1)
    eidx <- nplots*i
    if (eidx > length(ps))
      eidx <- length(ps)
    #print(sidx)
    #print(eidx)
    do.call(grid.arrange, c(ps[sidx:eidx], ncol=ncol, nrow=nrow))
    #grid.newpage()
  }  
  dev.off()
}

flags <- c(T,F)
gene_list <- c("GPC2","PRAME","FOXM1")

ps <- list()
ps_nocl <- list()
bCellline <- F
for (bCellline in flags) {
  for (gene in gene_list) {
    print(gene)
    if (bCellline) {
      RPKM.Data.Exhaustion   <- expressionTMM.RPKM.All %>% dplyr::filter(GeneName %in% c(gene)) %>% dplyr::arrange(GeneName)
    } else {
      RPKM.Data.Exhaustion   <- expressionTMM.RPKM %>% dplyr::filter(GeneName %in% c(gene)) %>% dplyr::arrange(GeneName)
    }
    Exhaustion.Transpose           <- as.data.frame(t(RPKM.Data.Exhaustion[,-c(1:7)]))
    colnames(Exhaustion.Transpose) <- RPKM.Data.Exhaustion$GeneName
    Exhaustion.Transpose <- Exhaustion.Transpose %>% tibble::rownames_to_column(var="Sample.Biowulf.ID.GeneExp")
    if (bCellline) {
      Exhaustion.Transpose$Sample.Biowulf.ID.GeneExp <- gsub("COG\\.N\\.", "COG-N-", Exhaustion.Transpose$Sample.Biowulf.ID.GeneExp)
      Exhaustion.Transpose$Sample.Biowulf.ID.GeneExp <- gsub("NB\\.", "NB-", Exhaustion.Transpose$Sample.Biowulf.ID.GeneExp)
      Exhaustion.Transpose$Sample.Biowulf.ID.GeneExp <- gsub("SMS\\.", "SMS-", Exhaustion.Transpose$Sample.Biowulf.ID.GeneExp)
      Exhaustion.Transpose$Sample.Biowulf.ID.GeneExp <- gsub("CHP\\.", "CHP-", Exhaustion.Transpose$Sample.Biowulf.ID.GeneExp)
      Exhaustion.Transpose$Sample.Biowulf.ID.GeneExp <- gsub("X6647", "6647", Exhaustion.Transpose$Sample.Biowulf.ID.GeneExp)
      Exhaustion.Transpose$Sample.Biowulf.ID.GeneExp <- gsub("X7556", "7556", Exhaustion.Transpose$Sample.Biowulf.ID.GeneExp)
      Exhaustion.Transpose$Sample.Biowulf.ID.GeneExp <- gsub("Felix\\.", "Felix-", Exhaustion.Transpose$Sample.Biowulf.ID.GeneExp)
      Exhaustion.Transpose$Sample.Biowulf.ID.GeneExp <- gsub("LA\\.N\\.", "LA-N-", Exhaustion.Transpose$Sample.Biowulf.ID.GeneExp)
      Exhaustion.Transpose$Sample.Biowulf.ID.GeneExp <- gsub("NBL\\.", "NBL-", Exhaustion.Transpose$Sample.Biowulf.ID.GeneExp)
      metaDataDFAll$FigS7 <- metaDataDFAll$DIAGNOSIS.Substatus.Tumor.Tissue
      metaDataDFAll$FigS7[which(metaDataDFAll$DIAGNOSIS.Substatus.Tumor.Tissue=="NS.Normal")] <- metaDataDFAll$DIAGNOSIS.Substatus.Normal.Tissue[which(metaDataDFAll$DIAGNOSIS.Substatus.Tumor.Tissue=="NS.Normal")]
      Exhaustion.Transpose.diag <- dplyr::inner_join(Exhaustion.Transpose, metaDataDFAll[,c("Sample.Biowulf.ID.GeneExp", 
                                                                                            "FigS7",
                                                                                            "Color.Jun",
                                                                                            "LIBRARY_TYPE")], 
                                                     by="Sample.Biowulf.ID.GeneExp") %>% dplyr::rename(Diagnosis=FigS7)
      
      Rm.Normal.Exhaustion.Transpose <- Exhaustion.Transpose.diag %>% filter(!grepl("NS.*", Exhaustion.Transpose.diag$Diagnosis))
      Normal.Exhaustion.Transpose <- Exhaustion.Transpose.diag %>% filter(grepl("NS.*", Exhaustion.Transpose.diag$Diagnosis))
      
      finalExhaustionMatrix.tumor <-  reshape2::melt(Rm.Normal.Exhaustion.Transpose[,-1], id.var = c("Color.Jun","Diagnosis","LIBRARY_TYPE"))
      Normal.Exhaustion.Transpose$Color.Jun <- "steelblue"
      finalExhaustionMatrix.normal <-  reshape2::melt(Normal.Exhaustion.Transpose[,-1], id.var = c("Color.Jun","Diagnosis","LIBRARY_TYPE"))
      finalExhaustionMatrix.tumor.order <- finalExhaustionMatrix.tumor %>% dplyr::arrange(Diagnosis)
      finalExhaustionMatrix.tidy.tumor.order <- finalExhaustionMatrix.tumor.order %>% dplyr::group_by(Diagnosis, variable) %>% 
        dplyr::mutate(Med=median(value)) %>% arrange(Diagnosis, variable, value) %>% 
        #arrange(desc(Med)) %>% 
        ungroup() %>% 
        mutate( Diagnosis.Marker = factor(paste(Diagnosis,variable,sep="."), levels= unique(paste(Diagnosis,variable,sep=".")),
                                          order = TRUE),
                LIBRARY_TYPE = factor(LIBRARY_TYPE)) %>% 
        arrange(Diagnosis.Marker)
      finalExhaustionMatrix.tidy.normal <- finalExhaustionMatrix.normal %>% dplyr::group_by(Diagnosis, variable) %>% 
        dplyr::mutate(Med=median(value)) %>% arrange(Diagnosis, variable) %>%
        #arrange(desc(Med)) %>% 
        ungroup() %>% 
        mutate( Diagnosis.Marker = factor(paste(Diagnosis,variable,sep="."), order = TRUE),
                LIBRARY_TYPE = factor(LIBRARY_TYPE, levels=c("Normal") )) %>% 
        arrange(Diagnosis.Marker)
    } else {
      Exhaustion.Transpose.diag <- dplyr::inner_join(Exhaustion.Transpose, rnaseqProject$validMetaDataDF[,c("Sample.Biowulf.ID.GeneExp",
                                                                                                            "Violin.normal", 
                                                                                                            "Color.Jun",  
                                                                                                            "LIBRARY_TYPE")], 
                                                     by="Sample.Biowulf.ID.GeneExp") %>% dplyr::rename(Diagnosis=Violin.normal) 
      Rm.Normal.Exhaustion.Transpose <- Exhaustion.Transpose.diag %>% filter(!grepl("NS.*", Exhaustion.Transpose.diag$Diagnosis))
      Normal.Exhaustion.Transpose <- Exhaustion.Transpose.diag %>% filter(grepl("NS.*", Exhaustion.Transpose.diag$Diagnosis))
      finalExhaustionMatrix.tumor <-  reshape2::melt(Rm.Normal.Exhaustion.Transpose[,-1], id.var = c("Color.Jun","Diagnosis","LIBRARY_TYPE"))
      Normal.Exhaustion.Transpose$Color.Jun <- "steelblue"
      finalExhaustionMatrix.normal <-  reshape2::melt(Normal.Exhaustion.Transpose[,-1], id.var = c("Color.Jun","Diagnosis","LIBRARY_TYPE"))
      
      
      finalExhaustionMatrix.tumor.order <- finalExhaustionMatrix.tumor %>% filter(grepl('Tumor', LIBRARY_TYPE))
      finalExhaustionMatrix.tumor.order$value <- as.numeric(as.character(finalExhaustionMatrix.tumor.order$value))
      finalExhaustionMatrix.tidy.tumor.order <- finalExhaustionMatrix.tumor.order %>% dplyr::group_by(Diagnosis, variable) %>% 
        dplyr::mutate(Med=median(value)) %>% arrange(Diagnosis, variable, value) %>% 
        arrange(desc(Med)) %>% 
        ungroup() %>% 
        mutate( Diagnosis.Marker = factor(paste(Diagnosis,variable,sep="."), levels= unique(paste(Diagnosis,variable,sep=".")),
                                          order = TRUE),
                LIBRARY_TYPE = factor(LIBRARY_TYPE, levels=c("Tumor") )) %>% 
        arrange(Diagnosis.Marker)
      
      # For Normal ordering
      normalS_order = c("NS.brain", "NS.heart", "NS.lung", "NS.liver", "NS.kidney", "NS.testis", "NS.ovary", "NS.other")
      normal_levels = unlist(lapply(levels(finalExhaustionMatrix.normal$variable), function(x) { paste(normalS_order,x,sep="." )} ) )
      finalExhaustionMatrix.tidy.normal <- finalExhaustionMatrix.normal %>% dplyr::group_by(Diagnosis, variable) %>% 
        dplyr::mutate(Med=median(value)) %>% arrange(Diagnosis, variable) %>%
        arrange(desc(Med)) %>% 
        ungroup() %>% 
        mutate( Diagnosis.Marker = factor(paste(Diagnosis,variable,sep="."), levels = normal_levels,order = TRUE),
                LIBRARY_TYPE = factor(LIBRARY_TYPE, levels=c("Normal") )) %>%         arrange(Diagnosis.Marker)
    }
    ## Merge both dataframes 
    finalExhaustionMatrix.tidy <- rbind(finalExhaustionMatrix.tidy.tumor.order, finalExhaustionMatrix.tidy.normal)
    #pdf(paste0("Figures/FPKM.v2.",gene,".Normal.violin.pdf"),height=10,width=20)
    #pdf("Figures/FigureS7.pdf")
    if (bCellline) {
      #finalExhaustionMatrix.tidy$Diagnosis.Marker <- gsub("NA\\.", "", finalExhaustionMatrix.tidy$Diagnosis.Marker)
      customColorsVector.tumor <- finalExhaustionMatrix.tidy.tumor.order %>% dplyr::group_by(Diagnosis) %>% dplyr::summarise(Color.Jun=max(Color.Jun))
      customColorsVector.tumor <- setNames(customColorsVector.tumor$Color.Jun, customColorsVector.tumor$Diagnosis)
      # Normal
      customColorsVector.normal <- setNames(rep("lightgrey", length(unique(as.character(finalExhaustionMatrix.tidy.normal$Diagnosis)))), 
                                            unique(as.character(finalExhaustionMatrix.tidy.normal$Diagnosis)))
      # merge
      customColorsVector <- c(customColorsVector.tumor, customColorsVector.normal); customColorsVector
      customColorsVector.dummy <- setNames(rep("white", length(unname(customColorsVector))),
                                           names(customColorsVector))
      p <- ggplot(finalExhaustionMatrix.tidy, aes(x=Diagnosis.Marker, y=value , fill=Diagnosis)) + 
        geom_violin(scale = "width",trim = T) + 
        scale_fill_manual(values = customColorsVector) + 
        stat_summary(fun=median.quartile,geom='point') +
        theme_bw() + 
        ylab( paste("log2(FPKM+1)") ) +
        xlab( "Diagnosis" ) +
        theme( title = element_text(size=13, face="bold")
               ,axis.title.y = element_text(size=13, face="bold")
               ,axis.text.x = element_text(size=8, face="bold", angle=90, hjust=0.95, vjust=0.2)
               ,axis.text.y = element_text(size=8, face="bold")
               ,panel.grid.major.x=element_blank()
               ,panel.grid.minor.x=element_blank()
               ,panel.spacing = unit(0, "cm")
               ,strip.text.y= element_blank()
               ,strip.text.x=element_text(size=13,face="bold")
               ,strip.background=element_blank()
               ,legend.position = "none"
        ) + facet_wrap( ~ variable , scales="free", nrow = 7) +
        scale_x_discrete(labels=setNames(as.character(finalExhaustionMatrix.tidy$Diagnosis), finalExhaustionMatrix.tidy$Diagnosis.Marker)) #+
      ps[[length(ps)+1]] <- p
      
    } else {
      ## Construct Color vector
      #customColorsVector <- setNames(unique(as.character(finalExhaustionMatrix$Color.Jun)), unique(as.character(finalExhaustionMatrix$Diagnosis)))
      # Tumor
      customColorsVector.tumor <- setNames(unique(as.character(finalExhaustionMatrix.tidy.tumor.order$Color.Jun)), 
                                           unique(as.character(finalExhaustionMatrix.tidy.tumor.order$Diagnosis)))
      # Normal
      customColorsVector.normal <- setNames(rep("lightgrey", 8), 
                                            unique(as.character(finalExhaustionMatrix.tidy.normal$Diagnosis)))
      # merge
      customColorsVector <- c(customColorsVector.tumor, customColorsVector.normal); customColorsVector
      customColorsVector.dummy <- setNames(rep("white", length(unname(customColorsVector))),
                                           names(customColorsVector))
      p <- ggplot(finalExhaustionMatrix.tidy, aes(x=Diagnosis.Marker, y=value , fill=Diagnosis)) + 
        ##ggplot(data, aes(x=Group, y=log2(ENSG00000182752))) + 
        geom_violin(scale = "width",trim = T) + 
        stat_summary(fun=median.quartile,geom='point') +
        #geom_point(data = test.dot.DF2, aes(x = Diagnosis.Marker, y = value), fill="white", alpha=0) +
        #geom_boxplot(width=0.2,notch = FALSE) + 
        #geom_jitter(width=0.1,color="lightgrey", ) +
        theme_bw() + 
        ylab( paste("log2(FPKM+1)") ) +
        #ylim(-2,4) +
        xlab( "Diagnosis" ) +
        #geom_hline(yintercept=0, size=0.1) + 
        scale_fill_manual(values = customColorsVector) + 
        theme( title = element_text(size=13, face="bold")
               ,axis.title.x = element_text(size=13, face="bold")
               ,axis.title.y = element_text(size=13, face="bold")
               ,axis.text.x = element_text(size=10, face="bold", angle=90, vjust=1)
               ,axis.text.y = element_text(size=10, face="bold")
               ,axis.ticks.x =element_blank()
               ,strip.text.y= element_blank()
               ,strip.text.x=element_text(size=13,face="bold")
               ,strip.background=element_blank()
               ,panel.grid.major.x=element_blank()
               ,panel.grid.minor.x=element_blank()
               ,panel.border = element_rect(colour = "black", fill=NA, size=0.0000000002, linetype = 2)
               ,panel.spacing = unit(0, "cm")
               ,strip.switch.pad.grid = unit(0, "cm")
               ,legend.position = "none"
        ) + facet_wrap( ~ variable , scales="free", nrow = 7) +
        scale_x_discrete(labels=setNames(as.character(finalExhaustionMatrix.tidy$Diagnosis), finalExhaustionMatrix.tidy$Diagnosis.Marker)) #+
        ps_nocl[[length(ps_nocl)+1]] <- p
    }
    #dev.off()
  }
}
plotMultiGG(ps, "output/Figures/FigureS7.pdf", nrow=3, ncol=1, height=12, width=8)
plotMultiGG(ps_nocl, "output/Figures/Figure4D.pdf", nrow=3, ncol=1, height=12, width=8)

#####################################################################################################

## TCR

metaData <- read.csv(meta_tcr_file, sep="\t")

emptyDF <- data.frame(count=c(0), freq=c(0), cdr3nt=c("NA"),cdr3aa=c("NF"),v=c("NF"),d=c("NF"),j=c("NF"),VEnd=c(0),DStart=c(0),
                      DEnd=c(0),JStart=c(0),SampleName=c(0))
emptyDFEntropy <- data.frame(VJcombo=c(), Counts =c(), Vcassette=c(), Jcassette=c(), aaCDR3_filtered = c(), ntCDR3= c())

emptyDFEntropyResults <- data.frame(FileName=c(), Hcdr3 =c(), Htot=c(), CLcdr3=c(), CLHvj= c(), CLtot= c(),
                                    Hcdr3_max=c(), Hvj_max =c(), Htot_max=c(), CLcdr3_max=c(), Num_CDR3= c(), Num_VJ= c(),
                                    Num_totCDR3 =c())

fileList <- list.files("input/MiXCR/cloneFiles.v2/")

AllClonesData             <- rbindlist( lapply(fileList, function(x){
  print(x)
  exomeData <- read.csv( paste("input/MiXCR/cloneFiles.v2/", x, sep=""), sep="\t", header = TRUE )
  if(nrow(exomeData)>0){
    exomeData$SampleName <- x
  } else {
    emptyDF$SampleName <- c(x)
    exomeData <- emptyDF
  }
  return(exomeData)
}) )

### Filter Clones by clone types ####
cloneObjIG                <- filterSpecificCloneTypes(cloneData = AllClonesData, cloneType = "IGH")
cloneObjIG.Expansion.GE3  <- cloneObjIG %>%  dplyr::filter(grepl("IGH",v) & count >= 3)
cloneObjTCR               <- filterSpecificCloneTypes(cloneData = AllClonesData, cloneType = "TRB")
cloneObjTCR.Expansion.GE3  <- cloneObjTCR %>% dplyr::filter(grepl("TRB",v) & count >= 3)

### Read normalised counts for each sample
readCounts <- readRDS("input/MiXCR/RNASeq.readcounts.rds")
readCountsSum <- apply(readCounts, 2, sum)
readCountsSum <- as.data.frame(readCountsSum)
readCountsSum.df <- readCountsSum %>% tibble::rownames_to_column(var="Sample.Biowulf.ID.GeneExp")

cloneType = "TRB"
#AllClonesEntropyData             <- sapply(fileList, makeEntropyInput,  cloneType=cloneType, inputDir="input/MiXCR/CloneFiles.v2/",
#                                           outputDir=paste0("input/MiXCR/CloneFilesEntropy.", cloneType, ".v2/") )
#fileList <- list.files("input/MiXCR/immunoseqv2/")
#ImmunoseqV2EntropyData             <- sapply(fileList, immunoseqv2 )

#################################################################################### For RNASEq ################################################################
### For now Select DF manually ####
#cloneType = "IGHClones"  ; countObj <- cloneObjIG %>% as.data.frame()
cloneType = "TRBClones"  ; countObj <- cloneObjTCR %>% as.data.frame()

### Attach metadata and generate countObj ####
countObj <- countObj %>% dplyr::rename(Sample.Data.ID=SampleName); 
countObj$Sample.Data.ID <- gsub("convert.|.clones.txt","", countObj$Sample.Data.ID)
#countObj$SAMPLE_ID <- gsub("-","_", countObj$SAMPLE_ID)
countObj.Annot <- dplyr::left_join(countObj, metaData, by="Sample.Data.ID") %>% 
  dplyr::select(one_of("count", "freq", "cdr3nt", "cdr3aa", "v", "d", "j", "VEnd", "DStart", "DEnd", "JStart", "Sample.Biowulf.ID.GeneExp", "Sample.ID.Alias",
                       "LIBRARY_TYPE","DIAGNOSIS.Substatus.Tumor.Normal.Tissue", "Color.Jun", "Patient.ID.updated" )) ; 
### Plot the clone expansion
countObj.Annot.NoCL <- countObj.Annot %>% filter(!grepl('CellLine',LIBRARY_TYPE)) %>% filter(!grepl('^NS', DIAGNOSIS.Substatus.Tumor.Normal.Tissue) )
countObj.Annot.NoCL <- countObj.Annot.NoCL %>% dplyr::rename(Diagnosis = DIAGNOSIS.Substatus.Tumor.Normal.Tissue)


countObj.Annot.complete <- countObj.Annot.NoCL[complete.cases(countObj.Annot.NoCL),]
## sanity check
dim(countObj)
dim(countObj.Annot)
dim(countObj.Annot.NoCL)
dim(countObj.Annot.complete)

countObj.Annot.NoCL.totalReads <- dplyr::left_join(countObj.Annot.complete, readCountsSum.df, by="Sample.Biowulf.ID.GeneExp")
countObj.Annot.NoCL.totalReads.complete <- countObj.Annot.NoCL.totalReads[complete.cases(countObj.Annot.NoCL.totalReads),]
dim(countObj.Annot.NoCL.totalReads.complete)

countObj.Annot.NoCL.totalReads$ReadsPerMillion <- ( countObj.Annot.NoCL.totalReads$count/countObj.Annot.NoCL.totalReads$readCountsSum)*1000000

######## Make Step Plots to show expansion ####
countObj.Annot.NoCL.totalReads <- countObj.Annot.NoCL.totalReads %>% dplyr::select(Sample.Biowulf.ID.GeneExp, 
                                                                                   Diagnosis, 
                                                                                   count, freq, readCountsSum, ReadsPerMillion,
                                                                                   Color.Jun)
#Color.Substatus)

#toPlotDF <- countObj.Annot.NoCL.totalReads %>% dplyr::mutate(ReadsPerMillion = if_else(ReadsPerMillion >= 2, 2, ReadsPerMillion))

toPlotDF <- countObj.Annot.NoCL.totalReads %>% dplyr::filter(count > 0) %>% 
  arrange(Sample.Biowulf.ID.GeneExp, count) %>%
  group_by(Sample.Biowulf.ID.GeneExp) %>% 
  mutate(rank = dense_rank( -count )) %>% 
  distinct()
# mutate(good_ranks = order(order(order_values, decreasing=TRUE)))
#View(toPlotDF)
#saveRDS(toPlotDF, )

val = c("NB.MYCN.NA", "ASPS", "HBL", "NB.Unknown", "RMS.FP", "RMS.FN", "NB.MYCN.A", "UDS", "OS", "EWS", "DSRCT", "SS", "CCSK", "ML", "WT", "YST", "Teratoma")
toPlotDF$Diagnosis <- factor(toPlotDF$Diagnosis,levels = val, ordered = TRUE)
toPlotDF <- toPlotDF %>% dplyr::arrange(Diagnosis)
customColorsVector <- setNames( unique(as.character(toPlotDF$Color.Jun)), unique(as.character(toPlotDF$Diagnosis)) )
# %<>% dplyr::mutate(Color.Jun = factor(Color.Jun, levels = unique(Color.Jun), ordered = TRUE))
### For RNASeq in Reverse
pdf("output/Figures/Figure3C.pdf", height = 15, width = 25)
#HC: we are missing some columns
#ggplot(toPlotDF[,c(1,3,5,6,7)]) +
ggplot(toPlotDF) +
  geom_step(aes(x = rank, y = ReadsPerMillion, group=Sample.Biowulf.ID.GeneExp, colour= as.character(toPlotDF$Diagnosis) ), 
            size=0.7 ) +
  scale_colour_manual(values=customColorsVector) +
  facet_wrap(~toPlotDF$Diagnosis) +
  theme_bw() +
  theme(legend.position="none") +
  theme(strip.text=element_text(size=16, face = "bold"),
        axis.text = element_text(size=14, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(colour = "skyblue", fill=NA, size=1)) +
  scale_x_continuous( breaks = seq(1,max(toPlotDF$rank),by=4) ) +
  coord_trans(y = "log2" ) +
  scale_y_continuous(minor_breaks = c(),
                     breaks = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 3, 4),
                     labels = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 3, 4) 
  ) 
dev.off()

dir.create("output/MiXCR", showWarnings = F)
dir.create("output/MiXCR/Results", showWarnings = F)
countObj.Annot.gb.Samples <- countObj.Annot.NoCL  %>% dplyr::group_by(Sample.Biowulf.ID.GeneExp) %>% 
  dplyr::summarise(
    TotalClones=n(),
    TotalCloneSum=sum(count),
    Diagnosis= paste(unique(Diagnosis), collapse = ',' )
    #Diagnosis= paste(unique(DIAGNOSIS.Substatus.Tumor.Normal.Tissue), collapse = ',' )
  )  %>% 
  dplyr::rename_(.dots=setNames(list("TotalClones"),c(cloneType))); dim(countObj.Annot.gb.Samples); head(countObj.Annot.gb.Samples); tbl_df(countObj.Annot.gb.Samples)
countObj.Annot.gb.Samples.Annotate  <- left_join(countObj.Annot.gb.Samples, metaData[,c("Sample.Biowulf.ID.GeneExp","Sample.ID.Alias","LIBRARY_TYPE","DIAGNOSIS.Substatus.Tumor.Normal.Tissue")], 
                                                 by="Sample.Biowulf.ID.GeneExp")
###Replace NA by 0
countObj.Annot.gb.Samples.Annotate[which(is.na(countObj.Annot.gb.Samples.Annotate$TotalCloneSum)), c(cloneType,"TotalCloneSum")] <- 0
countObj.Annot.gb.Samples.Annotate[which(countObj.Annot.gb.Samples.Annotate$TotalCloneSum == 0), "IGHClones"] <- 0

### Saving files
saveRDS(countObj.Annot.gb.Samples, paste("output/MiXCR/Results/countObj.Annot.gb.Samples.Annotate.v2",".", cloneType,".RDS" , sep="") )
write.table(countObj.Annot.gb.Samples, paste("output/MiXCR/Results/countObj.Annot.gb.Samples.v2.",".", cloneType,".txt" , sep=""), sep="\t", quote = F, row.names = F)
write.table(countObj.Annot.gb.Samples.Annotate, paste("output/MiXCR/Results/countObj.Annot.gb.Samples",".Annotate.v2.", cloneType,".txt" , sep=""), sep="\t", quote = F, row.names = F)

# Binding ImmuneScore with TCR ####
countObj.Annot.gb.Samples.Annotate.NoNS <- countObj.Annot.gb.Samples.Annotate %>% filter( ! LIBRARY_TYPE %in% c("Normal", "CellLine")) %>% dplyr::mutate(TotalCloneSum = log10(TotalCloneSum+1)) 
length(unique(countObj.Annot.gb.Samples.Annotate.NoNS$Diagnosis))

selectCol="TotalCloneSum" ; StatsFinalCol="Diagnosis" ; SampleNames <- "Sample.ID.Alias"
tcrcloneCountPre          <- countObj.Annot.gb.Samples.Annotate.NoNS %>% 
  dplyr::select_(.dots=c(paste0("selectCol"), paste0("StatsFinalCol"), paste0("SampleNames")))

tcrcloneCountPre.Diag   <- tcrcloneCountPre %>%  dplyr::rename_(.dots = setNames(list(SampleNames,StatsFinalCol),c("Samples","Diagnosis"))) 
ScoresPre               <- tcrcloneCountPre.Diag[,!(colnames(tcrcloneCountPre.Diag) %in% c("Samples")), drop=FALSE]
ScoresPre <- ScoresPre %>% dplyr::filter(!Diagnosis %in% c("Teratoma", "YST"))
orderOfFactor           <- as.character( unique(ScoresPre$Diagnosis) )
orderOfSignature        <- colnames(ScoresPre)[-ncol(ScoresPre)]
colList                 <- c(1:(ncol(ScoresPre)-1)) ; Scores <- ScoresPre

## Plot and Save ####
#HC: Figure 3A
customColors=as.data.frame(customColorsVector)
colnames(customColors) <- c("Color")
customColors$Diagnosis <- rownames(customColors)

pdf("output/Figures/Figure3A.pdf", width = 15, height = 10)
corUtilsFuncs$OneVariablePlotSort(colList, Scores=Scores, orderOfFactor, orderOfSignature, standardize =FALSE, yLab=cloneType, customColorDF=customColors, plotType = "StringBean", sizeOfDots = 1.8)
dev.off()

###### plot for percent vs cloneCopy ####

## Remove samples with no cdr3aa ####

countObj.Annot.NoNA <- countObj.Annot.NoCL.totalReads %>% dplyr::filter( count != 0); dim(countObj.Annot.NoNA)
countObj.Annot.PercentTCR <- countObj.Annot.NoNA  %>% dplyr::group_by(Sample.Biowulf.ID.GeneExp) %>% 
  dplyr::mutate( percentinSample = freq) %>% 
  #HC: no cdr3aa and id.updated
  dplyr::select(count, Diagnosis, percentinSample,ReadsPerMillion)
#dplyr::select(count, cdr3aa, Diagnosis, percentinSample,ReadsPerMillion,Patient.ID.updated)
dim(countObj.Annot.PercentTCR) ; #View(countObj.Annot.PercentTCR)

## Plot and Save ####

base_breaks <- function(n = 10){
  function(x) {
    axisTicks(log2(range(x, na.rm = TRUE)), log = TRUE, n = n)
  }
}
customColorsVector <- setNames( unique(as.character(countObj.Annot.NoNA$Color.Jun)), unique(as.character(countObj.Annot.NoNA$Diagnosis)))
#HC: jun color
toPlotDF$Diagnosis <- factor(toPlotDF$Diagnosis,levels = val, ordered = TRUE)
toPlotDF <- toPlotDF %>% dplyr::arrange(Diagnosis)
customColorsVector <- setNames( unique(as.character(toPlotDF$Color.Jun)), unique(as.character(toPlotDF$Diagnosis)) )

#HC: Figure3D
pdf("output/Figures/Figure3D.pdf")
pctPlot <- ggplot( data = countObj.Annot.PercentTCR, aes( ReadsPerMillion, percentinSample) ) + 
  geom_point(aes(colour = factor(Diagnosis)), size = 2.5) + 
  coord_trans(x="log10") +
  scale_colour_manual("Diagnosis", values=customColorsVector  )+
  scale_y_continuous( trans = log_trans(10), 
                      name = paste0("frequency of a TCRB clone"),
                      breaks = c(0.01,0.1,0.25,0.50,0.75,1),
                      labels = scales::percent
  ) +
  scale_x_continuous( #trans = log_trans(10), 
    name =  paste0(" Expression of each TCRB clone"),
    breaks = c(0.1,0.5,1,2,3,4)
  ) +
  theme_bw() +
  theme( panel.grid.major = element_line(colour = "grey50", size = 0.25), 
         panel.grid.minor = element_blank(),
         legend.position = "right")  #element_line(colour = "grey50", size = 0.25) ) + 
pctPlot
dev.off()

