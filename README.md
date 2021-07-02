# RNAseqLandscape
R code for paper "Immuno-Transcriptomic Profiling of Extracranial Pediatric Solid Malignancies"

## Installation

The easiest way to get this pipeline is to clone the repository.

SSH:
```
git clone git@github.com:CCRGeneticsBranch/RNAseqLandscape.git
```
HTTPS:
```
https://github.com/CCRGeneticsBranch/RNAseqLandscape.git
```
We need an extra data file (too big for GitHub). Please copy this to input/GeneRDS/:
```
cp /data/khanlab3/hsienchao/Landscape/RPKM_Data_Filt_Consolidated.GeneNames.all.TCGA.Khanlab.pc.log22019-03-19.rds RNAseqLandscape/input/GeneRDS/
```

## Requirements

We need the following R packages. Please make sure they are installed properly:

```
plyr
dplyr
data.table
ggplot2
tidyr
gridExtra
ggrepel
RColorBrewer
treemap
lazyeval
grid
ggridges
scales
randomcoloR
gplots
pheatmap
ape
amap
limma
edgeR
psych
tibble
Hmisc
```

## Run the script

Please change directory to the cloned folder:
```
cd RNAseqLandscape
```

Then run:
```
Rscripts src/expressionAnalysis.R
```

## Output

```
output/Figures: all the figures
output/FigureData: the data used to generate figures
output/DiffExpResults: the output of differential expression analyses
output/MiXCR: the output data for TCR analysis

