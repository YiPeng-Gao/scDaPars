# scDaPars: Dynamic Analysis of Alternative Polyadenylation from Single-Cell RNA-Seq
Yipeng Gao, Wei Li 2020-08-22

## Latest News
>2020-08-22: 
- Version 0.0.1 is released!

## Introduction

scDaPars is a bioinformatics algorithm to accurately quantify Alternative Polyadenylation (APA) events at both single-cell and single-gene resolution using standard scRNA-seq data. 

Step.1 scDaPars first takes scRNA-seq genome coverage data (bedgraph format) as input and forms a linear regression model to jointly infer the exact location of proximal poly(A) sites (Current Version of scDaPars do not support this function, raw PDUI values are calulated seperately using [DaPars2](https://github.com/3UTR/DaPars2/)).

Step.2 scDaPars constructs a nearest neighbor graph based on the sparse APA matrix generated in step.1 to identify a pool of candidate neighboring cells that have similar APA profiles.

Step.3 scDaPars uses a non-negative least square (NNLS) regression model to refine neighboring cells and impute PDUIs of dropout genes in each cell.

We welcome any suggestions on the package. For technical problems, please report to [Issues](https://github.com/YiPeng-Gao/scDaPars/issues). For suggestions and comments, please contact Yipeng (yipeng.gao@bcm.edu) or Dr. Wei Li (wei.li@uci.edu).
## Installation

The package is not on CRAN yet. For installation please use the following codes in ```R```

```
install.packages("devtools")
library(devtools)

install_github("YiPeng-Gao/scDaPars")
```
## Quick Start
The imputation steps of ```scDaPars``` takes the APA matrix from step.1 as input and in the simplest case, the imputation can be done with one single function ```scDaPars```:

```
scDaPars(raw_PDUI_file,         # full path of the raw PDUI matrix generated by step1 of scDaPars
         out_dir,               # full path of the output directory
         filter_gene_thre,      # the percent of cells a gene's APA must be detected in step.1
         filter_cell_thre)      # the percent of gene APA a cell must be detected in step.1
```
## Example
The dataset used in this example is a time-course scRNA-seq dataset containing 758 cells sequenced at 0, 12, 24, 36, 72 and 96 h of differentiation during human definitive endoderm (DE) emergence from [Chu et al.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1033-x) under GEO accession code GSE75748. After quality control, there is 739 cells remained for analysis.

>1. Genearting raw PDUI matrix
For generating raw PDUI matrix in step 1, we assume that the scRNA-seq data has been preprocessed, so that we have one wiggle file per cell. The raw PDUI files are then generated by [DaPars2](https://github.com/3UTR/DaPars2/). The raw PDUI matrix for this example "Dapars_hESC_combined_all_chromosome.txt" is included in the example folder. 

>2. install and load scDaPars R package
```
if(!require(devtools)) install.packages("devtools")
library(devtools)
```
```
## Loading required package: devtools

## Loading required package: usethis
```
```
devtools::install_github("YiPeng-Gao/scDaPars")
library(scDaPars)
```
```## Loading required package: penalized

## Loading required package: survival

## Welcome to penalized. For extended examples, see vignette("penalized").

## Loading required package: foreach

## Loading required package: RANN

## Loading required package: igraph

## 
## Attaching package: 'igraph'

## The following objects are masked from 'package:stats':
## 
##     decompose, spectrum

## The following object is masked from 'package:base':
## 
##     union
```
>3. Run scDaPars
```
scDaPars.res = scDaPars(raw_PDUI_file = "Dapars_hESC_combined_all_chromosome.txt", 
                        out_dir = "./scDaPars_result", 
                        filter_gene_thre = 0.2, 
                        filter_cell_thre = 0.1)
```
```
## [1] "Reading in raw PDUI matrix ..."
## [1] "number of genes in raw count matrix 13029"
## [1] "number of cells in raw count matrix 743"
## [1] "Reading finished!"
## [1] "Start Processing raw PDUI matrix"
## [1] "Pre-Processing finished"
## [1] "Start Imputation Steps ..."
## [1] "Find potential Neighboring Cells ..."
## [1] "Number of neighbors(clusters) is 4"
## [1] "Outliers is/are "
## [1] "Imputation with neighboring cells ..."
## [1] 100
## [1] 100
## [1] 100
## [1] 200
## [1] 300
## [1] "Imputation Steps Finished!"
## [1] "Writing imputed PDUI matrix ..."
```
```
head(scDaPars.res)[,1:6]
```
```
##                      SRR2978558  SRR2978559  SRR2978560  SRR2978561  SRR2978562
## ENST00000000412.3    1.00000000 0.804251535 0.857067085 1.000000000 0.706309880
## ENST00000002165.10_1 0.01711468 0.002728381 0.001751477 0.001465125 0.001576918
## ENST00000003100.13_3 0.05000000 0.607145404 0.496367748 0.733390448 0.611959391
## ENST00000005082.13_1 0.75459201 1.000000000 0.788302984 0.110000000 0.672660621
## ENST00000005257.7_2  0.15491832 0.000000000 0.289365021 0.285070704 0.224915218
## ENST00000005259.8_2  0.05000000 0.046576290 0.036624434 0.050765637 0.070246239
##                       SRR2978564
## ENST00000000412.3    0.710688649
## ENST00000002165.10_1 0.002705341
## ENST00000003100.13_3 0.523058864
## ENST00000005082.13_1 0.655225597
## ENST00000005257.7_2  1.000000000
## ENST00000005259.8_2  0.071559391
```
>4. Visualize scDaPars' results
Load cell type information for the example data
```
hESC_SRA = read.table("/Users/yipenggao/Documents/scDaPars/example/SraRunTable.txt", header = T, sep = ",", stringsAsFactors = F)
cell_type = hESC_SRA[which(hESC_SRA$Run %in% colnames(scDaPars.res)),]
cell_type = cell_type[match(colnames(scDaPars.res), cell_type$Run),]
head(cell_type)
```
```
##            Run Assay.Type AvgSpotLen  BioProject    BioSample Center.Name
## 933 SRR2978558    RNA-Seq         51 PRJNA305280 SAMN04322527         GEO
## 934 SRR2978559    RNA-Seq         51 PRJNA305280 SAMN04322528         GEO
## 935 SRR2978560    RNA-Seq         51 PRJNA305280 SAMN04322529         GEO
## 936 SRR2978561    RNA-Seq         51 PRJNA305280 SAMN04322530         GEO
## 937 SRR2978562    RNA-Seq         51 PRJNA305280 SAMN04322531         GEO
## 939 SRR2978564    RNA-Seq         51 PRJNA305280 SAMN04322533         GEO
##     Consent DATASTORE.filetype DATASTORE.provider
## 933  public                sra         gs,ncbi,s3
## 934  public                sra         gs,ncbi,s3
## 935  public                sra         gs,ncbi,s3
## 936  public                sra         gs,ncbi,s3
## 937  public                sra         gs,ncbi,s3
## 939  public                sra         gs,ncbi,s3
##                   DATASTORE.region Experiment facs_sorting GEO_Accession
## 933 gs.US,ncbi.public,s3.us-east-1 SRX1468145   not sorted    GSM1965870
## 934 gs.US,ncbi.public,s3.us-east-1 SRX1468146   not sorted    GSM1965871
## 935 gs.US,ncbi.public,s3.us-east-1 SRX1468147   not sorted    GSM1965872
## 936 gs.US,ncbi.public,s3.us-east-1 SRX1468148   not sorted    GSM1965873
## 937 gs.US,ncbi.public,s3.us-east-1 SRX1468149   not sorted    GSM1965874
## 939 gs.US,ncbi.public,s3.us-east-1 SRX1468151   not sorted    GSM1965876
##              Instrument LibraryLayout LibrarySelection  LibrarySource MBases
## 933 Illumina HiSeq 2500        SINGLE             cDNA TRANSCRIPTOMIC    105
## 934 Illumina HiSeq 2500        SINGLE             cDNA TRANSCRIPTOMIC    116
## 935 Illumina HiSeq 2500        SINGLE             cDNA TRANSCRIPTOMIC    141
## 936 Illumina HiSeq 2500        SINGLE             cDNA TRANSCRIPTOMIC     99
## 937 Illumina HiSeq 2500        SINGLE             cDNA TRANSCRIPTOMIC     90
## 939 Illumina HiSeq 2500        SINGLE             cDNA TRANSCRIPTOMIC    103
##     MBytes     Organism       passage Platform          ReleaseDate sample_acc
## 933     65 Homo sapiens passage 30-35 ILLUMINA 2016-07-22T00:00:00Z SRS1194214
## 934     72 Homo sapiens passage 30-35 ILLUMINA 2016-07-22T00:00:00Z SRS1194213
## 935     85 Homo sapiens passage 30-35 ILLUMINA 2016-07-22T00:00:00Z SRS1194212
## 936     61 Homo sapiens passage 30-35 ILLUMINA 2016-07-22T00:00:00Z SRS1194211
## 937     54 Homo sapiens passage 30-35 ILLUMINA 2016-07-22T00:00:00Z SRS1194210
## 939     62 Homo sapiens passage 30-35 ILLUMINA 2016-07-22T00:00:00Z SRS1194208
##     Sample.Name                          source_name SRA.Study
## 933  GSM1965870 H9 cells differentiated for 12 hours SRP067036
## 934  GSM1965871 H9 cells differentiated for 12 hours SRP067036
## 935  GSM1965872 H9 cells differentiated for 12 hours SRP067036
## 936  GSM1965873 H9 cells differentiated for 12 hours SRP067036
## 937  GSM1965874 H9 cells differentiated for 12 hours SRP067036
## 939  GSM1965876 H9 cells differentiated for 12 hours SRP067036
```
Perform UMAP analysis
```
scDaPars.res.umap = umap(t(scDaPars.res))
scDaPars.res.umap.data = data.frame(scDaPars.res.umap$layout)
colnames(scDaPars.res.umap.data) = c("Dim1", "Dim2")
scDaPars.res.umap.data$cellType = cell_type$source_name
```
Generate Scatter plot
```
ggplot(scDaPars.res.umap.data, aes(x=Dim1, y=Dim2, color = cellType)) +
  scale_color_manual(values = c("#7DD2D9", "#FFA500", "#e55b54", "#ad6c58", "#989797", "#166FD5")) +
  geom_point(size = 0.5) +
  theme_bw(base_size = 14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right", legend.title = element_blank())
```  
![](https://github.com/YiPeng-Gao/scDaPars/blob/master/Example/hESC_scDaPars_clustering.pdf)
