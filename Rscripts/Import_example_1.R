####### Data from liver met CRC #####
library(dplyr)
library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
library(sctransform)
library(hdf5r)
library(harmony)

#https://github.com/satijalab/seurat/wiki/Assay
#utils::methods(class = 'Assay')

# Or for the latest development version
devtools::install_github("const-ae/ggsignif")
library(ggplot2)
library(ggsignif)
library(readxl)
########
#GSE137720: each sample comes as group of 3 files: barcodes.gsv.gz, genes.csv.gz, matrix.mx.gz
#load each sample using read10x. NOTE: the files should not have sufix, so you need to place the 3 files in a subfolder with the name of the sample.
data_dirAMO <- 'RAW'
list.files(data_dirAMO) #This Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
LiverMET_AMO <- Read10X(data.dir = data_dirAMO)


#CreateSeuratObject and combine
#Create a Seurat object from each sample
LiverMET_AMO_obj <- CreateSeuratObject(counts = LiverMET_AMO, min.cells = 3)

head(LiverMET_AMO_obj@meta.data)
##################
#Remove dead cells
##################
LiverMET_AMO_obj[["percent.mt"]] <- PercentageFeatureSet(LiverMET_AMO_obj, pattern = "^MT-")
Idents(LiverMET_AMO_obj) <- "orig.ident"
VlnPlot(LiverMET_AMO_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
############
#SUBSET.DATA.DEAD.CELL.SINGLET.REMOVAL
############
##Not necessary for this object after cheching violin plot ####
#################################
#INSPECT.META.DATA.AND.SET.IDENTS
#################################
##NOTE: In this case genotype denotes sample numer from 1 to 5
head(LiverMET_AMO_obj@meta.data,100)
tail(LiverMET_AMO_obj@meta.data,100)
# Add barcode metadata: Extract genotype from barcode by strsplit on "-" and extracting second element of each resulting list 
sample <- sapply(strsplit(rownames(LiverMET_AMO_obj@meta.data), split="_"), "[[", 2)
# Add barcode metadata: Supply the genotype as additional metadata
LiverMET_AMO_obj <- AddMetaData(LiverMET_AMO_obj, metadata=data.frame(sample=sample, row.names=rownames(LiverMET_AMO_obj@meta.data)))
Idents(LiverMET_AMO_obj) <- "sample"

######origin###
# Add barcode metadata: Extract genotype from barcode by strsplit on "-" and extracting second element of each resulting list 
origin <- sapply(strsplit(rownames(LiverMET_AMO_obj@meta.data), split="_"), "[[", 3)
# Add barcode metadata: Supply the genotype as additional metadata
LiverMET_AMO_obj <- AddMetaData(LiverMET_AMO_obj, metadata=data.frame(origin=origin, row.names=rownames(LiverMET_AMO_obj@meta.data)))
Idents(LiverMET_AMO_obj) <- "origin"

############
#SCTransform
############
######Harmony############## skip PBMC
library(harmony)

Idents(LiverMET_AMO_obj) <- "origin"

LiverMET_AMO_CRC <- subset(LiverMET_AMO_obj, idents = c("CRC","LM"))

LiverMET_AMO_SCT1 <- SCTransform(LiverMET_AMO_CRC, vars.to.regress = "percent.mt", verbose = FALSE)
####################
#DIM.RED.AND.DIMPLOT
####################
LiverMET_AMO_SCT1 <- RunPCA(LiverMET_AMO_SCT1, verbose = T)

LiverMET_AMO_CRC_H<-RunHarmony(LiverMET_AMO_SCT1,group.by.vars = "sample", assay.use = "SCT")

LiverMET_AMO_CRC_H <- RunUMAP(LiverMET_AMO_CRC_H, reduction = "harmony",dims = 1:10, verbose = T)
LiverMET_AMO_CRC_H <- FindNeighbors(LiverMET_AMO_CRC_H,reduction = "harmony",dims = 1:10, verbose = T)
LiverMET_AMO_CRC_H <- FindClusters(LiverMET_AMO_CRC_H, verbose = T, resolution = 0.2)


Idents(LiverMET_AMO_CRC_H)<-"seurat_clusters"