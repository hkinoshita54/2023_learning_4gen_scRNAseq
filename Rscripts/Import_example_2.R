####################
#LOAD DATA TO SEURAT
####################
Regev_CRC<-Read10X_h5(filename = "H5_GSE178341/GSE178341_crc10x_full_c295v4_submit.h5", use.names = TRUE, unique.features = TRUE)

Regev_CRC_Seurat<-CreateSeuratObject(Regev_CRC, min.cells = 3)

head(Regev_CRC_Seurat@meta.data,1)
################
#ADDING METADATA
################
rownames<-rownames(Regev_CRC_Seurat@meta.data)
metadata<-read.csv(file="H5_GSE178341/GSE178341_metadata_AD.csv")
metadata1<-data.frame(metadata,row.names = rownames)
Regev_CRC_Seurat<-AddMetaData(Regev_CRC_Seurat, metadata1)
tail(Regev_CRC_Seurat@meta.data,1)
Idents(Regev_CRC_Seurat)<-"clTopLevel"
saveRDS(Regev_CRC_Seurat, file="RDSobjects/RDS_Regev_CRC_Seurat")
##################
#REMOVE.DEAD.CELLS-------------> It is CLEAN
##################
Regev_CRC_Seurat[["percent.mt"]] <- PercentageFeatureSet(Regev_CRC_Seurat, pattern = "^mt-")
Idents(Regev_CRC_Seurat) <- "orig.ident"
VlnPlot(Regev_CRC_Seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
plot1 <- FeatureScatter(Regev_CRC_Seurat, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.05)
plot2 <- FeatureScatter(Regev_CRC_Seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.05)
head(Regev_CRC_Seurat@meta.data,1)
