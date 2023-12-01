EpiCellType <- read.delim("Celltype_Epi.txt", header = TRUE, col.names = 1)
EMBO_CRCcombined_SCT <- AddMetaData(EMBO_CRCcombined_SCT,EpiCellType, col.name = "Celltype")



#####Checking dead cells###########
EMBO_CRCcombined2[["percent.mt"]] <- PercentageFeatureSet(EMBO_CRCcombined2, pattern = "^MT-")
Idents(EMBO_CRCcombined2) <- "orig.ident"
VlnPlot(EMBO_CRCcombined2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

##############
#DIM.RED.AND.DIMPLOT
####################
############
#SCTransform
############
EMBO_CRCcombined_SCT <- SCTransform(EMBO_CRCcombined2, verbose = TRUE)
EMBO_CRCcombined_SCT <- RunPCA(EMBO_CRCcombined_SCT, verbose = T)
ElbowPlot(EMBO_CRCcombined_SCT)
EMBO_CRCcombined_SCT <- RunUMAP(EMBO_CRCcombined_SCT, dims = 1:15, verbose = T)
EMBO_CRCcombined_SCT <- FindNeighbors(EMBO_CRCcombined_SCT, dims = 1:15, verbose = T)
EMBO_CRCcombined_SCT <- FindClusters(EMBO_CRCcombined_SCT, verbose = T, resolution = 0.5)

head(EMBO_CRCcombined_SCT@meta.data,1)
Idents(EMBO_CRCcombined_SCT) <- "Compartment"
DimPlot(EMBO_CRCcombined_SCT, label = T)
DimPlot(EMBO_CRCcombined_SCT, label = T)
DimPlot(EMBO_CRCcombined_SCT, label = T)

saveRDS(EMBO_CRCcombined_SCT,file="RDSobjects/EMBO_CRCcombined_SCT")
Compartment <- read.delim("Compartment.txt", header = TRUE, col.names = 1)
EMBO_CRCcombined_SCT <- AddMetaData(EMBO_CRCcombined_SCT,Compartment, col.name = "Compartment")

#####Create ident serrated and conventional and Mutation#######
Origin <- read.delim("Origin_P.txt", header = TRUE, col.names = 1)
EMBO_CRCcombined_SCT <- AddMetaData(EMBO_CRCcombined_SCT,Origin, col.name = "Origin")

Sample <- read.delim("SampleID_P.txt", header = TRUE, col.names = 1)
EMBO_CRCcombined_SCT <- AddMetaData(EMBO_CRCcombined_SCT,Sample, col.name = "SampleID")

head(EMBO_CRCcombined_SCT@meta.data,1)