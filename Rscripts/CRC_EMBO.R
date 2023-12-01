library(ggplot2)
library(Seurat)
library(ggsignif)

######
CRC_p007nframe.data <- read.delim("GSM5075659_p007n.tsv", header = TRUE, row.names = 1)
CRC_p007nframe.data<-CreateSeuratObject(CRC_p007nframe.data, project = "CRC", min.cells = 3)

CRC_p0026tframe.data <- read.delim("GSM5075683_p026t.tsv", header = TRUE, row.names = 1)
CRC_p0026tframe.data<-CreateSeuratObject(CRC_p0026tframe.data, project = "CRC", min.cells = 3)

CRC_p0025tframe.data <- read.delim("GSM5075682_p025t.tsv", header = TRUE, row.names = 1)
CRC_p0025tframe.data<-CreateSeuratObject(CRC_p0025tframe.data, project = "CRC", min.cells = 3)

CRC_p0025nframe.data <- read.delim("GSM5075681_p025n.tsv", header = TRUE, row.names = 1)
CRC_p0025nframe.data<-CreateSeuratObject(CRC_p0025nframe.data, project = "CRC", min.cells = 3)

CRC_p0021tframe.data <- read.delim("GSM5075680_p021t.tsv", header = TRUE, row.names = 1)
CRC_p0021tframe.data<-CreateSeuratObject(CRC_p0021tframe.data, project = "CRC", min.cells = 3)

CRC_p0021nframe.data <- read.delim("GSM5075679_p021n.tsv", header = TRUE, row.names = 1)
CRC_p0021nframe.data<-CreateSeuratObject(CRC_p0021nframe.data, project = "CRC", min.cells = 3)

CRC_p0020nframe.data <- read.delim("GSM5075677_p020n.tsv", header = TRUE, row.names = 1)
CRC_p0020nframe.data<-CreateSeuratObject(CRC_p0020nframe.data, project = "CRC", min.cells = 3)

CRC_p0020tframe.data <- read.delim("GSM5075678_p020t.tsv", header = TRUE, row.names = 1)
CRC_p0020tframe.data<-CreateSeuratObject(CRC_p0020tframe.data, project = "CRC", min.cells = 3)

CRC_p0017tframe.data <- read.delim("GSM5075676_p017t.tsv", header = TRUE, row.names = 1)
CRC_p0017tframe.data<-CreateSeuratObject(CRC_p0017tframe.data, project = "CRC", min.cells = 3)

CRC_p0017nframe.data <- read.delim("GSM5075675_p017n.tsv", header = TRUE, row.names = 1)
CRC_p0017nframe.data<-CreateSeuratObject(CRC_p0017nframe.data, project = "CRC", min.cells = 3)

CRC_p0016tframe.data <- read.delim("GSM5075674_p016t.tsv", header = TRUE, row.names = 1)
CRC_p0016tframe.data<-CreateSeuratObject(CRC_p0016tframe.data, project = "CRC", min.cells = 3)

CRC_p0016nframe.data <- read.delim("GSM5075673_p016n.tsv", header = TRUE, row.names = 1)
CRC_p0016nframe.data<-CreateSeuratObject(CRC_p0016nframe.data, project = "CRC", min.cells = 3)

CRC_p0014tframe.data <- read.delim("GSM5075672_p014t.tsv", header = TRUE, row.names = 1)
CRC_p0014tframe.data<-CreateSeuratObject(CRC_p0014tframe.data, project = "CRC", min.cells = 3)

CRC_p0014nframe.data <- read.delim("GSM5075671_p014n.tsv", header = TRUE, row.names = 1)
CRC_p0014nframe.data<-CreateSeuratObject(CRC_p0014nframe.data, project = "CRC", min.cells = 3)

CRC_p0013tframe.data <- read.delim("GSM5075670_p013t.tsv", header = TRUE, row.names = 1)
CRC_p0013tframe.data<-CreateSeuratObject(CRC_p0013tframe.data, project = "CRC", min.cells = 3)

CRC_p0013nframe.data <- read.delim("GSM5075669_p013n.tsv", header = TRUE, row.names = 1)
CRC_p0013nframe.data<-CreateSeuratObject(CRC_p0013nframe.data, project = "CRC", min.cells = 3)

CRC_p0012tframe.data <- read.delim("GSM5075668_p012t.tsv", header = TRUE, row.names = 1)
CRC_p0012tframe.data<-CreateSeuratObject(CRC_p0012tframe.data, project = "CRC", min.cells = 3)

CRC_p0012nframe.data <- read.delim("GSM5075667_p012n.tsv", header = TRUE, row.names = 1)
CRC_p0012nframe.data<-CreateSeuratObject(CRC_p0012nframe.data, project = "CRC", min.cells = 3)

CRC_p009t1frame.data <- read.delim("GSM5075665_p009t1.tsv", header = TRUE, row.names = 1)
CRC_p009t1frame.data<-CreateSeuratObject(CRC_p009t1frame.data, project = "CRC", min.cells = 3)

CRC_p009t2frame.data <- read.delim("GSM5075666_p009t2.tsv", header = TRUE, row.names = 1)
CRC_p009t2frame.data<-CreateSeuratObject(CRC_p009t2frame.data, project = "CRC", min.cells = 3)

CRC_p009n1frame.data <- read.delim("GSM5075663_p009n1.tsv", header = TRUE, row.names = 1)
CRC_p009n1frame.data<-CreateSeuratObject(CRC_p009n1frame.data, project = "CRC", min.cells = 3)

CRC_p009n2frame.data <- read.delim("GSM5075664_p009n2.tsv", header = TRUE, row.names = 1)
CRC_p009n2frame.data<-CreateSeuratObject(CRC_p009n2frame.data, project = "CRC", min.cells = 3)

CRC_p008tframe.data <- read.delim("GSM5075662_p008t.tsv", header = TRUE, row.names = 1)
CRC_p008tframe.data<-CreateSeuratObject(CRC_p008tframe.data, project = "CRC", min.cells = 3)

CRC_p008nframe.data <- read.delim("GSM5075661_p008n.tsv", header = TRUE, row.names = 1)
CRC_p008nframe.data<-CreateSeuratObject(CRC_p008nframe.data, project = "CRC", min.cells = 3)

CRC_p007tframe.data <- read.delim("GSM5075660_p007t.tsv", header = TRUE, row.names = 1)
CRC_p007tframe.data<-CreateSeuratObject(CRC_p007tframe.data, project = "CRC", min.cells = 3)



EMBO_CRCcombined <- merge(x= CRC_p007tframe.data, y = c(CRC_p008nframe.data, CRC_p008tframe.data, CRC_p009n2frame.data, CRC_p009n1frame.data,CRC_p009t2frame.data,CRC_p009t1frame.data,CRC_p0012nframe.data,CRC_p0012tframe.data,CRC_p0013nframe.data,CRC_p0013tframe.data,CRC_p0014nframe.data,CRC_p0014tframe.data,CRC_p0016tframe.data,CRC_p0016nframe.data, CRC_p0017tframe.data, CRC_p0017nframe.data,CRC_p0020tframe.data,CRC_p0020nframe.data,CRC_p0021tframe.data,CRC_p0021nframe.data,CRC_p0025tframe.data,CRC_p0025nframe.data,CRC_p0026tframe.data,CRC_p007nframe.data),project = "EMBO")



Celltype <- read.delim("Celltype.txt", header = TRUE, col.names = 1)

SampleID <- read.delim("sample_ID.txt", header = TRUE, col.names = 1)
EMBO_CRCcombined2 <- AddMetaData(EMBO_CRCcombined,Celltype, col.name = "X1")

EpiCellType <- read.delim("Celltype_Epi.txt", header = TRUE, col.names = 1)
EMBO_CRCcombined_SCT <- AddMetaData(EMBO_CRCcombined_SCT,EpiCellType, col.name = "Celltype")

head(EMBO_CRCcombined_SCT@meta.data,2)

saveRDS(EMBO_CRCcombined, file="RDSobjects/EMBO_CRCcombined")

saveRDS(EMBO_CRCcombined2, file="RDSobjects/EMBO_CRCcombined2")


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


########Set ident mutation### and Serrated
Idents(EMBO_CRCcombined_SCT) <- "SampleID"

KRASmut<-WhichCells(EMBO_CRCcombined_SCT, ident= c("p016t",	"p017t",	"p021t",	"p025t"))
EMBO_CRCcombined_SCT@meta.data$Mutation [KRASmut]
EMBO_CRCcombined_SCT<-SetIdent(EMBO_CRCcombined_SCT, cells = KRASmut, value = "KRASmut")
EMBO_CRCcombined_SCT$Mutation<-Idents (EMBO_CRCcombined_SCT)

RAS_RAF_WT<-WhichCells(EMBO_CRCcombined_SCT, ident= c("p008t","p009t1","p009t2","p012t","p013t"))
EMBO_CRCcombined_SCT@meta.data$Mutation [RAS_RAF_WT]
EMBO_CRCcombined_SCT<-SetIdent(EMBO_CRCcombined_SCT, cells = RAS_RAF_WT, value = "RAS_RAF_WT")
EMBO_CRCcombined_SCT$Mutation<-Idents (EMBO_CRCcombined_SCT)

BRAFmut<-WhichCells(EMBO_CRCcombined_SCT, ident= c("p007t","p014t","p020t","p026t"))
EMBO_CRCcombined_SCT@meta.data$Mutation [BRAFmut]
EMBO_CRCcombined_SCT<-SetIdent(EMBO_CRCcombined_SCT, cells = BRAFmut, value = "BRAFmut")
EMBO_CRCcombined_SCT$Mutation<-Idents (EMBO_CRCcombined_SCT)

Normal<-WhichCells(EMBO_CRCcombined_SCT, ident= c("p007n",
"p008n","p009n1","p009n2","p012n",
"p013n",
"p014n",
"p016n",
"p017n",
"p020n",
"p021n",
"p025n"))
EMBO_CRCcombined_SCT@meta.data$Mutation [Normal]
EMBO_CRCcombined_SCT<-SetIdent(EMBO_CRCcombined_SCT, cells = Normal, value = "Normal")
EMBO_CRCcombined_SCT$Mutation<-Idents (EMBO_CRCcombined_SCT)


Idents(EMBO_CRCcombined_SCT) <- "Compartment"
DimPlot(EMBO_CRCcombined_SCT)

FeaturePlot(EMBO_CRCcombined_SCT, "TNFSF15", cols = c("lightgrey","red","darkred"), split.by = "Mutation")

DimPlot(EMBO_CRCcombined_SCT,split.by = "Mutation")

######
########Set ident mutation### and Serrated
Idents(EMBO_CRCcombined_SCT) <- "SampleID"

Conventional<-WhichCells(EMBO_CRCcombined_SCT, ident= c("p016t",	"p017t",	"p021t",	"p025t","p009t1","p009t2","p013t"))
EMBO_CRCcombined_SCT@meta.data$Progression [Conventional]
EMBO_CRCcombined_SCT<-SetIdent(EMBO_CRCcombined_SCT, cells = Conventional, value = "Conventional")
EMBO_CRCcombined_SCT$Progression<-Idents (EMBO_CRCcombined_SCT)

Inflammation_associated<-WhichCells(EMBO_CRCcombined_SCT, ident= c("p008t"))
EMBO_CRCcombined_SCT@meta.data$Progression [Inflammation_associated]
EMBO_CRCcombined_SCT<-SetIdent(EMBO_CRCcombined_SCT, cells = Inflammation_associated, value = "Inflammation_associated")
EMBO_CRCcombined_SCT$Progression<-Idents (EMBO_CRCcombined_SCT)

Serrated<-WhichCells(EMBO_CRCcombined_SCT, ident= c("p007t","p020t","p026t"))
EMBO_CRCcombined_SCT@meta.data$Progression [Serrated]
EMBO_CRCcombined_SCT<-SetIdent(EMBO_CRCcombined_SCT, cells = Serrated, value = "Serrated")
EMBO_CRCcombined_SCT$Progression<-Idents (EMBO_CRCcombined_SCT)

NotAvailable<-WhichCells(EMBO_CRCcombined_SCT, ident= c("p012t","p014t"))
EMBO_CRCcombined_SCT@meta.data$Progression [NotAvailable]
EMBO_CRCcombined_SCT<-SetIdent(EMBO_CRCcombined_SCT, cells = NotAvailable, value = "NotAvailable")
EMBO_CRCcombined_SCT$Progression<-Idents (EMBO_CRCcombined_SCT)

Normal<-WhichCells(EMBO_CRCcombined_SCT, ident= c("p007n",
                                                  "p008n","p009n1","p009n2","p012n",
                                                  "p013n",
                                                  "p014n",
                                                  "p016n",
                                                  "p017n",
                                                  "p020n",
                                                  "p021n",
                                                  "p025n"))
EMBO_CRCcombined_SCT@meta.data$Progression [Normal]
EMBO_CRCcombined_SCT<-SetIdent(EMBO_CRCcombined_SCT, cells = Normal, value = "Normal")
EMBO_CRCcombined_SCT$Progression<-Idents (EMBO_CRCcombined_SCT)

#####
Idents(EMBO_CRCcombined_SCT) <- "X1"
DimPlot(EMBO_CRCcombined_SCT)
DimPlot(EMBO_CRCcombined_SCT,split.by = "Mutation")
head(EMBO_CRCcombined_SCT@meta.data,1)
######SUbset epithelial######

Idents(EMBO_CRCcombined_SCT) <- "Compartment"
DimPlot(EMBO_CRCcombined_SCT)

#####Read FSC signatures ####
FSC <- read.delim("./Signatures/FSC_HUMAN.txt")
EMBO_CRCcombined_SCT <- AddModuleScore(EMBO_CRCcombined_SCT, features = FSC, name = "FSC_")
FeaturePlot(All_tumor,"FSC_1", pt.size = 0.1, min.cutoff =0, cols = c("lightgrey","red","darkred"), label = F)



Idents(EMBO_CRCcombined_SCT) <- "Progression"

All_tumor <- subset(EMBO_CRCcombined_SCT, idents = c("NotAvailable","Serrated","Conventional","Inflammation_associated"))

Idents(All_tumor) <- "Compartment"
Idents(EMBO_CRCcombined_SCT) <- "Celltype"

DimPlot(EMBO_CRCcombined_SCT)
FeaturePlot(EMBO_CRCcombined_SCT,features = "TNFSF15", cols = c("lightgrey","red","darkred"))


Epithelium_EMBO <- subset(EMBO_CRCcombined_SCT, idents = c("Epithelial"))


Epithelium_EMBO  <- SCTransform(Epithelium_EMBO , assay="RNA", verbose = T)
Epithelium_EMBO  <- RunPCA(Epithelium_EMBO, verbose = T)
ElbowPlot(Epithelium_EMBO)
Epithelium_EMBO  <- RunUMAP(Epithelium_EMBO, dims = 1:15, verbose = T)
Epithelium_EMBO  <- FindNeighbors(Epithelium_EMBO , dims = 1:15, verbose = T)
Epithelium_EMBO  <- FindClusters(Epithelium_EMBO, verbose = T, resolution = 0.2)
DimPlot(Epithelium_EMBO, label = T, label.size = 3, pt.size = 0, split.by = "Progression")

Idents(Epithelium_EMBO) <- "Celltype"

DimPlot(Epithelium_EMBO, label = F, label.size = 3, pt.size = 0, split.by = "Progression")

FeaturePlot(Epithelium_EMBO, "TNFSF15", split.by = "Mutation")



######Counts create celltype.mutations######
####
####Counts#######
#Cell.counts
head(Epithelium_EMBO@meta.data,2)
#Add group_new and cluster info to meta.data
Idents(Epithelium_EMBO) <- "Celltype"
Epithelium_EMBO$Celltype.Mutation <- paste(Idents(Epithelium_EMBO), Epithelium_EMBO$Mutation, sep = "_")
head(Epithelium_EMBO@meta.data,3)

Idents(Epithelium_EMBO) <- "Celltype.Mutation"
counts.Epithelium_EMBO_byMutation <-table(Idents(Epithelium_EMBO))
write.table(as.matrix(counts.Epithelium_EMBO_byMutation),"./Results/Counts/counts.Epithelium_EMBO_byMutation.txt",sep="\t",col.names=T,row.names=T)



####
####Counts#######
#Cell.counts
head(Epithelium_EMBO@meta.data,2)
#Add group_new and cluster info to meta.data
Idents(Epithelium_EMBO) <- "Celltype"
Epithelium_EMBO$Celltype.Progression <- paste(Idents(Epithelium_EMBO), Epithelium_EMBO$Progression, sep = "_")
head(Epithelium_EMBO@meta.data,3)

Idents(Epithelium_EMBO) <- "Celltype.Progression"
counts.Epithelium_EMBO_byProgression <-table(Idents(Epithelium_EMBO))
write.table(as.matrix(counts.Epithelium_EMBO_byProgression),"./Results/Counts/counts.Epithelium_EMBO_byProgression.txt",sep="\t",col.names=T,row.names=T)



#####Try to classify tumor patients by aPKCs expression in epithelium low and high #################
#####To do GSVA#########
Idents(Epithelium_EMBO)<-"cluster.group"
avg.epithelium_LKO_DKO_Harmony<-AverageExpression(SEURAT_Intestine_UMAP_Epithelium_HarmonyClean1) 
write.table(as.matrix(avg.epithelium_LKO_DKO_Harmony$RNA), "./Results/GeneExpression///avg.epithelium_LKO_DKO_Harmony_RNA.txt", sep="\t",col.names = T, row.names = T)
write.table(as.matrix(avg.epithelium_LKO_DKO_Harmony$SCT), "./Results/GeneExpression//avg.epithelium_LKO_DKO_Harmony_SCT.txt", sep="\t",col.names = T, row.names = T)

Idents(SEURAT_Intestine_UMAP_Epithelium_HarmonyClean1)<-"seurat_clusters"

DimPlot(Epithelium_EMBO)
head(Epithelium_EMBO@meta.data,20)

##############
Cholesterol <- read.delim("./Cholesterol_HUMAN.txt")
Epithelium_EMBO <- AddModuleScore(Epithelium_EMBO, features = Cholesterol, name = "Cholesterol_")
FeaturePlot(Epithelium_EMBO,"Cholesterol_1", pt.size = 0.1, min.cutoff = -3, cols = c("lightgrey","red","darkred"), label = F, split.by = "Mutation")

Idents(Epithelium_EMBO) <- "Progression"

VlnPlot(Epithelium_EMBO,"Cholesterol_1", pt.size = 0.0,idents = c("Serrated","Conventional"))+geom_signif(comparisons = list(c("Serrated","Conventional")),test = "t.test",map_signif_level = TRUE, textsize =3) + ylim(NA,1) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)


#####Subset tumoral cells alone #####

Idents(Epithelium_EMBO) <- "Progression"
EpitheliumSerratedConventional_EMBO <- subset(Epithelium_EMBO, idents = c("Serrated","Conventional"))


Idents(EpitheliumSerratedConventional_EMBO) <- "Celltype"
EpitheliumSerratedConventionalTC_EMBO <- subset(EpitheliumSerratedConventional_EMBO, idents = c("TC3","TC1","TC2","TC4","Imma. Gob.","Stem/TA"))

EpitheliumSerratedConventionalTC_EMBO  <- SCTransform(EpitheliumSerratedConventionalTC_EMBO , assay="RNA", verbose = T)
EpitheliumSerratedConventionalTC_EMBO  <- RunPCA(EpitheliumSerratedConventionalTC_EMBO, verbose = T)
ElbowPlot(EpitheliumSerratedConventionalTC_EMBO)
EpitheliumSerratedConventionalTC_EMBO  <- RunUMAP(EpitheliumSerratedConventionalTC_EMBO, dims = 1:8, verbose = T)
EpitheliumSerratedConventionalTC_EMBO  <- FindNeighbors(EpitheliumSerratedConventionalTC_EMBO , dims = 1:8, verbose = T)
EpitheliumSerratedConventionalTC_EMBO  <- FindClusters(EpitheliumSerratedConventionalTC_EMBO, verbose = T, resolution = 0.2)
DimPlot(EpitheliumSerratedConventionalTC_EMBO, label = T, label.size = 3, pt.size = 0, split.by = "Progression")

Idents(EpitheliumSerratedConventionalTC_EMBO) <- "Celltype"

DimPlot(EpitheliumSerratedConventionalTC_EMBO, label = F, label.size = 3, pt.size = 0, split.by = "Progression")


Cholesterol <- read.delim("./Cholesterol_HUMAN.txt")
EpitheliumSerratedConventionalTC_EMBO <- AddModuleScore(EpitheliumSerratedConventionalTC_EMBO, features = Cholesterol, name = "Cholesterol_")
FeaturePlot(EpitheliumSerratedConventionalTC_EMBO,"Cholesterol_1", pt.size = 0.1, min.cutoff = -3, cols = c("lightgrey","red","darkred"), label = F, split.by = "Progression")

Idents(EpitheliumSerratedConventionalTC_EMBO) <- "Progression"

VlnPlot(EpitheliumSerratedConventionalTC_EMBO,"Cholesterol_1", pt.size = 0.0,idents = c("Serrated","Conventional"))+geom_signif(comparisons = list(c("Serrated","Conventional")),test = "t.test",map_signif_level = TRUE, textsize =3) + ylim(NA,1) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)

head(EpitheliumSerratedConventionalTC_EMBO@meta.data)

FeaturePlot(EpitheliumSerratedConventionalTC_EMBO,"Cholesterol_1", pt.size = 0.1, min.cutoff = -3, cols = c("lightgrey","red","darkred"), label = F, split.by = "SampleID")

saveRDS(EpitheliumSerratedConventionalTC_EMBO,file="RDSobjects/EpitheliumSerratedConventionalTC_EMBO")

########Subset by Mutation

Idents(Epithelium_EMBO) <- "Mutation"
DimPlot(Epithelium_EMBO)

EpitheliumMut_EMBO <- subset(Epithelium_EMBO, idents = c("BRAFmut","RAS_RAF_WT","KRASmut"))


Idents(EpitheliumMut_EMBO) <- "Celltype"
EpitheliumMut_EMBOTC_EMBO <- subset(EpitheliumMut_EMBO, idents = c("TC3","TC1","TC2","TC4","Imma. Gob.","Stem/TA"))

EpitheliumMut_EMBOTC_EMBO  <- SCTransform(EpitheliumMut_EMBOTC_EMBO , assay="RNA", verbose = T)
EpitheliumMut_EMBOTC_EMBO  <- RunPCA(EpitheliumMut_EMBOTC_EMBO, verbose = T)
ElbowPlot(EpitheliumMut_EMBOTC_EMBO)
EpitheliumMut_EMBOTC_EMBO  <- RunUMAP(EpitheliumMut_EMBOTC_EMBO, dims = 1:8, verbose = T)
EpitheliumMut_EMBOTC_EMBO  <- FindNeighbors(EpitheliumMut_EMBOTC_EMBO , dims = 1:8, verbose = T)
EpitheliumMut_EMBOTC_EMBO  <- FindClusters(EpitheliumMut_EMBOTC_EMBO, verbose = T, resolution = 0.2)
DimPlot(EpitheliumMut_EMBOTC_EMBO, label = T, label.size = 3, pt.size = 0, split.by = "Progression")

Idents(EpitheliumMut_EMBOTC_EMBO) <- "Celltype"

DimPlot(EpitheliumMut_EMBOTC_EMBO, label = F, label.size = 3, pt.size = 0, split.by = "Progression")

FSC <- read.delim("./Signatures/FSC_HUMAN.txt")
EpitheliumMut_EMBOTC_EMBO <- AddModuleScore(EpitheliumMut_EMBOTC_EMBO, features = FSC, name = "FSC_")
FeaturePlot(EpitheliumMut_EMBOTC_EMBO,"FSC_1", pt.size = 0.1, min.cutoff = 0, cols = c("lightgrey","red","darkred"), label = F, split.by = "Progression")

RSC <- read.delim("./Signatures/RSC_HUMAN.txt")
EpitheliumMut_EMBOTC_EMBO <- AddModuleScore(EpitheliumMut_EMBOTC_EMBO, features = RSC, name = "RSC_")
FeaturePlot(EpitheliumMut_EMBOTC_EMBO,"RSC_1", pt.size = 0.1, min.cutoff = 0, cols = c("lightgrey","red","darkred"), label = F, split.by = "Progression")

TGC <- read.delim("./Signatures/TGC_HUMAN.txt")
EpitheliumMut_EMBOTC_EMBO <- AddModuleScore(EpitheliumMut_EMBOTC_EMBO, features = TGC, name = "TGC_")
FeaturePlot(EpitheliumMut_EMBOTC_EMBO,"TGC_1", pt.size = 0.1, min.cutoff = 0, cols = c("lightgrey","red","darkred"), label = F, split.by = "Progression")

TCC <- read.delim("./Signatures/TCC_HUMAN.txt")
EpitheliumMut_EMBOTC_EMBO <- AddModuleScore(EpitheliumMut_EMBOTC_EMBO, features = TCC, name = "TCC_")
FeaturePlot(EpitheliumMut_EMBOTC_EMBO,"TCC_1", pt.size = 0.1, min.cutoff = 0.5, cols = c("lightgrey","red","darkred"), label = F, split.by = "Progression")

FeaturePlot(EpitheliumMut_EMBOTC_EMBO,"SOX2", pt.size = 0.1, min.cutoff = 0.0, cols = c("lightgrey","red","darkred"), label = F, split.by = "Mutation")
DimPlot(EpitheliumMut_EMBOTC_EMBO, split.by = "Progression")

VlnPlot(EpitheliumMut_EMBOTC_EMBO, "SOX2")

Idents(EpitheliumMut_EMBOTC_EMBO) <- "Mutation"

VlnPlot(EpitheliumMut_EMBOTC_EMBO,"Cholesterol_1", pt.size = 0.0,idents = c("BRAFmut","RAS_RAF_WT","KRASmut"))+geom_signif(comparisons = list(c("BRAFmut","RAS_RAF_WT")),test = "t.test",map_signif_level = TRUE, textsize =3) + ylim(NA,1) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)

#####Try to classify tumor patients by aPKCs expression in epithelium low and high #################
#####To do GSVA#########
head(EpitheliumMut_EMBOTC_EMBO@meta.data,10)

Idents(EpitheliumMut_EMBOTC_EMBO)<-"SampleID"
avg.EpitheliumMut_EMBOTC_EMBO<-AverageExpression(EpitheliumMut_EMBOTC_EMBO) 
write.table(as.matrix(avg.EpitheliumMut_EMBOTC_EMBO$RNA), "./Results/avg.EpitheliumMut_EMBOTC_EMBO_RNA.txt", sep="\t",col.names = T, row.names = T)
write.table(as.matrix(avg.EpitheliumMut_EMBOTC_EMBO$SCT), "./Results//avg.EpitheliumMut_EMBOTC_EMBO_SCT.txt", sep="\t",col.names = T, row.names = T)

Idents(avg.EpitheliumMut_EMBOTC_EMBO)<-"seurat_clusters"


Idents(EpitheliumMut_EMBOTC_EMBO)<-"SampleID"

aPKClow<-WhichCells(EpitheliumMut_EMBOTC_EMBO, ident= c("p021t","p025t"))
EpitheliumMut_EMBOTC_EMBO@meta.data$aPKC [aPKClow]
EpitheliumMut_EMBOTC_EMBO<-SetIdent(EpitheliumMut_EMBOTC_EMBO, cells = aPKClow, value = "aPKClow")
EpitheliumMut_EMBOTC_EMBO$aPKC<-Idents (EpitheliumMut_EMBOTC_EMBO)
aPKChigh<-WhichCells (EpitheliumMut_EMBOTC_EMBO, ident= c("p008t","p013t","p017t"))
EpitheliumMut_EMBOTC_EMBO@meta.data$aPKC [aPKChigh]
EpitheliumMut_EMBOTC_EMBO<-SetIdent(EpitheliumMut_EMBOTC_EMBO, cells = aPKChigh, value = "aPKChigh")
EpitheliumMut_EMBOTC_EMBO$aPKC<-Idents (EpitheliumMut_EMBOTC_EMBO)

Idents(EpitheliumMut_EMBOTC_EMBO)<-"aPKC"

EpitheliumMut_EMBOTC_EMBO_aPKC <- subset(EpitheliumMut_EMBOTC_EMBO, idents = c("aPKClow","aPKChigh"))


Cholesterol <- read.delim("./Cholesterol_HUMAN.txt")
EpitheliumMut_EMBOTC_EMBO_aPKC <- AddModuleScore(EpitheliumMut_EMBOTC_EMBO_aPKC, features = Cholesterol, name = "Cholesterol_")
FeaturePlot(EpitheliumMut_EMBOTC_EMBO_aPKC,"Cholesterol_1", pt.size = 0.1, min.cutoff = -3, cols = c("lightgrey","red","darkred"), label = F, split.by = "Mutation")

FeaturePlot(EpitheliumMut_EMBOTC_EMBO_aPKC,"Cholesterol_1", pt.size = 0.1, min.cutoff = -3, cols = c("lightgrey","red","darkred"), label = F, split.by = "aPKC")


VlnPlot(EpitheliumMut_EMBOTC_EMBO,"Cholesterol_1", pt.size = 0.0,idents = c("aPKClow","aPKChigh"))+geom_signif(comparisons = list(c("aPKClow","aPKChigh")),test = "t.test",map_signif_level = TRUE, textsize =3) + ylim(NA,1) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)

save.image()

DimPlot(EpitheliumMut_EMBOTC_EMBO)

######TO perform harmony for tumor epihtelial cells ######
######Harmony##############
library(harmony)

EpitheliumMut_EMBOTC_EMBO_H<-RunHarmony(EpitheliumMut_EMBOTC_EMBO,group.by.vars = "SampleID", assay.use = "SCT")

EpitheliumMut_EMBOTC_EMBO_H <- RunUMAP(EpitheliumMut_EMBOTC_EMBO_H, reduction = "harmony",dims = 1:5, verbose = T)
EpitheliumMut_EMBOTC_EMBO_H <- FindNeighbors(EpitheliumMut_EMBOTC_EMBO_H,reduction = "harmony",dims = 1:5, verbose = T)
EpitheliumMut_EMBOTC_EMBO_H <- FindClusters(EpitheliumMut_EMBOTC_EMBO_H, verbose = T, resolution = 0.2)

Idents(EpitheliumMut_EMBOTC_EMBO_H)<-"Celltype"

Idents(EpitheliumMut_EMBOTC_EMBO_H)<-"SampleID"

######Which cells AMO####

DimPlot(EpitheliumMut_EMBOTC_EMBO_H, cols =c("darksalmon","purple","chocolate4", "#00B8E7","darksalmon","olivedrab3"))

DimPlot(EpitheliumMut_EMBOTC_EMBO_H)

FSC <- read.delim("./Signatures/FSC_HUMAN.txt")
EpitheliumMut_EMBOTC_EMBO_H <- AddModuleScore(EpitheliumMut_EMBOTC_EMBO_H, features = FSC, name = "FSC_")
FeaturePlot(EpitheliumMut_EMBOTC_EMBO,"FSC_1", pt.size = 0.1, min.cutoff = 0, cols = c("lightgrey","red","darkred"), label = F, split.by = "Progression")

RSC <- read.delim("./Signatures/RSC_HUMAN.txt")
EpitheliumMut_EMBOTC_EMBO <- AddModuleScore(EpitheliumMut_EMBOTC_EMBO, features = RSC, name = "RSC_")
FeaturePlot(EpitheliumMut_EMBOTC_EMBO,"RSC_1", pt.size = 0.1, min.cutoff = 0, cols = c("lightgrey","red","darkred"), label = F, split.by = "Progression")

TGC <- read.delim("./Signatures/TGC_HUMAN.txt")
EpitheliumMut_EMBOTC_EMBO <- AddModuleScore(EpitheliumMut_EMBOTC_EMBO, features = TGC, name = "TGC_")
FeaturePlot(EpitheliumMut_EMBOTC_EMBO,"TGC_1", pt.size = 0.1, min.cutoff = 0, cols = c("lightgrey","red","darkred"), label = F, split.by = "Progression")

TCC <- read.delim("./Signatures/TCC_HUMAN.txt")
EpitheliumMut_EMBOTC_EMBO <- AddModuleScore(EpitheliumMut_EMBOTC_EMBO, features = TCC, name = "TCC_")
FeaturePlot(EpitheliumMut_EMBOTC_EMBO,"TCC_1", pt.size = 0.1, min.cutoff = 0.5, cols = c("lightgrey","red","darkred"), label = F, split.by = "Progression")

Cholesterol <- read.delim("./Signatures/Cholesterol_HUMAN.txt")
EpitheliumMut_EMBOTC_EMBO <- AddModuleScore(EpitheliumMut_EMBOTC_EMBO, features = Cholesterol, name = "Cholesterol_")
FeaturePlot(EpitheliumMut_EMBOTC_EMBO,"Cholesterol_1", pt.size = 0.1, min.cutoff = 0.0, cols = c("lightgrey","red","darkred"), label = F, split.by = "Progression")

Idents(EpitheliumMut_EMBOTC_EMBO_H) <- "Progression"


VlnPlot(EpitheliumMut_EMBOTC_EMBO_H,"Cholesterol_1", pt.size = 0.0,idents = c("Conventional","Serrated"))+geom_signif(comparisons = list(c("Conventional","Serrated")),test = "t.test",map_signif_level = TRUE, textsize =3) + ylim(NA,1) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)



DotPlot(EpitheliumMut_EMBOTC_EMBO_H, features = c("TCC_1","RSC_1","FSC_1","TGC_1","Cholesterol_1"), dot.scale = 8, 
        col.min = 0,
        col.max = NA,
) + scale_colour_gradient2(low = "grey", mid = "red", high = "black",midpoint = 1.2)+RotatedAxis()

DotPlot(EpitheliumMut_EMBOTC_EMBO_H, features = c("RANBP1","MKI67","HMGB2","BIRC5","CLU","CD44","WFDC2","SOX4","ANXA1","ANXA2","MUC17","PRAP1","EMP1","ITLN1","SPINK4","TFF3","CLCA1"), dot.scale = 8, 
        col.min = 0,
        col.max = NA,
) + scale_colour_gradient2(low = "grey", mid = "red", high = "black",midpoint = 1.2)+RotatedAxis()


Idents(EpitheliumMut_EMBOTC_EMBO_H) <- "Celltype"

VlnPlot(EpitheliumMut_EMBOTC_EMBO_H,"Cholesterol_1", pt.size = 0.0,idents = c("TC1","Stem/TA","TC2","TC4","Imma.Gob.","TC3"))+geom_signif(comparisons = list(c("TC4","TC1")),test = "t.test",map_signif_level = TRUE, textsize =3) + ylim(NA,1) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)


FeaturePlot(EpitheliumMut_EMBOTC_EMBO_H, features = c("FSC_1","Cholesterol_1"),min.cutoff = -0.5, blend = T,cols=c("red","green"))

FeaturePlot(EpitheliumMut_EMBOTC_EMBO_H, features = c("LGR5"),cols=c("lightgrey","red","darkred"), split.by = "Progression")
FeaturePlot(EpitheliumMut_EMBOTC_EMBO_H, features = c("OLFM4"),cols=c("lightgrey","red","darkred"), split.by = "Progression")
FeaturePlot(EpitheliumMut_EMBOTC_EMBO_H, features = c("PRKCI"),cols=c("lightgrey","red","darkred"), split.by = "Progression")

VlnPlot(EpitheliumMut_EMBOTC_EMBO_H,"PRKCI")


FeaturePlot(EpitheliumMut_EMBOTC_EMBO_H, features = c("LGR5"),cols=c("lightgrey","red","darkred"), split.by = "Progression")

VlnPlot(EpitheliumMut_EMBOTC_EMBO_H,"LGR5")

########
####Set ident Olfm4 positive celltype group####################
Idents(EpitheliumMut_EMBOTC_EMBO_H)<-"Celltype"

Lhigh<-WhichCells(EpitheliumMut_EMBOTC_EMBO_H,slot = "counts", expression = PRKCI>0.99)
EpitheliumMut_EMBOTC_EMBO_H@meta.data$PRKCI_status<-"Lhigh"
EpitheliumMut_EMBOTC_EMBO_H<-SetIdent(EpitheliumMut_EMBOTC_EMBO_H, cells = Lhigh, value ="Lhigh" )
EpitheliumMut_EMBOTC_EMBO_H$PRKCI_status<-Idents(EpitheliumMut_EMBOTC_EMBO_H)
####Lgr5 neg
Llow<-WhichCells(EpitheliumMut_EMBOTC_EMBO_H,slot = "counts", expression = PRKCI<0.99)
EpitheliumMut_EMBOTC_EMBO_H@meta.data$PRKCI_status<-"Llow"
EpitheliumMut_EMBOTC_EMBO_H<-SetIdent(EpitheliumMut_EMBOTC_EMBO_H, cells = Llow, value ="Llow" )
EpitheliumMut_EMBOTC_EMBO_H$PRKCI_status<-Idents(EpitheliumMut_EMBOTC_EMBO_H)

######
Idents(EpitheliumMut_EMBOTC_EMBO_H)<-"PRKCI_status"

DimPlot(EpitheliumMut_EMBOTC_EMBO_H)

VlnPlot(EpitheliumMut_EMBOTC_EMBO_H,"OLFM4", pt.size = 0.0,idents = c("Llow","Lhigh"))+geom_signif(comparisons = list(c("Llow","Lhigh")),test = "t.test",map_signif_level = TRUE, textsize =3) + ylim(NA,6) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)
VlnPlot(EpitheliumMut_EMBOTC_EMBO_H,"LGR5", pt.size = 0.0,idents = c("Llow","Lhigh"))+geom_signif(comparisons = list(c("Llow","Lhigh")),test = "t.test",map_signif_level = TRUE, textsize =3) + ylim(NA,3) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)


FeaturePlot(EpitheliumMut_EMBOTC_EMBO_H, features = c("PRKCI"),cols=c("lightgrey","darkred"), split.by = "PRKCI_status")
VlnPlot(EpitheliumMut_EMBOTC_EMBO_H, "PRKCI")

FeaturePlot(EpitheliumMut_EMBOTC_EMBO_H, features = c("OLFM4"),cols=c("lightgrey","red","darkred"), split.by = "PRKCI_status")

DimPlot(EpitheliumMut_EMBOTC_EMBO_H)

######
YapWang_HUMAN<- read.delim("./Signatures/YAPwang_Human.txt")
EpitheliumMut_EMBOTC_EMBO_H <- AddModuleScore(EpitheliumMut_EMBOTC_EMBO_H, features = YapWang_HUMAN, name = "YapWang_HUMAN_")

WntVander_HUMAN<- read.delim("./Signatures/WNT_FlierHuman.txt")
EpitheliumMut_EMBOTC_EMBO_H <- AddModuleScore(EpitheliumMut_EMBOTC_EMBO_H, features = WntVander_HUMAN, name = "WntVander_HUMAN_")

iCMS3signature<- read.delim("./Signatures/iCMS3_Human.txt")
EpitheliumMut_EMBOTC_EMBO_H <- AddModuleScore(EpitheliumMut_EMBOTC_EMBO_H, features = iCMS3signature, name = "iCMS3signature_")

iCMS2signature<- read.delim("./Signatures/iCMS2_Human.txt")
EpitheliumMut_EMBOTC_EMBO_H <- AddModuleScore(EpitheliumMut_EMBOTC_EMBO_H, features = iCMS2signature, name = "iCMS2signature_")

MerlosStem_HUMAN<- read.delim("./Signatures/MerlosStem_Human.txt")
EpitheliumMut_EMBOTC_EMBO_H <- AddModuleScore(EpitheliumMut_EMBOTC_EMBO_H, features = MerlosStem_HUMAN, name = "MerlosStem_HUMAN_")

KleemanStem_HUMAN<- read.delim("./Signatures/Kleeman_Stem_HUMAN.txt")
EpitheliumMut_EMBOTC_EMBO_H <- AddModuleScore(EpitheliumMut_EMBOTC_EMBO_H, features = KleemanStem_HUMAN, name = "KleemanStem_HUMAN_")

JDP2_HUMAN<- read.delim("./Signatures/JDP2_Human.txt")
EpitheliumMut_EMBOTC_EMBO_H <- AddModuleScore(EpitheliumMut_EMBOTC_EMBO_H, features = JDP2_HUMAN, name = "JDP2_HUMAN_")

JunD_HUMAN<- read.delim("./Signatures/JunD_HUMAN.txt")
EpitheliumMut_EMBOTC_EMBO_H <- AddModuleScore(EpitheliumMut_EMBOTC_EMBO_H, features = JunD_HUMAN, name = "JunD_HUMAN_")

JunB_HUMAN<- read.delim("./Signatures/JunB_HUMAN.txt")
EpitheliumMut_EMBOTC_EMBO_H <- AddModuleScore(EpitheliumMut_EMBOTC_EMBO_H, features = JunB_HUMAN, name = "JunB_HUMAN_")

YAPGrego_HUMAN<- read.delim("./Signatures/YAPGrego_Human.txt")
EpitheliumMut_EMBOTC_EMBO_H <- AddModuleScore(EpitheliumMut_EMBOTC_EMBO_H, features = YAPGrego_HUMAN, name = "YAPGrego_HUMAN_")

YAPCommon_HUMAN<- read.delim("./Signatures/YAP_commonH.txt")
EpitheliumMut_EMBOTC_EMBO_H <- AddModuleScore(EpitheliumMut_EMBOTC_EMBO_H, features = YAPCommon_HUMAN, name = "YAPCommon_HUMAN_")

AP1Common_HUMAN<- read.delim("./Signatures/AP1_commonH.txt")
EpitheliumMut_EMBOTC_EMBO_H <- AddModuleScore(EpitheliumMut_EMBOTC_EMBO_H, features = AP1Common_HUMAN, name = "AP1Common_HUMAN_")


DotPlot(EpitheliumMut_EMBOTC_EMBO_H, features = c("YapWang_HUMAN_1","WntVander_HUMAN_1","iCMS3signature_1","iCMS2signature_1","MerlosStem_HUMAN_1","KleemanStem_HUMAN_1","JDP2_HUMAN_1","JunD_HUMAN_1","JunB_HUMAN_1"), dot.scale = 8, 
        col.min = 0,
        col.max = 1.5,
) + scale_colour_gradient2(low = "grey", mid = "red", high = "black",midpoint = 0.75)+RotatedAxis()


FeatureScatter(object = EpitheliumMut_EMBOTC_EMBO_H, feature1 = 'YapWang_HUMAN_1', feature2 = 'JunD_HUMAN_1', cols=c("darksalmon","purple","chocolate4", "#00B8E7","darksalmon","olivedrab3"))
FeatureScatter(object = EpitheliumMut_EMBOTC_EMBO_H, feature1 = 'YapWang_HUMAN_1', feature2 = 'WntVander_HUMAN_1', cols=c("darksalmon","purple","chocolate4", "#00B8E7","darksalmon","olivedrab3"))
FeatureScatter(object = EpitheliumMut_EMBOTC_EMBO_H, feature1 = 'KleemanStem_HUMAN_1', feature2 = 'JunD_HUMAN_1', cols=c("darksalmon","purple","chocolate4", "#00B8E7","darksalmon","olivedrab3"))

####generate patient ID aPKC low ######
#######Average expression ######
Idents(EpitheliumMut_EMBOTC_EMBO_H)<-"SampleID"
avg.EpitheliumMut_EMBOTC_EMBO_H<-AverageExpression(EpitheliumMut_EMBOTC_EMBO_H) 
write.table(as.matrix(avg.EpitheliumMut_EMBOTC_EMBO_H$RNA), "./Results/GeneExpression///avg.EpitheliumMut_EMBOTC_EMBO_H_RNA.txt", sep="\t",col.names = T, row.names = T)
write.table(as.matrix(avg.EpitheliumMut_EMBOTC_EMBO_H$SCT), "./Results/GeneExpression//avg.EpitheliumMut_EMBOTC_EMBO_H_SCT.txt", sep="\t",col.names = T, row.names = T)

#############
#########generate ident #### aPKC
Idents(EpitheliumMut_EMBOTC_EMBO_H)<-"SampleID"
aPKChigh<-WhichCells(EpitheliumMut_EMBOTC_EMBO_H, ident= c("p008t","p013t","p017t"))
EpitheliumMut_EMBOTC_EMBO_H@meta.data$aPKC [aPKChigh]
EpitheliumMut_EMBOTC_EMBO_H<-SetIdent(EpitheliumMut_EMBOTC_EMBO_H, cells = aPKChigh, value = "aPKChigh")
EpitheliumMut_EMBOTC_EMBO_H$aPKC<-Idents (EpitheliumMut_EMBOTC_EMBO_H)

aPKClow<-WhichCells (EpitheliumMut_EMBOTC_EMBO_H, ident= c("p026t","p014t","p025t","p021t"))
EpitheliumMut_EMBOTC_EMBO_H@meta.data$aPKC [aPKClow]
EpitheliumMut_EMBOTC_EMBO_H<-SetIdent(EpitheliumMut_EMBOTC_EMBO_H, cells = aPKClow, value = "aPKClow")
EpitheliumMut_EMBOTC_EMBO_H$aPKC<-Idents (EpitheliumMut_EMBOTC_EMBO_H)

aPKCmedium<-WhichCells (EpitheliumMut_EMBOTC_EMBO_H, ident= c("p009t1","p009t2","p020t","p012t","p016t","p007t"))
EpitheliumMut_EMBOTC_EMBO_H@meta.data$aPKC [aPKCmedium]
EpitheliumMut_EMBOTC_EMBO_H<-SetIdent(EpitheliumMut_EMBOTC_EMBO_H, cells = aPKCmedium, value = "aPKCmedium")
EpitheliumMut_EMBOTC_EMBO_H$aPKC<-Idents (EpitheliumMut_EMBOTC_EMBO_H)


Idents(EpitheliumMut_EMBOTC_EMBO_H)<-"aPKC"

FeaturePlot(EpitheliumMut_EMBOTC_EMBO_H,"OLFM4", split.by ="aPKC", cols = c("lightgrey","red","darkred"), min.cutoff = 0.8)


Idents(EpitheliumMut_EMBOTC_EMBO_H) <- "aPKC"
EpitheliumMut_EMBOTC_EMBO_H_aPKC <- subset(EpitheliumMut_EMBOTC_EMBO_H, idents = c("aPKClow","aPKChigh"))

DotPlot(EpitheliumMut_EMBOTC_EMBO_H_aPKC, features = c("OLFM4","LGR5","AXIN2","ASCL2","CLU","CD44","ANXA1","ANXA10","MUC17"), dot.scale = 8, 
        col.min = 0,
        col.max = 1,
) + scale_colour_gradient2(low = "grey", mid = "red", high = "black",midpoint = 0.5)+RotatedAxis()


DotPlot(EpitheliumMut_EMBOTC_EMBO_H_aPKC, features = c("YapWang_HUMAN_1","WntVander_HUMAN_1","iCMS3signature_1","iCMS2signature_1","MerlosStem_HUMAN_1","KleemanStem_HUMAN_1","JDP2_HUMAN_1","JunD_HUMAN_1","JunB_HUMAN_1"), dot.scale = 8, 
        col.min = 0,
        col.max = 1,
) + scale_colour_gradient2(low = "grey", mid = "red", high = "black",midpoint = 0.5)+RotatedAxis()





VlnPlot(EpitheliumMut_EMBOTC_EMBO_H,"OLFM4",pt.size = 0.0,idents = c("aPKChigh","aPKClow"))+geom_signif(comparisons = list(c("aPKChigh","aPKClow")),test = "t.test",map_signif_level = TRUE, textsize =3) + ylim(NA,6) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)
VlnPlot(EpitheliumMut_EMBOTC_EMBO_H,"LGR5",pt.size = 0.0,idents = c("aPKChigh","aPKClow"))+geom_signif(comparisons = list(c("aPKChigh","aPKClow")),test = "t.test",map_signif_level = TRUE, textsize =3) + ylim(NA,3) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)
VlnPlot(EpitheliumMut_EMBOTC_EMBO_H,"CLU",pt.size = 0.0,idents = c("aPKChigh","aPKClow"))+geom_signif(comparisons = list(c("aPKChigh","aPKClow")),test = "t.test",map_signif_level = TRUE, textsize =3) + ylim(NA,3) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)
VlnPlot(EpitheliumMut_EMBOTC_EMBO_H,"ANXA1",pt.size = 0.0,idents = c("aPKChigh","aPKClow"))+geom_signif(comparisons = list(c("aPKChigh","aPKClow")),test = "t.test",map_signif_level = TRUE, textsize =3) + ylim(NA,3) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)
VlnPlot(EpitheliumMut_EMBOTC_EMBO_H,"CD44",pt.size = 0.0,idents = c("aPKChigh","aPKClow"))+geom_signif(comparisons = list(c("aPKChigh","aPKClow")),test = "t.test",map_signif_level = TRUE, textsize =3) + ylim(NA,3) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)





head(EpitheliumMut_EMBOTC_EMBO_H@meta.data)




#Cell.counts
head(EpitheliumMut_EMBOTC_EMBO_H@meta.data,2)
#Add group_new and cluster info to meta.data
Idents(EpitheliumMut_EMBOTC_EMBO_H) <- "Celltype"
EpitheliumMut_EMBOTC_EMBO_H$Celltype.aPKC <- paste(Idents(EpitheliumMut_EMBOTC_EMBO_H), EpitheliumMut_EMBOTC_EMBO_H$aPKC, sep = "_")
head(EpitheliumMut_EMBOTC_EMBO_H@meta.data,3)

Idents(EpitheliumMut_EMBOTC_EMBO_H) <- "Celltype.aPKC"
counts.EpitheliumMut_EMBOTC_EMBO_H <-table(Idents(EpitheliumMut_EMBOTC_EMBO_H))
write.table(as.matrix(counts.EpitheliumMut_EMBOTC_EMBO_H),"./Results/Counts/counts.EpitheliumMut_EMBOTC_EMBO_H.txt",sep="\t",col.names=T,row.names=T)


######Checking signaling Stem cell paper ##########
Idents(EpitheliumMut_EMBOTC_EMBO_H) <- "Celltype"

TRSC<-WhichCells(EpitheliumMut_EMBOTC_EMBO_H, ident= c("TC3","Stem/TA"))
EpitheliumMut_EMBOTC_EMBO_H@meta.data$phenotype [TRSC]
EpitheliumMut_EMBOTC_EMBO_H<-SetIdent(EpitheliumMut_EMBOTC_EMBO_H, cells = TRSC, value = "TRSC")
EpitheliumMut_EMBOTC_EMBO_H$phenotype<-Idents (EpitheliumMut_EMBOTC_EMBO_H)
TFMC<-WhichCells (EpitheliumMut_EMBOTC_EMBO_H, ident= c("TC4"))
EpitheliumMut_EMBOTC_EMBO_H@meta.data$phenotype [TFMC]
EpitheliumMut_EMBOTC_EMBO_H<-SetIdent(EpitheliumMut_EMBOTC_EMBO_H, cells = TFMC, value = "TFMC")
EpitheliumMut_EMBOTC_EMBO_H$phenotype<-Idents (EpitheliumMut_EMBOTC_EMBO_H)
TcTA<-WhichCells (EpitheliumMut_EMBOTC_EMBO_H, ident= c("TC1"))
EpitheliumMut_EMBOTC_EMBO_H@meta.data$phenotype [TcTA]
EpitheliumMut_EMBOTC_EMBO_H<-SetIdent(EpitheliumMut_EMBOTC_EMBO_H, cells = TcTA, value = "TcTA")
EpitheliumMut_EMBOTC_EMBO_H$phenotype<-Idents (EpitheliumMut_EMBOTC_EMBO_H)
TGC<-WhichCells (EpitheliumMut_EMBOTC_EMBO_H, ident= c("Imma. Gob."))
EpitheliumMut_EMBOTC_EMBO_H@meta.data$phenotype [TGC]
EpitheliumMut_EMBOTC_EMBO_H<-SetIdent(EpitheliumMut_EMBOTC_EMBO_H, cells = TGC, value = "TGC")
EpitheliumMut_EMBOTC_EMBO_H$phenotype<-Idents (EpitheliumMut_EMBOTC_EMBO_H)
TC2<-WhichCells (EpitheliumMut_EMBOTC_EMBO_H, ident= c("TC2"))
EpitheliumMut_EMBOTC_EMBO_H@meta.data$phenotype [TC2]
EpitheliumMut_EMBOTC_EMBO_H<-SetIdent(EpitheliumMut_EMBOTC_EMBO_H, cells = TC2, value = "TC2")
EpitheliumMut_EMBOTC_EMBO_H$phenotype<-Idents (EpitheliumMut_EMBOTC_EMBO_H)

Idents(EpitheliumMut_EMBOTC_EMBO_H) <- "phenotype"

my_levels <- c("TcTA","TRSC","TFMC","TGC","TC2")
EpitheliumMut_EMBOTC_EMBO_H@meta.data$phenotype <- factor(EpitheliumMut_EMBOTC_EMBO_H@meta.data$phenotype, levels = my_levels)
Idents(EpitheliumMut_EMBOTC_EMBO_H) <- "phenotype"


DotPlot(EpitheliumMut_EMBOTC_EMBO_H, features = c("YapWang_HUMAN_1","YAPGrego_HUMAN_1","JDP2_HUMAN_1","JunD_HUMAN_1","JunB_HUMAN_1"), dot.scale = 8, 
        col.min = 0,
        col.max = NA,
) + scale_colour_gradient2(low = "grey", mid = "red", high = "darkred",midpoint = 0.85)+RotatedAxis()

VlnPlot(EpitheliumMut_EMBOTC_EMBO_H,"JDP2_HUMAN_1", pt.size = 0.0,idents = c("TC2","TGC","TFMC","TRSC","TcTA"))+geom_signif(comparisons = list(c("TRSC","TFMC")),test = "t.test",map_signif_level = TRUE, textsize =3) + ylim(NA,1) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)
VlnPlot(EpitheliumMut_EMBOTC_EMBO_H,"JunD_HUMAN_1", pt.size = 0.0,idents = c("TC2","TGC","TFMC","TRSC","TcTA"))+geom_signif(comparisons = list(c("TRSC","TFMC")),test = "t.test",map_signif_level = TRUE, textsize =3) + ylim(NA,1) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)
VlnPlot(EpitheliumMut_EMBOTC_EMBO_H,"JunB_HUMAN_1", pt.size = 0.0,idents = c("TC2","TGC","TFMC","TRSC","TcTA"))+geom_signif(comparisons = list(c("TRSC","TFMC")),test = "t.test",map_signif_level = TRUE, textsize =3) + ylim(NA,1) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)
VlnPlot(EpitheliumMut_EMBOTC_EMBO_H,"YapWang_HUMAN_1", pt.size = 0.0,idents = c("TC2","TGC","TFMC","TRSC","TcTA"))+geom_signif(comparisons = list(c("TRSC","TFMC")),test = "t.test",map_signif_level = TRUE, textsize =3) + ylim(NA,0.75) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)
VlnPlot(EpitheliumMut_EMBOTC_EMBO_H,"YAPGrego_HUMAN_1", pt.size = 0.0,idents = c("TC2","TGC","TFMC","TRSC","TcTA"))+geom_signif(comparisons = list(c("TRSC","TFMC")),test = "t.test",map_signif_level = TRUE, textsize =3) + ylim(NA,0.5) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)
VlnPlot(EpitheliumMut_EMBOTC_EMBO_H,"YAPCommon_HUMAN_1", pt.size = 0.0,idents = c("TC2","TGC","TFMC","TRSC","TcTA"))+geom_signif(comparisons = list(c("TRSC","TFMC")),test = "t.test",map_signif_level = TRUE, textsize =3) + ylim(NA,0.5) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)
VlnPlot(EpitheliumMut_EMBOTC_EMBO_H,"AP1Common_HUMAN_1", pt.size = 0.0,idents = c("TC2","TGC","TFMC","TRSC","TcTA"))+geom_signif(comparisons = list(c("TRSC","TFMC")),test = "t.test",map_signif_level = TRUE, textsize =3) + ylim(NA,0.5) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)


######
####Subset double positive cells####

DoublePos<-WhichCells(EpitheliumMut_EMBOTC_EMBO_H,slot = "counts", expression = AP1Common_HUMAN_1>0.05&YAPCommon_HUMAN_1>0.05)
EpitheliumMut_EMBOTC_EMBO_H@meta.data$YAPAP1<-"DoublePos"
EpitheliumMut_EMBOTC_EMBO_H<-SetIdent(EpitheliumMut_EMBOTC_EMBO_H, cells = DoublePos, value ="DoublePos")
EpitheliumMut_EMBOTC_EMBO_H$YAPAP1<-Idents(EpitheliumMut_EMBOTC_EMBO_H)

YAPhighAP1low<-WhichCells(EpitheliumMut_EMBOTC_EMBO_H,slot = "counts", expression = AP1Common_HUMAN_1<0.05&YAPCommon_HUMAN_1>0.05)
EpitheliumMut_EMBOTC_EMBO_H@meta.data$YAPAP1<-"YAPhighAP1low"
EpitheliumMut_EMBOTC_EMBO_H<-SetIdent(EpitheliumMut_EMBOTC_EMBO_H, cells = YAPhighAP1low, value ="YAPhighAP1low")
EpitheliumMut_EMBOTC_EMBO_H$YAPAP1<-Idents(EpitheliumMut_EMBOTC_EMBO_H)


YAPlowAP1high<-WhichCells(EpitheliumMut_EMBOTC_EMBO_H,slot = "counts", expression = AP1Common_HUMAN_1>0.05&YAPCommon_HUMAN_1<0.05)
EpitheliumMut_EMBOTC_EMBO_H@meta.data$YAPAP1<-"YAPlowAP1high"
EpitheliumMut_EMBOTC_EMBO_H<-SetIdent(EpitheliumMut_EMBOTC_EMBO_H, cells = YAPlowAP1high, value ="YAPlowAP1high")
EpitheliumMut_EMBOTC_EMBO_H$YAPAP1<-Idents(EpitheliumMut_EMBOTC_EMBO_H)


DoubleNeg<-WhichCells(EpitheliumMut_EMBOTC_EMBO_H,slot = "counts", expression = AP1Common_HUMAN_1<0.05&YAPCommon_HUMAN_1<0.05)
EpitheliumMut_EMBOTC_EMBO_H@meta.data$YAPAP1<-"DoubleNeg"
EpitheliumMut_EMBOTC_EMBO_H<-SetIdent(EpitheliumMut_EMBOTC_EMBO_H, cells = DoubleNeg, value ="DoubleNeg")
EpitheliumMut_EMBOTC_EMBO_H$YAPAP1<-Idents(EpitheliumMut_EMBOTC_EMBO_H)


Idents(EpitheliumMut_EMBOTC_EMBO_H) <- "YAPAP1"
DimPlot(EpitheliumMut_EMBOTC_EMBO_H, cols = c("lightgrey","darkgoldenrod1","cadetblue3","tomato2"))

FeaturePlot(EpitheliumMut_EMBOTC_EMBO_H, features = c("AP1Common_HUMAN_1","YAPCommon_HUMAN_1"),min.cutoff = -0.5, blend = T,cols=c("red","green"))

#####Counts proportion double pos####
#Add group_new and cluster info to meta.data
Idents(EpitheliumMut_EMBOTC_EMBO_H) <- "phenotype"
EpitheliumMut_EMBOTC_EMBO_H$phenotype.YAPAP1 <- paste(Idents(EpitheliumMut_EMBOTC_EMBO_H), EpitheliumMut_EMBOTC_EMBO_H$YAPAP1, sep = "_")
head(EpitheliumMut_EMBOTC_EMBO_H@meta.data,3)

Idents(EpitheliumMut_EMBOTC_EMBO_H) <- "phenotype.YAPAP1"
counts.YAP_JNK_Epithelium_EMBO <-table(Idents(EpitheliumMut_EMBOTC_EMBO_H))
write.table(as.matrix(counts.YAP_JNK_Epithelium_EMBO),"./Results/Counts/counts.YAP_JNK_Epithelium_EMBO.txt",sep="\t",col.names=T,row.names=T)

#########For Stem cell paper, just non-secretory populations###
Idents(EpitheliumMut_EMBOTC_EMBO_H) <- "phenotype"

EpitheliumMut_EMBOTC_EMBO_H_NonSecretory <- subset(EpitheliumMut_EMBOTC_EMBO_H, idents = c("TFMC","TRSC","TcTA"))
library(harmony)

EpitheliumMut_EMBOTC_EMBO_H_NonSecretory<-RunHarmony(EpitheliumMut_EMBOTC_EMBO_H_NonSecretory,group.by.vars = "SampleID", assay.use = "SCT")

EpitheliumMut_EMBOTC_EMBO_H_NonSecretory <- RunUMAP(EpitheliumMut_EMBOTC_EMBO_H_NonSecretory, reduction = "harmony",dims = 1:5, verbose = T)
EpitheliumMut_EMBOTC_EMBO_H_NonSecretory <- FindNeighbors(EpitheliumMut_EMBOTC_EMBO_H_NonSecretory,reduction = "harmony",dims = 1:5, verbose = T)
EpitheliumMut_EMBOTC_EMBO_H_NonSecretory <- FindClusters(EpitheliumMut_EMBOTC_EMBO_H_NonSecretory, verbose = T, resolution = 0.2)

Idents(EpitheliumMut_EMBOTC_EMBO_H_NonSecretory)<-"phenotype"

DimPlot(EpitheliumMut_EMBOTC_EMBO_H_NonSecretory)

DimPlot(EpitheliumMut_EMBOTC_EMBO_H_NonSecretory, cols =c("olivedrab3","darksalmon","chocolate4"))


DotPlot(EpitheliumMut_EMBOTC_EMBO_H_NonSecretory, features = c("YapWang_HUMAN_1","YAPGrego_HUMAN_1","JDP2_HUMAN_1","JunD_HUMAN_1","JunB_HUMAN_1"), dot.scale = 8, 
        col.min = 0,
        col.max = NA,
) + scale_colour_gradient2(low = "grey", mid = "red", high = "darkred",midpoint = 0.6)+RotatedAxis()

#####
VlnPlot(EpitheliumMut_EMBOTC_EMBO_H_NonSecretory,"YAPCommon_HUMAN_1", pt.size = 0.0,idents = c("TFMC","TRSC","TcTA"))+geom_signif(comparisons = list(c("TRSC","TFMC")),test = "t.test",map_signif_level = TRUE, textsize =3) + ylim(NA,0.5) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)
VlnPlot(EpitheliumMut_EMBOTC_EMBO_H_NonSecretory,"AP1Common_HUMAN_1", pt.size = 0.0,idents = c("TFMC","TRSC","TcTA"))+geom_signif(comparisons = list(c("TRSC","TcTA")),test = "t.test",map_signif_level = TRUE, textsize =3) + ylim(NA,0.5) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)


Idents(EpitheliumMut_EMBOTC_EMBO_H_NonSecretory) <- "YAPAP1"
DimPlot(EpitheliumMut_EMBOTC_EMBO_H_NonSecretory, cols = c("lightgrey","darkgoldenrod1","cadetblue3","tomato2"))



Idents(EpitheliumMut_EMBOTC_EMBO_H_NonSecretory)<-"phenotype"

VlnPlot(EpitheliumMut_EMBOTC_EMBO_H_NonSecretory,"Cholesterol_1", pt.size = 0.0,idents = c("TFMC","TRSC","TcTA"))+geom_signif(comparisons = list(c("TRSC","TFMC")),test = "t.test",map_signif_level = TRUE, textsize =3) + ylim(NA,1) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)

FeaturePlot(EpitheliumMut_EMBOTC_EMBO_H_NonSecretory,"Cholesterol_1", pt.size = 0.1, min.cutoff = -0.1, cols = c("lightgrey","red","darkred"), label = F)

#########
ReactomeActivationHuman<- read.delim("./Signatures/ReactomeActivationHuman.txt")
EpitheliumMut_EMBOTC_EMBO_H_NonSecretory <- AddModuleScore(EpitheliumMut_EMBOTC_EMBO_H_NonSecretory, features = ReactomeActivationHuman, name = "ReactomeActivationHuman_")
FeaturePlot(EpitheliumMut_EMBOTC_EMBO_H_NonSecretory,"ReactomeActivationHuman_1", pt.size = 0.1, min.cutoff = -0.2, cols = c("lightgrey","red","darkred"), label = F)


HortonHuman<- read.delim("./Signatures/Horton_HUMAN.txt")
EpitheliumMut_EMBOTC_EMBO_H_NonSecretory <- AddModuleScore(EpitheliumMut_EMBOTC_EMBO_H_NonSecretory, features = HortonHuman, name = "HortonHuman_")
FeaturePlot(EpitheliumMut_EMBOTC_EMBO_H_NonSecretory,"HortonHuman_1", pt.size = 0.1, min.cutoff = -0.2, cols = c("lightgrey","red","darkred"), label = F)




VlnPlot(EpitheliumMut_EMBOTC_EMBO_H_NonSecretory,cols = c("grey75","indianred2","indianred4"),"HortonHuman_1", pt.size = 0.0,idents = c("TcTA","TRSC","TFMC"))+geom_signif(comparisons = list(c("TFMC","TRSC")),map_signif_level = TRUE, textsize =3) + ylim(NA,2) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)




########################
####Subset cholesterol high cells####

CholesPos<-WhichCells(EpitheliumMut_EMBOTC_EMBO_H_NonSecretory,slot = "counts", expression =Cholesterol_1>0.1)
EpitheliumMut_EMBOTC_EMBO_H_NonSecretory@meta.data$Choles<-"CholesPos"
EpitheliumMut_EMBOTC_EMBO_H_NonSecretory<-SetIdent(EpitheliumMut_EMBOTC_EMBO_H_NonSecretory, cells = CholesPos, value ="CholesPos")
EpitheliumMut_EMBOTC_EMBO_H_NonSecretory$Choles<-Idents(EpitheliumMut_EMBOTC_EMBO_H_NonSecretory)


CholesNeg<-WhichCells(EpitheliumMut_EMBOTC_EMBO_H_NonSecretory,slot = "counts", expression = Cholesterol_1<0.1)
EpitheliumMut_EMBOTC_EMBO_H_NonSecretory@meta.data$Choles<-"CholesNeg"
EpitheliumMut_EMBOTC_EMBO_H_NonSecretory<-SetIdent(EpitheliumMut_EMBOTC_EMBO_H_NonSecretory, cells = CholesNeg, value ="CholesNeg")
EpitheliumMut_EMBOTC_EMBO_H_NonSecretory$Choles<-Idents(EpitheliumMut_EMBOTC_EMBO_H_NonSecretory)

Idents(EpitheliumMut_EMBOTC_EMBO_H_NonSecretory) <- "Choles"
DimPlot(EpitheliumMut_EMBOTC_EMBO_H_NonSecretory, cols = c("lightgrey","tomato2"))



#####Counts proportion double pos####
#Add group_new and cluster info to meta.data
Idents(EpitheliumMut_EMBOTC_EMBO_H_NonSecretory) <- "phenotype"
EpitheliumMut_EMBOTC_EMBO_H_NonSecretory$phenotype.Choles <- paste(Idents(EpitheliumMut_EMBOTC_EMBO_H_NonSecretory), EpitheliumMut_EMBOTC_EMBO_H_NonSecretory$Choles, sep = "_")
head(EpitheliumMut_EMBOTC_EMBO_H_NonSecretory@meta.data,3)

Idents(EpitheliumMut_EMBOTC_EMBO_H_NonSecretory) <- "phenotype.Choles"
counts.Choles_Epithelium_Koreans <-table(Idents(EpitheliumMut_EMBOTC_EMBO_H_NonSecretory))
write.table(as.matrix(counts.Choles_Epithelium_Koreans),"./Results/Counts/counts.Choles_Epithelium_EMBO.txt",sep="\t",col.names=T,row.names=T)






######Peter checkins=g SOX2 and SOX2 targets####

SOX2targets <- read.delim("./Signatures/Sox2_RP_HUMAN.txt")
EpitheliumMut_EMBOTC_EMBO_H <- AddModuleScore(EpitheliumMut_EMBOTC_EMBO_H, features = SOX2targets, name = "SOX2targets_")
FeaturePlot(EpitheliumMut_EMBOTC_EMBO_H,"SOX2targets_1", pt.size = 0.1, min.cutoff = 0.0, cols = c("lightgrey","red","darkred"), label = F, split.by = "Progression")



Idents(EpitheliumMut_EMBOTC_EMBO_H) <- "Celltype.Progression"

DotPlot(EpitheliumMut_EMBOTC_EMBO_H, features = c("PRKCI","PRKCZ","SOX2","MUC5AC","ANXA10","ANXA1","CLU","OLFM4","SOX2targets_1","FSC_1"),dot.scale = 8, 
        col.min = 0,
        col.max = 1)+scale_colour_gradient2(low = "grey", mid = "red", high = "black",midpoint = 0.6)+RotatedAxis()

Idents(EpitheliumMut_EMBOTC_EMBO_H) <- "Celltype.Progression"

Idents(EpitheliumMut_EMBOTC_EMBO_H) <- "Celltype"

EpitheliumMut_TC4 <- subset(EpitheliumMut_EMBOTC_EMBO_H, idents = c("TC4"))

Idents(EpitheliumMut_TC4) <- "Celltype.Progression"
head(EpitheliumMut_TC4@meta.data,10)

VlnPlot(EpitheliumMut_TC4,"SOX2", pt.size = 0.0,idents = c("TC4_Conventional","TC4_Serrated"))+geom_signif(comparisons = list(c("TC4_Conventional","TC4_Serrated")),test = "t.test",map_signif_level = TRUE, textsize =3) + ylim(NA,2) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)


######
save.image()
