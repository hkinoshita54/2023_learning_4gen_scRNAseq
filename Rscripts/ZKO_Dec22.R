######
library(Seurat)
library(DropletUtils)
library(ggplot2)
library(DoubletFinder)
library(SingleR)
library(dplyr)
library(celldex)
library(knitr)
library(RColorBrewer)



#####Try to merge Seurat objects WT and DKOn (from WCM integration) and DKO tumor from Cancer cell#############

###Fibroblast and tumors First combination Fibroblast DKO veh and Tumor DKO veh#########



Epithelium_WT_DKOn.aggr <-readRDS(file = "./Results/Combined_Epithelium_WCM_SD_WT_DKO")
Epithelium_PEGPH20.aggr <-readRDS(file = "./Results/RDSobjects/Epithelium_DKO_HARM_UMAP_MT12_Clean4")


DimPlot(Epithelium_WT_DKOn.aggr, split.by = "group.batch")
DimPlot(Epithelium_PEGPH20.aggr, split.by = "group")

Idents(Epithelium_PEGPH20.aggr) <- "group"

Epithelium_T.aggr <- subset(Epithelium_PEGPH20.aggr, idents = c("DKOtVEH"))

###########

####Merge seurat objetcs


Epithelium_T.aggr$dataset <- 'SD'
Epithelium_WT_DKOn.aggr$dataset <- 'WTDKO'


# merge all datasets, adding a cell ID to make sure cell names are unique
Combined_WT_DKOn_DKONT <- merge(
  x = Epithelium_T.aggr,
  y = list(Epithelium_WT_DKOn.aggr),
  add.cell.ids = c("SD", "WTDKO")
)

tail(Combined_WT_DKOn_DKONT@meta.data, 10)

#####
Combined_WT_DKOn_DKONT  <- SCTransform(Combined_WT_DKOn_DKONT , assay="RNA", vars.to.regress = "percent.mt", verbose = T)
Combined_WT_DKOn_DKONT <- RunPCA(Combined_WT_DKOn_DKONT, verbose = T)


#c("dataset", "donor", "batch_id"))
#####Harmony#######
Combined_WT_DKOn_DKONT_Harmony<-RunHarmony(Combined_WT_DKOn_DKONT,group.by.vars =c("Batch"), assay.use = "SCT")

###########
ElbowPlot(Combined_WT_DKOn_DKONT_Harmony)
Combined_WT_DKOn_DKONT_Harmony <- RunUMAP(Combined_WT_DKOn_DKONT_Harmony, reduction = "harmony",dims = 1:11, verbose = T)
Combined_WT_DKOn_DKONT_Harmony <- FindNeighbors(Combined_WT_DKOn_DKONT_Harmony,reduction = "harmony",dims = 1:11, verbose = T)
Combined_WT_DKOn_DKONT_Harmony <- FindClusters(Combined_WT_DKOn_DKONT_Harmony, verbose = T, resolution = 0.8)

VlnPlot(Combined_WT_DKOn_DKONT_Harmony,"nFeature_RNA")

Epithelium_DKO_HARM_UMAP_MT12 <- subset(x = Combined_WT_DKOn_DKONT_Harmony, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 12 )


#Cell.counts
head(Epithelium_DKO_HARM_UMAP_MT12@meta.data,2)
#Add group_new and cluster info to meta.data
Idents(Epithelium_DKO_HARM_UMAP_MT12) <- "celltype"
Epithelium_DKO_HARM_UMAP_MT12$celltype.group <- paste(Idents(Epithelium_DKO_HARM_UMAP_MT12), Epithelium_DKO_HARM_UMAP_MT12$group, sep = "_")
head(Epithelium_DKO_HARM_UMAP_MT12@meta.data,3)

Idents(Epithelium_DKO_HARM_UMAP_MT12) <- "celltype.group"


Epithelium_DKO_HARM_UMAP_MT12_FINAL <- subset(Epithelium_DKO_HARM_UMAP_MT12, idents = c("Goblets_WT","Goblets_DKO","Tuft_WT","Tuft_DKO","TA_WT","TA_DKO","Cycling_TA_WT","Cycling_TA_DKO","Enterocytes_WT","Enterocytes_DKO","Paneth_WT","Paneth_DKO","EE_WT","EE_DKO","TumorCycling_DKOtVEH","TumorGoblet_DKOtVEH","FetalStem_DKOtVEH","RevStemCell_DKOtVEH","Stem_WT","Stem_DKO","SecretoryProgenitors_WT","SecretoryProgenitors_DKO"))

#####
###########
ElbowPlot(Epithelium_DKO_HARM_UMAP_MT12_FINAL)
Epithelium_DKO_HARM_UMAP_MT12_FINAL <- RunUMAP(Epithelium_DKO_HARM_UMAP_MT12_FINAL, reduction = "harmony",dims = 1:12, verbose = T)
Epithelium_DKO_HARM_UMAP_MT12_FINAL <- FindNeighbors(Epithelium_DKO_HARM_UMAP_MT12_FINAL,reduction = "harmony",dims = 1:12, verbose = T)
Epithelium_DKO_HARM_UMAP_MT12_FINAL <- FindClusters(Epithelium_DKO_HARM_UMAP_MT12_FINAL, verbose = T, resolution = 0.3)

Idents(Epithelium_DKO_HARM_UMAP_MT12_FINAL) <- "celltype"


DimPlot(Epithelium_DKO_HARM_UMAP_MT12_FINAL, label = T, split.by = "group")

my_levels <- c("WT","DKO","DKOtVEH")
Epithelium_DKO_HARM_UMAP_MT12_FINAL@meta.data$group <- factor(Epithelium_DKO_HARM_UMAP_MT12_FINAL@meta.data$group, levels = my_levels)

#####Correct this because they do not want to show tumor cluster because they don not understand the concept of mergins seurat objects

Idents(Epithelium_DKO_HARM_UMAP_MT12_FINAL)<-"celltype"
Stem<-WhichCells(Epithelium_DKO_HARM_UMAP_MT12_FINAL, ident= c("Stem"))
Epithelium_DKO_HARM_UMAP_MT12_FINAL@meta.data$celltypebroad [Stem]
Epithelium_DKO_HARM_UMAP_MT12_FINAL<-SetIdent(Epithelium_DKO_HARM_UMAP_MT12_FINAL, cells = Stem, value = "Stem")
Epithelium_DKO_HARM_UMAP_MT12_FINAL$celltypebroad<-Idents (Epithelium_DKO_HARM_UMAP_MT12_FINAL)
Tumor<-WhichCells (Epithelium_DKO_HARM_UMAP_MT12_FINAL, ident= c("TumorCycling","TumorGoblet","FetalStem","RevStemCell"))
Epithelium_DKO_HARM_UMAP_MT12_FINAL@meta.data$celltypebroad [Tumor]
Epithelium_DKO_HARM_UMAP_MT12_FINAL<-SetIdent(Epithelium_DKO_HARM_UMAP_MT12_FINAL, cells = Tumor, value = "Tumor")
Epithelium_DKO_HARM_UMAP_MT12_FINAL$celltypebroad<-Idents (Epithelium_DKO_HARM_UMAP_MT12_FINAL)
Goblets<-WhichCells (Epithelium_DKO_HARM_UMAP_MT12_FINAL, ident= c("Goblets"))
Epithelium_DKO_HARM_UMAP_MT12_FINAL@meta.data$celltypebroad [Goblets]
Epithelium_DKO_HARM_UMAP_MT12_FINAL<-SetIdent(Epithelium_DKO_HARM_UMAP_MT12_FINAL, cells = Goblets, value = "Goblets")
Epithelium_DKO_HARM_UMAP_MT12_FINAL$celltypebroad<-Idents (Epithelium_DKO_HARM_UMAP_MT12_FINAL)
TA<-WhichCells (Epithelium_DKO_HARM_UMAP_MT12_FINAL, ident= c("TA"))
Epithelium_DKO_HARM_UMAP_MT12_FINAL@meta.data$celltypebroad [TA]
Epithelium_DKO_HARM_UMAP_MT12_FINAL<-SetIdent(Epithelium_DKO_HARM_UMAP_MT12_FINAL, cells = TA, value = "TA")
Epithelium_DKO_HARM_UMAP_MT12_FINAL$celltypebroad<-Idents (Epithelium_DKO_HARM_UMAP_MT12_FINAL)
Enterocytes<-WhichCells (Epithelium_DKO_HARM_UMAP_MT12_FINAL, ident= c("Enterocytes"))
Epithelium_DKO_HARM_UMAP_MT12_FINAL@meta.data$celltypebroad [Enterocytes]
Epithelium_DKO_HARM_UMAP_MT12_FINAL<-SetIdent(Epithelium_DKO_HARM_UMAP_MT12_FINAL, cells = Enterocytes, value = "Enterocytes")
Epithelium_DKO_HARM_UMAP_MT12_FINAL$celltypebroad<-Idents (Epithelium_DKO_HARM_UMAP_MT12_FINAL)
EE<-WhichCells (Epithelium_DKO_HARM_UMAP_MT12_FINAL, ident= c("EE"))
Epithelium_DKO_HARM_UMAP_MT12_FINAL@meta.data$celltypebroad [EE]
Epithelium_DKO_HARM_UMAP_MT12_FINAL<-SetIdent(Epithelium_DKO_HARM_UMAP_MT12_FINAL, cells = EE, value = "EE")
Epithelium_DKO_HARM_UMAP_MT12_FINAL$celltypebroad<-Idents (Epithelium_DKO_HARM_UMAP_MT12_FINAL)
Tuft<-WhichCells (Epithelium_DKO_HARM_UMAP_MT12_FINAL, ident= c("Tuft"))
Epithelium_DKO_HARM_UMAP_MT12_FINAL@meta.data$celltypebroad [Tuft]
Epithelium_DKO_HARM_UMAP_MT12_FINAL<-SetIdent(Epithelium_DKO_HARM_UMAP_MT12_FINAL, cells = Tuft, value = "Tuft")
Epithelium_DKO_HARM_UMAP_MT12_FINAL$celltypebroad<-Idents (Epithelium_DKO_HARM_UMAP_MT12_FINAL)
Paneth<-WhichCells (Epithelium_DKO_HARM_UMAP_MT12_FINAL, ident= c("Paneth"))
Epithelium_DKO_HARM_UMAP_MT12_FINAL@meta.data$celltypebroad [Paneth]
Epithelium_DKO_HARM_UMAP_MT12_FINAL<-SetIdent(Epithelium_DKO_HARM_UMAP_MT12_FINAL, cells = Paneth, value = "Paneth")
Epithelium_DKO_HARM_UMAP_MT12_FINAL$celltypebroad<-Idents (Epithelium_DKO_HARM_UMAP_MT12_FINAL)
Cycling_TA<-WhichCells (Epithelium_DKO_HARM_UMAP_MT12_FINAL, ident= c("Cycling_TA"))
Epithelium_DKO_HARM_UMAP_MT12_FINAL@meta.data$celltypebroad [Cycling_TA]
Epithelium_DKO_HARM_UMAP_MT12_FINAL<-SetIdent(Epithelium_DKO_HARM_UMAP_MT12_FINAL, cells = Cycling_TA, value = "Cycling_TA")
Epithelium_DKO_HARM_UMAP_MT12_FINAL$celltypebroad<-Idents (Epithelium_DKO_HARM_UMAP_MT12_FINAL)
SecretoryProgenitors<-WhichCells (Epithelium_DKO_HARM_UMAP_MT12_FINAL, ident= c("SecretoryProgenitors"))
Epithelium_DKO_HARM_UMAP_MT12_FINAL@meta.data$celltypebroad [SecretoryProgenitors]
Epithelium_DKO_HARM_UMAP_MT12_FINAL<-SetIdent(Epithelium_DKO_HARM_UMAP_MT12_FINAL, cells = SecretoryProgenitors, value = "SecretoryProgenitors")
Epithelium_DKO_HARM_UMAP_MT12_FINAL$celltypebroad<-Idents (Epithelium_DKO_HARM_UMAP_MT12_FINAL)

Idents(Epithelium_DKO_HARM_UMAP_MT12_FINAL) <- "celltypebroad"

DimPlot(Epithelium_DKO_HARM_UMAP_MT12_FINAL,label=F,cols=c("pink","Khaki2","hotpink","turquoise","skyblue3","tomato2","darkgoldenrod2","brown","indianred4","palegreen4"), split.by = "group")
DimPlot(Epithelium_DKO_HARM_UMAP_MT12_FINAL,label=F,cols=c("pink","Khaki2","hotpink","turquoise","skyblue3","tomato2","darkgoldenrod2","brown","indianred4","palegreen4"))


DimPlot(Epithelium_DKO_HARM_UMAP_MT12_FINAL,label=F,cols=c("chocolate4","darkorchid4","olivedrab3","darksalmon","tomato2","turquoise","palegreen4","brown","Khaki2","hotpink","pink","darkgoldenrod2","skyblue3"), split.by = "group")
DimPlot(Epithelium_DKO_HARM_UMAP_MT12_FINAL,label=F,cols=c("palegreen4","Khaki2","darkgoldenrod2","pink","tomato2","brown","skyblue3","turquoise","hotpink","olivedrab3","darksalmon","chocolate4","darkorchid4"))

saveRDS(Epithelium_DKO_HARM_UMAP_MT12_FINAL,file="Results/WT_DKO_Tumor_Dec22")


Idents(Epithelium_DKO_HARM_UMAP_MT12_FINAL) <- "group"
FeaturePlot(Epithelium_DKO_HARM_UMAP_MT12_FINAL, "Lgr5",cols=c("lightgrey","red","darkred"),split.by = "group", min.cutoff = 0.35)
FeaturePlot(Epithelium_DKO_HARM_UMAP_MT12_FINAL, "Olfm4",cols=c("lightgrey","red","darkred"),split.by = "group", min.cutoff = 1.5)
FeaturePlot(Epithelium_DKO_HARM_UMAP_MT12_FINAL, "Olfm4",cols=c("lightgrey","red","darkred"),split.by = "group", min.cutoff = 1.5)
FeaturePlot(Epithelium_DKO_HARM_UMAP_MT12_FINAL, "Clu",cols=c("lightgrey","red","darkred"),split.by = "group", min.cutoff = 0)
FeaturePlot(Epithelium_DKO_HARM_UMAP_MT12_FINAL, "Anxa10",cols=c("lightgrey","red","darkred"),split.by = "group", min.cutoff = 0)

###Percentage positive cells ######
######Percentage positive cells for each cell type#############
####Set ident Olfm4 positive celltype group####################
Idents(Epithelium_DKO_HARM_UMAP_MT12_FINAL)<-"seurat_clusters"

Lgr5pos<-WhichCells(Epithelium_DKO_HARM_UMAP_MT12_FINAL,slot = "counts", expression = Lgr5>0.5)
Epithelium_DKO_HARM_UMAP_MT12_FINAL@meta.data$Lgr5_status<-"Lgr5pos"
Epithelium_DKO_HARM_UMAP_MT12_FINAL<-SetIdent(Epithelium_DKO_HARM_UMAP_MT12_FINAL, cells = Lgr5pos, value ="Lgr5pos" )
Epithelium_DKO_HARM_UMAP_MT12_FINAL$Lgr5_status<-Idents(Epithelium_DKO_HARM_UMAP_MT12_FINAL)
####Lgr5 neg
Lgr5neg<-WhichCells(Epithelium_DKO_HARM_UMAP_MT12_FINAL,slot = "counts", expression = Lgr5<0.5)
Epithelium_DKO_HARM_UMAP_MT12_FINAL@meta.data$Lgr5_status<-"Lgr5neg"
Epithelium_DKO_HARM_UMAP_MT12_FINAL<-SetIdent(Epithelium_DKO_HARM_UMAP_MT12_FINAL, cells = Lgr5neg, value ="Lgr5neg" )
Epithelium_DKO_HARM_UMAP_MT12_FINAL$Lgr5_status<-Idents(Epithelium_DKO_HARM_UMAP_MT12_FINAL)

######
Idents(Epithelium_DKO_HARM_UMAP_MT12_FINAL)<-"Lgr5_status"
DimPlot(Epithelium_DKO_HARM_UMAP_MT12_FINAL)

######Percentage positive cells for each cell type#############
####Set ident Lgr5 positive celltype group####################
Idents(Epithelium_DKO_HARM_UMAP_MT12_FINAL)<-"seurat_clusters"

Olfm4pos<-WhichCells(Epithelium_DKO_HARM_UMAP_MT12_FINAL,slot = "counts", expression = Olfm4>0.5)
Epithelium_DKO_HARM_UMAP_MT12_FINAL@meta.data$Olfm4_status<-"Olfm4pos"
Epithelium_DKO_HARM_UMAP_MT12_FINAL<-SetIdent(Epithelium_DKO_HARM_UMAP_MT12_FINAL, cells = Olfm4pos, value ="Olfm4pos" )
Epithelium_DKO_HARM_UMAP_MT12_FINAL$Olfm4_status<-Idents(Epithelium_DKO_HARM_UMAP_MT12_FINAL)
####Olfm4 neg
Olfm4neg<-WhichCells(Epithelium_DKO_HARM_UMAP_MT12_FINAL,slot = "counts", expression = Olfm4<0.5)
Epithelium_DKO_HARM_UMAP_MT12_FINAL@meta.data$Olfm4_status<-"Olfm4neg"
Epithelium_DKO_HARM_UMAP_MT12_FINAL<-SetIdent(Epithelium_DKO_HARM_UMAP_MT12_FINAL, cells = Olfm4neg, value ="Olfm4neg" )
Epithelium_DKO_HARM_UMAP_MT12_FINAL$Olfm4_status<-Idents(Epithelium_DKO_HARM_UMAP_MT12_FINAL)

################Set ident sample.cell_type plus gene positive####
#Add cell_type and cluster info to meta.data
Idents(Epithelium_DKO_HARM_UMAP_MT12_FINAL) <- "celltype.group"
head(Epithelium_DKO_HARM_UMAP_MT12_FINAL@meta.data,2)
Epithelium_DKO_HARM_UMAP_MT12_FINAL$Olfm4_status.celltype.group <- paste(Idents(Epithelium_DKO_HARM_UMAP_MT12_FINAL), Epithelium_DKO_HARM_UMAP_MT12_FINAL$Olfm4_status, sep = "_")
head(Epithelium_DKO_HARM_UMAP_MT12_FINAL@meta.data,3)

###############
################Set ident sample.cell_type plus gene positive####
#Add cell_type and cluster info to meta.data
Idents(Epithelium_DKO_HARM_UMAP_MT12_FINAL) <- "celltype.group"
head(Epithelium_DKO_HARM_UMAP_MT12_FINAL@meta.data,2)
Epithelium_DKO_HARM_UMAP_MT12_FINAL$Lgr5_status.celltype.group <- paste(Idents(Epithelium_DKO_HARM_UMAP_MT12_FINAL), Epithelium_DKO_HARM_UMAP_MT12_FINAL$Lgr5_status, sep = "_")
head(Epithelium_DKO_HARM_UMAP_MT12_FINAL@meta.data,3)

########

Idents(Epithelium_DKO_HARM_UMAP_MT12_FINAL) <- "Lgr5_status.celltype.group"
counts.by.CellType_LGR5.FIG <-table(Idents(Epithelium_DKO_HARM_UMAP_MT12_FINAL))
write.table(as.matrix(counts.by.CellType_LGR5.FIG),"./Results/Counts/counts.by.CellType_LGR5.txt",sep="\t",col.names=T,row.names=T)

Idents(Epithelium_DKO_HARM_UMAP_MT12_FINAL) <- "Olfm4_status.celltype.group"
counts.by.CellType_OLFM4.FIG <-table(Idents(Epithelium_DKO_HARM_UMAP_MT12_FINAL))
write.table(as.matrix(counts.by.CellType_OLFM4.FIG),"./Results/Counts/counts.by.CellType_OLFM4.txt",sep="\t",col.names=T,row.names=T)

###To use later #####
####Subset for Dotplot
Idents(Epithelium_DKO_HARM_UMAP_MT12_FINAL) <- "celltype"

Progenitors_DKO_HARM_UMAP_MT12_FINAL <- subset(Epithelium_DKO_HARM_UMAP_MT12_FINAL, idents = c("Stem","TumorCycling","TumorGoblet","FetalStem","RevStemCell"))

Idents(Progenitors_DKO_HARM_UMAP_MT12_FINAL)<-"celltype"
Stem<-WhichCells(Progenitors_DKO_HARM_UMAP_MT12_FINAL, ident= c("Stem"))
Progenitors_DKO_HARM_UMAP_MT12_FINAL@meta.data$celltypebroad [Stem]
Progenitors_DKO_HARM_UMAP_MT12_FINAL<-SetIdent(Progenitors_DKO_HARM_UMAP_MT12_FINAL, cells = Stem, value = "Stem")
Progenitors_DKO_HARM_UMAP_MT12_FINAL$celltypebroad<-Idents (Progenitors_DKO_HARM_UMAP_MT12_FINAL)
Tumor<-WhichCells (Progenitors_DKO_HARM_UMAP_MT12_FINAL, ident= c("TumorCycling","TumorGoblet","FetalStem","RevStemCell"))
Progenitors_DKO_HARM_UMAP_MT12_FINAL@meta.data$celltypebroad [Tumor]
Progenitors_DKO_HARM_UMAP_MT12_FINAL<-SetIdent(Progenitors_DKO_HARM_UMAP_MT12_FINAL, cells = Tumor, value = "Tumor")
Progenitors_DKO_HARM_UMAP_MT12_FINAL$celltypebroad<-Idents (Progenitors_DKO_HARM_UMAP_MT12_FINAL)

#####
Idents(Progenitors_DKO_HARM_UMAP_MT12_FINAL) <- "celltypebroad"
head(Progenitors_DKO_HARM_UMAP_MT12_FINAL@meta.data,2)
Progenitors_DKO_HARM_UMAP_MT12_FINAL$celltypebroad.group <- paste(Idents(Progenitors_DKO_HARM_UMAP_MT12_FINAL), Progenitors_DKO_HARM_UMAP_MT12_FINAL$group, sep = "_")
head(Progenitors_DKO_HARM_UMAP_MT12_FINAL@meta.data,3)


####
Idents(Progenitors_DKO_HARM_UMAP_MT12_FINAL) <- "celltypebroad.group"

########
DotPlot(Progenitors_DKO_HARM_UMAP_MT12_FINAL, features = c("Olfm4","Lgr5","Ascl2","Axin2","Gkn3"), dot.scale = 8, 
        col.min = 0,
        col.max = 1,
) + scale_colour_gradient2(low = "lightgrey", mid = "red", high = "darkred",midpoint = 0.5)+RotatedAxis()

#####
DotPlot(Progenitors_DKO_HARM_UMAP_MT12_FINAL, features = c("Munoz_1","RegevStemTF_1","RegevStem_1","ISCI_1","ISCII_1","ISCIII_1"), dot.scale = 8, 
        col.min = 0,
        col.max = 1,
) + scale_colour_gradient2(low = "lightgrey", mid = "red", high = "darkred",midpoint =0.5)+RotatedAxis()


#########

Idents(Progenitors_DKO_HARM_UMAP_MT12_FINAL)<-"celltypebroad.group"

my_levels <- c("Stem_WT","Stem_DKO","Tumor_DKOtVEH")
Progenitors_DKO_HARM_UMAP_MT12_FINAL@meta.data$celltypebroad.group <- factor(Progenitors_DKO_HARM_UMAP_MT12_FINAL@meta.data$celltypebroad.group, levels = my_levels)

VlnPlot(Progenitors_DKO_HARM_UMAP_MT12_FINAL,cols = c("bisque","indianred2","indianred4"),"Lgr5", pt.size = 0.0,idents = c("Stem_WT","Stem_DKO","Tumor_DKOtVEH"))+geom_signif(comparisons = list(c("Stem_WT","Stem_DKO")),map_signif_level = TRUE, textsize =3) + ylim(NA,4) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)
VlnPlot(Progenitors_DKO_HARM_UMAP_MT12_FINAL,cols = c("bisque","indianred2","indianred4"),"Olfm4", pt.size = 0.0,idents = c("Stem_WT","Stem_DKO","Tumor_DKOtVEH"))+geom_signif(comparisons = list(c("Stem_WT","Stem_DKO")),map_signif_level = TRUE, textsize =3) + ylim(NA,6) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)
VlnPlot(Progenitors_DKO_HARM_UMAP_MT12_FINAL,cols = c("bisque","indianred2","indianred4"),"CBC_1", pt.size = 0.0,idents = c("Stem_WT","Stem_DKO","Tumor_DKOtVEH"))+geom_signif(comparisons = list(c("Stem_WT","Stem_DKO")),map_signif_level = TRUE, textsize =3) + ylim(NA,3) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)
VlnPlot(Progenitors_DKO_HARM_UMAP_MT12_FINAL,cols = c("bisque","indianred2","indianred4"),"Metapalasia_1", pt.size = 0.0,idents = c("Stem_WT","Stem_DKO","Tumor_DKOtVEH"))+geom_signif(comparisons = list(c("Stem_WT","Stem_DKO")),map_signif_level = TRUE, textsize =3) + ylim(NA,1) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)
VlnPlot(Progenitors_DKO_HARM_UMAP_MT12_FINAL,cols = c("bisque","indianred2","indianred4"),"RSC_1", pt.size = 0.0,idents = c("Stem_WT","Stem_DKO","Tumor_DKOtVEH"))+geom_signif(comparisons = list(c("Stem_WT","Stem_DKO")),map_signif_level = TRUE, textsize =3) + ylim(NA,3) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)
VlnPlot(Progenitors_DKO_HARM_UMAP_MT12_FINAL,cols = c("bisque","indianred2","indianred4"),"Fetal_1", pt.size = 0.0,idents = c("Stem_WT","Stem_DKO","Tumor_DKOtVEH"))+geom_signif(comparisons = list(c("Stem_WT","Stem_DKO")),map_signif_level = TRUE, textsize =3) + ylim(NA,0.5) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)
VlnPlot(Progenitors_DKO_HARM_UMAP_MT12_FINAL,cols = c("bisque","indianred2","indianred4"),"Munoz_1", pt.size = 0.0,idents = c("Stem_WT","Stem_DKO","Tumor_DKOtVEH"))+geom_signif(comparisons = list(c("Stem_WT","Stem_DKO")),map_signif_level = TRUE, textsize =3) + ylim(NA,0.5) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)
VlnPlot(Progenitors_DKO_HARM_UMAP_MT12_FINAL,cols = c("bisque","indianred2","indianred4"),"RegevStemTF_1", pt.size = 0.0,idents = c("Stem_WT","Stem_DKO","Tumor_DKOtVEH"))+geom_signif(comparisons = list(c("Stem_WT","Stem_DKO")),map_signif_level = TRUE, textsize =3) + ylim(NA,0.5) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)
VlnPlot(Progenitors_DKO_HARM_UMAP_MT12_FINAL,cols = c("bisque","indianred2","indianred4"),"RegevStem_1", pt.size = 0.0,idents = c("Stem_WT","Stem_DKO","Tumor_DKOtVEH"))+geom_signif(comparisons = list(c("Stem_WT","Stem_DKO")),map_signif_level = TRUE, textsize =3) + ylim(NA,0.75) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)
VlnPlot(Progenitors_DKO_HARM_UMAP_MT12_FINAL,cols = c("bisque","indianred2","indianred4"),"ISCI_1", pt.size = 0.0,idents = c("Stem_WT","Stem_DKO","Tumor_DKOtVEH"))+geom_signif(comparisons = list(c("Stem_WT","Stem_DKO")),map_signif_level = TRUE, textsize =3) + ylim(NA,1) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)
VlnPlot(Progenitors_DKO_HARM_UMAP_MT12_FINAL,cols = c("bisque","indianred2","indianred4"),"ISCII_1", pt.size = 0.0,idents = c("Stem_WT","Stem_DKO","Tumor_DKOtVEH"))+geom_signif(comparisons = list(c("Stem_WT","Stem_DKO")),map_signif_level = TRUE, textsize =3) + ylim(NA,2) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)
VlnPlot(Progenitors_DKO_HARM_UMAP_MT12_FINAL,cols = c("bisque","indianred2","indianred4"),"ISCIII_1", pt.size = 0.0,idents = c("Stem_WT","Stem_DKO","Tumor_DKOtVEH"))+geom_signif(comparisons = list(c("Stem_WT","Stem_DKO")),map_signif_level = TRUE, textsize =3) + ylim(NA,0.5) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)



#####Counts #######
Idents(Epithelium_DKO_HARM_UMAP_MT12_FINAL) <- "celltype.group"
counts.by.CellType.group <-table(Idents(Epithelium_DKO_HARM_UMAP_MT12_FINAL))
write.table(as.matrix(counts.by.CellType.group),"./Results/Counts/counts.by.CellType.group.txt",sep="\t",col.names=T,row.names=T)


CBC <- c("Olfm4","Lgr5","Ascl2","Axin2","Gkn3")
Epithelium_DKO_HARM_UMAP_MT12_FINAL <- AddModuleScore(Epithelium_DKO_HARM_UMAP_MT12_FINAL, features = list(CBC), name="CBC_")
Idents(Epithelium_DKO_HARM_UMAP_MT12_FINAL)<-"seurat_clusters"
FeaturePlot(Epithelium_DKO_HARM_UMAP_MT12_FINAL, "CBC_1", pt.size = -0.1, min.cutoff = 0, cols = c("lightgrey","red","darkred"), label = F, split.by = "group")
FeaturePlot(Epithelium_DKO_HARM_UMAP_MT12_FINAL, "CBC_1", pt.size = -0.1, min.cutoff = 0, cols = c("lightgrey","red","darkred"), label = F)


Metapalasia <- c("Tff2","Muc5ac","Anxa1","Ahr","Aqp5","Pdx1","Il18","Fscn1","Relb","Anxa10","Tacstd2","Msln","Muc3","Gsdmd","Mdk","Reg4")
Epithelium_DKO_HARM_UMAP_MT12_FINAL <- AddModuleScore(Epithelium_DKO_HARM_UMAP_MT12_FINAL, features = list(Metapalasia), name="Metapalasia_")
Idents(Epithelium_DKO_HARM_UMAP_MT12_FINAL)<-"seurat_clusters"
FeaturePlot(Epithelium_DKO_HARM_UMAP_MT12_FINAL,"Metapalasia_1", pt.size = 0.1, min.cutoff = 0, cols = c("lightgrey","red","darkred"), label = F, split.by = "group")

Fetal <- read_xlsx("./Results/Signatures//FetalMustata.xlsx")
Epithelium_DKO_HARM_UMAP_MT12_FINAL <- AddModuleScore(Epithelium_DKO_HARM_UMAP_MT12_FINAL, features = Fetal, name = "Fetal_")
FeaturePlot(Epithelium_DKO_HARM_UMAP_MT12_FINAL,"Fetal_1", pt.size = 0.1, min.cutoff = 0, cols = c("lightgrey","red","darkred"), label = F, split.by = "group")
FeaturePlot(Epithelium_DKO_HARM_UMAP_MT12_FINAL,"Fetal_1", pt.size = 0.1, min.cutoff = 0, cols = c("lightgrey","red","darkred"), label = F)


ISCI <- read.table("./Results/Signatures/ISC1_signature.txt", header = TRUE)
Epithelium_DKO_HARM_UMAP_MT12_FINAL <- AddModuleScore(Epithelium_DKO_HARM_UMAP_MT12_FINAL, features = ISCI, name = "ISCI_")
FeaturePlot(Epithelium_DKO_HARM_UMAP_MT12_FINAL, "ISCI_1", pt.size = 0, min.cutoff = 0,cols = c("lightgrey","red","darkred"), label = T, split.by = "group")

ISCII <- read.table("./Results/Signatures/ISC2_signature.txt", header = TRUE)
Epithelium_DKO_HARM_UMAP_MT12_FINAL <- AddModuleScore(Epithelium_DKO_HARM_UMAP_MT12_FINAL, features = ISCII, name = "ISCII_")
FeaturePlot(Epithelium_DKO_HARM_UMAP_MT12_FINAL, "ISCII_1", pt.size = 0, min.cutoff = 0,cols = c("lightgrey","red","darkred"), label = T, split.by = "group")

ISCIII <- read.table("./Results/Signatures/ISC3_signature.txt", header = TRUE)
Epithelium_DKO_HARM_UMAP_MT12_FINAL <- AddModuleScore(Epithelium_DKO_HARM_UMAP_MT12_FINAL, features = ISCIII, name = "ISCIII_")
FeaturePlot(Epithelium_DKO_HARM_UMAP_MT12_FINAL, "ISCIII_1", pt.size = 0, min.cutoff = 0,cols = c("lightgrey","red","darkred"), label = T, split.by = "group")

RegevStemTF <- read.table("./Results/Signatures/Regev_Lgr5.txt", header = TRUE)
Epithelium_DKO_HARM_UMAP_MT12_FINAL <- AddModuleScore(Epithelium_DKO_HARM_UMAP_MT12_FINAL, features = RegevStemTF, name = "RegevStemTF_")
FeaturePlot(Epithelium_DKO_HARM_UMAP_MT12_FINAL, "RegevStemTF_1", pt.size = 0, min.cutoff = 0,cols = c("lightgrey","red","darkred"), label = T, split.by = "group")

Munoz <- read.table("./Results/Signatures/Munoz_Lgr5.txt", header = TRUE)
Epithelium_DKO_HARM_UMAP_MT12_FINAL <- AddModuleScore(Epithelium_DKO_HARM_UMAP_MT12_FINAL, features = Munoz, name = "Munoz_")
FeaturePlot(Epithelium_DKO_HARM_UMAP_MT12_FINAL, "Munoz_1", pt.size = 0, min.cutoff = 0,cols = c("lightgrey","red","darkred"), label = T, split.by = "group")

KimLgr5 <- read.table("./Results/Signatures/Kim_Lgr5.txt", header = TRUE)
Epithelium_DKO_HARM_UMAP_MT12_FINAL <- AddModuleScore(Epithelium_DKO_HARM_UMAP_MT12_FINAL, features = KimLgr5, name = "KimLgr5_")
FeaturePlot(Epithelium_DKO_HARM_UMAP_MT12_FINAL, "KimLgr5_1", pt.size = 0, min.cutoff = 0,cols = c("lightgrey","red","darkred"), label = T, split.by = "group")

RSC <- c("Basp1","Clu","Cd44","Cxadr")
Epithelium_DKO_HARM_UMAP_MT12_FINAL <- AddModuleScore(Epithelium_DKO_HARM_UMAP_MT12_FINAL, features = list(RSC), name="RSC_")
FeaturePlot(Epithelium_DKO_HARM_UMAP_MT12_FINAL,"RSC_1", pt.size = 0.1, min.cutoff = 0, cols = c("lightgrey","red","darkred"), label = F, split.by = "group")
FeaturePlot(Epithelium_DKO_HARM_UMAP_MT12_FINAL,"RSC_1", pt.size = 0.1, min.cutoff = 0, cols = c("lightgrey","red","darkred"), label = F)

RegevStem <- read.delim("./Results/Signatures/Stem_REgev_Survey.txt", header = TRUE)
Progenitors_DKO_HARM_UMAP_MT12_FINAL <- AddModuleScore(Progenitors_DKO_HARM_UMAP_MT12_FINAL, features = RegevStem, name = "RegevStem_")
FeaturePlot(Progenitors_DKO_HARM_UMAP_MT12_FINAL, "RegevStem_1", pt.size = 0, min.cutoff = 0,cols = c("lightgrey","red","darkred"), label = T, split.by = "group")


#####Features################
######
TA <- c("Dmbt1", "Cdk2ap2", "Rbp7", "S100a16","Gjb3","Fgfbp1","Nrarp","Pycard")
Epithelium_DKO_HARM_UMAP_MT12_FINAL <- AddModuleScore(Epithelium_DKO_HARM_UMAP_MT12_FINAL, features = list(TA), name="TA_")
FeaturePlot(Epithelium_DKO_HARM_UMAP_MT12_FINAL, "TA_1", pt.size = 0.1, min.cutoff = -0.25, cols = c("lightgrey","red","darkred"), label = F)
######
Cycling <- c("Top2a", "Mki67", "Pcna", "Birc5","Pcna","Hmgb2")
Epithelium_DKO_HARM_UMAP_MT12_FINAL <- AddModuleScore(Epithelium_DKO_HARM_UMAP_MT12_FINAL, features = list(Cycling), name="Cycling_")
FeaturePlot(Epithelium_DKO_HARM_UMAP_MT12_FINAL, "Cycling_1", pt.size = 0.1, min.cutoff = -0.5, cols = c("lightgrey","red","darkred"), label = F)

#######
Enterocytes <- c("Krt20", "Vil1", "Alpi", "Fabp1","Apoa1","Apoa4","Reg3a","Reg3g")
Epithelium_DKO_HARM_UMAP_MT12_FINAL <- AddModuleScore(Epithelium_DKO_HARM_UMAP_MT12_FINAL, features = list(Enterocytes), name="Enterocytes_")
FeaturePlot(Epithelium_DKO_HARM_UMAP_MT12_FINAL, "Enterocytes_1", pt.size = 0.1, min.cutoff = 0, cols = c("lightgrey","red","darkred"), label = F)

#######
Paneth <- c("Lyz1", "Defa17", "Defa22", "Defa24","Ang4")
Epithelium_DKO_HARM_UMAP_MT12_FINAL <- AddModuleScore(Epithelium_DKO_HARM_UMAP_MT12_FINAL, features = list(Paneth), name="Paneth_")
FeaturePlot(Epithelium_DKO_HARM_UMAP_MT12_FINAL, "Paneth_1", pt.size = 0.1, min.cutoff = 0, cols = c("lightgrey","red","darkred"), label = F)
#######
Goblets <- c("Agr2", "Muc2", "Spink4", "Tff3")
Epithelium_DKO_HARM_UMAP_MT12_FINAL <- AddModuleScore(Epithelium_DKO_HARM_UMAP_MT12_FINAL, features = list(Goblets), name="Goblets_")
FeaturePlot(Epithelium_DKO_HARM_UMAP_MT12_FINAL, "Goblets_1", pt.size = 0.1, min.cutoff = 1, cols = c("lightgrey","red","darkred"), label = F)
EE <- c("Chga", "Chgb", "Tak1", "Tph1")
Epithelium_DKO_HARM_UMAP_MT12_FINAL <- AddModuleScore(Epithelium_DKO_HARM_UMAP_MT12_FINAL, features = list(EE), name="EE_")
FeaturePlot(Epithelium_DKO_HARM_UMAP_MT12_FINAL, "EE_1", pt.size = 0.1, min.cutoff = -0.1, cols = c("lightgrey","red","darkred"), label = F)
Tuft <- c("Dclk1", "Avil", "Gfi1b", "Trpm5")
Epithelium_DKO_HARM_UMAP_MT12_FINAL <- AddModuleScore(Epithelium_DKO_HARM_UMAP_MT12_FINAL, features = list(Tuft), name="Tuft_")
FeaturePlot(Epithelium_DKO_HARM_UMAP_MT12_FINAL, "Tuft_1", pt.size = 0.1, min.cutoff = 0, cols = c("lightgrey","red","darkred"), label = F)

#####
######
Idents(Epithelium_DKO_HARM_UMAP_MT12_FINAL)<-"group"
DimPlot(Epithelium_DKO_HARM_UMAP_MT12_FINAL, cols=c("bisque","indianred2","indianred4"))

Idents(Epithelium_DKO_HARM_UMAP_MT12_FINAL)<-"celltype"


DotPlot(Epithelium_DKO_HARM_UMAP_MT12_FINAL, features = c("Lgr5","Olfm4","Birc5","Top2a","Dmbt1","Pycard","Vil1","Alpi","Muc2","Zg16","Neurog3","Chga","Dclk1","Avil","Lyz1","Defa24","Clu","Cd44","Anxa10","Ly6a"), dot.scale = 8, 
        col.min = 0,
        col.max = 3,
) + scale_colour_gradient2(low = "lightgrey", mid = "red", high = "darkred",midpoint = 1.5)+RotatedAxis() 


#####My_levels#####
Idents(Epithelium_DKO_HARM_UMAP_MT12_FINAL)<-"celltype"

my_levels <- c("TumorGoblet","FetalStem","RevStemCell","TumorCycling","Paneth","Tuft","EE","Goblets","Enterocytes","SecretoryProgenitors","TA","Cycling_TA","Stem")
Epithelium_DKO_HARM_UMAP_MT12_FINAL@meta.data$celltype <- factor(Epithelium_DKO_HARM_UMAP_MT12_FINAL@meta.data$celltype, levels = my_levels)


######Epithelium####

FeaturePlot(Epithelium_DKO_HARM_UMAP_MT12_FINAL, "Dcn", cols = c("Lightgrey","red","darkred"), min.cutoff = 0)

############################################
##########Biomarkers################
####BIOMARKERS#####
Idents(Epithelium_DKO_HARM_UMAP_MT12_FINAL) <- "celltype"
Epithelium_DKO_HARM_UMAP_MT12_FINAL.biomarkers <- FindAllMarkers(Epithelium_DKO_HARM_UMAP_MT12_FINAL, only.pos = FALSE, min.pct = 0.2, logfc.threshold = 0.3)
write.table(as.matrix(Epithelium_DKO_HARM_UMAP_MT12_FINAL.biomarkers),"./Results/Biomarkers/Epithelium_DKO_HARM_UMAP_MT12_FINAL.biomarkers.txt")

top5_Tumor <- Epithelium_DKO_HARM_UMAP_MT12_FINAL.biomarkers %>% group_by(cluster) %>% top_n(5, wt=avg_log2FC)


#####My_levels#####
Idents(Epithelium_DKO_HARM_UMAP_MT12_FINAL)<-"celltype"

my_levels <- c("Stem","Cycling_TA","TA","SecretoryProgenitors","Enterocytes","Goblets","EE","Tuft","Paneth","TumorCycling","RevStemCell","FetalStem","TumorGoblet")
Epithelium_DKO_HARM_UMAP_MT12_FINAL@meta.data$celltype <- factor(Epithelium_DKO_HARM_UMAP_MT12_FINAL@meta.data$celltype, levels = my_levels)

#####My_levels#####
Idents(Epithelium_DKO_HARM_UMAP_MT12_FINAL)<-"celltypebroad"

my_levels <- c("Stem","Cycling_TA","TA","SecretoryProgenitors","Enterocytes","Goblets","EE","Tuft","Paneth","Tumor")
Epithelium_DKO_HARM_UMAP_MT12_FINAL@meta.data$celltypebroad <- factor(Epithelium_DKO_HARM_UMAP_MT12_FINAL@meta.data$celltypebroad, levels = my_levels)

#####My_levels#####
Idents(Epithelium_DKO_HARM_UMAP_MT12_FINAL)<-"celltypebroad"

my_levels <- c("Tumor","Paneth","Tuft","EE","Goblets","Enterocytes","SecretoryProgenitors","TA","Cycling_TA","Stem")
Epithelium_DKO_HARM_UMAP_MT12_FINAL@meta.data$celltypebroad <- factor(Epithelium_DKO_HARM_UMAP_MT12_FINAL@meta.data$celltypebroad, levels = my_levels)

####BIOMARKERS#####
Idents(Epithelium_DKO_HARM_UMAP_MT12_FINAL) <- "celltypebroad"
BRaodEpithelium_DKO_HARM_UMAP_MT12_FINAL.biomarkers <- FindAllMarkers(Epithelium_DKO_HARM_UMAP_MT12_FINAL, only.pos = FALSE, min.pct = 0.2, logfc.threshold = 0.3)
write.table(as.matrix(BRaodEpithelium_DKO_HARM_UMAP_MT12_FINAL.biomarkers),"./Results/Biomarkers/BRaodEpithelium_DKO_HARM_UMAP_MT12_FINAL.biomarkers.txt")

top5_Tumor <- BRaodEpithelium_DKO_HARM_UMAP_MT12_FINAL.biomarkers %>% group_by(cluster) %>% top_n(5, wt=avg_log2FC)

DoHeatmap(Epithelium_DKO_HARM_UMAP_MT12_FINAL, features =top5_Tumor$gene, angle = 0,draw.lines=T) + scale_fill_gradientn(colors = c("lightblue", "white", "darkred")) 

saveRDS(Epithelium_DKO_HARM_UMAP_MT12_FINAL,file="Results/WT_DKO_DKOt_dec22_HARM_UMAP")


####Now work on WT-ZKO-LKO and DKO object #####
###Set same parameters

Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_QC <-readRDS(file = "./Results/WT_ZKO_LKO_DKO_HARM_UMAP")


Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12 <- subset(x = Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_QC, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 12 )

DimPlot(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12)

Idents(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12)<-"celltype"


DimPlot(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12,label=F,cols=c("chocolate4","darkorchid4","olivedrab3","darksalmon","tomato2","turquoise","palegreen4","brown","Khaki2","hotpink","pink","darkgoldenrod2","skyblue3"), split.by = "group")
DimPlot(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12,label=F,cols=c("pink","skyblue3","palegreen4","hotpink","turquoise","darkgoldenrod2","tomato2","brown","Khaki2"))


Idents(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12)<-"group"


VlnPlot(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12,cols = c("bisque","skyblue3","aquamarine3","indianred2"),"Cdx2", pt.size = 0.0,idents = c("WT","ZKO","LKO", "DKO"))+geom_signif(comparisons = list(c("WT","ZKO")),map_signif_level = TRUE, textsize =3) + ylim(NA,4) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)


#Counts#
#Cell.counts
head(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12@meta.data,2)
#Add group_new and cluster info to meta.data
Idents(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12) <- "celltype"
Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12$celltype.group <- paste(Idents(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12), Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12$group, sep = "_")
head(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12@meta.data,3)

Idents(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12) <- "celltype.group"
counts.Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12 <-table(Idents(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12))
write.table(as.matrix(counts.Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12),"./Results/Counts/counts.Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12.txt",sep="\t",col.names=T,row.names=T)



##############
########
Idents(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12) <- "celltype.group"

DotPlot(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, features = c("Olfm4","Lgr5","Ascl2","Axin2","Gkn3"), dot.scale = 8, 
        col.min = 0,
        col.max = 1,
) + scale_colour_gradient2(low = "lightgrey", mid = "red", high = "darkred",midpoint = 0.5)+RotatedAxis()


Idents(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12) <- "celltype"

Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12 <- subset(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, idents = c("Stem"))


Idents(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12) <- "group"


VlnPlot(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12,cols = c("bisque","skyblue3","aquamarine3","indianred2"),"Olfm4", pt.size = 0.0,idents = c("WT","ZKO","LKO", "DKO"))+geom_signif(comparisons = list(c("WT","ZKO")),map_signif_level = TRUE, textsize =3) + ylim(NA,6) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)
VlnPlot(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12,cols = c("bisque","skyblue3","aquamarine3","indianred2"),"Lgr5", pt.size = 0.0,idents = c("WT","ZKO","LKO", "DKO"))+geom_signif(comparisons = list(c("WT","ZKO")),map_signif_level = TRUE, textsize =3) + ylim(NA,3) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)
VlnPlot(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12,cols = c("bisque","skyblue3","aquamarine3","indianred2"),"Axin2", pt.size = 0.0,idents = c("WT","ZKO","LKO", "DKO"))+geom_signif(comparisons = list(c("WT","DKO")),map_signif_level = TRUE, textsize =3) + ylim(NA,3) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)
VlnPlot(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12,cols = c("bisque","skyblue3","aquamarine3","indianred2"),"Ascl2", pt.size = 0.0,idents = c("WT","ZKO","LKO", "DKO"))+geom_signif(comparisons = list(c("WT","LKO")),map_signif_level = TRUE, textsize =3) + ylim(NA,3) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)
VlnPlot(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12,cols = c("bisque","skyblue3","aquamarine3","indianred2"),"Gkn3", pt.size = 0.0,idents = c("WT","ZKO","LKO", "DKO"))+geom_signif(comparisons = list(c("WT","ZKO")),map_signif_level = TRUE, textsize =4) + ylim(NA,4) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)
VlnPlot(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12,cols = c("bisque","skyblue3","aquamarine3","indianred2"),"Clu", pt.size = 0.0,idents = c("WT","ZKO","LKO", "DKO"))+geom_signif(comparisons = list(c("WT","ZKO")),map_signif_level = TRUE, textsize =4) + ylim(NA,4) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)

FeaturePlot(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, "Lgr5", pt.size = 0.1, min.cutoff = 0.2, cols = c("lightgrey","red","darkred"), label = F, split.by = "group")
FeaturePlot(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, "Olfm4", pt.size = 0.1, min.cutoff = 1.5, cols = c("lightgrey","red","darkred"), label = F, split.by = "group")

##add signatures###
ISCI <- read.table("./Results/Signatures/ISC1_signature.txt", header = TRUE)
Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12 <- AddModuleScore(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, features = ISCI, name = "ISCI_")
FeaturePlot(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, "ISCI_1", pt.size = 0, min.cutoff = 0,cols = c("lightgrey","red","darkred"), label = T, split.by = "group")

ISCII <- read.table("./Results/Signatures/ISC2_signature.txt", header = TRUE)
Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12 <- AddModuleScore(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, features = ISCII, name = "ISCII_")
FeaturePlot(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, "ISCII_1", pt.size = 0, min.cutoff = 0,cols = c("lightgrey","red","darkred"), label = T, split.by = "group")

ISCIII <- read.table("./Results/Signatures/ISC3_signature.txt", header = TRUE)
Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12 <- AddModuleScore(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, features = ISCIII, name = "ISCIII_")
FeaturePlot(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, "ISCIII_1", pt.size = 0, min.cutoff = 0,cols = c("lightgrey","red","darkred"), label = T, split.by = "group")

RegevStemTF <- read.table("./Results/Signatures/Regev_Lgr5.txt", header = TRUE)
Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12 <- AddModuleScore(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, features = RegevStemTF, name = "RegevStemTF_")
FeaturePlot(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, "RegevStemTF_1", pt.size = 0, min.cutoff = 0,cols = c("lightgrey","red","darkred"), label = T, split.by = "group")

Munoz <- read.table("./Results/Signatures/Munoz_Lgr5.txt", header = TRUE)
Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12 <- AddModuleScore(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, features = Munoz, name = "Munoz_")
FeaturePlot(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, "Munoz_1", pt.size = 0, min.cutoff = 0,cols = c("lightgrey","red","darkred"), label = T, split.by = "group")

KimLgr5 <- read.table("./Results/Signatures/Kim_Lgr5.txt", header = TRUE)
Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12 <- AddModuleScore(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, features = KimLgr5, name = "KimLgr5_")
FeaturePlot(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, "KimLgr5_1", pt.size = 0, min.cutoff = 0,cols = c("lightgrey","red","darkred"), label = T, split.by = "group")

MerlosLgr5 <- read.table("./Results/Signatures/MErlos_LGR5.txt", header = TRUE)
Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12 <- AddModuleScore(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, features = MerlosLgr5, name = "MerlosLgr5_")
FeaturePlot(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, "MerlosLgr5_1", pt.size = 0, min.cutoff = 0,cols = c("lightgrey","red","darkred"), label = T, split.by = "group")


DotPlot(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, features = c("Munoz_1","RegevStemTF_1","RegevStem_1","ISCI_1","ISCII_1","ISCIII_1"), dot.scale = 8, 
        col.min = 0,
        col.max = 1,
) + scale_colour_gradient2(low = "lightgrey", mid = "red", high = "darkred",midpoint =0.5)+RotatedAxis()

####Violin####
Idents(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12) <- "celltype.group"

VlnPlot(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12,cols = c("bisque","skyblue3","aquamarine3","indianred2"),"Munoz_1", pt.size = 0.0,idents = c("Stem_WT","Stem_DKO","Stem_LKO","Stem_ZKO"))+geom_signif(comparisons = list(c("Stem_WT","Stem_ZKO")),map_signif_level = TRUE, textsize =3) + ylim(NA,0.5) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)
VlnPlot(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12,cols = c("bisque","skyblue3","aquamarine3","indianred2"),"RegevStemTF_1", pt.size = 0.0,idents = c("Stem_WT","Stem_DKO","Stem_LKO","Stem_ZKO"))+geom_signif(comparisons = list(c("Stem_WT","Stem_ZKO")),map_signif_level = TRUE, textsize =3) + ylim(NA,0.5) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)
VlnPlot(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12,cols = c("bisque","skyblue3","aquamarine3","indianred2"),"RegevStem_1", pt.size = 0.0,idents = c("Stem_WT","Stem_DKO","Stem_LKO","Stem_ZKO"))+geom_signif(comparisons = list(c("Stem_WT","Stem_LKO")),map_signif_level = TRUE, textsize =3) + ylim(NA,0.5) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)
VlnPlot(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12,cols = c("bisque","skyblue3","aquamarine3","indianred2"),"ISCI_1", pt.size = 0.0,idents = c("Stem_WT","Stem_DKO","Stem_LKO","Stem_ZKO"))+geom_signif(comparisons = list(c("Stem_WT","Stem_ZKO")),map_signif_level = TRUE, textsize =3) + ylim(NA,1) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)
VlnPlot(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12,cols = c("bisque","skyblue3","aquamarine3","indianred2"),"ISCII_1", pt.size = 0.0,idents = c("Stem_WT","Stem_DKO","Stem_LKO","Stem_ZKO"))+geom_signif(comparisons = list(c("Stem_WT","Stem_ZKO")),map_signif_level = TRUE, textsize =3) + ylim(NA,2) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)
VlnPlot(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12,cols = c("bisque","skyblue3","aquamarine3","indianred2"),"ISCIII_1", pt.size = 0.0,idents = c("Stem_WT","Stem_DKO","Stem_LKO","Stem_ZKO"))+geom_signif(comparisons = list(c("Stem_WT","Stem_ZKO")),map_signif_level = TRUE, textsize =3) + ylim(NA,0.5) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)

Idents(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12)<-"celltype.group"

my_levels <- c("Stem_WT","Stem_ZKO","Stem_LKO","Stem_DKO")
Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12@meta.data$celltype.group <- factor(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12@meta.data$celltype.group, levels = my_levels)

###Signatures celltype#####
#####Features################
CBC <- c("Olfm4","Lgr5","Ascl2","Axin2","Gkn3")
Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12 <- AddModuleScore(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, features = list(CBC), name="CBC_")
Idents(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12)<-"seurat_clusters"
FeaturePlot(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, "CBC_1", pt.size = -0.1, min.cutoff = 0, cols = c("lightgrey","red","darkred"), label = F, split.by = "group")
FeaturePlot(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, "CBC_1", pt.size = -0.1, min.cutoff = 0, cols = c("lightgrey","red","darkred"), label = F)
######
TA <- c("Dmbt1", "Cdk2ap2", "Rbp7", "S100a16","Gjb3","Fgfbp1","Nrarp","Pycard")
Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12 <- AddModuleScore(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, features = list(TA), name="TA_")
FeaturePlot(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, "TA_1", pt.size = 0.1, min.cutoff = -0.25, cols = c("lightgrey","red","darkred"), label = F)
######
Cycling <- c("Top2a", "Mki67", "Pcna", "Birc5","Pcna","Hmgb2")
Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12 <- AddModuleScore(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, features = list(Cycling), name="Cycling_")
FeaturePlot(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, "Cycling_1", pt.size = 0.1, min.cutoff = -0.5, cols = c("lightgrey","red","darkred"), label = F)

#######
Enterocytes <- c("Krt20", "Vil1", "Alpi", "Fabp1","Apoa1","Apoa4","Reg3a","Reg3g")
Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12 <- AddModuleScore(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, features = list(Enterocytes), name="Enterocytes_")
FeaturePlot(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, "Enterocytes_1", pt.size = 0.1, min.cutoff = 0, cols = c("lightgrey","red","darkred"), label = F)

#######
Paneth <- c("Lyz1", "Defa17", "Defa22", "Defa24","Ang4")
Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12 <- AddModuleScore(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, features = list(Paneth), name="Paneth_")
FeaturePlot(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, "Paneth_1", pt.size = 0.1, min.cutoff = 0, cols = c("lightgrey","red","darkred"), label = F)
#######
Goblets <- c("Agr2", "Muc2", "Spink4", "Tff3")
Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12 <- AddModuleScore(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, features = list(Goblets), name="Goblets_")
FeaturePlot(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, "Goblets_1", pt.size = 0.1, min.cutoff = 1, cols = c("lightgrey","red","darkred"), label = F)
EE <- c("Chga", "Chgb", "Tak1", "Tph1")
Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12 <- AddModuleScore(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, features = list(EE), name="EE_")
FeaturePlot(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, "EE_1", pt.size = 0.1, min.cutoff = -0.1, cols = c("lightgrey","red","darkred"), label = F)
Tuft <- c("Dclk1", "Avil", "Gfi1b", "Trpm5")
Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12 <- AddModuleScore(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, features = list(Tuft), name="Tuft_")
FeaturePlot(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, "Tuft_1", pt.size = 0.1, min.cutoff = 0, cols = c("lightgrey","red","darkred"), label = F)


Idents(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12)<-"celltype"

my_levels <- c("Stem","Cycling_TA","TA","SecretoryProgenitors","Enterocytes","Goblets","EE","Tuft","Paneth")
Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12@meta.data$celltype <- factor(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12@meta.data$celltype, levels = my_levels)

my_levels <- c("Paneth","Tuft","EE","Goblets","Enterocytes","SecretoryProgenitors","TA","Cycling_TA","Stem")
Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12@meta.data$celltype <- factor(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12@meta.data$celltype, levels = my_levels)


###Dotplot####
DotPlot(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, features = c("Lgr5","Olfm4","Birc5","Top2a","Dmbt1","Pycard","Vil1","Alpi","Muc2","Zg16","Neurog3","Chga","Dclk1","Avil","Lyz1","Defa24"), dot.scale = 8, 
        col.min = 0,
        col.max = 2,
) + scale_colour_gradient2(low = "lightgrey", mid = "red", high = "darkred",midpoint = 1)+RotatedAxis() 


Idents(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12) <- "celltype"
Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12.biomarkers <- FindAllMarkers(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, only.pos = FALSE, min.pct = 0.2, logfc.threshold = 0.3)
write.table(as.matrix(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12.biomarkers),"./Results/Biomarkers/Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12.biomarkers.txt")

top5_Tumor <- Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12.biomarkers %>% group_by(cluster) %>% top_n(5, wt=avg_log2FC)
DoHeatmap(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, features =top5_Tumor$gene, angle = 0,draw.lines=T) + scale_fill_gradientn(colors = c("lightblue", "white", "darkred")) 


#####JNK signatures#####

OzanneAp1 <- read.delim("./Results/Signatures/Ap1_Ozanne_mouse.txt", header = TRUE)
Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12 <- AddModuleScore(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, features = OzanneAp1, name = "OzanneAp1_")
FeaturePlot(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, "OzanneAp1_1", pt.size = 0, min.cutoff = 0,cols = c("lightgrey","red","darkred"), label = T, split.by = "group")

MathewsAP1 <- read.delim("./Results/Signatures/Mattews_AP1_mouse.txt", header = TRUE)
Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12 <- AddModuleScore(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, features = MathewsAP1, name = "MathewsAP1_")
FeaturePlot(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, "MathewsAP1_1", pt.size = 0, min.cutoff = 0,cols = c("lightgrey","red","darkred"), label = T, split.by = "group")

JNKscore <- read.delim("./Results/Signatures/JKN_EMBO_Mouse.txt", header = TRUE)
Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12 <- AddModuleScore(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, features = JNKscore, name = "JNKscore_")
FeaturePlot(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, "JNKscore_1", pt.size = 0, min.cutoff = 0,cols = c("lightgrey","red","darkred"), label = T, split.by = "group")

JNKhigh <- read.delim("./Results/Signatures/JNK_high_mouse.txt", header = TRUE)
Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12 <- AddModuleScore(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, features = JNKhigh, name = "JNKhigh_")
FeaturePlot(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, "JNKhigh_1", pt.size = 0, min.cutoff = 0,cols = c("lightgrey","red","darkred"), label = T, split.by = "group")

CjunTargets <- read.delim("./Results/Signatures/cJUNtargets_mouse.txt", header = TRUE)
Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12 <- AddModuleScore(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, features = CjunTargets, name = "CjunTargets_")
FeaturePlot(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, "CjunTargets_1", pt.size = 0, min.cutoff = 0,cols = c("lightgrey","red","darkred"), label = T, split.by = "group")

RegJNK <- read.delim("./Results/Signatures/Reg_JUNKactivity_mouse.txt", header = TRUE)
Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12 <- AddModuleScore(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, features = RegJNK, name = "RegJNK_")
FeaturePlot(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, "RegJNK_1", pt.size = 0, min.cutoff = 0,cols = c("lightgrey","red","darkred"), label = T, split.by = "group")

JNKcascade <- read.delim("./Results/Signatures/GO_JNKcascade_mouse.txt", header = TRUE)
Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12 <- AddModuleScore(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, features = JNKcascade, name = "JNKcascade_")
FeaturePlot(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, "JNKcascade_1", pt.size = 0, min.cutoff = 0,cols = c("lightgrey","red","darkred"), label = T, split.by = "group")

JNKposcascade <- read.delim("./Results/Signatures/GO_posREGJNKcascade_mouse.txt", header = TRUE)
Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12 <- AddModuleScore(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, features = JNKposcascade, name = "JNKposcascade_")
FeaturePlot(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, "JNKposcascade_1", pt.size = 0, min.cutoff = 0,cols = c("lightgrey","red","darkred"), label = T, split.by = "group")

AP1 <- read.delim("./Results/Signatures/PID_AP1_Mouse.txt", header = TRUE)
Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12 <- AddModuleScore(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, features = AP1, name = "AP1_")
FeaturePlot(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, "AP1_1", pt.size = 0, min.cutoff = 0,cols = c("lightgrey","red","darkred"), label = T, split.by = "group")

RegevStem <- read.delim("./Results/Signatures/Stem_REgev_Survey.txt", header = TRUE)
Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12 <- AddModuleScore(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, features = RegevStem, name = "RegevStem_")
FeaturePlot(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, "RegevStem_1", pt.size = 0, min.cutoff = 0,cols = c("lightgrey","red","darkred"), label = T, split.by = "group")


DotPlot(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, features = c("OzanneAp1_1","MathewsAP1_1","JNKscore_1","JNKhigh_1","CjunTargets_1","RegJNK_1","JNKcascade_1","JNKposcascade_1","AP1_1"), dot.scale = 8, 
        col.min = 0,
        col.max = 0.5,
) + scale_colour_gradient2(low = "lightgrey", mid = "red", high = "darkred",midpoint = 0.25)+RotatedAxis() 


#####
#####JNK signatures#####

OzanneAp1 <- read.delim("./Results/Signatures/Ap1_Ozanne_mouse.txt", header = TRUE)
Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12 <- AddModuleScore(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, features = OzanneAp1, name = "OzanneAp1_")
FeaturePlot(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, "OzanneAp1_1", pt.size = 0, min.cutoff = 0,cols = c("lightgrey","red","darkred"), label = F)

MathewsAP1 <- read.delim("./Results/Signatures/Mattews_AP1_mouse.txt", header = TRUE)
Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12 <- AddModuleScore(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, features = MathewsAP1, name = "MathewsAP1_")
FeaturePlot(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, "MathewsAP1_1", pt.size = 0, min.cutoff = 0,cols = c("lightgrey","red","darkred"), label = T, split.by = "group")

JNKscore <- read.delim("./Results/Signatures/JKN_EMBO_Mouse.txt", header = TRUE)
Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12 <- AddModuleScore(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, features = JNKscore, name = "JNKscore_")
FeaturePlot(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, "JNKscore_1", pt.size = 0, min.cutoff = 0,cols = c("lightgrey","red","darkred"), label = T, split.by = "group")

JNKhigh <- read.delim("./Results/Signatures/JNK_high_mouse.txt", header = TRUE)
Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12 <- AddModuleScore(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, features = JNKhigh, name = "JNKhigh_")
FeaturePlot(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, "JNKhigh_1", pt.size = 0, min.cutoff = 0,cols = c("lightgrey","red","darkred"), label = T, split.by = "group")

CjunTargets <- read.delim("./Results/Signatures/cJUNtargets_mouse.txt", header = TRUE)
Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12 <- AddModuleScore(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, features = CjunTargets, name = "CjunTargets_")
FeaturePlot(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, "CjunTargets_1", pt.size = 0, min.cutoff = 0,cols = c("lightgrey","red","darkred"), label = T, split.by = "group")

RegJNK <- read.delim("./Results/Signatures/Reg_JUNKactivity_mouse.txt", header = TRUE)
Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12 <- AddModuleScore(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, features = RegJNK, name = "RegJNK_")
FeaturePlot(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, "RegJNK_1", pt.size = 0, min.cutoff = 0,cols = c("lightgrey","red","darkred"), label = T, split.by = "group")

JNKcascade <- read.delim("./Results/Signatures/GO_JNKcascade_mouse.txt", header = TRUE)
Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12 <- AddModuleScore(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, features = JNKcascade, name = "JNKcascade_")
FeaturePlot(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, "JNKcascade_1", pt.size = 0, min.cutoff = 0,cols = c("lightgrey","red","darkred"), label = T, split.by = "group")

JNKposcascade <- read.delim("./Results/Signatures/GO_posREGJNKcascade_mouse.txt", header = TRUE)
Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12 <- AddModuleScore(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, features = JNKposcascade, name = "JNKposcascade_")
FeaturePlot(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, "JNKposcascade_1", pt.size = 0, min.cutoff = 0,cols = c("lightgrey","red","darkred"), label = T, split.by = "group")

AP1 <- read.delim("./Results/Signatures/PID_AP1_Mouse.txt", header = TRUE)
Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12 <- AddModuleScore(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, features = AP1, name = "AP1_")
FeaturePlot(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, "AP1_1", pt.size = 0, min.cutoff = 0,cols = c("lightgrey","red","darkred"), label = F)

######
saveRDS(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12,file="Results/WT_ZKO-LKO_DKO_dec22_HARM_UMAP")

####

Idents(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12) <- "group"


DotPlot(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, features = c("OzanneAp1_1","MathewsAP1_1","JNKscore_1","JNKhigh_1","CjunTargets_1","RegJNK_1","JNKcascade_1","JNKposcascade_1","AP1_1"), dot.scale = 8, 
        col.min = 0,
        col.max = 2,
) + scale_colour_gradient2(low = "lightgrey", mid = "red", high = "darkred",midpoint = 1)+RotatedAxis() 


######
DotPlot(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, features = c("OzanneAp1_1","JNKscore_1"), dot.scale = 8, 
        col.min = 0,
        col.max = 1,
) + scale_colour_gradient2(low = "lightgrey", mid = "red", high = "darkred",midpoint = 0.5)+RotatedAxis() 

#####
VlnPlot(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12,cols = c("bisque","skyblue3","aquamarine3","indianred2"),"OzanneAp1_1", pt.size = 0.0,idents = c("Stem_WT","Stem_DKO","Stem_LKO","Stem_ZKO"))+geom_signif(comparisons = list(c("Stem_WT","Stem_DKO")),map_signif_level = TRUE, textsize =3) + ylim(NA,0.5) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)
VlnPlot(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12,cols = c("bisque","skyblue3","aquamarine3","indianred2"),"JNKscore_1", pt.size = 0.0,idents = c("Stem_WT","Stem_DKO","Stem_LKO","Stem_ZKO"))+geom_signif(comparisons = list(c("Stem_ZKO","Stem_LKO")),map_signif_level = TRUE, textsize =3) + ylim(NA,0.2) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)
VlnPlot(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12,cols = c("bisque","skyblue3","aquamarine3","indianred2"),"JNKhigh_1", pt.size = 0.0,idents = c("Stem_WT","Stem_DKO","Stem_LKO","Stem_ZKO"))+geom_signif(comparisons = list(c("Stem_WT","Stem_LKO")),map_signif_level = TRUE, textsize =3) + ylim(NA,0.3) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)

####
VlnPlot(CyclingTA_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12,cols = c("bisque","skyblue3","aquamarine3","indianred2"),"OzanneAp1_1", pt.size = 0.0,idents = c("Cycling_TA_WT","Cycling_TA_DKO","Cycling_TA_ZKO","Cycling_TA_LKO"))+geom_signif(comparisons = list(c("Cycling_TA_WT","Cycling_TA_DKO")),map_signif_level = TRUE, textsize =3) + ylim(NA,1) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)
VlnPlot(CyclingTA_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12,cols = c("bisque","skyblue3","aquamarine3","indianred2"),"JNKscore_1", pt.size = 0.0,idents = c("Cycling_TA_WT","Cycling_TA_DKO","Cycling_TA_ZKO","Cycling_TA_LKO"))+geom_signif(comparisons = list(c("Cycling_TA_WT","Cycling_TA_DKO")),map_signif_level = TRUE, textsize =3) + ylim(NA,0.2) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)
VlnPlot(CyclingTA_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12,cols = c("bisque","skyblue3","aquamarine3","indianred2"),"JNKhigh_1", pt.size = 0.0,idents = c("Cycling_TA_WT","Cycling_TA_DKO","Cycling_TA_ZKO","Cycling_TA_LKO"))+geom_signif(comparisons = list(c("Cycling_TA_WT","Cycling_TA_DKO")),map_signif_level = TRUE, textsize =3) + ylim(NA,0.5) + stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)


######
Idents(CyclingTA_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12) <- "celltype.group"

CyclingTA_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12 <- subset(Combined_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, idents = c("Cycling_TA"))
###

my_levels <- c("Cycling_TA_WT","Cycling_TA_ZKO","Cycling_TA_LKO","Cycling_TA_DKO")
CyclingTA_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12@meta.data$celltype.group <- factor(CyclingTA_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12@meta.data$celltype.group, levels = my_levels)


########
OzanneAp1 <- read.delim("./Results/Signatures/Ap1_Ozanne_mouse.txt", header = TRUE)
CyclingTA_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12 <- AddModuleScore(CyclingTA_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, features = OzanneAp1, name = "OzanneAp1_")
FeaturePlot(CyclingTA_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, "OzanneAp1_1", pt.size = 0, min.cutoff = 0,cols = c("lightgrey","red","darkred"), label = T, split.by = "group")

MathewsAP1 <- read.delim("./Results/Signatures/Mattews_AP1_mouse.txt", header = TRUE)
CyclingTA_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12 <- AddModuleScore(CyclingTA_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, features = MathewsAP1, name = "MathewsAP1_")
FeaturePlot(CyclingTA_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, "MathewsAP1_1", pt.size = 0, min.cutoff = 0,cols = c("lightgrey","red","darkred"), label = T, split.by = "group")

JNKscore <- read.delim("./Results/Signatures/JKN_EMBO_Mouse.txt", header = TRUE)
CyclingTA_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12 <- AddModuleScore(CyclingTA_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, features = JNKscore, name = "JNKscore_")
FeaturePlot(CyclingTA_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, "JNKscore_1", pt.size = 0, min.cutoff = 0,cols = c("lightgrey","red","darkred"), label = T, split.by = "group")

JNKhigh <- read.delim("./Results/Signatures/JNK_high_mouse.txt", header = TRUE)
CyclingTA_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12 <- AddModuleScore(CyclingTA_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, features = JNKhigh, name = "JNKhigh_")
FeaturePlot(CyclingTA_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, "JNKhigh_1", pt.size = 0, min.cutoff = 0,cols = c("lightgrey","red","darkred"), label = T, split.by = "group")

CjunTargets <- read.delim("./Results/Signatures/cJUNtargets_mouse.txt", header = TRUE)
CyclingTA_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12 <- AddModuleScore(CyclingTA_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, features = CjunTargets, name = "CjunTargets_")
FeaturePlot(CyclingTA_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, "CjunTargets_1", pt.size = 0, min.cutoff = 0,cols = c("lightgrey","red","darkred"), label = T, split.by = "group")

RegJNK <- read.delim("./Results/Signatures/Reg_JUNKactivity_mouse.txt", header = TRUE)
CyclingTA_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12 <- AddModuleScore(CyclingTA_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, features = RegJNK, name = "RegJNK_")
FeaturePlot(CyclingTA_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, "RegJNK_1", pt.size = 0, min.cutoff = 0,cols = c("lightgrey","red","darkred"), label = T, split.by = "group")

JNKcascade <- read.delim("./Results/Signatures/GO_JNKcascade_mouse.txt", header = TRUE)
CyclingTA_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12 <- AddModuleScore(CyclingTA_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, features = JNKcascade, name = "JNKcascade_")
FeaturePlot(CyclingTA_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, "JNKcascade_1", pt.size = 0, min.cutoff = 0,cols = c("lightgrey","red","darkred"), label = T, split.by = "group")

JNKposcascade <- read.delim("./Results/Signatures/GO_posREGJNKcascade_mouse.txt", header = TRUE)
CyclingTA_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12 <- AddModuleScore(CyclingTA_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, features = JNKposcascade, name = "JNKposcascade_")
FeaturePlot(CyclingTA_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, "JNKposcascade_1", pt.size = 0, min.cutoff = 0,cols = c("lightgrey","red","darkred"), label = T, split.by = "group")

AP1 <- read.delim("./Results/Signatures/PID_AP1_Mouse.txt", header = TRUE)
CyclingTA_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12 <- AddModuleScore(CyclingTA_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, features = AP1, name = "AP1_")
FeaturePlot(CyclingTA_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, "AP1_1", pt.size = 0, min.cutoff = 0,cols = c("lightgrey","red","darkred"), label = T, split.by = "group")

Idents(CyclingTA_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12) <- "celltype.group"

DotPlot(CyclingTA_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, features = c("OzanneAp1_1","MathewsAP1_1","JNKscore_1","JNKhigh_1","CjunTargets_1","RegJNK_1","JNKcascade_1","JNKposcascade_1","AP1_1"), dot.scale = 8, 
        col.min = 0,
        col.max = 1,
) + scale_colour_gradient2(low = "lightgrey", mid = "red", high = "darkred",midpoint = 0.5)+RotatedAxis() 


#######To do trajectories using Cytotrace############
#####Only Stem cells#####
####Write counts table#######
counts_matrix = GetAssayData(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, slot="counts")
# write the output table in txt format
write.table(round(counts_matrix, digits=3), file='sc.10xStemWTZKOLKODKO.txt', quote=T, sep="\t") 

#######Export metadata#######
#####Check metadata###
head(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12@meta.data,1)
####Run library
library(data.table)

meta.data.cluster <- as.data.frame(as.matrix(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12@meta.data))
fwrite(x = meta.data.cluster, row.names = TRUE, file = "outfileStemWTZKOLKODKO.csv")

Idents(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12)<-"celltype.group"
Cytotrace_name <- read.delim("./Results/Cytotrace/CytoTRACE_results_Dec22.txt", header = TRUE, col.names = 1)
Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12 <- AddMetaData(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12,Cytotrace_name, col.name = "Cytotrace_name")

head(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12@meta.data, 5)


library(RColorBrewer)
cols <- brewer.pal(100, "Spectral")

FeaturePlot(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12,"Cytotrace_name", split.by = "group") & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))

Idents(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, "celltype.group")


RidgePlot(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, "Cytotrace_name", sort = T)



#####GSEA for Stem cells 4 GT ####

#######
Idents(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12)<-"group"
Stem_ZKO_vsStem_WT <-FindMarkers(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, ident.1 = c("ZKO"), ident.2 = c("WT"), verbose = T, min.pct = 0.2)
write.table(as.matrix(Stem_ZKO_vsStem_WT),"./Results/DEG//Stem_ZKO_vsStem_WT.txt", sep ="\t", col.names = T,row.names = T)

######
Idents(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12)<-"group"
Stem_LKO_vsStem_WT <-FindMarkers(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, ident.1 = c("LKO"), ident.2 = c("WT"), verbose = T, min.pct = 0.2)
write.table(as.matrix(Stem_LKO_vsStem_WT),"./Results/DEG//Stem_LKO_vsStem_WT.txt", sep ="\t", col.names = T,row.names = T)

######
Idents(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12)<-"group"
Stem_DKO_vsStem_WT <-FindMarkers(Progenitors_Epithelium_WCM2_WCM3_SD_HarmonyCLean_MT12, ident.1 = c("DKO"), ident.2 = c("WT"), verbose = T, min.pct = 0.2)
write.table(as.matrix(Stem_DKO_vsStem_WT),"./Results/DEG//Stem_DKO_vsStem_WT.txt", sep ="\t", col.names = T,row.names = T)



save.image()

















####Disregards from here #####




#####Mito percentage is 15%

FeaturePlot(Combined_WT_DKOn_DKONT_Harmony, "Epcam", cols = c("Lightgrey","red","darkred"), split.by = "group", min.cutoff = 2.5)
VlnPlot(Combined_WT_DKOn_DKONT_Harmony,"percent.mt")

Idents(Combined_WT_DKOn_DKONT_Harmony) <- "celltype"


#####Levels ####

Idents(Combined_WT_DKOn_DKONT_Harmony)<-"group"

my_levels <- c("WT","DKOn","DKOtVEH")
Combined_WT_DKOn_DKONT_Harmony@meta.data$group <- factor(Combined_WT_DKOn_DKONT_Harmony@meta.data$group, levels = my_levels)


DotPlot(Combined_WT_DKOn_DKONT_Harmony, features = c("Birc5","Top2a","Muc2","Tff3","Vil1","Alpi","Dmbt1","Pycard","Dclk1","Avil","Lyz1","Defa24","Olfm4","Lgr5","Neurog3","Chga"), dot.scale = 8, 
        col.min = 0,
        col.max = 1,
) + scale_colour_gradient2(low = "lightgrey", mid = "red", high = "darkred",midpoint = 0.5)+RotatedAxis()

FeaturePlot(Combined_WT_DKOn_DKONT_Harmony_NoTumor, "Muc2", pt.size = 0.1, min.cutoff = 0, cols = c("lightgrey","red","darkred"), label = T, split.by = "group")
FeaturePlot(Combined_WT_DKOn_DKONT_Harmony_NoTumor, "Cdx2", pt.size = 0.1, min.cutoff = 0, cols = c("lightgrey","red","darkred"), label = F, split.by = "group")
FeaturePlot(Combined_WT_DKOn_DKONT_Harmony_NoTumor, "Sox9", pt.size = 0.1, min.cutoff = 0, cols = c("lightgrey","red","darkred"), label = F, split.by = "group")
FeaturePlot(Combined_WT_DKOn_DKONT_Harmony_NoTumor, "Msh6", pt.size = 0.1, min.cutoff = 0, cols = c("lightgrey","red","darkred"), label = F, split.by = "group")
FeaturePlot(Combined_WT_DKOn_DKONT_Harmony_NoTumor, "Olfm4", pt.size = 0.1, min.cutoff = 0, cols = c("lightgrey","red","darkred"), label = F, split.by = "group")
FeaturePlot(Combined_WT_DKOn_DKONT_Harmony_NoTumor, "Muc5ac", pt.size = 0.1, min.cutoff = 0, cols = c("lightgrey","red","darkred"), label = F, split.by = "group")


#####Group####
Idents(Combined_WT_DKOn_DKONT_Harmony_NoTumor) <- "celltype"
Combined_WT_DKOn_DKONT_Harmony_NoTumor$celltype.group <- paste(Idents(Combined_WT_DKOn_DKONT_Harmony_NoTumor), Combined_WT_DKOn_DKONT_Harmony_NoTumor$group, sep = "_")
head(Combined_WT_DKOn_DKONT_Harmony_NoTumor@meta.data,3)

Idents(Combined_WT_DKOn_DKONT_Harmony_NoTumor) <- "celltype.group"

Combined_WT_DKOn_DKONT_Harmony_NoTumor_Absortive <- subset(Combined_WT_DKOn_DKONT_Harmony_NoTumor, idents = c("Stem_WT","Cycling_TA_WT","CyclingTA_DKOn","CyclingTA_DKOtVEH","TA_WT","TA_DKOn","TA_DKOtVEH","Enterocytes_WT","Enterocytes_DKOn","Enterocytes_DKOtVEH"))

Idents(Combined_WT_DKOn_DKONT_Harmony_NoTumor_Absortive) <- "celltype.group"

DimPlot(Combined_WT_DKOn_DKONT_Harmony_NoTumor_Absortive)

my_levels <- c("Stem_WT","Cycling_TA_WT","CyclingTA_DKOn","CyclingTA_DKOtVEH","TA_WT","TA_DKOn","TA_DKOtVEH","Enterocytes_WT","Enterocytes_DKOn","Enterocytes_DKOtVEH")
Combined_WT_DKOn_DKONT_Harmony_NoTumor_Absortive@meta.data$celltype.group <- factor(Combined_WT_DKOn_DKONT_Harmony_NoTumor_Absortive@meta.data$celltype.group, levels = my_levels)

DotPlot(Combined_WT_DKOn_DKONT_Harmony_NoTumor_Absortive, features = c("Anxa1","Anxa10","Tff2","Pdx1","Fscn1","Tacstd2","Aqp5","Mdk","Clu","Ly6a"),dot.scale = 8, 
        col.min = 0,
        col.max = 1,
) + scale_colour_gradient2(low = "lightgrey", mid = "red", high = "darkred",midpoint = 0.5)+RotatedAxis()


DotPlot(Combined_WT_DKOn_DKONT_Harmony_NoTumor_Absortive, features = c("Cdx2","Sox9","Olfm4","Lgr5","Ascl2","Axin2","Cd44"),dot.scale = 8, 
        col.min = 0,
        col.max = 1,
) + scale_colour_gradient2(low = "lightgrey", mid = "red", high = "darkred",midpoint = 0.5)+RotatedAxis()


YAP <- read_xlsx("./Results/Signatures/YAP_Gregorieff.xlsx")
Combined_WT_DKOn_DKONT_Harmony_NoTumor_Absortive <- AddModuleScore(Combined_WT_DKOn_DKONT_Harmony_NoTumor_Absortive, features = YAP, name = "YAP_")
FeaturePlot(Combined_WT_DKOn_DKONT_Harmony_NoTumor_Absortive,"YAP_1", pt.size = 0.1, min.cutoff = -0.1, cols = c("lightgrey","red","darkred"), label = F, split.by = "group")

#Counts#
#Cell.counts
head(Combined_WT_DKOn_DKONT_Harmony_NoTumor@meta.data,2)


Idents(Combined_WT_DKOn_DKONT_Harmony_NoTumor) <- "celltype.group"
counts.Combined_WT_DKOn_DKONT_Harmony_NoTumor <-table(Idents(Combined_WT_DKOn_DKONT_Harmony_NoTumor))
write.table(as.matrix(counts.Combined_WT_DKOn_DKONT_Harmony_NoTumor),"./Results/Counts/counts.Combined_WT_DKOn_DKONT_Harmony_NoTumor.txt",sep="\t",col.names=T,row.names=T)


FeaturePlot(Combined_WT_DKOn_DKONT_Harmony_NoTumor, "Clu", cols = c("lightgrey","red","darkred"), split.by = "group")

save.image()


