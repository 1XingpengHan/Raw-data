

if (T) {
  dir.create("scripts")
  dir.create("results")
  dir.create("files")
  dir.create("figures")
  dir.create("origin_datas/GEO",recursive = T)
  # dir.create("origin_datas/TCGA")
}
library(stringr)
library(tidydr)
library(openxlsx)
library(data.table)
library(reshape2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(clusterProfiler)
library(pheatmap)
library(ComplexHeatmap)
library(GSVA)
library(GSEABase)
library(fgsea)
library(corrplot)
library(colorspace)
library(survival)
library(survminer)
library(maftools)
library(vegan)
library(forcats)
library(ggpubr)
library(ggsci)
library(ggplot2)
library(rstatix)
library(ggstatsplot)
library(ggcor)
library(ggstance)
options(stringsAsFactors = F)
source('/home/pub252/projects/codes/mg_base.R')

library(Seurat)
# library(SeuratObject)
library(Matrix)
library(dplyr)
library(ggplot2)
library(magrittr)
library(gtools)
library(stringr)
library(tidyverse)
library(patchwork)
library(data.table)
library(RColorBrewer)
library(ggpubr)



#01.单细胞图谱###########
dir.create('results/01.cell_annotation')
###样本信息
data.sample=read.delim('origin_datas/GEO/GSE182434_cell_annotation.txt.gz')
head(data.sample)
colnames(data.sample)[3]='Samples'
table(data.sample$Tissue)
table(data.sample$TumorNormal)
data.sample=data.sample[data.sample$Tissue!='FL',]
table(data.sample$TumorNormal,data.sample$Tissue)

###单细胞count表达矩阵
data.count=fread('origin_datas/GEO/GSM4569783_Sepsis1_processed_UMI.txt.gz')
dim(data.count)
data.count=as.data.frame(data.count)
rownames(data.count)=data.count[,1]
length(rownames(data.count))
length(unique(rownames(data.count)))
data.count=data.count[,-1]
data.count[1:5,1:5]


data.count=data.count[,data.sample$ID]
dim(data.count)


##1.1 
library(Seurat)
setwd('origin_datas/GEO/')
scRNAlist <- list()
file <- list.files()
dir <- paste0("./",file)
samples_name = file
for(i in 1:length(file) ) {
  
  counts <- Read10X(data.dir=dir[i])
  # counts <- counts$`Gene Expression`
  scRNAlist[[i]] <- CreateSeuratObject(counts, project=samples_name[i], min.cells=3, min.features=200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]], add.cell.id=samples_name[i])  
  
  #
  if(T){ 
    scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^MT-")
  }
  #计
  if(T){
    scRNAlist[[i]][["percent.rb"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^RP[SL]")
  }
  #
  if(T){
    HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
    HB.genes <- CaseMatch(HB.genes, rownames(scRNAlist[[i]]))
    scRNAlist[[i]][["percent.HB"]]=PercentageFeatureSet(scRNAlist[[i]],features=HB.genes)
    
  }
}

### 
names(scRNAlist) <- samples_name
# 
scRNA <- merge(scRNAlist[[1]], scRNAlist[2:length(scRNAlist)])
sce <- scRNA
sce@meta.data$Samples <- substr(sce@meta.data$orig.ident,7,9)
sce@meta.data$Samples <- gsub("nor", "control", sce@meta.data$Samples, ignore.case = FALSE, perl = FALSE,
                              fixed = FALSE, useBytes = FALSE)
sce@meta.data$Samples <- gsub("ade", "LUAD", sce@meta.data$Samples, ignore.case = FALSE, perl = FALSE,
                              fixed = FALSE, useBytes = FALSE)
unique(sce@meta.data$Samples)
rm(datalist)
#
raw_meta=sce@meta.data
raw_count <- table(raw_meta$orig.ident)
raw_count
sum(raw_count)#  12814
pearplot_befor<-VlnPlot(sce,group.by ='orig.ident',
                        features = c("nFeature_RNA", "nCount_RNA","percent.mt"),
                        pt.size = 0,
                        ncol = 3)
pearplot_befor
sce=subset(sce, subset=nFeature_RNA>200 & nFeature_RNA<6000 & percent.mt<10)

######
sce <- NormalizeData(sce)
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
sce <- ScaleData(sce, features = rownames(sce))


#
sce <- RunPCA(sce, features = VariableFeatures(sce))
colnames(sce@meta.data)
##
library(harmony)
my_harmony_embeddings <- HarmonyMatrix(
  data_mat  = as.matrix(sce@reductions$pca@cell.embeddings),
  meta_data = sce@meta.data,
  vars_use  = 'orig.ident',
  do_pca = FALSE)

rownames(my_harmony_embeddings) <- rownames(sce@reductions$pca@cell.embeddings)
sce[["harmony"]] <- CreateDimReducObject(embeddings = my_harmony_embeddings, key = "harmony_", assay = DefaultAssay(sce))                              
#
ElbowPlot(sce,ndims = 50)+theme(text = element_text(family = 'Times',size = 12))
dev.off()
##1.2 
sce <- RunUMAP(sce, dims=1:20, reduction="harmony")
####
after_batch=DimPlot(sce,group.by='Samples',reduction="umap",label = F,pt.size = 0.2)+
  theme_dr(xlength = 0.3, ylength = 0.3,arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(text = element_text(family = 'Times',size = 12),panel.grid = element_blank())
after_batch


library(clustree)
sce <- FindNeighbors(sce, dims = 1:20, reduction="harmony")
#
sce <- FindClusters(object = sce,resolution = .1)
saveRDS(sce, file = "./clean_data.Rds")
DefaultAssay(sce) <- "RNA"
colnames(sce@meta.data)
length(table(sce@meta.data$seurat_clusters))
# sce=subset(sce, subset=seurat_clusters%in% c(0:7,9))
my.cols=brewer.pal(12,"Set3")[-c(2,9)]

seurat_clusters_umap=DimPlot(sce,group.by='seurat_clusters',reduction="umap",label = T,pt.size = 0.2,cols = my.cols)+
  theme_dr(xlength = 0.3, ylength = 0.3,arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(text = element_text(family = 'Times',size = 12),panel.grid = element_blank())
seurat_clusters_umap+
  DimPlot(sce,group.by='CellType',reduction="umap",label = T,pt.size = 0.2,cols = my.cols)+
  theme_dr(xlength = 0.3, ylength = 0.3,arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(text = element_text(family = 'Times',size = 12),panel.grid = element_blank())

table(sce$seurat_clusters,sce$CellType)

plotlength=1
p1 <- VlnPlot(sce, features = "nFeature_RNA", group.by='orig.ident', pt.size = 0) & NoLegend() & labs(x = NULL)
ggsave(p1, file = "./filter_data_VlnPlot_nFeature_RNA_nopoint.pdf", width = plotlength*length(levels(sce)), height = 6)
p1 <- VlnPlot(sce, features = "nCount_RNA", group.by='orig.ident', pt.size = 0) & NoLegend() & labs(x = NULL)
ggsave(p1, file = "./filter_data_VlnPlot_nCount_RNA_nopoint.pdf", width = plotlength*length(levels(sce)), height = 6)
p1 <- VlnPlot(sce, features = "percent.mt", group.by='orig.ident', pt.size = 0) & NoLegend() & labs(x = NULL)
ggsave(p1, file = "./filter_data_VlnPlot_percent.mt_nopoint.pdf", width = plotlength*length(levels(sce)), height = 6)
p1 <- VlnPlot(sce, features = "percent.rb", group.by='orig.ident', pt.size = 0) & NoLegend() & labs(x = NULL)
ggsave(p1, file = "./filter_data_VlnPlot_percent.rb_nopoint.pdf", width = plotlength*length(levels(sce)), height = 6)
##1.3 
# #
Logfc = 0.25
#
Minpct = 0.25
DefaultAssay(sce) <- "RNA"
colnames(sce@meta.data)
Idents(sce)<-'seurat_clusters'

sce.markers <- FindAllMarkers(object = sce,logfc.threshold = Logfc, min.pct = Minpct,only.pos = T)
head(sce.markers)

sce.markers["pct.diff"]=sce.markers$pct.1-sce.markers$pct.2
sce.markers <- sce.markers[sce.markers$p_val_adj<0.05,]
table(sce.markers$cluster)
length(unique(sce.markers$gene))#3845
head(sce.markers)
table(sce.markers$cluster)
write.csv(sce.markers,'seurat_clusters_degs.csv')
# 
# 
Top5 <- sce.markers %>% group_by(cluster) %>% slice_max(n =5, order_by = avg_logFC)
length(Top5$gene)
length(unique(Top5$gene))
###
DotPlot(object = sce, features = unique(Top5$gene),
        cols=c("snow", "blue"),scale = T,col.min = 0)+
  RotatedAxis()+ ggtitle("Marker Genes")+
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab('')+ylab('')+coord_flip()
table(sce$seurat_clusters,sce$CellType)
#

marker_gene <- unique(c("CD3D","CD4","CD8A","CCL5","CCR7","FOXP3","GNLY","IL2RA","IL7R","NCAM1","NKG7"))
marker_gene <- unique(c("CDH5","CLDN5","ENG","VWF","EPCAM","KRT8","KRT10","NOTCH3","RGS5",
                        "CD74","CD14","LYZ","C1QC","FCGR3A","CD1C","FCER1A","LRG1","KIT",
                        "TPSAB1","MS4A2","GYPA"))
FeaturePlot(sce, features = marker_gene, reduction = "umap", raster = FALSE, pt.size = 0.25, order = T,
            cols = c("lightgrey", "red"))
VlnPlot(sce, features=c("MFSD2A"), pt.size=0)+NoLegend()+theme(axis.title.x=element_blank())

marker <- data.frame(cluster = 0:8,cell = 0:8)

marker[marker$cluster %in% c(0),2] <- 'Dendritic cells'  ##LAMP3,MFSD2A
marker[marker$cluster %in% c(1),2] <- 'Endothelial cells'  ##CLDN5,FCN3
marker[marker$cluster %in% c(2),2] <- 'T cells'  ##CD3E,CD3D
marker[marker$cluster %in% c(3),2] <- 'Cancer stem cell'  ##CEACAM5，ANPEP
marker[marker$cluster %in% c(4),2] <- 'Fibroblast cells'  ##DCN,LUM
marker[marker$cluster %in% c(5),2] <- 'Epithelial cells'  ##SCEL，TNNC1
marker[marker$cluster %in% c(6),2] <- 'Myeloid cells'  ##MS4A4A,SPI1
marker[marker$cluster %in% c(7),2] <- 'Mast cells'  ##MS4A2,TPSAB1
marker[marker$cluster %in% c(8),2] <- 'Ciliated cell'  ##MORN5,ENKUR
# marker[marker$cluster %in% c(9),2] <- ''
marker

# sce=subset(sce, subset=seurat_clusters%in% c(0:5,7:11))
sce@meta.data$cell_type <- sapply(sce@meta.data$seurat_clusters,function(x){marker[x,2]})
# my.cols=brewer.pal(12,"Set3")[-9]
# cell_type_umap=DimPlot(sce,group.by='cell_type',reduction="umap",label = F,pt.size = 0.2,cols =my.cols)+
#   theme_dr(xlength = 0.3, ylength = 0.3,arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
#   theme(panel.grid = element_blank(),text = element_text(family = 'Times',size = 12))
# cell_type_umap=LabelClusters(cell_type_umap,id = 'cell_type',family='Times')
# cell_type_umap
custome_theme_1=theme(axis.text.y=element_text(family="Times",face="plain")
                      ,axis.text.x=element_text(family="Times",face="plain")
                      ,plot.title = element_text(hjust = 0.5,family="Times",face="plain")
                      ,axis.title.x=element_text(family="Times",face="plain")
                      ,axis.title.y=element_text(family="Times",face="plain")
                      ,legend.title = element_text(family="Times",face="plain")
                      ,legend.text = element_text(family="Times",face="plain"))
#
sc_umap = DimPlot(sce,cols=my.cols,group.by='cell_type',
                  reduction="umap",
                  label = F, 
                  pt.size = 0.2,
                  label.size = 5) +
  theme(axis.line = element_line(size=0.1, colour = "black"), 
        #axis.text = element_blank(), 
        #axis.title = element_blank(),
        axis.ticks = element_blank()
  ) +ggtitle('')
sc_umap
# LabelClusters(sc_umap,id = 'cell_type',family = 'Times',size = 6,fontface = 'bold',color = 'red')
sc_umap_Sample = DimPlot(sce,cols=my.cols,group.by='orig.ident',
                         reduction="umap",
                         label = "F", 
                         pt.size = 0.2,
                         label.size = 0) +
  theme(axis.line = element_line(size=0.1, colour = "black"), 
        #axis.text = element_blank(), 
        #axis.title = element_blank(),
        axis.ticks = element_blank()
  ) +ggtitle('')

sc_umap_Sample
p1=sc_umap+custome_theme_1
p2=sc_umap_Sample+custome_theme_1
p1=ggarrange(AugmentPlot(p1,dpi = 300), ncol = 1, nrow = 1,
             legend ='right',legend.grob=get_legend(p1),
             font.label = list(size = 12, face = "bold",family ='Times')) # 
ggsave(filename = 'sc_umap.pdf',plot = p1,he=8,width = 9)

p2=ggarrange(AugmentPlot(p2,dpi = 300), ncol = 1, nrow = 1,
             legend ='right',legend.grob=get_legend(p2),
             font.label = list(size = 12, face = "bold",family ='Times')) # 
ggsave(filename = 'sample_umap.pdf',plot = p2,he=8,width = 9)

library(ggpubr)
fig1AB=ggarrange(p2,p1,
                 ncol = 2, nrow = 1,
                 labels = c("A","B"), # 添加标签
                 font.label = list(size = 12, face = "bold",family ='Times')) # 
fig1AB
ggsave(filename = 'Fig1AB.png',plot = fig1AB,he=8,width = 18,dpi = 300)
ggsave(filename = 'Fig1AB.pdf',plot = fig1AB,he=8,width = 18)

table(sce$cell_type)

# 0：CD4+T    	IL7R,LTB
# 1：CD8+T		GZMK，CST7，CD8A，IFNG
# 3: Tregs		FOXP3,CTLA4,TIGIT,ICOS
# 2,7: B cell		SPIB,CD79A,MS4A1,JCHAIN
# 4:Myeloid cell	CST3,LYZ,C1QB
# 5: NK cell		XCL1,GNLY,XCL2,KLRD1

marker_gene=c("MS4A4A","SPI1",
              "MS4A2","TPSAB1",
              "CEACAM5","ANPEP",
              "CD3E","CD3D",
              "MORN5","ENKUR",
              "LAMP3","MFSD2A",
              "DCN","LUM",
              "SCEL","TNNC1",
              "CLDN5","FCN3")
Idents(sce)='cell_type'
marker.dot=DotPlot(object = sce, features = marker_gene,cols=c("#80B1D3", "#D95F02"),scale = T,col.min = 0)+
  RotatedAxis()+ ggtitle("Marker Genes")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.y = element_text(size=12),text = element_text(family = 'Times',size = 12)) +
  xlab('')+ylab('')#+coord_flip()
marker.dot
ggsave('marker.dot.pdf',marker.dot,height = 5,width = 8)
VlnPlot(sce, features=marker_gene, pt.size=0, ncol=5)+NoLegend()+theme(axis.title.x=element_blank())
ggsave('violin_Plot.pdf',height = 15,width =15)

saveRDS(sce,file = 'sce.rds')

##1.4 
#####
bar <-  as.data.frame(with(sce@meta.data, table(Samples, cell_type)))
ggplot(data=bar, aes(x=Samples, y=Freq, fill=cell_type))+ 
  geom_bar(stat="identity", position=position_fill())+
  scale_fill_manual(values=my.cols)+theme_classic()+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(face="bold", size=12, angle=0, hjust=0.5), axis.text.y=element_text(face="bold", size=12), legend.text=element_text(face="bold", size=12), legend.title=element_blank(), legend.position="right")
ggsave("cell_type_number.pdf",  width=6, height=6)
##########  
# Type_label = c("PBMMC", "HHD")
# bar$type = factor(bar$type, levels=Type_label)
bar = bar %>% group_by(Samples) %>% mutate(percent=100*Freq/sum(Freq))

ggplot(data=bar, aes(x=cell_type, y=percent, fill=Samples,label = sprintf("%.2f", percent)))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_manual(values=c("#4E79A7","#FF9D9A"))+theme_classic()+
  ggtitle("Percent(%)")+
  geom_text(position = position_dodge(width = 0.9), vjust = -0.5, size = 4)+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(angle=30, hjust=1, size=12, face="bold"), axis.text.y=element_text(face="bold", size=12), legend.text=element_text(face="bold", size=12), legend.title=element_blank())
ggsave("barplot_pair_number.pdf",  width=10, height=5)

# head(sce@meta.data)
# table(sce@meta.data$COO)
# meta.data=sce@meta.data
# head(meta.data)
# meta.data=meta.data[meta.data$COO!='N/A',]
# table(meta.data$COO)
# bar = meta.data %>% group_by(COO, cell_type) %>% count()
# Type_label = unique(meta.data$COO)
# bar$COO = factor(bar$COO, levels=Type_label)
# bar = bar %>% group_by(COO) %>% mutate(percent=100*n/sum(n))
# 
# p=ggplot(data=bar, aes(x=cell_type, y=percent, fill=COO,label = sprintf("%.2f", percent)))+
#   geom_bar(stat="identity", position=position_dodge())+
#   scale_fill_manual(values=c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#CCEBC5"))+theme_classic()+
#   ggtitle("Percent(%)")+xlab('')+
#   geom_text(position = position_dodge(width = 0.9), vjust = -0.5, size = 4)+
#   theme(axis.text.x=element_text(angle=30, hjust=1, size=12, face="bold"),
#         legend.text=element_text(family = 'Times', size=12),
#         legend.title=element_blank(), text = element_text(family = 'Times'))
# p
# 
# ggsave('results/02.T_NK_cell/cell_barplot.pdf',p,height = 5,width = 8)


