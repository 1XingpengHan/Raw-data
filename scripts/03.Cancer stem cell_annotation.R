
dir.create('results/03.T_annotation')

####
library(Seurat)
sce <- readRDS("/home/pub252/users/wangtingyu/20240407_LUAD/results/01.cell_annotation/sce.rds")
sce <- subset(sce, cell_type=="Cancer stem cell")
saveRDS(sce,file = 'Cancer stem cell.rds')
setwd("/home/pub252/users/wangtingyu/20240407_LUAD/results/03.T_annotation/")
#####################################################################################
library(harmony)
s.genes = Seurat::cc.genes.updated.2019$s.genes
g2m.genes = Seurat::cc.genes.updated.2019$g2m.genes
sce = CellCycleScoring(sce, s.features=s.genes, g2m.features=g2m.genes, set.ident=TRUE)
sce = SCTransform(sce, vars.to.regress=c("nCount_RNA", "percent.mt", "percent.rb", "S.Score", "G2M.Score"), verbose=FALSE)
sce = RunPCA(sce, verbose=FALSE)
my_harmony_embeddings <- HarmonyMatrix(
  data_mat  = as.matrix(sce@reductions$pca@cell.embeddings),
  meta_data = sce@meta.data,
  vars_use  = 'Samples',
  do_pca = FALSE)

rownames(my_harmony_embeddings) <- rownames(sce@reductions$pca@cell.embeddings)
sce[["harmony"]] <- CreateDimReducObject(embeddings = my_harmony_embeddings, key = "harmony_", assay = DefaultAssay(sce))
ElbowPlot(sce, ndims=50)
sce = FindNeighbors(sce, dims=1:20, reduction="harmony")
sce = RunUMAP(sce, dims=1:20, reduction="harmony")
p = DimPlot(sce, reduction="umap", group.by="Samples", pt.size=3)
p
theme(legend.position="right", plot.title=element_blank())+scale_color_manual(values=colors)
ggsave("UMAP_Sample.pdf", p, height=6, width=6)
###################################### 
mydata = FindClusters(sce, resolution=.3)
UMAPPlot(mydata, pt.size=1, label=T, label.size=5)+NoLegend()
markers = FindAllMarkers(mydata, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
markers["pct.diff"]=markers$pct.1-markers$pct.2
write.csv(markers, "./Cancer stem cell_subcluster_DEG.csv")
saveRDS(mydata,"T_subcluster.rds")
mydata <- readRDS("mydata.rds")
############## ##################
VlnPlot(mydata, features=c("CD3D"), pt.size=0)+NoLegend()+theme(axis.title.x=element_blank())

marker_gene <- unique(c("LAMP3",
                        "CCR7",
                        "FSCN1",
                        "CCL19"
                        
                        
))
FeaturePlot(mydata, features = marker_gene, reduction = "umap", raster = FALSE, pt.size = 0.25, order = T,
            cols = c("lightgrey", "red"))
# 0: CHI3L2+ T cells
# 1: GRASP+ T cells
marker <- data.frame(cluster = 0:2,cell = 0:2)
marker[marker$cluster %in% c(0),2] <- 'Paneth cell'  ##GPR171,DNAJA1
marker[marker$cluster %in% c(1),2] <- 'LGR5+ stem cell'  ##PLAC8,PLBD1
marker[marker$cluster %in% c(2),2] <- 'Basal cell'  ##CEBPD,TCF7L2
# marker[marker$cluster %in% c(3),2] <- 'NK cells'  ##NKG7,GNLY
# marker[marker$cluster %in% c(4),2] <- 'Megakaryocyte'  ##GP9,ITGA2B
# marker[marker$cluster %in% c(5),2] <- 'Mast cells'  ##GATA2

table(mydata$seurat_clusters)



mydata = subset(mydata, seurat_clusters %in% c(0,1,2,4,6,7,8,9))
mydata@meta.data$cell_type <- sapply(mydata@meta.data$seurat_clusters,function(x){marker[x,2]})
# mydata@meta.data$seurat_clusters = droplevels(mydata@meta.data$seurat_clusters, exclude=c(3))
# ################
# cell_label = c('CHI3L2+ T cells','GRASP+ T cells')
# names(cell_label) <- levels(mydata)
# mydata <- RenameIdents(mydata, cell_label)
mydata[["cell_type"]] = Idents(mydata)
# colors=c4a('brewer.accent',10)
my.cols=brewer.pal(12,"Set3")[-c(2,9)]
custome_theme_1=theme(axis.text.y=element_text(family="Times",face="plain")
                      ,axis.text.x=element_text(family="Times",face="plain")
                      ,plot.title = element_text(hjust = 0.5,family="Times",face="plain")
                      ,axis.title.x=element_text(family="Times",face="plain")
                      ,axis.title.y=element_text(family="Times",face="plain")
                      ,legend.title = element_text(family="Times",face="plain")
                      ,legend.text = element_text(family="Times",face="plain"))
p = DimPlot(mydata,cols=my.cols,group.by='cell_type',
            reduction="umap",
            label = F, 
            pt.size = 2,
            label.size = 5) +
  theme(axis.line = element_line(size=0.1, colour = "black"), 
        #axis.text = element_blank(), 
        #axis.title = element_blank(),
        axis.ticks = element_blank()
  ) +ggtitle('')
p
p1=p+custome_theme_1
p1=ggarrange(AugmentPlot(p1,dpi = 300), ncol = 1, nrow = 1,
             legend ='right',legend.grob=get_legend(p1),
             font.label = list(size = 12, face = "bold",family ='Times')) # 
p1
ggsave("UMAP_Cancer stem cell_subcluster.pdf", p, width=8, height=6)
genes = c("DSC3","PKP1","ARSE","GJB1","CXCL3","PLCG2")
DotPlot(mydata, features=genes,group.by = 'cell_type', cols=c("lightgray", "red"))+theme_bw()+
  theme(text=element_text(family="Times"), axis.text.x=element_text(angle=30, hjust=1, size=12, face="bold"), 
        axis.text.y=element_text(face="bold", size=12), axis.title.x=element_blank(), axis.title.y=element_blank())+
  coord_flip()
ggsave("dotplot_gene_marker.pdf",height = 6,width = 8)
##
# genes = c("LGALS2","PGLYRP1","PLBD1","GATA2","GP9","TCF7L2")
VlnPlot(mydata, features=genes, group.by = 'cell_type',pt.size=0)+NoLegend()+theme(axis.title.x=element_blank())
ggsave("Cancer stem cell_vln_gene_marker.pdf",height =8,width = 6)
#####
bar <-  as.data.frame(with(mydata@meta.data, table(Samples, cell_type)))
ggplot(data=bar, aes(x=Samples, y=Freq, fill=cell_type))+ 
  geom_bar(stat="identity", position=position_fill())+
  scale_fill_manual(values=my.cols)+theme_classic()+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(face="bold", size=12, angle=0, hjust=0.5), 
        axis.text.y=element_text(face="bold", size=12), legend.text=element_text(face="bold", size=12), legend.title=element_blank(), legend.position="right")
# ggsave(cell_type_number.pdf",  width=6, height=6)
#################################################################### 
bar = bar %>% group_by(Samples) %>% mutate(percent=100*Freq/sum(Freq))

ggplot(data=bar, aes(x=cell_type, y=percent, fill=Samples,label = sprintf("%.2f", percent)))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_manual(values=c("#4E79A7","#FF9D9A"))+theme_classic()+
  ggtitle("Percent(%)")+
  geom_text(position = position_dodge(width = 0.9), vjust = -0.5, size = 4)+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(angle=30, hjust=1, size=12, face="bold"), 
        axis.text.y=element_text(face="bold", size=12), legend.text=element_text(face="bold", size=12), legend.title=element_blank())
ggsave("Cancer stem cell_barplot_pair_number.pdf",  width=5, height=3)

