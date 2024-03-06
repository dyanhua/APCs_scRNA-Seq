library(Seurat)
library(ggplot2)
library(dplyr) 
library(magrittr)
library(harmony)
options(stringsAsFactors = F)
set.seed(42)

#### Load all scRNA-seq data
Allsamples_ident=readRDS("Allsamples_ident.rds.gz")
####cluster and markers, fig1a, fig1b, sfig1b
cluster_colors=c("#DC050C","#e1d77e","#ffbc40","brown","#1965B0","#50aabe","#98fb98","#977af2","#30d5c8","#00ffff")
names(cluster_colors)=rev(c("Mast cells","B cells","Myeloid cells","T/NKs","Adipocytes","Eodothelial cells","Mesothelium cells","APCs"))
DimPlot(Allsamples_ident, label = F,cols =cluster_colors,repel = T,pt.size = 0.001)
DimPlot(Allsamples_ident, label = F,cols =cluster_colors,split.by = "Group",repel = T,pt.size = 0.001)
markers=c("TPSB2","TPSAB1","JCHAIN","CD79A","CD14","CD68","NKG7","CD3D","PLIN1","ADIPOQ","PECAM1","CLDN5","MSLN","KRT18","APOD","PDGFRA")
DotPlot(Allsamples_ident, features = rev(markers), cols =c("lightblue","red"),dot.scale = 8,scale.by = "radius") + RotatedAxis()+
  scale_color_gradient2(low = "#1861b8", mid = "white", high =  "#dd5035",midpoint = 0.75) 
library(RColorBrewer)
FeaturePlot(Allsamples_ident,features = c("PDGFRA","ADIPOQ","VWF","KRT18","TPSAB1","CD68","CD3D","CD79A"),min.cutoff = 0.2,order = T,cols = c("grey","orange","red"),ncol = 8)

####statics, sfig1e
sname=Allsamples_ident
ClusterFreq <- sname@meta.data[,c("Celltype","orig.ident")] %>% table %>%
  data.frame() %>% set_colnames(c("Celltype","Sample","Per"))
ClusterPer <- ClusterFreq %>% tidyr::spread(Celltype,Per)
ClusterPer[,2:ncol(ClusterPer)] <- t(apply(ClusterPer[,2:ncol(ClusterPer)],1,function(x){x/sum(x)*100}))
ClusterPer <- ClusterPer %>% tidyr::gather(key=Celltype,value=Per,-Sample)
ClusterPer$Celltype=factor(ClusterPer$Celltype,levels =c("Mast cells","B cells","Myeloid cells","T/NKs","Adipocytes","Eodothelial cells","Mesothelium cells","APCs"))
ggplot(data = ClusterPer, aes(x = Sample, fill=Celltype,y=Per)) + 
  geom_bar(stat = "identity",width = 0.6)+
  scale_y_continuous(expand=c(0,0))+labs(y="percentage(%)")+coord_flip()+
  scale_fill_manual(values=cluster_colors)+
  theme(axis.text=element_text(color = "black"),
        panel.background = element_blank(),panel.grid=element_blank(),
        legend.title=element_blank(), axis.text.y=element_text(color = "black"),
        axis.title.x = element_blank(),axis.line = element_line(color="black"),legend.text = element_text(size=10))

####sfig1f
ClusterPer$Group=gsub("OT.*","Obese T2D",gsub("ON.*","Obese",gsub("H.*","Lean control",ClusterPer$Sample)))
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum <- ddply(data, groupnames, .fun=summary_func, varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}
ClusterPerdf <- data_summary(ClusterPer, varname="Per", groupnames=c("Group", "Celltype"))
ggplot(data = ClusterPerdf,  aes(x=Celltype, y=Per, fill=Group)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_errorbar(aes(ymin=Per -sd, ymax=Per +sd), width=.2,position=position_dodge(.9))+
  geom_jitter(data = ClusterPer,  aes(x=Celltype, y=Per, fill=Group),width =0.3,shape = 21,size=2)+
  scale_fill_manual(values = c("#92ca68","#c0ba73","#f08a7a"))

#### immune clusters define
immune=subset(Allsamples_ident,idents=c("Myeloid cells","T/NKs"))
immune  <- NormalizeData(immune) %>% FindVariableFeatures(nfeatures = 2000)
VGENES=VariableFeatures(immune)
VGENES=setdiff(VGENES,VGENES[grep("^MT|^RPL|^RPS|^HSP|^HBB|^HBA|^HBG|^HBM|^HBE|^HBZ",VGENES)])
immune <- ScaleData(immune,features =VGENES,vars.to.regress =c("percent.mito","nCount_RNA","percent.HB")) %>% RunPCA(npcs = 30, verbose = FALSE)%>% 
immune_harmony <- immune %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(immune_harmony, 'harmony')
ElbowPlot(immune_harmony)
immune_harmony <- immune_harmony %>% 
  RunUMAP(reduction = "harmony", dims = 1:10) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:10)
immune_harmony <- FindClusters(immune_harmony, resolution =0.5, algorithm = 1)
immune_ident <- RenameIdents(immune_harmony, `4` = "Naive T",`1` = "Tcm",`0` = "Tem",`3` = "NKT",`8` = "NK",`10` = "Profi",
                           `2` = "Macrophage",`6` = "Macrophage",`9` = "Macrophage",`5` = "Monocyte",`7` = "cDC2",`11` = "cDC1")
immune_ident$subtypes=immune_ident@active.ident

####sfig1c
immune_color=c("#ffd700","#fa5c00","#a0522d","#4b0082","#800080","#ffc0cb","#008b8b","#6495ed","#add8e6","#0027fa")
DimPlot(immune_ident,label = T,order = F,pt.size = 0.001,cols =immune_color)

#####markers sfig1d
immuneexp <- log1p(AverageExpression(immune_ident, verbose = FALSE)$RNA)
features1=rev(c("CD3D","CD3G","CD8A","CD4","CCR7","LEF1","SELL","TCF7","GPR183","CD40LG","IL7R","GZMK","GZMA","CST7","FGFBP2","GNLY","GZMB","PRF1","NKG7","KLRD1","XCL1","XCL2","MKI67","TOP2A"))
pheatmap(na.omit(immuneexp[features1,1:6]),cluster_cols = F,cluster_rows  = F,scale = "row",colorRampPalette (colors = c ("#496e9c","#70a8cc","#c7e6fd","white","#ffd7db","#f7a28a","#cd5c5c")) (100))
features2=rev(c("CD14","FCGR3A","C1QC","C1QB","FOLR2","PLTP","S100A8","S100A9","VCAN","FCN1","FCER1A","CD1C","CD1E","CLEC9A","IDO1","XCR1"))
pheatmap(na.omit(immuneexp[features2,c(7:10)]),cluster_cols = F,cluster_rows  = F,scale = "row",colorRampPalette (colors = c ("#496e9c","#70a8cc","#c7e6fd","white","#ffd7db","#f7a28a","#cd5c5c")) (100))

#### APCs clusters define
APCs=subset(Allsamples_ident,idents=c("APCs"))
APCs  <- NormalizeData(APCs) %>% FindVariableFeatures(selection.method = "vst",nfeatures = 2000)%>% ScaleData() %>% RunPCA(npcs = 30, verbose = FALSE)
APCs_harmony <- APCs %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(APCs_harmony, 'harmony')
ElbowPlot(APCs_harmony)
APCs_harmony <- APCs_harmony %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20)
APCs_harmony <- FindClusters(APCs_harmony, resolution =0.3, algorithm = 1)
APCs_ident <- RenameIdents(APCs_harmony, `2` = "CD55+ APCs",`3` = "CD9+ APCs", `0` = "ICAM1+ APCs",`1` = "CD142+ APCs")
APCs_ident$subtypes=APCs_ident@active.ident

#### fig1c, fig1d, sfig2a
apc_colors=c("#ffb6c1","#c71585","#8fbc8f","#add8e6")
DimPlot(APCs_ident,cols = apc_colors,pt.size = 0.1)
DimPlot(APCs_ident,split.by = "Sample",cols = apc_colors,pt.size = 0.1)
DimPlot(APCs_ident,pt.size = 0.001,group.by = "Sample",order = T)
markers=c("PDGFRA","CD55","CD9","ICAM1","F3","THY1","PI16","COL14A1","NABP1","CXCL12","CD34","SLPI","RARRES1","ABHD5","ABCA8","CD44","SEMA3C","SFRP4","ATF3","SVEP1")
VlnPlot(APCs_ident,group.by = "celltype",features = markers,ncol =5,cols =apc_colors,pt.size = 0)

####statics fig1e, fig1f,
sname=APCs_ident
ClusterFreq <- sname@meta.data[,c("subtypes","orig.ident")] %>% table %>%
  data.frame() %>% set_colnames(c("Celltype","Sample","Per"))
ClusterPer <- ClusterFreq %>% tidyr::spread(Celltype,Per)
ClusterPer[,2:ncol(ClusterPer)] <- t(apply(ClusterPer[,2:ncol(ClusterPer)],1,function(x){x/sum(x)*100}))
ClusterPer <- ClusterPer %>% tidyr::gather(key=Celltype,value=Per,-Sample)
ggplot(data = ClusterPer, aes(x = Sample, fill=Celltype,y=Per)) + #log2(as.numeric(CopyNumber)))
  geom_bar(stat = "identity",width = 0.6)+
  scale_y_continuous(expand=c(0,0))+labs(y="percentage(%)")+coord_flip()+
  scale_fill_manual(values=apc_colors)+
  theme(axis.text=element_text(color = "black"),
        panel.background = element_blank(),panel.grid=element_blank(),
        legend.title=element_blank(), axis.text.y=element_text(color = "black"),
        axis.title.x = element_blank(),axis.line = element_line(color="black"),legend.text = element_text(size=10))
ClusterPer$Group=gsub("OT.*","OT",gsub("ON.*","ON",gsub("H.*","HC",ClusterPer$Sample)))
ClusterPerdf <- data_summary(ClusterPer, varname="Per", groupnames=c("Group", "Celltype"))
head(ClusterPerdf)
ggplot(data = ClusterPerdf,  aes(x=Celltype, y=Per, fill=Group)) +
  geom_bar(stat="identity", position=position_dodge())+geom_errorbar(aes(ymin=Per -sd, ymax=Per +sd), width=.2,position=position_dodge(.9))+
  geom_jitter(data = ClusterPer,  aes(x=Celltype, y=Per, fill=Group),width =0.2,shape = 21,size=2.5)

#### fig1g, sfig2e, sfig2f, sfig2d
######Schwalie et al,Nature paper;
naturegenes=read.xlsx(xlsxFile =  "Supplementary_Table_7.xlsx",sheet = "Supplementary_Table_7")
APCs_ident@meta.data$Nature_G1<- apply(APCs_ident@assays$RNA@data[rownames(APCs_ident@assays$RNA@data) %in% toupper(naturegenes$G1[1:20]),],2,mean)
APCs_ident@meta.data$Nature_G2<- apply(APCs_ident@assays$RNA@data[rownames(APCs_ident@assays$RNA@data) %in% toupper(naturegenes$G2[1:20]),],2,mean)
APCs_ident@meta.data$Nature_G3<- apply(APCs_ident@assays$RNA@data[rownames(APCs_ident@assays$RNA@data) %in% toupper(naturegenes$G3[1:20]),],2,mean)
APCs_ident@meta.data$Nature_G4<- apply(APCs_ident@assays$RNA@data[rownames(APCs_ident@assays$RNA@data) %in% toupper(naturegenes$G4[1:20]),],2,mean)
VlnPlot(APCs_ident,features = c("Nature_G1","Nature_G2","Nature_G3","Nature_G4"),pt.size = 0,cols = apc_colors,ncol = 4)

######fig1g-science papers
science_g1 <-c("DPP4","PI16","CADM3","ALDH1A3","CD55","SEMA3C","PCOLCE2","WNT2","LIMCH1","IGFBP5","PCSK6","BMP7","RRAS2","AIF1L","DACT2")
science_g2 <-c("DLK1","PREF1","PPARG","COL15A1","GSC","NTF3","VCAM1","LPL","ADAMTS4","COL27A1","SPARCL1","BRNPER","SPON1","LIGP1","PDE3A","THSD7A","CASP4","CXCL9","GGT5","ENPEP","CLDN15","TRIM25","SDC1")
science_g3 <-c("F3","CLEC11A","IFITM1","MGP","ZIM1","GDF10","CYGB","IGFBP4","BGN","MGP","MFAP4","SPOCK2","GAS6","PLTP","BACE2","VIT","RBP1","ACTN1","DPEP1","VSTM4","NR2F2","EPHA3","C2","MEOX2","FMO2","ANGPT2","ARHGDIB","FSTL3","RASF11A","NDUFA4L2","PRDM6","MSX1","NTRK3","NAP1L5","FKBP1B","KCNS3","HSPB2")
APCs_ident@meta.data$g1_Dpp4<- apply(APCs_ident@assays$RNA@data[rownames(APCs_ident@assays$RNA@data) %in% science_g1,],2,mean)
APCs_identt@meta.data$g2_Icam1<- apply(APCs_ident@assays$RNA@data[rownames(APCs_ident@assays$RNA@data) %in% science_g2,],2,mean)
APCs_ident@meta.data$g3_Cd142<- apply(APCs_ident@assays$RNA@data[rownames(APCs_ident@assays$RNA@data) %in% science_g3,],2,mean)
VlnPlot(APCs_ident,features = c("Science_g1","Science_g2","Science_g3"),pt.size = 0,cols = apc_colors,ncol = 3)

#####human nature papers
humannaturegenes=read.xlsx(xlsxFile = "41586_2022_4518_MOESM4_ESM.xlsx",sheet =  "ASPC markers",colNames = T)
APCs_ident@meta.data$hASP1<- apply(APCs_ident@assays$RNA@data[rownames(APCs_ident@assays$RNA@data) %in% humannaturegenes[humannaturegenes$X7=="hASPC1",8][1:20],],2,mean)
APCs_ident@meta.data$hASP2<- apply(APCs_ident@assays$RNA@data[rownames(APCs_ident@assays$RNA@data) %in% humannaturegenes[humannaturegenes$X7=="hASPC2",8][1:20],],2,mean)
APCs_ident@meta.data$hASP3<- apply(APCs_ident@assays$RNA@data[rownames(APCs_ident@assays$RNA@data) %in% humannaturegenes[humannaturegenes$X7=="hASPC3",8][1:20],],2,mean)
APCs_ident@meta.data$hASP4<- apply(APCs_ident@assays$RNA@data[rownames(APCs_ident@assays$RNA@data) %in% humannaturegenes[humannaturegenes$X7=="hASPC4",8][1:20],],2,mean)
APCs_ident@meta.data$hASP5<- apply(APCs_ident@assays$RNA@data[rownames(APCs_ident@assays$RNA@data) %in% humannaturegenes[humannaturegenes$X7=="hASPC5",8][1:20],],2,mean)
APCs_ident@meta.data$hASP6<- apply(APCs_ident@assays$RNA@data[rownames(APCs_ident@assays$RNA@data) %in% humannaturegenes[humannaturegenes$X7=="hASPC6",8][1:20],],2,mean)
VlnPlot(APCs_ident,features = c("HumanNature_c1","HumanNature_c2","HumanNature_c3","HumanNature_c4","HumanNature_c5","HumanNature_c6"),pt.size = 0,cols = apc_colors,ncol = 6)

#### fig1h
library(slingshot)
sds_apcs <- slingshot(Embeddings(APCs_ident, "umap"), clusterLabels =APCs_ident@active.ident, 
                      start.clus = "CD55+ APCs", stretch = 0)
cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}
#We need color palettes for both cell types and Leiden clusters. These would be the same colors seen in the Seurat plots.
cell_colors_clust <- cell_pal(APCs_ident@active.ident,colorRampPalette(c("#ffb6c1","#c71585","#8fbc8f","#add8e6")) ) 
par(mar=c(5.1, 4.1, 4.1, 14), xpd=TRUE)
plot(reducedDim(sds_apcs), col = cell_colors_clust, pch = 16, cex = 0.5)
lines(sds_apcs, lwd = 1, col = 'black')
legend("right",inset = c(-0.7, 0),legend = unique(APCs_ident@active.ident), col = unique(cell_colors_clust), pch = 16, box.lwd = 0)
dev.off()

####psudotime
nc <- 3
pt <- slingPseudotime(sds_apcs)
nms <- colnames(pt)
nr <- ceiling(length(nms)/nc)
pal <- viridis(100, end = 0.95)
colors=paste(pal[cut(pt[,1], breaks = 100)],pal[cut(pt[,2], breaks = 100)],sep = "_") ##route1 & 2 combined
colors=gsub("_.*","",gsub("NA_","",gsub("_NA","",colors)))
plot_col <- function (nlev, col) {
  plot(x = c(1, (nlev + 1)), y = c(1, 2),
       xaxs = 'i', yaxs = 'i', xaxt = 'n', yaxt = 'n', type = 'n', ann = F)
  rect(xleft = 1:nlev,
       ybottom = rep(0, (nlev + 1)),
       xright = 2:(nlev + 1),
       ytop = rep(2, (nlev + 1)),
       col = col, border = 1)
  par(new = T)
  plot(x = c(1, (nlev + 1)), y = c(1, 2),
       xaxs = 'i', yaxs = 'i', xaxt = 'n', yaxt = 'n', type = 'n', ann = F)
}
library(s2dverification)
getPalette <- function(n) {
  # n should less than (length(colors)-1)*numOfColor
  colors <- c("#440154FF", "#472F7DFF", "#38578CFF", "#277E8EFF", "#1FA088FF","#43BF70FF","#96D83FFF","#DEE318FF")
  numOfColor <- 500
  resCols <- c()
  for (i in 1:(length(colors))-1) {
    tmpCols <- colorRampPalette(c(colors[i], colors[i+1]), numOfColor)
    resCols <- c(resCols, tmpCols(numOfColor))
  }
  colIndex <- seq(1, length(resCols), length.out=n)
  return(resCols[colIndex])
}
colorMin <- 0
colorMax <- 14
vals <- round(seq(colorMin, colorMax, length.out=50), 2)
pdf("slingshot_colorbar.pdf", width = 8, height = 10)
ColorBar(vals,
         color_fun=getPalette,
         label_scale=6,
         tick_scale=6,
         extra_margin=c(1, 0, 1, 0),
         label_digits=3)
dev.off()

#### fig1i
sname=APCs_ident
sname=newimport(sname)
sname <- estimateSizeFactors(sname)
sname <- estimateDispersions(sname)
sname <- detectGenes(sname, min_expr = 0.1)
head(pData(sname))
DEgenes=FindAllMarkers(subset(sname,downsample=2000), only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.15)
sname_monocle_seurat <- setOrderingFilter(sname, unique(DEgenes$gene))
plot_ordering_genes(sname_monocle_seurat)
sname_monocle_seurat <- reduceDimension(sname_monocle_seurat, max_components = 2, reduction_method = "DDRTree")
sname_monocle_seurat <- orderCells(sname_monocle_seurat)#root_state = 7
plot_cell_trajectory(apc1_monocle_clusters, color_by = "celltype",show_branch_points = T,cell_size = 0.5)+scale_color_manual(values=apc_colors)





