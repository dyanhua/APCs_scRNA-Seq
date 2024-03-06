library(WGCNA)
library ( pheatmap )
#####load DEGs exp matrix
data=APCs_exp[up_genes,]
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
datExpr0=as.data.frame(t(data));
names(datExpr0)=rownames(data)
gene_id=rownames(data)
dim(datExpr0)
##1. cluster and qc
sampleTree = hclust(dist(datExpr0), method = "average")
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
abline(h = 80000, col = "red");
clust = cutreeStatic(sampleTree, cutHeight = 80000, minSize = 10)
table(clust)

##2. choose power
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)
# Plot the results:
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

##3, construct gene net
net = blockwiseModules(datExpr0, power = 14, maxBlockSize = 20000,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.3,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "de-1-SUM-TOM",
                       verbose = 3)
netLabels = net$colors
netColors = labels2colors(net$colors)
netMEs = net$MEs;
geneTree = net$dendrograms[[1]];
plotDendroAndColors(geneTree, netColors,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
save(netMEs, netLabels, netColors, geneTree,
     file = "de-net-networkConstruction-auto.RData")  ##save data
load("de-net-networkConstruction-auto.RData")

##4, module link to trait
allTraits=read.table("trait.txt",header = T)
names(allTraits)
# they match those of datExpr0:
name=rownames(datExpr0)
traitRows = match(name, allTraits$name)
datTraits = allTraits[traitRows,-1]
rownames(datTraits) = allTraits[traitRows, 1]
# show that row names agree
table(rownames(datTraits)==rownames(datExpr0))
# Choose a module assignment
moduleColorschoose=netColors
# Define numbers of genes and samples
nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)
# Recalculate MEs with color labels
#moduleColorschoose=gsub("green","red",moduleColorschoose)
MEs0 = moduleEigengenes(datExpr0,moduleColorschoose)$eigengenes
#MEs0 =MEs0[,-4]
MEschoose = orderMEs(MEs0)
modTraitCor = cor(MEschoose, datTraits, use = "p")
modTraitP = corPvalueStudent(modTraitCor, nSamples)
#Since we have a moderately large number of modules and traits,
#a suitable graphical representation will help in reading
#the table. We color code each association by the correlation value:
# Will display correlations and their p-values
textMatrix = paste(signif(modTraitCor, 2), "\n(",signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
dev.off()
# Display the correlation values within a heatmap plot
#### fig2i
p=labeledHeatmap(Matrix = modTraitCor, xLabels = names(datTraits),
                 yLabels = names(MEschoose), ySymbols = names(MEschoose),
                 colorLabels =FALSE,colors=blueWhiteRed(50),textMatrix=textMatrix,
                 setStdMargins = FALSE, cex.text = 0.6, zlim = c(-1,1),
                 main = paste("Module-trait relationships"))

##5, module genes expression pattern and functions
#### modules expression pattern
datExpr_scale1 = apply(datExpr0,2,scale) 
colnames(datExpr_scale1) = colnames(datExpr_scale1)
datExpr_scale_cols1 = numbers2colors(datExpr_scale1,colors = blueWhiteRed(100))
scale_cols1 = t(datExpr_scale_cols1)
colnames(scale_cols1) = rownames(datExpr0) 
plotDendroAndColors(geneTree, cbind(netColors,scale_cols1),
                    c("module", colnames(scale_cols1)),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
for (i in c("blue","brown","green","red","turquoise","yellow")) {
  exp=datExpr0[,netColors==i]
  pheatmap(exp,cluster_cols = TRUE, cluster_rows = FALSE,scale= "column")
  myheatmap <- pheatmap(exp,cluster_cols = TRUE, cluster_rows = FALSE,scale= "column")
  exp=t(rbind(exp,myheatmap$tree_col$order+100))
  exp=t(exp[order(exp[,13]),1:12])
  assign(paste(i,"_exp_order",sep=""),exp)
}

#### module function
for (i in c("blue","brown","green","red","turquoise","yellow")) {
  id = gene_id[netColors==i]
  id_gene <- bitr(id,fromType ="SYMBOL" ,toType = "ENTREZID",OrgDb =org.Hs.eg.db,drop = TRUE)
  id_ego <- enrichGO(id_gene[,2], 'org.Hs.eg.db', ont="BP", pvalueCutoff=0.01,qvalueCutoff = 0.05,readable = TRUE)
  filename=paste(i,"_ego_bp.csv",sep="")
  write.csv(summary(id_ego),paste(i,"_ego_bp.csv",sep=""),row.names =F)
  id_ekegg= enrichKEGG(id_gene[,2], organism = "human",  pvalueCutoff = 0.05,
                       pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 1)
  id_ekegg=setReadable(id_ekegg,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
  write.csv(summary(id_ekegg),paste(i,"_ekegg.csv",sep=""),row.names = F)
}

####fig2h
exp=cbind(turquoise_exp_order,green_exp_order,brown_exp_order,yellow_exp_order,blue_exp_order,red_exp_order)
pheatmap(exp,cluster_cols = F, cluster_rows = FALSE,scale= "column",color = blueWhiteRed(100),show_colnames = F)

####sfig3c
library(ComplexHeatmap)
data=datExpr0[,netColors %in% c("brown","yellow")]
names=c("ZFAND5","KLF4","ZNF331","ADAMTS1","SEMA4A","MMP19","POSTN","MYC","NR4A1","CSRNP1","NNMT","TIPARP","INTS6","CHSY1","PGAP1","GADD45B","SERTAD1","ANK2","RGS2","BGN")
for (i in 1:nrow(data)) data[i, ] <- scale(log(unlist(data[i, ]) + 1, 2))
data <- as.matrix(data)
samples <- rep(c('Lean control', 'Obese', 'Obese T2D'), c(4,4,4))  
class <- rep(c('CD55+ APCs', 'CD9+ APCs', 'ICAM1+ APCs',"CD142+ APCs"), c(3,3,3,3))  
RColorBrewer::display.brewer.all()
library(RColorBrewer)
heat <- Heatmap(data, cluster_rows = T,cluster_columns = F,
                col = colorRampPalette(colors = rev(brewer.pal(9,"RdBu")))(11), 
                heatmap_legend_param = list(grid_height = unit(10,'mm')),
                show_row_names = FALSE,
                top_annotation = HeatmapAnnotation(Group = samples[c(1,5,9,2,6,10,3,7,11,4,8,12)], 
                                                   simple_anno_size = unit(2, 'mm'), 
                                                   col = list(Group = c('Lean control' = '#FF9289', 'Obese' = '#00DAE0', 'Obese T2D' = 'pink')),
                                                   show_annotation_name = FALSE), 
                column_names_gp = gpar(fontsize = 10), row_names_gp = gpar(fontsize = 6))
p<- heat + rowAnnotation(link = anno_mark(at = which(rownames(data) %in% c(names)), 
                                          labels = c(names), labels_gp = gpar(fontsize = 6)))

