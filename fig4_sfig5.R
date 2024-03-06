####fig4a
mynet <- read.delim("/cellphonedb_out/count_network.txt", check.names = FALSE)
net<- graph_from_data_frame(mynet)
E(net)$width  <- E(net)$count/50  # 边点权重（粗细）
length(unique(mynet$SOURCE)) 
net2 <- net 
for (i in 1:length(unique(mynet$SOURCE))){
  E(net)[map(unique(mynet$SOURCE),function(x) {
    get.edge.ids(net,vp = c(unique(mynet$SOURCE)[i],x))
  })%>% unlist()]$color <- allcolour[i]
  plot(net, edge.arrow.size=0.1, 
       edge.curved=0.01,
       vertex.color=allcolour,
       vertex.frame.color="#555555",
       vertex.label.color="black",
       layout = layout_with_gem,
       vertex.label.cex=0.7) 
}
####fig4c
dedata=read.csv("Disease_Control_degs.csv")
library(ComplexHeatmap)

for (i in 1:nrow(dedata)) dedata[i, ] <- scale(log(unlist(dedata[i, ]) + 1, 2))
dedata <- as.matrix(dedata)
samples <- rep(c('Disease group', 'Control group'), c(3, 3))    
heat <- Heatmap(dedata, 
                col = colorRampPalette(colors = rev(brewer.pal(9,"RdBu")))(11), #定义热图由低值到高值的渐变颜色
                heatmap_legend_param = list(grid_height = unit(10,'mm')),  #图例高度设置
                show_row_names = FALSE,  #不展示基因名称
                top_annotation = HeatmapAnnotation(Group = samples, 
                                                   simple_anno_size = unit(2, 'mm'), 
                                                   col = list(Group = c('Disease group' = '#FF9289', 'Control group' = '#00DAE0')),  #定义样本分组的颜色
                                                   show_annotation_name = FALSE), 
                column_names_gp = gpar(fontsize = 10), row_names_gp = gpar(fontsize = 6))

names=c("IQGAP1","MTHFD1","ILF2","NEO1","VEGFD","WDR1","ECM2","GOT2","SERPINF1","CCDC80","TFPI","CALML3","COPS4","C3","SELENBP1")
heat + rowAnnotation(link = anno_mark(at = which(rownames(dedata) %in% c(names)), 
                                      labels = c(names), labels_gp = gpar(fontsize = 10)))

####fig4f, sfig5a
choose<-c("ANXA1_FPR1","CSF1_CSF1R","CXCL12_CXCR4","FGF2_FGFR1","FN1_ITGA4_ITGB1","FN1_ITGA5_ITGB1","GAS6_AXL","GAS6_MERTK","GRN_SORT1","MDK_LRP1","MDK_NCL","MIF_CD74_CXCR4","RARRES2_CMKLR1")
##also plot all pairs with P<0.01
####cellphonedb results
cellphoneout='/APCs_Myeloid/'
mypvals <- read.delim(paste0(cellphoneout,"pvalues.txt"), check.names = FALSE)
mymeans <- read.delim(paste0(cellphoneout,"means.txt"), check.names = FALSE)

mymeans %>%  dplyr::filter(interacting_pair %in% choose)%>%
  dplyr::select("interacting_pair",starts_with("CD9+ APCs"),ends_with("CD9+ APCs"))  %>%  
  reshape2::melt() -> meansdf
colnames(meansdf)<- c("interacting_pair","CC","means")
meansdf %>% dplyr::filter(CC %in% c("CD9+ APCs|Macrophage","CD9+ APCs|Monocyte")) -> meansdf

mypvals %>% dplyr::filter(interacting_pair %in% choose)%>%
  dplyr::select("interacting_pair",starts_with("CD9+ APCs"),ends_with("CD9+ APCs"))%>%  
  reshape2::melt()-> pvalsdf
colnames(pvalsdf)<- c("interacting_pair","CC","pvals")
pvalsdf %>% dplyr::filter(CC %in% c("CD9+ APCs|Macrophage","CD9+ APCs|Monocyte"))-> pvalsdf
pvalsdf$joinlab<- paste0(pvalsdf$interacting_pair,"_",pvalsdf$CC)
meansdf$joinlab<- paste0(meansdf$interacting_pair,"_",meansdf$CC)
pldf <- merge(pvalsdf,meansdf,by = "joinlab")
summary((filter(pldf,means >1))$means)

pldf%>% #filter(pvals<=0.01) %>% 
  ggplot(aes(CC.x,interacting_pair.x) )+ 
  geom_point(aes(color=means,size=-log10(pvals+0.0001)) ) +
  scale_size_continuous(range = c(1,3))+
  scale_color_gradientn(colours =brewer.pal(7,'YlOrRd'))+theme_bw()+ 
  theme(axis.text.x = element_text(angle = -45,hjust = -0.1,vjust = 0.8)) 

####cellchat results
cellchat=readRDS("cellchat_APCs_Myeloid_int.rds.gz")
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("FGF","CCL","CXCL","MIF","CSF","MK","PTN","ANNEXIN","GAS","GRN","CHEMERIN","FN1"),geneLR.return = F)
pairLR.use %>% dplyr::filter(interaction_name %in% choose)-> pairLR.use
netVisual_bubble(cellchat, sources.use = c(6), targets.use = c(3,4), pairLR.use = pairLR.use, remove.isolate = TRUE)+scale_color_gradientn(colours =brewer.pal(7,'YlOrRd'))
netVisual_bubble(cellchat, sources.use = c(6), targets.use = c(3,4), remove.isolate = TRUE)+scale_color_gradientn(colours =brewer.pal(7,'YlOrRd'))

