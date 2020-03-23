setwd("E:/HMF3A_reprogramming/R_analysis/Day10_oct4_nanog/combined_data/")

spheroid_combined<-read.csv("spheroid_combined.csv",stringsAsFactors = F)
spheroid_nanog<-read.csv("spheroid_2d_int_NANOG.csv",stringsAsFactors = F)
spheroid_oct4<-read.csv("spheroid_2d_int_OCT4.csv",stringsAsFactors = F)
spheroid_actin<-read.csv("spheroid_2d_int_ACTIN.csv",stringsAsFactors = F)

#obtaining the clusters_first time
{
  dird<-"E:/HMF3A_reprogramming/R_analysis/Day10_oct4_nanog/D10_Spheroid_clustering/"
  dir.create(dird)
  setwd(dird)
  
  parameters_clustering_data<-c( "Vol..unit.","Surf..unit." ,"Spher..unit.","Feret..unit.","Ell_MajRad","Ell_Elon",
                                 "Moment1","Moment2","Moment3","Moment4","Moment5","DCMin..unit.","DCMax..unit.","DCMean..unit.","DCSD..unit.",
                                 "Area","Circ.","Feret","MinFeret","AR")
  labelling_geometrical_data<-c("Volume", "Surface Area","Sphericity","Feret (3D)", "Ellipsoid R1","Ellipsoid Ellongation",
                                "Moment1","Moment2","Moment3","Moment4","Moment5", "DC_Min.","DC_Max","DC_Mean.","DC_SD",
                                "Proj. Area", "Proj. Circularity","Proj. Feret","Proj. Min Feret","Proj. A.R")
  
  library(corrplot)
  require(plyr)
  library(dendextend)
  library(gplots)
  
  #get the data for clustering
  data<-spheroid_combined[!duplicated(spheroid_combined$Label),]
  as.character(data$Label)->rownames(data)
  data_susbet<-data[,which(colnames(data) %in% parameters_clustering_data)]
  data_susbet<-na.omit(data_susbet)
  data_susbet<-as.data.frame(scale(data_susbet,center = TRUE, scale = TRUE))

  #Remove highly correlated features
  cor_mat<-cor(data_susbet)
  cor_mat[lower.tri(cor_mat)] <- 0
  diag(cor_mat) <- 0
  filtered_data_susbet<-data_susbet[,!apply(cor_mat,2,function(x) any(abs(x) > 0.8))]
  cor_mat_filtered<-cor(filtered_data_susbet)
  
  #cluster and get braches st colors
  distmat<-as.dist(1-cor(t(filtered_data_susbet), method="pearson"))
  rowclusters = hclust(distmat, method="average")
  mycl <- cutree(rowclusters, h=max(rowclusters$height/1.5))
  levels(as.factor(mycl))

  clusterCols <- c("blue","forestgreen","gold","orchid")
  dend1 <- as.dendrogram(rowclusters)
  dend1 <- color_branches(dend1, k = 4, col = clusterCols)
  col_labels <- get_leaves_branches_col(dend1)
  col_labels <- col_labels[order(order.dendrogram(dend1))]
  
  clusterCols <- c("forestgreen","gold","blue","orchid")
  

  png(filename="Heatmap_clusters_all_samples.png", units="in",width=4, height=4 , pointsize=7, res=1200)
  par(font.axis = 2,font.lab=2, font=2,mar=c(10,4,1,1))
  data_heat<-heatmap.2(as.matrix(data_susbet),scale="none",Rowv = dend1, Colv=T, dendrogram = "row", 
                       density.info = "none",trace="none",col=redblue(100),labRow = NA, key.xlab="Norm. Values",
                       keysize =1.5,key.title=" ",cexCol = 0.7)
  dev.off()
  
  #finding the clusters
  data_heat_cl<-as.data.frame((t(data_heat$carpet)))
  data_heat_cl$Label<-row.names(data_heat_cl)
  mycl<-as.data.frame(mycl)
  mycl$Label<-rownames(mycl)
  
  comb<-merge(data,mycl,by="Label")
  rownames(comb)<-comb$Label
  cluster1<-subset(comb, comb$mycl==1)
  cluster2<-subset(comb, comb$mycl==2)
  cluster3<-subset(comb, comb$mycl==3)
  cluster4<-subset(comb, comb$mycl==4)

  
}

#obtaining the cluster_second time
{
  #get the data for clustering
  data<-spheroid_combined[!duplicated(spheroid_combined$Label),]
  as.character(data$Label)->rownames(data)
  data_susbet<-data[,which(colnames(data) %in% parameters_clustering_data)]
  data_susbet<-data_susbet[which(! rownames(data_susbet) %in% cluster4$Label),]
  
  data_susbet<-na.omit(data_susbet)
  data_susbet<-as.data.frame(scale(data_susbet,center = TRUE, scale = TRUE))
  
  #Remove highly correlated features
  cor_mat<-cor(data_susbet)
  cor_mat[lower.tri(cor_mat)] <- 0
  diag(cor_mat) <- 0
  filtered_data_susbet<-data_susbet[,!apply(cor_mat,2,function(x) any(abs(x) > 0.8))]
  cor_mat_filtered<-cor(filtered_data_susbet)
  
  #cluster and get braches st colors
  distmat<-as.dist(1-cor(t(filtered_data_susbet), method="pearson"))
  rowclusters = hclust(distmat, method="average")
  mycl <- cutree(rowclusters, h=max(rowclusters$height/1.5))
  levels(as.factor(mycl))
  
  clusterCols <- c("forestgreen","gold","orchid")
  dend1 <- as.dendrogram(rowclusters)
  dend1 <- color_branches(dend1, k = 3, col = clusterCols)
  col_labels <- get_leaves_branches_col(dend1)
  col_labels <- col_labels[order(order.dendrogram(dend1))]
  
  clusterCols <- c("orchid","forestgreen","gold")
  
  png(filename="features.png", units="in",width=4, height=4 , pointsize=7, res=1200)
  corrplot(cor_mat,method = "color", order = "hclust", tl.col = "black")
  dev.off()
  
  png(filename="features_filtered.png", units="in",width=4, height=4 , pointsize=7, res=1200)
  corrplot(cor_mat_filtered,method = "color", order = "hclust", tl.col = "black")
  dev.off()
  
  png(filename="filtered_Heatmap_clusters.png", units="in",width=4, height=4 , pointsize=7, res=1200)
  par(font.axis = 2,font.lab=2, font=2,mar=c(10,4,1,1))
  data_heat<-heatmap.2(as.matrix(filtered_data_susbet),scale="none",Rowv = dend1, Colv=FALSE, dendrogram = "row", 
                       density.info = "none",trace="none",col=redblue(100),labRow = NA, key.xlab="Norm. Values",
                       keysize =1.5,key.title=" ",cexCol = 0.8)
  dev.off()
  
  png(filename="Heatmap_clusters.png", units="in",width=4, height=4 , pointsize=7, res=1200)
  par(font.axis = 2,font.lab=2, font=2,mar=c(10,4,1,1))
  data_heat<-heatmap.2(as.matrix(data_susbet),scale="none",Rowv = dend1, Colv=T, dendrogram = "row", 
                       density.info = "none",trace="none",col=redblue(100),labRow = NA, key.xlab="Norm. Values",
                       keysize =1.5,key.title=" ",cexCol = 0.7)
  dev.off()
  
  #finding the clusters
  data_heat_cl<-as.data.frame((t(data_heat$carpet)))
  data_heat_cl$Label<-row.names(data_heat_cl)
  mycl<-as.data.frame(mycl)
  mycl$Label<-rownames(mycl)
  
  comb<-merge(data,mycl,by="Label")
  rownames(comb)<-comb$Label
  cluster1<-subset(comb, comb$mycl==1)
  cluster2<-subset(comb, comb$mycl==2)
  cluster3<-subset(comb, comb$mycl==3)
  cluster4<-subset(comb, comb$mycl==4)
  
  #plotting the cluster values
  for ( i in 1:length(parameters_clustering_data)){
    index<-which(colnames(comb)==parameters_clustering_data[i])
    
    filename_png<-paste(parameters_clustering_data[i],"box_black.png",sep="_")
    png(filename=filename_png, units="in",width=2, height=2, pointsize=7,  res=1200)
    par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2,bg='black',fg='white',col.lab='white',col.axis='white')
    boxplot(cluster1[,index],cluster2[,index],cluster3[,index], las=1, lty=1,cex.axis=0.7,
            names=c("Cluster\n1","Cluster\n2","Cluster\n3"),col=(clusterCols),ylab=labelling_geometrical_data[i])
    dev.off()
    filename_png<-paste(parameters_clustering_data[i],"box_white.png",sep="_")
    png(filename=filename_png, units="in",width=2, height=2, pointsize=7,  res=1200)
    par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2,bg='white',fg='black',col.lab='black',col.axis='black')
    boxplot(cluster1[,index],cluster2[,index],cluster3[,index], las=1, lty=1,cex.axis=0.7,
            names=c("Cluster\n1","Cluster\n2","Cluster\n3"),col=(clusterCols),ylab=labelling_geometrical_data[i])
    dev.off()
  }
  
}
#Reprogramming_states
{
  c1<-subset(spheroid_nanog,spheroid_nanog$Label %in% cluster1$Label)$Median
  c2<-subset(spheroid_nanog,spheroid_nanog$Label %in% cluster2$Label)$Median
  c3<-subset(spheroid_nanog,spheroid_nanog$Label %in% cluster3$Label)$Median

  
  filename_png<-paste("nanog_clusters","box_black.png",sep="_")
  png(filename=filename_png, units="in",width=2, height=2, pointsize=7,  res=1200)
  par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2,bg='black',fg='white',col.lab='white',col.axis='white')
  boxplot(c1,c2,c3, las=1, lty=1,cex.axis=0.7,
          names=c("Cluster\n1","Cluster\n2","Cluster\n3"),col=(clusterCols),ylab="Spheroid Mean Nanog levels")
  dev.off()
  
  filename_png<-paste("nanog_clusters","box_white.png",sep="_")
  png(filename=filename_png, units="in",width=2, height=2, pointsize=7,  res=1200)
  par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2,bg='white',fg='black',col.lab='black',col.axis='black')
  boxplot(c1,c2,c3, las=1, lty=1,cex.axis=0.7,
          names=c("Cluster\n1","Cluster\n2","Cluster\n3"),col=(clusterCols),ylab="Spheroid Mean Nanog levels")
  dev.off()
  
  c1<-subset(spheroid_actin,spheroid_actin$Label %in% cluster1$Label)$Median
  c2<-subset(spheroid_actin,spheroid_actin$Label %in% cluster2$Label)$Median
  c3<-subset(spheroid_actin,spheroid_actin$Label %in% cluster3$Label)$Median

  
  filename_png<-paste("actin_clusters","box_black.png",sep="_")
  png(filename=filename_png, units="in",width=2, height=2, pointsize=7,  res=1200)
  par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2,bg='black',fg='white',col.lab='white',col.axis='white')
  boxplot(c1,c2,c3, las=1, lty=1,cex.axis=0.7,
          names=c("Cluster\n1","Cluster\n2","Cluster\n3"),col=(clusterCols),ylab="Spheroid Mean Actin levels")
  dev.off()
  
  filename_png<-paste("actin_clusters","box_white.png",sep="_")
  png(filename=filename_png, units="in",width=2, height=2, pointsize=7,  res=1200)
  par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2,bg='white',fg='black',col.lab='black',col.axis='black')
  boxplot(c1,c2,c3, las=1, lty=1,cex.axis=0.7,
          names=c("Cluster\n1","Cluster\n2","Cluster\n3"),col=(clusterCols),ylab="Spheroid Mean Actin levels")
  dev.off()
  
  
  c1<-subset(spheroid_oct4,spheroid_oct4$Label %in% cluster1$Label)$Median
  c2<-subset(spheroid_oct4,spheroid_oct4$Label %in% cluster2$Label)$Median
  c3<-subset(spheroid_oct4,spheroid_oct4$Label %in% cluster3$Label)$Median

  
  filename_png<-paste("oct4_clusters","box_black.png",sep="_")
  png(filename=filename_png, units="in",width=2, height=2, pointsize=7,  res=1200)
  par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2,bg='black',fg='white',col.lab='white',col.axis='white')
  boxplot(c1,c2,c3, las=1, lty=1,cex.axis=0.7,
          names=c("Cluster\n1","Cluster\n2","Cluster\n3"),col=(clusterCols),ylab="Spheroid Mean Oct4 levels")
  dev.off()
  
  filename_png<-paste("Oct4_clusters","box_white.png",sep="_")
  png(filename=filename_png, units="in",width=2, height=2, pointsize=7,  res=1200)
  par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2,bg='white',fg='black',col.lab='black',col.axis='black')
  boxplot(c1,c2,c3, las=1, lty=1,cex.axis=0.7,
          names=c("Cluster\n1","Cluster\n2","Cluster\n3"),col=(clusterCols),ylab="Spheroid Mean Oct4 levels")
  dev.off()
  
  
}

#visualising the clusters
#PCA
{
  comb_subset<-comb[,which(colnames(comb) %in% colnames(filtered_data_susbet))]
  comb_subset<-as.data.frame(scale(comb_subset,center = TRUE, scale = TRUE))
  
  nuc_pca <- princomp(na.omit(comb_subset),  cor = TRUE, scores = TRUE)
  r<-nuc_pca$scores
  
  png(filename="pca_variance_vs_PCs_white.png", units="in", width=2, height=2 , pointsize=5, res=1200)
  par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2,bg='white',fg='black',col.lab='black',col.axis='black' )
  plot(nuc_pca, type="l", las=1, main="", pch=19, las=1)
  box()
  dev.off()
  
  pca_scores<-as.data.frame(r)
  d1<-subset(pca_scores,rownames(pca_scores)%in% cluster1$Label)
  d2<-subset(pca_scores,rownames(pca_scores)%in% cluster2$Label)
  d3<-subset(pca_scores,rownames(pca_scores)%in% cluster3$Label)

  png(filename="pca_samples1_2.png", units="in", width=2, height=2 , pointsize=7, res=1200)
  par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2)
  plot(pca_scores$Comp.1,pca_scores$Comp.2, las=1, col="white", xlab="PC1",ylab="PC2")
  points(d3$Comp.1,d3$Comp.2,col=clusterCols[3],pch=19)
  points(d2$Comp.1,d2$Comp.2,col=clusterCols[2],pch=19)
  points(d1$Comp.1,d1$Comp.2,col=clusterCols[1], pch=19)
  dev.off()
  
  png(filename="pca_samples1_3.png", units="in", width=2, height=2 , pointsize=7, res=1200)
  par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2)
  plot(pca_scores$Comp.1,pca_scores$Comp.3, las=1, col="white")
  points(d3$Comp.1,d3$Comp.3,col=clusterCols[3],pch=19)
  points(d2$Comp.1,d2$Comp.3,col=clusterCols[2],pch=19)
  points(d1$Comp.1,d1$Comp.3,col=clusterCols[1], pch=19)
  dev.off()
  
  png(filename="pca_samples2_3.png", units="in", width=2, height=2 , pointsize=7, res=1200)
  par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2)
  plot(pca_scores$Comp.2,pca_scores$Comp.3, las=1, col="white")
  points(d3$Comp.2,d3$Comp.3,col=clusterCols[3],pch=19)
  points(d2$Comp.2,d2$Comp.3,col=clusterCols[2],pch=19)
  points(d1$Comp.2,d1$Comp.3,col=clusterCols[1], pch=19)
  dev.off()
  
  
}

#Umap
{
  library(uwot)
  data<-comb
  data_subset<-data[,which(colnames(data) %in% parameters_clustering_data)]
  
  comb_subset<-as.data.frame(scale(comb_subset,center = TRUE, scale = TRUE))
  
  upmap_spheroids <- umap(comb_subset,n_neighbors =10,min_dist = 2, spread =3)
  
  upmap_spheroids <- data.frame(UMAP1 = upmap_spheroids[, 1],UMAP2 = upmap_spheroids[, 2],Label = rownames(comb_subset))
  d1<-subset(upmap_spheroids,upmap_spheroids$Label %in% cluster1$Label)
  d2<-subset(upmap_spheroids,upmap_spheroids$Label %in% cluster2$Label)
  d3<-subset(upmap_spheroids,upmap_spheroids$Label %in% cluster3$Label)

  png(filename="umap_cluster.png", units="in", width=2, height=2 , pointsize=7, res=1200)
  par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2)
  plot(upmap_spheroids$UMAP1,upmap_spheroids$UMAP2, las=1, col="white", xlab="Comp1",ylab="Comp2")
  points(d3$UMAP1,d3$UMAP2,col=clusterCols[3],pch=19)
  points(d2$UMAP1,d2$UMAP2,col=clusterCols[2],pch=19)
  points(d1$UMAP1,d1$UMAP2,col=clusterCols[1], pch=19)
  dev.off()
}

#tsne
{
  library(Rtsne)
  tsne_spheroids<- Rtsne(comb_subset,pca = T, perplexity = 15, theta = 0.5)
  tsne_spheroids <- data.frame( TSNE1 = tsne_spheroids$Y[, 1],TSNE2 = tsne_spheroids$Y[, 2],Label = rownames(comb_subset))
  
  d1<-subset(tsne_spheroids,tsne_spheroids$Label %in% cluster1$Label)
  d2<-subset(tsne_spheroids,tsne_spheroids$Label %in% cluster2$Label)
  d3<-subset(tsne_spheroids,tsne_spheroids$Label %in% cluster3$Label)

  png(filename="tsne_cluster.png", units="in", width=2, height=2 , pointsize=7, res=1200)
  par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2)
  plot(tsne_spheroids$TSNE1,tsne_spheroids$TSNE2, las=1, col="white", xlab="Comp1",ylab="Comp2")
  points(d3$TSNE1,d3$TSNE2,col=clusterCols[3],pch=19)
  points(d2$TSNE1,d2$TSNE2,col=clusterCols[2],pch=19)
  points(d1$TSNE1,d1$TSNE2,col=clusterCols[1], pch=19)
  dev.off()
}

#reprogramming_vsspheroid_properties
{
  T1<-subset(spheroid_nanog,spheroid_nanog$trial=="T1")
  T2<-subset(spheroid_nanog,spheroid_nanog$trial=="T2")
  T3<-subset(spheroid_nanog,spheroid_nanog$trial=="T3")
  T1$Median<-T1$Median/max(T1$Median)
  T2$Median<-T2$Median/max(T2$Median)
  T3$Median<-T3$Median/max(T3$Median)
  spheroid_nanog<-rbind(T1,T2,T3)
  
  T1<-subset(spheroid_oct4,spheroid_oct4$trial=="T1")
  T2<-subset(spheroid_oct4,spheroid_oct4$trial=="T2")
  T3<-subset(spheroid_oct4,spheroid_oct4$trial=="T3")
  T1$Median<-T1$Median/max(T1$Median)
  T2$Median<-T2$Median/max(T2$Median)
  T3$Median<-T3$Median/max(T3$Median)
  spheroid_oct4<-rbind(T1,T2,T3)
  
  T1<-subset(spheroid_actin,spheroid_actin$trial=="T1")
  T2<-subset(spheroid_actin,spheroid_actin$trial=="T2")
  T3<-subset(spheroid_actin,spheroid_actin$trial=="T3")
  T1$Median<-T1$Median/max(T1$Median)
  T2$Median<-T2$Median/max(T2$Median)
  T3$Median<-T3$Median/max(T3$Median)
  spheroid_actin<-rbind(T1,T2,T3)
  
  #nanog 
  {
    merged_data<-merge(cbind(Label=spheroid_combined$Label,Volume=spheroid_combined$Vol..unit., Circ=spheroid_combined$Circ.,
                             Spherecity=spheroid_combined$Spher..unit.,Ellongation=spheroid_combined$Ell_Elon, 
                             DCmin=spheroid_combined$DCMin..unit.),
                       cbind(Label=spheroid_nanog$Label,Median=spheroid_nanog$Median),by="Label")
    
    png(filename="Nanog_spheroid_volume.png", units="in", width=2, height=2 , pointsize=7, res=1200)
    par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2)
    plot(as.numeric(as.character(merged_data$Volume)),as.numeric(as.character(merged_data$Median)),
         las=1,pch=18,xlab="Spheroid Volume",ylab="Spheroid Mean Nanog")   
    dev.off()
    
    png(filename="Nanog_spheroid_circ.png", units="in", width=2, height=2 , pointsize=7, res=1200)
    par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2)
    plot(as.numeric(as.character(merged_data$Circ)),as.numeric(as.character(merged_data$Median)),
         las=1,pch=18,xlab="Spheroid Circularity",ylab="Spheroid Mean Nanog")              
    dev.off()
    
    png(filename="Nanog_spheroid_circ.png", units="in", width=2, height=2 , pointsize=7, res=1200)
    par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2)
    plot(as.numeric(as.character(merged_data$Circ)),as.numeric(as.character(merged_data$Median)),
         las=1,pch=18,xlab="Spheroid Circularity",ylab="Spheroid Mean Nanog")              
    dev.off()
    
    png(filename="Nanog_spheroid_sphericity.png", units="in", width=2, height=2 , pointsize=7, res=1200)
    par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2)
    plot(as.numeric(as.character(merged_data$Spherecity)),as.numeric(as.character(merged_data$Median)),
         las=1,pch=18,xlab="Spheroid Spherecity",ylab="Spheroid Mean Nanog")              
    dev.off()
    
    png(filename="Nanog_spheroid_DCmin.png", units="in", width=2, height=2 , pointsize=7, res=1200)
    par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2)
    plot(as.numeric(as.character(merged_data$DCmin)),as.numeric(as.character(merged_data$Median)),
         las=1,pch=18,xlab="Spheroid DC min",ylab="Spheroid Mean Nanog")              
    dev.off()
    
    png(filename="Nanog_spheroid_ellongation.png", units="in", width=2, height=2 , pointsize=7, res=1200)
    par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2)
    plot(as.numeric(as.character(merged_data$Ellongation)),as.numeric(as.character(merged_data$Median)),
         las=1,pch=18,xlab="Spheroid Ellongation",ylab="Spheroid Mean Nanog")              
    dev.off()
  }
  
  #oct4 
  {
    merged_data<-merge(cbind(Label=spheroid_combined$Label,Volume=spheroid_combined$Vol..unit., Circ=spheroid_combined$Circ.,
                             Spherecity=spheroid_combined$Spher..unit.,Ellongation=spheroid_combined$Ell_Elon, 
                             DCmin=spheroid_combined$DCMin..unit.),
                       cbind(Label=spheroid_oct4$Label,Median=spheroid_oct4$Median),by="Label")
    
    png(filename="Oct4_spheroid_volume.png", units="in", width=2, height=2 , pointsize=7, res=1200)
    par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2)
    plot(as.numeric(as.character(merged_data$Volume)),as.numeric(as.character(merged_data$Median)),
         las=1,pch=18,xlab="Spheroid Volume",ylab="Spheroid Mean Oct4")              
    dev.off()
    
    png(filename="Oct4_spheroid_circ.png", units="in", width=2, height=2 , pointsize=7, res=1200)
    par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2)
    plot(as.numeric(as.character(merged_data$Circ)),as.numeric(as.character(merged_data$Median)),
         las=1,pch=18,xlab="Spheroid Circularity",ylab="Spheroid Mean Oct4")              
    dev.off()
    
    png(filename="Oct4_spheroid_sphericity.png", units="in", width=2, height=2 , pointsize=7, res=1200)
    par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2)
    plot(as.numeric(as.character(merged_data$Spherecity)),as.numeric(as.character(merged_data$Median)),
         las=1,pch=18,xlab="Spheroid Spherecity",ylab="Spheroid Mean Oct4")              
    dev.off()
    
    png(filename="Oct4_spheroid_DCmin.png", units="in", width=2, height=2 , pointsize=7, res=1200)
    par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2)
    plot(as.numeric(as.character(merged_data$DCmin)),as.numeric(as.character(merged_data$Median)),
         las=1,pch=18,xlab="Spheroid DC min",ylab="Spheroid Mean Oct4")              
    dev.off()
    
    png(filename="Oct4_spheroid_ellongation.png", units="in", width=2, height=2 , pointsize=7, res=1200)
    par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2)
    plot(as.numeric(as.character(merged_data$Ellongation)),as.numeric(as.character(merged_data$Median)),
         las=1,pch=18,xlab="Spheroid Ellongation",ylab="Spheroid Mean Oct4")              
    dev.off()
  }
  
  #actin 
  {
    merged_data<-merge(cbind(Label=spheroid_combined$Label,Volume=spheroid_combined$Vol..unit., Circ=spheroid_combined$Circ.,
                             Spherecity=spheroid_combined$Spher..unit.,Ellongation=spheroid_combined$Ell_Elon, 
                             DCmin=spheroid_combined$DCMin..unit.),
                       cbind(Label=spheroid_actin$Label,Median=spheroid_actin$Median),by="Label")
    
    png(filename="actin_spheroid_volume.png", units="in", width=2, height=2 , pointsize=7, res=1200)
    par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2)
    plot(as.numeric(as.character(merged_data$Volume)),as.numeric(as.character(merged_data$Median)),
         las=1,pch=18,xlab="Spheroid Volume",ylab="Spheroid Mean Actin")              
    dev.off()
    
    png(filename="actin_spheroid_circ.png", units="in", width=2, height=2 , pointsize=7, res=1200)
    par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2)
    plot(as.numeric(as.character(merged_data$Circ)),as.numeric(as.character(merged_data$Median)),
         las=1,pch=18,xlab="Spheroid Circularity",ylab="Spheroid Mean Actin")              
    dev.off()
    
    png(filename="actin_spheroid_sphericity.png", units="in", width=2, height=2 , pointsize=7, res=1200)
    par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2)
    plot(as.numeric(as.character(merged_data$Spherecity)),as.numeric(as.character(merged_data$Median)),
         las=1,pch=18,xlab="Spheroid Spherecity",ylab="Spheroid Mean Actin")              
    dev.off()
    
    png(filename="actin_spheroid_DCmin.png", units="in", width=2, height=2 , pointsize=7, res=1200)
    par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2)
    plot(as.numeric(as.character(merged_data$DCmin)),as.numeric(as.character(merged_data$Median)),
         las=1,pch=18,xlab="Spheroid DC min",ylab="Spheroid Mean Actin")              
    dev.off()
    
    png(filename="actin_spheroid_ellongation.png", units="in", width=2, height=2 , pointsize=7, res=1200)
    par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2)
    plot(as.numeric(as.character(merged_data$Ellongation)),as.numeric(as.character(merged_data$Median)),
         las=1,pch=18,xlab="Spheroid Ellongation",ylab="Spheroid Mean Actin")              
    dev.off()
  }
  
}
