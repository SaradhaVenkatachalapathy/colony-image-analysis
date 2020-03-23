#libbraries
{library(gplots)
  library(RColorBrewer)
  library(corrplot)
  require(car)
  require(sp)
  require(rgeos)
  require(corrplot)
  require(plyr)
  library(RColorBrewer)
  library(MASS)
  library(Rtsne)
  library(e1071)
  library(gplots)
  library(dendextend)
  parameters_geometrical_data<-c( "Vol..unit.","Spher..unit.","Feret..unit.","Ell_Elon","Ell_Flatness",
                                  "Moment1","Moment2","Moment3","Moment4","Moment5","DCMin..unit.","DCMax..unit.","DCMean..unit.","DCSD..unit.",
                                  "Area","Perim.","Major","Minor","Circ.","Feret","MinFeret","AR","Round")
  
  labelling_geometrical_data<-c("Volume", "Sphericity","Feret (3D)","Ellipsoid Ellongation",
                                "Ellipsoid Flatness", "Moment1","Moment2","Moment3","Moment4","Moment5", "DC_Min.","DC_Max","DC_Mean.","DC_SD",
                                "Proj. Area", "Proj. Perimeter","Proj. Major Axis", "Proj. Minor Axis", "Proj. Circularity","Proj. Feret","Proj. Min Feret", 
                                "Proj. A.R","Proj. Roundness")
  
  }
# combbined_all_spheroids
{
  subfolders<-c("actin_pmlc_oct4/","vimentin_ecad_oct4/","ki67/","H3k9ac/")
  path=paste("E:/HMF3A_reprogramming/R_analysis/",subfolders[1],"combined_data/spheroid_combined.csv",sep="")
  spheroid_combined<-read.csv(path,header=T, stringsAsFactors = F)
  spheroid_combined$meg_exp<-subfolders[1]
  for (i in 2:length(subfolders)){
    path=paste("E:/HMF3A_reprogramming/R_analysis/",subfolders[i],"combined_data/spheroid_combined.csv",sep="")
    temp<-read.csv(path,header=T, stringsAsFactors = F)
    temp$meg_exp<-subfolders[i]
    spheroid_combined<-rbind(spheroid_combined,temp)
    rm(temp)
    
  }
  
  subfolders<-c("actin_pmlc_oct4/","none")
  path=paste("E:/HMF3A_reprogramming/R_analysis/",subfolders[1],"combined_data/spheroid_2d_int_oct4.csv",sep="")
  spheroid_oct4_combined<-read.csv(path,header=T, stringsAsFactors = F)
  spheroid_oct4_combined$meg_exp<-subfolders[1]
  
  subfolders<-c("actin_pmlc_oct4/","none")
  path=paste("E:/HMF3A_reprogramming/R_analysis/",subfolders[1],"combined_data/spheroid_int_actin.csv",sep="")
  spheroid_actin_combined<-read.csv(path,header=T, stringsAsFactors = F)
  spheroid_actin_combined$meg_exp<-subfolders[1]
  
  subfolders<-c("actin_pmlc_oct4/","none")
  path=paste("E:/HMF3A_reprogramming/R_analysis/",subfolders[1],"combined_data/cellular_int_pmlc.csv",sep="")
  spheroid_pmlc_combined<-read.csv(path,header=T, stringsAsFactors = F)
  spheroid_pmlc_combined$meg_exp<-subfolders[1]
  
  subfolders<-c("vimentin_ecad_oct4/","none")
  path=paste("E:/HMF3A_reprogramming/R_analysis/",subfolders[1],"combined_data/spheroid_int_VIMENTIN.csv",sep="")
  spheroid_vimentin_combined<-read.csv(path,header=T, stringsAsFactors = F)
  spheroid_vimentin_combined$meg_exp<-subfolders[1]
  
  subfolders<-c("vimentin_ecad_oct4/","none")
  path=paste("E:/HMF3A_reprogramming/R_analysis/",subfolders[1],"combined_data/spheroid_int_ECAD.csv",sep="")
  spheroid_ecad_combined<-read.csv(path,header=T, stringsAsFactors = F)
  spheroid_ecad_combined$meg_exp<-subfolders[1]

  subfolders<-c("H3k9ac/","none")
  path=paste("E:/HMF3A_reprogramming/R_analysis/",subfolders[1],"combined_data/spheroid_2d_int_H3K9AC.csv",sep="")
  spheroid_h3k9ac_combined<-read.csv(path,header=T, stringsAsFactors = F)
  spheroid_h3k9ac_combined$meg_exp<-subfolders[1]
}

dird<-"E:/HMF3A_reprogramming/R_analysis/Spheroid_clustering/"
dir.create(dird)
setwd(dird)
#obtaining and visualising the clusters
{
  
  data<-spheroid_combined
  data$X<-rownames(data)
  temp<-data[,which(colnames(data) %in% parameters_geometrical_data)]
  temp<-na.omit(temp)
  temp<-as.data.frame(scale(temp,center = TRUE, scale = TRUE))

  cor_mat<-cor(temp)
  png(filename="features.png", units="in",width=4, height=4 , pointsize=7, res=1200)
  corrplot(cor_mat,method = "color", order = "hclust", tl.col = "black")
  dev.off()
  
  cor_mat[lower.tri(cor_mat)] <- 0
  diag(cor_mat) <- 0
  filtered_temp<-temp[,!apply(cor_mat,2,function(x) any(x > 0.80))]
  cor_mat<-cor(filtered_temp)
  png(filename="features_filtered.png", units="in",width=4, height=4 , pointsize=7, res=1200)
  corrplot(cor_mat,method = "color", order = "hclust", tl.col = "black")
  dev.off()
  
  distmat<-as.dist(1-cor(t(filtered_temp), method="spearman"))
  rowclusters = hclust(distmat, method="average")
  mycl <- cutree(rowclusters, h=max(rowclusters$height/1.7))
  clusterCols <- c("blue","forestgreen","gold","orchid")
  dend1 <- as.dendrogram(rowclusters)
  dend1 <- color_branches(dend1, k = 4, col = clusterCols)
  col_labels <- get_leaves_branches_col(dend1)
  col_labels <- col_labels[order(order.dendrogram(dend1))]
  Colors=gray.colors(15, start = 1, end = 0.2)
  
  png(filename="3d_Heatmap.png_clusters.png", units="in",width=4, height=4 , pointsize=7, res=1200)
  par(font.axis = 2,font.lab=2, font=2,mar=c(10,4,1,1))
  data_heat<-heatmap.2(as.matrix(filtered_temp),scale="none",Rowv = dend1, Colv=FALSE, dendrogram = "row", 
                       density.info = "none",trace="none",col=redblue(100),labRow = NA, key.xlab="Intensity Fraction" ,keysize =1.5,key.title=" ",
                       distfun=function(x) as.dist(1-cor(t(x), method="pearson")), hclustfun=function(x) hclust(x, method="average"),cexCol = 0.7)
  dev.off()
  
  data_heat_cl<-as.data.frame((t(data_heat$carpet)))
  data_heat_cl$X<-row.names(data_heat_cl)
  mycl<-as.data.frame(mycl)
  mycl$X<-row.names(mycl)

  comb<-merge(data,mycl,by="X")
  
  cluster1<-subset(comb, comb$mycl==1)
  cluster2<-subset(comb, comb$mycl==2)
  cluster3<-subset(comb, comb$mycl==3)
  cluster4<-subset(comb, comb$mycl==4)
  
  library(uwot)
  data<-comb
  
  rownames(data)<-comb$X
  
  temp<-data[,which(colnames(data) %in% parameters_geometrical_data)]
  temp<-na.omit(temp)
  temp<-as.data.frame(scale(temp,center = TRUE, scale = TRUE))
  
  upmap_spheroids <- umap(temp,
                          n_neighbors = 10,
                          min_dist = 3, spread = 5
  )
  
  upmap_spheroids <- data.frame(
    UMAP1 = upmap_spheroids[, 1],
    UMAP2 = upmap_spheroids[, 2],
    X = rownames(temp)
  )
  
  d1<-subset(upmap_spheroids,upmap_spheroids$X %in% cluster1$X)
  d2<-subset(upmap_spheroids,upmap_spheroids$X %in% cluster2$X)
  d3<-subset(upmap_spheroids,upmap_spheroids$X %in% cluster3$X)
  d4<-subset(upmap_spheroids,upmap_spheroids$X %in% cluster4$X)
  
  png(filename="umap_cluster.png", units="in", width=2, height=2 , pointsize=5, res=1200)
  par(font.axis = 2,font.lab=2,mar=c(4,4,3,1), font=2)
  plot(upmap_spheroids$UMAP1,upmap_spheroids$UMAP2, las=1, col="white", xlab="Comp1",ylab="Comp2")
  points(d4$UMAP1,d4$UMAP2,col=clusterCols[4],pch=19)
  points(d3$UMAP1,d3$UMAP2,col=clusterCols[3],pch=19)
  points(d2$UMAP1,d2$UMAP2,col=clusterCols[2],pch=19)
  points(d1$UMAP1,d1$UMAP2,col=clusterCols[1], pch=19)
  dev.off()
  
}

#PCA

{
  
  
  data<-comb
  rownames(data)<-comb$X
  D0_lab<-(subset(data,data$sample.x=="D0"))$X
  D2_lab<-(subset(data,data$sample.x=="D2"))$X
  D4_lab<-(subset(data,data$sample.x=="D4"))$X
  D6_lab<-(subset(data,data$sample.x=="D6"))$X
  D8_lab<-(subset(data,data$sample.x=="D8"))$X
  
  temp<-data[,which(colnames(data) %in% parameters_geometrical_data)]
  temp<-na.omit(temp)
  temp<-as.data.frame(scale(temp,center = TRUE, scale = TRUE))
  
  nuc_pca <- princomp(na.omit(filtered_temp),  cor = TRUE, scores = TRUE)
  r<-nuc_pca$scores
  png(filename="pca_variance_vs_PCs.png", units="in", width=2, height=2 , pointsize=5, res=1200)
  par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2)
  plot(nuc_pca, type="l", las=1, main="", pch=19, las=1)
  box()
  dev.off()
  
  s<-abs(with(nuc_pca, unclass(loadings)))
  
  png(filename="pc1_loadings.png", units="in", width=2, height=2 , pointsize=5, res=1200)
  so<-s[order(s[,1]),1]
  par(font.axis = 2,font.lab=2,mar=c(5,12,3,1), font=2)
  barplot(so, horiz=T, las=1, xlab="loading on pC1")
  box()
  dev.off()
  
  png(filename="pc2_loadings.png", units="in", width=2, height=2 , pointsize=5, res=1200)
  so<-s[order(s[,2]),2]
  par(font.axis = 2,font.lab=2,mar=c(5,12,3,1), font=2)
  barplot(so, horiz=T, las=1, xlab="loading on pC1")
  box()
  dev.off()
  
  pca_scores<-as.data.frame(r)
  d1<-subset(pca_scores,rownames(pca_scores)%in% D0_lab)
  d2<-subset(pca_scores,rownames(pca_scores)%in% D2_lab)
  d3<-subset(pca_scores,rownames(pca_scores)%in% D4_lab)
  d4<-subset(pca_scores,rownames(pca_scores)%in% D6_lab)
  d5<-subset(pca_scores,rownames(pca_scores)%in% D8_lab)
  
  cols_time<-colorRampPalette(c("red","yellow"))( 5) 
  
  png(filename="pca_samples.png", units="in", width=2, height=2 , pointsize=5, res=1200)
  par(font.axis = 2,font.lab=2,mar=c(4,4,3,1), font=2)
  plot(pca_scores$Comp.1,pca_scores$Comp.2, las=1, col="white")
  points(d5$Comp.1,d5$Comp.2,col=cols_time[5],pch=19)
  points(d4$Comp.1,d4$Comp.2,col=cols_time[4],pch=19)
  points(d3$Comp.1,d3$Comp.2,col=cols_time[3],pch=19)
  points(d2$Comp.1,d2$Comp.2,col=cols_time[2],pch=19)
  points(d1$Comp.1,d1$Comp.3,col=cols_time[1], pch=19)
  dev.off()
  
  png(filename="pca_samplesendpoints.png", units="in", width=2, height=2 , pointsize=5, res=1200)
  par(font.axis = 2,font.lab=2,mar=c(4,4,3,1), font=2)
  plot(pca_scores$Comp.1,pca_scores$Comp.2, las=1, col="white")
  points(d5$Comp.1,d5$Comp.2,col=cols_time[5],pch=19)
  points(d1$Comp.1,d1$Comp.3,col=cols_time[1], pch=19)
  dev.off()
  
  d1<-subset(pca_scores,rownames(pca_scores)%in% cluster1$X)
  d2<-subset(pca_scores,rownames(pca_scores)%in% cluster2$X)
  d3<-subset(pca_scores,rownames(pca_scores)%in% cluster3$X)
  d4<-subset(pca_scores,rownames(pca_scores)%in% cluster4$X)
  
  png(filename="pca_cluster.png", units="in", width=2, height=2 , pointsize=5, res=1200)
  par(font.axis = 2,font.lab=2,mar=c(4,4,3,1), font=2)
  plot(pca_scores$Comp.1,pca_scores$Comp.2, las=1, col="white")
  points(d4$Comp.1,d4$Comp.2,col=clusterCols[4],pch=19)
  points(d3$Comp.1,d3$Comp.2,col=clusterCols[3],pch=19)
  points(d2$Comp.1,d2$Comp.2,col=clusterCols[2],pch=19)
  points(d1$Comp.1,d1$Comp.3,col=clusterCols[1], pch=19)
  dev.off()
  
  
}

clusterCols <- c("blue","orchid","gold","forestgreen")
#Defining clusters
{
  {
    D0_c1<-nrow(subset(cluster1,cluster1$sample=="D0"))/nrow(cluster1)
    D0_c2<-nrow(subset(cluster2,cluster2$sample=="D0"))/nrow(cluster2)
    D0_c3<-nrow(subset(cluster3,cluster3$sample=="D0"))/nrow(cluster3)
    D0_c4<-nrow(subset(cluster4,cluster4$sample=="D0"))/nrow(cluster4)
    c(D0_c1,D0_c2,D0_c3,D0_c4)
    D2_c1<-nrow(subset(cluster1,cluster1$sample=="D2"))/nrow(cluster1)
    D2_c2<-nrow(subset(cluster2,cluster2$sample=="D2"))/nrow(cluster2)
    D2_c3<-nrow(subset(cluster3,cluster3$sample=="D2"))/nrow(cluster3)
    D2_c4<-nrow(subset(cluster4,cluster4$sample=="D2"))/nrow(cluster4)
    c(D2_c1,D2_c2,D2_c3,D2_c4)
    D4_c1<-nrow(subset(cluster1,cluster1$sample=="D4"))/nrow(cluster1)
    D4_c2<-nrow(subset(cluster2,cluster2$sample=="D4"))/nrow(cluster2)
    D4_c3<-nrow(subset(cluster3,cluster3$sample=="D4"))/nrow(cluster3)
    D4_c4<-nrow(subset(cluster4,cluster4$sample=="D4"))/nrow(cluster4)
    c(D4_c1,D4_c2,D4_c3,D4_c4)
    D6_c1<-nrow(subset(cluster1,cluster1$sample=="D6"))/nrow(cluster1)
    D6_c2<-nrow(subset(cluster2,cluster2$sample=="D6"))/nrow(cluster2)
    D6_c3<-nrow(subset(cluster3,cluster3$sample=="D6"))/nrow(cluster3)
    D6_c4<-nrow(subset(cluster4,cluster4$sample=="D6"))/nrow(cluster4)
    c(D6_c1,D6_c2,D6_c3,D6_c4)
    D8_c1<-nrow(subset(cluster1,cluster1$sample=="D8"))/nrow(cluster1)
    D8_c2<-nrow(subset(cluster2,cluster2$sample=="D8"))/nrow(cluster2)
    D8_c3<-nrow(subset(cluster3,cluster3$sample=="D8"))/nrow(cluster3)
    D8_c4<-nrow(subset(cluster4,cluster4$sample=="D8"))/nrow(cluster4)
    c(D8_c1,D8_c2,D8_c3,D8_c4)
    
    clucter_timeline<-matrix(nrow=4,ncol=5)
    clucter_timeline[1,]<-c(D0_c1,D2_c1,D4_c1,D6_c1,D8_c1)
    clucter_timeline[2,]<-c(D0_c2,D2_c2,D4_c2,D6_c2,D8_c2)
    clucter_timeline[3,]<-c(D0_c3,D2_c3,D4_c3,D6_c3,D8_c3)
    clucter_timeline[4,]<-c(D0_c4,D2_c4,D4_c4,D6_c4,D8_c4)
    colnames(clucter_timeline)<-c("D0","D2","D4","D6","D8")
    rownames(clucter_timeline)<-c("c1","c2","c3","c4")
    
  }
  
  coltimeline<-colorRampPalette(c("red","yellow"))( 5) 
  
  png(filename="Timelineclusters.png", units="in",width=2, height=2 , pointsize=6, res=1200)
  par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2,bg='white',fg='black',col.lab='black',col.axis='black')
  barplot(t(clucter_timeline),col=coltimeline,las=1,ylab="%Cells",legend = colnames(clucter_timeline),xlim=c(0,7))
  box()
  dev.off()
  
  png(filename="Cluster_volume.png", units="in",width=2, height=2 , pointsize=6, res=1200)
  par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2,bg='white',fg='black',col.lab='black',col.axis='black')
  boxplot(cluster1$Vol..unit.,cluster2$Vol..unit.,cluster3$Vol..unit.,cluster4$Vol..unit.,col=clusterCols,las=1,
          ylab="Spheroid Volume",names=c("c1","c2","c3","c4"))
  dev.off()
  
  png(filename="Cluster_sphericiy.png", units="in",width=2, height=2 , pointsize=6, res=1200)
  par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2,bg='white',fg='black',col.lab='black',col.axis='black')
  boxplot(cluster1$Spher..unit.,cluster2$Spher..unit.,cluster3$Spher..unit.,cluster4$Spher..unit.,col=clusterCols,las=1,
          ylab="Spheroid Sphericity",names=c("c1","c2","c3","c4"))
  dev.off()
  
  png(filename="Cluster_Ell_Flatness.png", units="in",width=2, height=2 , pointsize=6, res=1200)
  par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2,bg='white',fg='black',col.lab='black',col.axis='black')
  boxplot(cluster1$Ell_Flatness,cluster2$Ell_Flatness,cluster3$Ell_Flatness,cluster4$Ell_Flatness,col=clusterCols,las=1,
          ylab="Spheroid Ell_Flatness",names=c("c1","c2","c3","c4"))
  dev.off()
  
  png(filename="Cluster_Feret.png", units="in",width=2, height=2 , pointsize=6, res=1200)
  par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2,bg='white',fg='black',col.lab='black',col.axis='black')
  boxplot(cluster1$Feret,cluster2$Feret,cluster3$Feret,cluster4$Feret,col=clusterCols,las=1,
          ylab="Spheroid Feret",names=c("c1","c2","c3","c4"))
  dev.off()
}

spheroid_oct4_combined$Label<-substring(spheroid_oct4_combined$Label,5, nchar(spheroid_oct4_combined$Label)-4)

c1_oct4<-subset(spheroid_oct4_combined,spheroid_oct4_combined$Label %in% cluster1$Label)
c2_oct4<-subset(spheroid_oct4_combined,spheroid_oct4_combined$Label %in% cluster2$Label)
c3_oct4<-subset(spheroid_oct4_combined,spheroid_oct4_combined$Label %in% cluster3$Label)
c4_oct4<-subset(spheroid_oct4_combined,spheroid_oct4_combined$Label %in% cluster4$Label)

png(filename="Cluster_oct4.png", units="in",width=2, height=2 , pointsize=6, res=1200)
par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2,bg='white',fg='black',col.lab='black',col.axis='black')
boxplot(c1_oct4$Mean,c2_oct4$Mean,c3_oct4$Mean,c4_oct4$Mean,col=clusterCols,las=1,
        ylab="Spheroid OCt4 Levels",names=c("c1","c2","c3","c4"))
dev.off()


merged_data_all<-merge(spheroid_oct4_combined, spheroid_combined,by="Label")
merged_data<-subset(merged_data_all,merged_data_all$meg_exp.x=="actin_pmlc_oct4/")
png(filename="Exp1.png", units="in",width=2, height=2 , pointsize=6, res=1200)
D8<-subset(merged_data,merged_data$sample.x=="D8")
D6<-subset(merged_data,merged_data$sample.x=="D6")
D4<-subset(merged_data,merged_data$sample.x=="D4")
D2<-subset(merged_data,merged_data$sample.x=="D2")
D0<-subset(merged_data,merged_data$sample.x=="D0")
plot(merged_data$Vol..unit.~merged_data$Mean.x,pch=18,col=NA, xlab="Mean Spheroid Oct4",ylab="Spheroid Volume",las=1)
points(D8$Vol..unit.~D8$Mean.x,col=coltimeline[5],pch=19)
points(D6$Vol..unit.~D6$Mean.x,col=coltimeline[4],pch=19)
points(D4$Vol..unit.~D4$Mean.x,col=coltimeline[3],pch=19)
points(D2$Vol..unit.~D2$Mean.x,col=coltimeline[2],pch=19)
points(D0$Vol..unit.~D0$Mean.x,col=coltimeline[1],pch=19)
dev.off()