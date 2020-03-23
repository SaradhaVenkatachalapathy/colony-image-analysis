#required functions
{
  require(plyr)
  library(MASS)
  library(corrplot)
  require(sp)
  require(rgeos)
  require(car)
  plotat <- function(RANGE) {
    if(length(RANGE) != 2) stop("RANGE argument must have a length of 2")
    if(RANGE[1] > RANGE[2]) stop("First element in RANGE must be smaller then second element")
    prettyres <- pretty(sprintf("%.2f",RANGE[1]):sprintf("%.2f",RANGE[2]), 7)
    while((min(prettyres) < RANGE[1]) == FALSE) {
      prdiff <- prettyres[2] - prettyres[1]
      prettyres[length(prettyres) + 1] <- prettyres[1] - prdiff
      prettyres <- sort(prettyres)
    } 
    while((max(prettyres) > RANGE[2]) == FALSE) {
      prdiff <- prettyres[2] - prettyres[1]
      prettyres[length(prettyres) + 1] <- prettyres[length(prettyres)] + prdiff
      prettyres <- sort(prettyres)    
    }   
    plotticks <- as.numeric(sprintf("%.2f",prettyres))
    plotticks
  }
  
  ## ellipseplot function
  ellipseplot <- function(x, y, factr, 
                          elev=0.95, # Ellipse probability level
                          pcol=NULL, # manual addition of colors, must meet length of factors
                          cexsize=1, # point size
                          ppch=21, # Point type, must meet length of factors
                          pbgcol=TRUE,
                          axissize=1, 
                          linewidth=1, 
                          font=1) {
    
    ## Set factor levels
    if(is.factor(factr)) {
      f <- factr
    } 
    else {
      f <- factor(factr, levels=unique(as.character(factr)))
    }
    intfactr <- as.integer(f) # Set integer vector that matches factor levels
    
    # Checking to make sure length of ppch equals number of factor levels
    if((length(ppch) > 1 & length(unique(intfactr)) != length(ppch))) stop("Can only increase point shape if equal to factor levels")
    
    ## Get data for ellipses
    edf <- data.frame(LV1 = x, LV2=y, factr = f) # create data frame with data and factor
    ellipses <- dlply(edf, .(factr), function(x) {
      LV1 <- x[,1]
      LV2 <- x[,2]
      dataEllipse(LV1, LV2, levels=elev, robust=TRUE, draw=FALSE) # Get confidence ellipse points from dataEllipse() function by factor level
    })
    ## Get range of x and y data
    xrange <- plotat(range(c(as.vector(sapply(ellipses, function(x) x[,1])), min(x), max(x))))
    yrange <- plotat(range(c(as.vector(sapply(ellipses, function(x) x[,2])), min(y), max(y))))
    
    
    ## Set colors for plots
    if(is.null(pcol) != TRUE) { # If colors are supplied by user
      ptcol <- pcol
      pgcol <- paste(pcol, "7e", sep="") # adds opaqueness
    } else { # Default
      pgcol <- c("#e41a1c7e","#377eb87e","#4daf4a7e","#984ea37e","#807f7d7e") # Defaults at 5 colors
      ptcol <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#807f7d") # For opaqueness
    }
    # Plotting graphic
    plot(x,y, type="n", xlab="", ylab="", main="", xlim=range(xrange), ylim=range(yrange), axes=FALSE)
    axis(1, at=xrange, labels=xrange, cex.axis=axissize,lwd=linewidth, font=font)
    axis(2, las=2, cex.axis=axissize,lwd=linewidth, font=font)
    box(lwd=linewidth, font=font)
    abline(h=0, v=0, col="gray", lty=2) # Adds lines at 0
    legpch <- c() # vector to collect legend pch data
    legcol <- c() # vector to collect legend col data
    
    ## Adds points, ellipse, and determines color specifications for legend 
    if(pbgcol==TRUE)  {
      for(i in 1:length(unique(intfactr))){
        points(x[intfactr==i], y[intfactr==i], pch=ppch[i], col=ptcol[i], bg=ptcol[i],cex=cexsize)
        polygon(ellipses[[i]], col=pgcol[i], border=ptcol[i])
        legpch[i] <- ppch[i]
        legcol[i] <- ptcol[i]
      }
    } else {
      for(i in 1:length(unique(intfactr))){
        points(x[intfactr==i], y[intfactr==i], pch=ppch[i], col="black", bg=ptcol[i],cex=cexsize)
        polygon(ellipses[[i]], col=pgcol[i], border=ptcol[i])
        legpch[i] <- ppch[i]
        legcol[i] <- ptcol[i]       
      }
    }
    
  }     
  
  ## Function for creating a SpatialPolygons object from data.frame of coords
  xy2SP <- function(xy, ID=NULL) {
    if(is.null(ID)) ID <- sample(1e12, size=1)
    SpatialPolygons(list(Polygons(list(Polygon(xy)), ID=ID)),
                    proj4string=CRS("+proj=merc"))
  }
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  
  normalize_2d_mean_trials <-function(data){
    
    if(length(levels(as.factor(data$trial))==3)){
      T1<-subset(data,data$trial=="T1")
      T1$Mean<-range01(T1$Mean)
      T1$IntDen<-range01(T1$IntDen)
      T1$Mode<-range01(T1$Mode)
      T1$Min<-range01(T1$Min)
      T1$Max<-range01(T1$Max)
      T1$Median<-range01(T1$Median)
      T1$RawIntDen<-range01(T1$RawIntDen)
      
      T2<-subset(data,data$trial=="T2")
      T2$Mean<-range01(T2$Mean)
      T2$IntDen<-range01(T2$IntDen)
      T2$Mode<-range01(T2$Mode)
      T2$Min<-range01(T2$Min)
      T2$Max<-range01(T2$Max)
      T2$Median<-range01(T2$Median)
      T2$RawIntDen<-range01(T2$RawIntDen)
      
      T3<-subset(data,data$trial=="T3")
      T3$Mean<-range01(T3$Mean)
      T3$IntDen<-range01(T3$IntDen)
      T3$Mode<-range01(T3$Mode)
      T3$Min<-range01(T3$Min)
      T3$Max<-range01(T3$Max)
      T3$Median<-range01(T3$Median)
      T3$RawIntDen<-range01(T3$RawIntDen)
      
      return(rbind(T1,T2,T3))
    }
    
    else if(length(levels(as.factor(data$trial))==4)){
      T1<-subset(data,data$trial=="T1")
      T1$Mean<-range01(T1$Mean)
      T1$IntDen<-range01(T1$IntDen)
      T1$Mode<-range01(T1$Mode)
      T1$Min<-range01(T1$Min)
      T1$Max<-range01(T1$Max)
      T1$Median<-range01(T1$Median)
      T1$RawIntDen<-range01(T1$RawIntDen)
      
      T2<-subset(data,data$trial=="T2")
      T2$Mean<-range01(T2$Mean)
      T2$IntDen<-range01(T2$IntDen)
      T2$Mode<-range01(T2$Mode)
      T2$Min<-range01(T2$Min)
      T2$Max<-range01(T2$Max)
      T2$Median<-range01(T2$Median)
      T2$RawIntDen<-range01(T2$RawIntDen)
      
      T3<-subset(data,data$trial=="T3")
      T3$Mean<-range01(T3$Mean)
      T3$IntDen<-range01(T3$IntDen)
      T3$Mode<-range01(T3$Mode)
      T3$Min<-range01(T3$Min)
      T3$Max<-range01(T3$Max)
      T3$Median<-range01(T3$Median)
      T3$RawIntDen<-range01(T3$RawIntDen)
      
      T4<-subset(data,data$trial=="T4")
      T4$Mean<-range01(T4$Mean)
      T4$IntDen<-range01(T4$IntDen)
      T4$Mode<-range01(T4$Mode)
      T4$Min<-range01(T4$Min)
      T4$Max<-range01(T4$Max)
      T4$Median<-range01(T4$Median)
      T4$RawIntDen<-range01(T4$RawIntDen)
      
      
      return(rbind(T1,T2,T3,T4))
    }
  }
  normalize_mean_trials <-function(data){
    
    if(length(levels(as.factor(data$trial))==3)){
      T1<-subset(data,data$trial=="T1")
      T1$Mean<-range01(T1$Mean)
      T1$IntDen<-range01(T1$IntDen)
      T1$Mode<-range01(T1$Mode)
      T1$Min<-range01(T1$Min)
      T1$Max<-range01(T1$Max)
      
      T2<-subset(data,data$trial=="T2")
      T2$Mean<-range01(T2$Mean)
      T2$IntDen<-range01(T2$IntDen)
      T2$Mode<-range01(T2$Mode)
      T2$Min<-range01(T2$Min)
      T2$Max<-range01(T2$Max)
      
      T3<-subset(data,data$trial=="T3")
      T3$Mean<-range01(T3$Mean)
      T3$IntDen<-range01(T3$IntDen)
      T3$Mode<-range01(T3$Mode)
      T3$Min<-range01(T3$Min)
      T3$Max<-range01(T3$Max)
      
      return(rbind(T1,T2,T3))
    }
    
    else if(length(levels(as.factor(data$trial))==4)){
      T1<-subset(data,data$trial=="T1")
      T1$Mean<-range01(T1$Mean)
      T1$IntDen<-range01(T1$IntDen)
      T1$Mode<-range01(T1$Mode)
      T1$Min<-range01(T1$Min)
      T1$Max<-range01(T1$Max)
      
      T2<-subset(data,data$trial=="T2")
      T2$Mean<-range01(T2$Mean)
      T2$IntDen<-range01(T2$IntDen)
      T2$Mode<-range01(T2$Mode)
      T2$Min<-range01(T2$Min)
      T2$Max<-range01(T2$Max)
      
      T3<-subset(data,data$trial=="T3")
      T3$Mean<-range01(T3$Mean)
      T3$IntDen<-range01(T3$IntDen)
      T3$Mode<-range01(T3$Mode)
      T3$Min<-range01(T3$Min)
      T3$Max<-range01(T3$Max)
      
      T4<-subset(data,data$trial=="T4")
      T4$Mean<-range01(T4$Mean)
      T4$IntDen<-range01(T4$IntDen)
      T4$Mode<-range01(T4$Mode)
      T4$Min<-range01(T4$Min)
      T4$Max<-range01(T4$Max)
      
      
      return(rbind(T1,T2,T3,T4))
    }
  }
  calclda <- function(variables,loadings)
  {
    # find the number of samples in the data set
    as.data.frame(variables)
    numsamples <- nrow(variables)
    # make a vector to store the discriminant function
    ld <- numeric(numsamples)
    # find the number of variables
    numvariables <- length(variables)
    # calculate the value of the discriminant function for each sample
    for (i in 1:numsamples)
    {
      valuei <- 0
      for (j in 1:numvariables)
      {
        valueij <- variables[i,j]
        loadingj <- loadings[j]
        valuei <- valuei + (valueij * loadingj)
      }
      ld[i] <- valuei
    }
    # standardise the discriminant function so that its mean value is 0:
    ld <- as.data.frame(scale(ld, center=TRUE, scale=FALSE))
    ld <- ld[[1]]
    return(ld)
  }
  
}
#read nucleus combined data from 5 mega sets 
{
  setwd("~/Desktop/Reprogramming/actin_pmlc_oct4/combined_data/")
  
  nucleus_combined_set1 <- read.csv("nucleus_combined.csv",stringsAsFactors = F)
  nucleus_combined_set1<-nucleus_combined_set1[!(duplicated(nucleus_combined_set1$Label)),]
  nucleus_combined_set1$meg_trial<-paste(nucleus_combined_set1$trial,"set1",sep="_")
  
  setwd("~/Desktop/Reprogramming/vimentin_ecad_oct4/combined_data/")
  
  nucleus_combined_set2 <- read.csv("nucleus_combined.csv",stringsAsFactors = F)
  nucleus_combined_set2<-nucleus_combined_set2[!(duplicated(nucleus_combined_set2$Label)),]
  nucleus_combined_set2$meg_trial<-paste(nucleus_combined_set2$trial,"set2",sep="_")
  
  setwd("~/Desktop/Reprogramming/ki67/combined_data/")
  
  nucleus_combined_set3 <- read.csv("nucleus_combined.csv",stringsAsFactors = F)
  nucleus_combined_set3<-nucleus_combined_set3[!(duplicated(nucleus_combined_set3$Label)),]
  nucleus_combined_set3<-nucleus_combined_set3[,which(colnames(nucleus_combined_set3) %in% colnames(nucleus_combined_set2))]
  nucleus_combined_set3$meg_trial<-paste(nucleus_combined_set3$trial,"set3",sep="_")
  
  
  setwd("~/Desktop/Reprogramming/H3k9ac/combined_data/")
  
  nucleus_combined_set4 <- read.csv("nucleus_combined.csv",stringsAsFactors = F)
  nucleus_combined_set4<-nucleus_combined_set4[!(duplicated(nucleus_combined_set4$Label)),]
  nucleus_combined_set4$meg_trial<-paste(nucleus_combined_set4$trial,"set4",sep="_")
  
  setwd("~/Desktop/Reprogramming/actin_laminac_oct4/combined_data/")
  
  nucleus_combined_set5 <- read.csv("nucleus_combined.csv",stringsAsFactors = F)
  nucleus_combined_set5<-nucleus_combined_set5[!(duplicated(nucleus_combined_set5$Label)),]
  nucleus_combined_set5$meg_trial<-paste(nucleus_combined_set5$trial,"set5",sep="_")
  
  nucleus_combined<-rbind(nucleus_combined_set1,nucleus_combined_set2,nucleus_combined_set3,nucleus_combined_set4,nucleus_combined_set5)
  nucleus_combined<-nucleus_combined[!(duplicated(nucleus_combined$Label)),]
  rownames(nucleus_combined)<-nucleus_combined$Label
  
  rm(nucleus_combined_set1,nucleus_combined_set2,nucleus_combined_set3,nucleus_combined_set4,nucleus_combined_set5)
  
}
#nuc_int_data
{
  
  setwd("~/Desktop/Reprogramming/actin_pmlc_oct4/combined_data/")
  
  actin <- read.csv("cellular_int_actin.csv",stringsAsFactors = F)
  actin<-actin[!(duplicated(actin$Label)),]
  actin<-normalize_mean_trials(actin)
  
  oct4_set1 <- read.csv("nuclear_2d_int_oct4.csv",stringsAsFactors = F)
  oct4_set1<-oct4_set1[!(duplicated(oct4_set1$Label)),]
  oct4_set1<-normalize_2d_mean_trials(oct4_set1)
  
  pmlc <- read.csv("cellular_int_pmlc.csv",stringsAsFactors = F)
  pmlc<-pmlc[!(duplicated(pmlc$Label)),]
  pmlc<-normalize_mean_trials(pmlc)
  
  setwd("~/Desktop/Reprogramming/ki67/combined_data/")
  ki67 <- read.csv("nuclear_int_KI67.csv",stringsAsFactors = F)
  ki67<-ki67[!(duplicated(ki67$Label)),]
  ki67<-normalize_mean_trials(ki67)
  
  setwd("~/Desktop/Reprogramming/H3k9ac/combined_data/")
  h3k9ac <- read.csv("nuclear_2d_int_H3K9AC.csv",stringsAsFactors = F)
  h3k9ac<-h3k9ac[!(duplicated(h3k9ac$Label)),]
  h3k9ac<-normalize_2d_mean_trials(h3k9ac)
  
  setwd("~/Desktop/Reprogramming/actin_laminac_oct4/combined_data/")
  laminac <- read.csv("cellular_int_LAMINAC.csv",stringsAsFactors = F)
  laminac<-laminac[!(duplicated(laminac$Label)),]
  laminac<-normalize_mean_trials(laminac)
  
  oct4_set2 <- read.csv("nuclear_2d_int_oct4.csv",stringsAsFactors = F)
  oct4_set2<-oct4_set2[!(duplicated(oct4_set2$Label)),]
  oct4_set2<-normalize_2d_mean_trials(oct4_set2)
  
  setwd("~/Desktop/Reprogramming/vimentin_ecad_oct4/combined_data/")
  
  vimentin <- read.csv("cellular_int_VIMENTIN.csv",stringsAsFactors = F)
  vimentin<-vimentin[!(duplicated(vimentin$Label)),]
  vimentin<-normalize_mean_trials(vimentin)
  
  ecad <- read.csv("cellular_int_ECAD.csv",stringsAsFactors = F)
  ecad<-ecad[!(duplicated(ecad$Label)),]
  ecad<-normalize_mean_trials(ecad)
  
  oct4_set3 <- read.csv("nuclear_2d_int_oct4.csv",stringsAsFactors = F)
  oct4_set3<-oct4_set3[!(duplicated(oct4_set3$Label)),]
  oct4_set3<-normalize_2d_mean_trials(oct4_set3)
  
  oct4<-rbind(oct4_set1,oct4_set2,oct4_set3)
  
  nuc_int_data<-list(actin=actin,oct4=oct4,pmlc=pmlc,
                     vimentin=vimentin,ecad=ecad,
                     ki67=ki67,h3k9ac=h3k9ac,
                     laminac=laminac)
  rm(actin,oct4,pmlc,vimentin,ecad,ki67,h3k9ac,laminac,oct4_set1,oct4_set2,oct4_set3)
  
}

#parameters
{
  
  
  parameters_2d_int_data<-c("Mean","StdDev","Mode","Min","Max","IntDen","Median","Skew","Kurt","RawIntDen")
  labelling_2d_int_data<-paste("Proj.Spheroid Level",parameters_2d_int_data,sep="\n")
  nuc_labelling_2d_int_data<-paste("Proj.Nuclear Level",parameters_2d_int_data,sep="\n")
  
  parameters_int_data<-c("AtCenter","IntDen","Min","Max","Mean","Sigma","Mode","Mode.NonZero")
  labelling_int_data<-paste("Spheroid Level",parameters_int_data,sep="\n")
  nuclabelling_int_data<-paste("Nuclear Level",parameters_int_data,sep="\n")
  celllabelling_int_data<-paste("Cellular Level",parameters_int_data,sep="\n")
  
  parameters_geo_texture_data<-c( "RatioVolbox","Vol..unit.","Surf..unit." ,"Comp..unit.","Spher..unit.",
                                  "Feret..unit.","Ell_MajRad","Ell_Elon","Ell_Flatness","RatioVolEllipsoid",
                                  "Moment1","Moment2","Moment3","Moment4","Moment5","DCMin..unit.",
                                  "DCMax..unit.","DCMean..unit.","DCSD..unit.",
                                  "Area","Perim.","Major","Minor","Circ.","Feret","MinFeret","AR","Round",
                                  
                                  "HCcontent","ECcontent","HCvolume","ECvolume","HC_EC_content","HC_EC_volume",
                                  "EFD1","EFD2","EFD3","EFD4","EFD5","EFD6","EFD7","EFD8","EFD9","EFD10",                                  
                                  
                                  "Deg_0_Angular_Second_Moment._step5","Deg_0_Contrast._step5","Deg_0_Correlation._step5",
                                  "Deg_0_Inverse_Difference_Moment._step5","Deg_0_Entropy._step5",
                                  "Deg_90_Angular_Second_Moment._step5","Deg_90_Contrast._step5","Deg_90_Correlation._step5",
                                  "Deg_90_Inverse_Difference_Moment._step5","Deg_90_Entropy._step5",
                                  "Deg_180_Angular_Second_Moment._step5","Deg_180_Contrast._step5", "Deg_180_Correlation._step5",
                                  "Deg_180_Inverse_Difference_Moment._step5","Deg_180_Entropy._step5",                   
                                  "Deg_270_Angular_Second_Moment._step5","Deg_270_Contrast._step5","Deg_270_Correlation._step5",
                                  "Deg_270_Inverse_Difference_Moment._step5","Deg_270_Entropy._step5",
                                  
                                  "Deg_0_Angular_Second_Moment._step15","Deg_0_Contrast._step15","Deg_0_Correlation._step15",
                                  "Deg_0_Inverse_Difference_Moment._step15" ,"Deg_0_Entropy._step15",
                                  "Deg_90_Angular_Second_Moment._step15","Deg_90_Contrast._step15","Deg_90_Correlation._step15",
                                  "Deg_90_Inverse_Difference_Moment._step15","Deg_90_Entropy._step15",                    
                                  "Deg_180_Angular_Second_Moment._step15","Deg_180_Contrast._step15","Deg_180_Correlation._step15",              
                                  "Deg_180_Inverse_Difference_Moment._step15","Deg_180_Entropy._step15",
                                  "Deg_270_Angular_Second_Moment._step15","Deg_270_Contrast._step15","Deg_270_Correlation._step15",
                                  "Deg_270_Inverse_Difference_Moment._step15","Deg_270_Entropy._step15",
                                  
                                  "Deg_0_Angular_Second_Moment._step25","Deg_0_Contrast._step25","Deg_0_Correlation._step25",
                                  "Deg_0_Inverse_Difference_Moment._step25","Deg_0_Entropy._step25",
                                  "Deg_90_Angular_Second_Moment._step25","Deg_90_Contrast._step25","Deg_90_Correlation._step25",               
                                  "Deg_90_Inverse_Difference_Moment._step25","Deg_90_Entropy._step25",
                                  "Deg_180_Angular_Second_Moment._step25","Deg_180_Contrast._step25","Deg_180_Correlation._step25",
                                  "Deg_180_Inverse_Difference_Moment._step25","Deg_180_Entropy._step25",
                                  "Deg_270_Angular_Second_Moment._step25","Deg_270_Contrast._step25",               
                                  "Deg_270_Correlation._step25","Deg_270_Inverse_Difference_Moment._step25","Deg_270_Entropy._step25",
                                  
                                  "Deg_0_Angular_Second_Moment._step35","Deg_0_Contrast._step35","Deg_0_Correlation._step35",                
                                  "Deg_0_Inverse_Difference_Moment._step35","Deg_0_Entropy._step35",
                                  "Deg_90_Angular_Second_Moment._step35","Deg_90_Contrast._step35","Deg_90_Correlation._step35",  
                                  "Deg_90_Inverse_Difference_Moment._step35","Deg_90_Entropy._step35",
                                  "Deg_180_Angular_Second_Moment._step35","Deg_180_Contrast._step35","Deg_180_Correlation._step35",              
                                  "Deg_180_Inverse_Difference_Moment._step35","Deg_180_Entropy._step35",
                                  "Deg_270_Angular_Second_Moment._step35","Deg_270_Contrast._step35","Deg_270_Correlation._step35",
                                  "Deg_270_Inverse_Difference_Moment._step35","Deg_270_Entropy._step35",
                                  
                                  "Deg_0_Angular_Second_Moment._step45","Deg_0_Contrast._step45","Deg_0_Correlation._step45",
                                  "Deg_0_Inverse_Difference_Moment._step45","Deg_0_Entropy._step45",
                                  "Deg_90_Angular_Second_Moment._step45","Deg_90_Contrast._step45","Deg_90_Correlation._step45",
                                  "Deg_90_Inverse_Difference_Moment._step45","Deg_90_Entropy._step45",                   
                                  "Deg_180_Angular_Second_Moment._step45","Deg_180_Contrast._step45","Deg_180_Correlation._step45",               
                                  "Deg_180_Inverse_Difference_Moment._step45","Deg_180_Entropy._step45",
                                  "Deg_270_Angular_Second_Moment._step45","Deg_270_Contrast._step45","Deg_270_Correlation._step45",
                                  "Deg_270_Inverse_Difference_Moment._step45","Deg_270_Entropy._step45",
                                  
                                  "Deg_0_Angular_Second_Moment._step100","Deg_0_Contrast._step100","Deg_0_Correlation._step100",
                                  "Deg_0_Inverse_Difference_Moment._step100","Deg_0_Entropy._step100",                                     
                                  "Deg_90_Angular_Second_Moment._step100","Deg_90_Contrast._step100","Deg_90_Correlation._step100",
                                  "Deg_90_Inverse_Difference_Moment._step100" ,"Deg_90_Entropy._step100",
                                  "Deg_180_Angular_Second_Moment._step100","Deg_180_Contrast._step100","Deg_180_Correlation._step100",
                                  "Deg_180_Inverse_Difference_Moment._step100","Deg_180_Entropy._step100",               
                                  "Deg_270_Angular_Second_Moment._step100","Deg_270_Contrast._step100","Deg_270_Correlation._step100",            
                                  "Deg_270_Inverse_Difference_Moment._step100","Deg_270_Entropy._step100")
  
  labelling_geo_texture_data<-c("Bounding Box Ratio","Volume", "Surface Area","Compaction","Sphericity",
                                "Feret (3D)", "Ellipsoid R1","Ellipsoid Ellongation","Ellipsoid Flatness", "Ellipsoid Ratio",
                                "Moment1","Moment2","Moment3","Moment4","Moment5", 
                                "DC_Min.","DC_Max","DC_Mean.","DC_SD",
                                "Proj. Area", "Proj. Perimeter","Proj. Major Axis", "Proj. Minor Axis", "Proj. Circularity",
                                "Proj. Feret","Proj. Min Feret","Proj. A.R","Proj. Roundness","Proj. Solidity",
                                
                                "HCcontent","ECcontent","HCvolume","ECvolume","HC_EC_content","HC_EC_volume",
                                "EFD1","EFD2","EFD3","EFD4","EFD5","EFD6","EFD7","EFD8","EFD9","EFD10",                                  
                                
                                "Deg_0_Angular_Second_Moment _step5","Deg_0_Contrast _step5","Deg_0_Correlation _step5",
                                "Deg_0_Inverse_Difference_Moment _step5","Deg_0_Entropy _step5",
                                "Deg_90_Angular_Second_Moment _step5","Deg_90_Contrast _step5","Deg_90_Correlation _step5",
                                "Deg_90_Inverse_Difference_Moment _step5","Deg_90_Entropy _step5",
                                "Deg_180_Angular_Second_Moment _step5","Deg_180_Contrast _step5", "Deg_180_Correlation _step5",
                                "Deg_180_Inverse_Difference_Moment _step5","Deg_180_Entropy _step5",                   
                                "Deg_270_Angular_Second_Moment _step5","Deg_270_Contrast _step5","Deg_270_Correlation _step5",
                                "Deg_270_Inverse_Difference_Moment _step5","Deg_270_Entropy _step5",
                                
                                "Deg_0_Angular_Second_Moment _step15","Deg_0_Contrast _step15","Deg_0_Correlation _step15",
                                "Deg_0_Inverse_Difference_Moment _step15" ,"Deg_0_Entropy _step15",
                                "Deg_90_Angular_Second_Moment _step15","Deg_90_Contrast _step15","Deg_90_Correlation _step15",
                                "Deg_90_Inverse_Difference_Moment _step15","Deg_90_Entropy _step15",                   
                                "Deg_180_Angular_Second_Moment _step15","Deg_180_Contrast _step15","Deg_180_Correlation _step15",              
                                "Deg_180_Inverse_Difference_Moment _step15","Deg_180_Entropy _step15",
                                "Deg_270_Angular_Second_Moment _step15","Deg_270_Contrast _step15","Deg_270_Correlation _step15",
                                "Deg_270_Inverse_Difference_Moment _step15","Deg_270_Entropy _step15",
                                
                                "Deg_0_Angular_Second_Moment _step25","Deg_0_Contrast _step25","Deg_0_Correlation _step25",
                                "Deg_0_Inverse_Difference_Moment _step25","Deg_0_Entropy _step25",
                                "Deg_90_Angular_Second_Moment _step25","Deg_90_Contrast _step25","Deg_90_Correlation _step25",               
                                "Deg_90_Inverse_Difference_Moment _step25","Deg_90_Entropy _step25",
                                "Deg_180_Angular_Second_Moment _step25","Deg_180_Contrast _step25","Deg_180_Correlation _step25",
                                "Deg_180_Inverse_Difference_Moment _step25","Deg_180_Entropy _step25",
                                "Deg_270_Angular_Second_Moment _step25","Deg_270_Contrast _step25",               
                                "Deg_270_Correlation _step25","Deg_270_Inverse_Difference_Moment _step25","Deg_270_Entropy _step25",
                                
                                "Deg_0_Angular_Second_Moment _step35","Deg_0_Contrast _step35","Deg_0_Correlation _step35",                
                                "Deg_0_Inverse_Difference_Moment _step35","Deg_0_Entropy _step35",
                                "Deg_90_Angular_Second_Moment _step35","Deg_90_Contrast _step35","Deg_90_Correlation _step35",  
                                "Deg_90_Inverse_Difference_Moment _step35","Deg_90_Entropy _step35",
                                "Deg_180_Angular_Second_Moment _step35","Deg_180_Contrast _step35","Deg_180_Correlation _step35",              
                                "Deg_180_Inverse_Difference_Moment _step35","Deg_180_Entropy _step35",
                                "Deg_270_Angular_Second_Moment _step35","Deg_270_Contrast _step35","Deg_270_Correlation _step35",
                                "Deg_270_Inverse_Difference_Moment _step35","Deg_270_Entropy _step35",
                                
                                "Deg_0_Angular_Second_Moment _step45","Deg_0_Contrast _step45","Deg_0_Correlation _step45",
                                "Deg_0_Inverse_Difference_Moment _step45","Deg_0_Entropy _step45",
                                "Deg_90_Angular_Second_Moment _step45","Deg_90_Contrast _step45","Deg_90_Correlation _step45",
                                "Deg_90_Inverse_Difference_Moment _step45","Deg_90_Entropy _step45",                   
                                "Deg_180_Angular_Second_Moment _step45","Deg_180_Contrast _step45","Deg_180_Correlation _step45",               
                                "Deg_180_Inverse_Difference_Moment _step45","Deg_180_Entropy _step45",
                                "Deg_270_Angular_Second_Moment _step45","Deg_270_Contrast _step45","Deg_270_Correlation _step45",
                                "Deg_270_Inverse_Difference_Moment _step45","Deg_270_Entropy _step45",
                                
                                "Deg_0_Angular_Second_Moment _step100","Deg_0_Contrast _step100","Deg_0_Correlation _step100",
                                "Deg_0_Inverse_Difference_Moment _step100","Deg_0_Entropy _step100",                                     
                                "Deg_90_Angular_Second_Moment _step100","Deg_90_Contrast _step100","Deg_90_Correlation _step100",
                                "Deg_90_Inverse_Difference_Moment _step100" ,"Deg_90_Entropy _step100",
                                "Deg_180_Angular_Second_Moment _step100","Deg_180_Contrast _step100","Deg_180_Correlation _step100",
                                "Deg_180_Inverse_Difference_Moment _step100","Deg_180_Entropy _step100",               
                                "Deg_270_Angular_Second_Moment _step100","Deg_270_Contrast _step100","Deg_270_Correlation _step100",            
                                "Deg_270_Inverse_Difference_Moment _step100","Deg_270_Entropy _step100")
  
}
#subset labels
{
  label_T1_set1<-subset(nucleus_combined,nucleus_combined$meg_trial=="T1_set1")$Label
  label_T2_set1<-subset(nucleus_combined,nucleus_combined$meg_trial=="T2_set1")$Label
  label_T3_set1<-subset(nucleus_combined,nucleus_combined$meg_trial=="T3_set1")$Label
  label_T1_set2<-subset(nucleus_combined,nucleus_combined$meg_trial=="T1_set2")$Label
  label_T2_set2<-subset(nucleus_combined,nucleus_combined$meg_trial=="T2_set2")$Label
  label_T3_set2<-subset(nucleus_combined,nucleus_combined$meg_trial=="T3_set2")$Label
  label_T1_set3<-subset(nucleus_combined,nucleus_combined$meg_trial=="T1_set3")$Label
  label_T2_set3<-subset(nucleus_combined,nucleus_combined$meg_trial=="T2_set3")$Label
  label_T3_set3<-subset(nucleus_combined,nucleus_combined$meg_trial=="T3_set3")$Label
  label_T1_set4<-subset(nucleus_combined,nucleus_combined$meg_trial=="T1_set4")$Label
  label_T2_set4<-subset(nucleus_combined,nucleus_combined$meg_trial=="T2_set4")$Label
  label_T3_set4<-subset(nucleus_combined,nucleus_combined$meg_trial=="T3_set4")$Label
  label_T1_set5<-subset(nucleus_combined,nucleus_combined$meg_trial=="T1_set5")$Label
  label_T2_set5<-subset(nucleus_combined,nucleus_combined$meg_trial=="T2_set5")$Label
  label_T3_set5<-subset(nucleus_combined,nucleus_combined$meg_trial=="T3_set5")$Label
  label_T4_set5<-subset(nucleus_combined,nucleus_combined$meg_trial=="T4_set5")$Label
  
  
  D0_lab<-subset(nucleus_combined,nucleus_combined$sample=="D0")$Label
  D2_lab<-subset(nucleus_combined,nucleus_combined$sample=="D2")$Label
  D4_lab<-subset(nucleus_combined,nucleus_combined$sample=="D4")$Label
  D6_lab<-subset(nucleus_combined,nucleus_combined$sample=="D6")$Label
  D8_lab<-subset(nucleus_combined,nucleus_combined$sample=="D8")$Label
  
}
cols_timeline<-colorRampPalette(c("red","yellow"))( 5) 
cols_trials<-colorRampPalette(c("gray","black")) (length(levels(as.factor(nucleus_combined$meg_trial))))

dir.create("/Users/saradhavenkatachalapathy/Desktop/reprogramming/timeline_nucleus_analysis_scaled/")
setwd("/Users/saradhavenkatachalapathy/Desktop/reprogramming/timeline_nucleus_analysis_scaled/")
dird<-"/Users/saradhavenkatachalapathy/Desktop/reprogramming/timeline_nucleus_analysis_scaled/"

dir.create("/Users/saradhavenkatachalapathy/Desktop/reprogramming/timeline_nucleus_analysis_scaled/equal_sampling/")
setwd("/Users/saradhavenkatachalapathy/Desktop/reprogramming/timeline_nucleus_analysis_scaled/equal_sampling/")
dird<-"/Users/saradhavenkatachalapathy/Desktop/reprogramming/timeline_nucleus_analysis_scaled/equal_sampling/"
#LDA with time
{
  
  #pca scaled acc to exp
  {
    nucleus_subset<-nucleus_combined[,which(colnames(nucleus_combined) %in% parameters_geo_texture_data)]
    
    scaled_nucleus_subset<-subset(nucleus_subset,rownames(nucleus_subset) %in% label_T1_set1)
    scaled_nucleus_subset<-scale(scaled_nucleus_subset,center = T, scale = T)
    temp2<-subset(nucleus_subset,rownames(nucleus_subset) %in% label_T2_set1)
    temp2<-scale(temp2,center = T, scale = T)
    scaled_nucleus_subset<-rbind(scaled_nucleus_subset,temp2)
    temp2<-subset(nucleus_subset,rownames(nucleus_subset) %in% label_T3_set1)
    temp2<-scale(temp2,center = T, scale = T)
    scaled_nucleus_subset<-rbind(scaled_nucleus_subset,temp2)
    temp2<-subset(nucleus_subset,rownames(nucleus_subset) %in% label_T1_set2)
    temp2<-scale(temp2,center = T, scale = T)
    scaled_nucleus_subset<-rbind(scaled_nucleus_subset,temp2)
    temp2<-subset(nucleus_subset,rownames(nucleus_subset) %in% label_T2_set2)
    temp2<-scale(temp2,center = T, scale = T)
    scaled_nucleus_subset<-rbind(scaled_nucleus_subset,temp2)
    temp2<-subset(nucleus_subset,rownames(nucleus_subset) %in% label_T3_set2)
    temp2<-scale(temp2,center = T, scale = T)
    scaled_nucleus_subset<-rbind(scaled_nucleus_subset,temp2)
    temp2<-subset(nucleus_subset,rownames(nucleus_subset) %in% label_T1_set3)
    temp2<-scale(temp2,center = T, scale = T)
    scaled_nucleus_subset<-rbind(scaled_nucleus_subset,temp2)
    temp2<-subset(nucleus_subset,rownames(nucleus_subset) %in% label_T2_set3)
    temp2<-scale(temp2,center = T, scale = T)
    scaled_nucleus_subset<-rbind(scaled_nucleus_subset,temp2)
    temp2<-subset(nucleus_subset,rownames(nucleus_subset) %in% label_T3_set3)
    temp2<-scale(temp2,center = T, scale = T)
    scaled_nucleus_subset<-rbind(scaled_nucleus_subset,temp2)
    temp2<-subset(nucleus_subset,rownames(nucleus_subset) %in% label_T1_set4)
    temp2<-scale(temp2,center = T, scale = T)
    scaled_nucleus_subset<-rbind(scaled_nucleus_subset,temp2)
    temp2<-subset(nucleus_subset,rownames(nucleus_subset) %in% label_T2_set4)
    temp2<-scale(temp2,center = T, scale = T)
    scaled_nucleus_subset<-rbind(scaled_nucleus_subset,temp2)
    temp2<-subset(nucleus_subset,rownames(nucleus_subset) %in% label_T3_set4)
    temp2<-scale(temp2,center = T, scale = T)
    scaled_nucleus_subset<-rbind(scaled_nucleus_subset,temp2)
    temp2<-subset(nucleus_subset,rownames(nucleus_subset) %in% label_T1_set5)
    temp2<-scale(temp2,center = T, scale = T)
    scaled_nucleus_subset<-rbind(scaled_nucleus_subset,temp2)
    temp2<-subset(nucleus_subset,rownames(nucleus_subset) %in% label_T2_set5)
    temp2<-scale(temp2,center = T, scale = T)
    scaled_nucleus_subset<-rbind(scaled_nucleus_subset,temp2)
    temp2<-subset(nucleus_subset,rownames(nucleus_subset) %in% label_T3_set5)
    temp2<-scale(temp2,center = T, scale = T)
    scaled_nucleus_subset<-rbind(scaled_nucleus_subset,temp2)
    temp2<-subset(nucleus_subset,rownames(nucleus_subset) %in% label_T4_set5)
    temp2<-scale(temp2,center = T, scale = T)
    scaled_nucleus_subset<-rbind(scaled_nucleus_subset,temp2)
    
    nucleus_subset<-as.matrix(na.omit((scaled_nucleus_subset)))
    
    re<-as.data.frame(nucleus_subset)
    d1<-subset(re, rownames(re) %in% D0_lab)
    d2<-subset(re, rownames(re) %in% D2_lab)
    d3<-subset(re, rownames(re) %in% D4_lab)
    d4<-subset(re, rownames(re) %in% D6_lab)
    d5<-subset(re, rownames(re) %in% D8_lab)
    
    d1$samplem<-"D0"
    d2$samplem<-"D2"
    d3$samplem<-"D4"
    d4$samplem<-"D6"
    d5$samplem<-"D8"
    
    re<-rbind(d1, d2)
    re<-rbind(re, d3)
    re<-rbind(re, d4)
    re<-rbind(re, d5)
    
    re$Label<-row.names(re)
    lengtha<-min(nrow(d1),nrow(d2),nrow(d3),nrow(d4),nrow(d5))
    
    trainingData<-ddply(re,.(samplem),function(x) x[sample(nrow(x),round(lengtha,0)),])
    
    rownames(trainingData)<-trainingData$Label
    nucleus_subset<-trainingData[,which(colnames(trainingData) %in% parameters_geo_texture_data)]
    nucleus_pca <- princomp(nucleus_subset,  cor = TRUE, scores = TRUE)
    pca_scores<-nucleus_pca$scores
    
    a<-subset(pca_scores,rownames(pca_scores)%in% D0_lab)
    b<-subset(pca_scores,rownames(pca_scores)%in% D2_lab)
    c<-subset(pca_scores,rownames(pca_scores)%in% D4_lab)
    d<-subset(pca_scores,rownames(pca_scores)%in% D6_lab)
    e<-subset(pca_scores,rownames(pca_scores)%in% D8_lab)
    
    png(filename="pca_samples.png", units="in", width=2, height=2 , pointsize=7, res=1200)
    par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2)
    plot(pca_scores[,1],pca_scores[,2],col=NA,xlab="PC1",ylab="PC2", las=1)
    points(a[,1],a[,2],col=cols_timeline[1],pch=19,cex=0.3)
    points(b[,1],b[,2],col=cols_timeline[2],pch=19,cex=0.3)
    points(c[,1],c[,2],col=cols_timeline[3],pch=19,cex=0.3)
    points(d[,1],d[,2],col=cols_timeline[4],pch=19,cex=0.3)
    points(e[,1],e[,2],col=cols_timeline[5],pch=19,cex=0.3)
    dev.off()
    
    png(filename="pca_samplesendpoints.png", units="in", width=2, height=2 , pointsize=7, res=1200)
    par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2)
    plot(pca_scores[,1],pca_scores[,2],col=NA,xlab="PC1",ylab="PC2", las=1)
    points(e[,1],e[,2],col=cols_timeline[5],pch=19,cex=0.3)
    points(a[,1],a[,2],col=cols_timeline[1],pch=19,cex=0.3)
    dev.off()
    
    rm(a,b,c,d,e,d1,d2,d3,d4,d5,nucleus_subset,re,scaled_nucleus_subset,trainingData,temp2)
    
  }
  
  #LDA
  {
    re<-as.data.frame(pca_scores)
    d1<-subset(re, rownames(re) %in% D0_lab)
    d2<-subset(re, rownames(re) %in% D2_lab)
    d3<-subset(re, rownames(re) %in% D4_lab)
    d4<-subset(re, rownames(re) %in% D6_lab)
    d5<-subset(re, rownames(re) %in% D8_lab)
    d1$samplem<-"D0"
    d2$samplem<-"D2"
    d3$samplem<-"D4"
    d4$samplem<-"D6"
    d5$samplem<-"D8"
    
    re<-rbind(d1, d2)
    re<-rbind(re, d3)
    re<-rbind(re, d4)
    re<-rbind(re, d5)
    
    lda_a <- lda(samplem ~. , data = re)
    
    lda_val<-as.data.frame(matrix(nrow=nrow(re),ncol=3))
    lda_val[,1]<-calclda(re[,1:3], lda_a$scaling[,1])
    lda_val[,2]<-calclda(re[,1:3], lda_a$scaling[,2])
    lda_val[,3]<-calclda(re[,1:3], lda_a$scaling[,3])
    lda_val[,4]<-calclda(re[,1:3], lda_a$scaling[,4])
    lda_val[,5]<-re$samplem
    colorset = cols_timeline[as.factor(lda_val[,5])]
    
    png(filename="LDA_whole_1_2.png", units="in",width=2, height=2 , pointsize=7, res=1200)
    par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2)
    plot(lda_val[,1],lda_val[,2],col=colorset,pch=19, ylab="LD2",xlab="LD1",las=1,cex=0.2)
    dev.off()
    
    lda_predict <- predict(lda_a, newdata = re)
    lda_pred_class_train<-as.data.frame(matrix(nrow=nrow(re),ncol=2))
    lda_pred_class_train[,1]<- lda_predict$class
    lda_pred_class_train[,2]<- re$samplem
    colnames(lda_pred_class_train)<-c("predicted","actual")
    x<-table(lda_pred_class_train)  
    
    png(filename="LDA_whole_pred.png", units="in",width=2, height=2 , pointsize=5, res=1200)
    par(font.axis = 1,font.lab=1,mar=c(1,1,1,2), font=2)
    corrplot(x,tl.col = "black",method = "color",is.corr=F, cl.pos = "b",tl.cex = 0.8, cl.ratio = 0.3)
    dev.off()
    write.csv(x,file="confusion_matrix_samples_lda.csv")
    
    #ellipse plot
    {
      png(filename="lda_biplot_1_2_ellipses_samples.png", units="in",width=2, height=2 , pointsize=7, res=1200)
      par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2)
      ellipseplot(lda_val[,1],                          # data for x-axis
                  lda_val[,2],                          # data for y-axis
                  as.factor(lda_val[,5]),                    # factor with classes
                  elev=0.95,
                  pcol=cols_timeline,    # colors for plotting (must match # of factors)
                  pbgcol=TRUE,                           # point borders black?
                  cexsize=0.1,                            # size of points 
                  ppch=rep(21,5),                          # shape of points (must match # of factors)
                  axissize=0.8,                           # Set axis text size
                  linewidth=1,                          # Set axis line size
                  font=1
      )
      title(xlab="LD1",    # % variance explained on PC1
            ylab="LD2",    # % variance explained on PC2 
            main=" ",                  # Title
            cex.lab=1,                    # size of label text
            cex.main=1                    # size of title text
      )
      dev.off()
      
      png(filename="lda_biplot_1_3_ellipses_samples.png", units="in",width=2, height=2 , pointsize=7, res=1200)
      par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2)
      ellipseplot(lda_val[,1],                          # data for x-axis
                  lda_val[,3],                          # data for y-axis
                  as.factor(lda_val[,5]),                    # factor with classes
                  elev=0.95,
                  pcol=cols_timeline,    # colors for plotting (must match # of factors)
                  pbgcol=TRUE,                           # point borders black?
                  cexsize=0.1,                            # size of points 
                  ppch=rep(21,5),                          # shape of points (must match # of factors)
                  axissize=0.8,                           # Set axis text size
                  linewidth=1,                          # Set axis line size
                  font=1
      )
      title(xlab="LD1",    # % variance explained on PC1
            ylab="LD3",    # % variance explained on PC2 
            main=" ",                  # Title
            cex.lab=1,                    # size of label text
            cex.main=1                    # size of title text
      )
      dev.off()
      
      png(filename="lda_biplot_1_4_ellipses_samples.png", units="in",width=2, height=2 , pointsize=7, res=1200)
      par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2)
      ellipseplot(lda_val[,1],                          # data for x-axis
                  lda_val[,4],                          # data for y-axis
                  as.factor(lda_val[,5]),                    # factor with classes
                  elev=0.95,
                  pcol=cols_timeline,    # colors for plotting (must match # of factors)
                  pbgcol=TRUE,                           # point borders black?
                  cexsize=0.1,                            # size of points 
                  ppch=rep(21,5),                          # shape of points (must match # of factors)
                  axissize=0.8,                           # Set axis text size
                  linewidth=1,                          # Set axis line size
                  font=1
      )
      title(xlab="LD1",    # % variance explained on PC1
            ylab="LD4",    # % variance explained on PC2 
            main=" ",                  # Title
            cex.lab=1,                    # size of label text
            cex.main=1                    # size of title text
      )
      dev.off()
      
      png(filename="lda_biplot_2_3_ellipses_samples.png", units="in",width=2, height=2 , pointsize=7, res=1200)
      par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2)
      ellipseplot(lda_val[,2],                          # data for x-axis
                  lda_val[,3],                          # data for y-axis
                  as.factor(lda_val[,5]),                    # factor with classes
                  elev=0.95,
                  pcol=cols_timeline,    # colors for plotting (must match # of factors)
                  pbgcol=TRUE,                           # point borders black?
                  cexsize=0.1,                            # size of points 
                  ppch=rep(21,5),                          # shape of points (must match # of factors)
                  axissize=0.8,                           # Set axis text size
                  linewidth=1,                          # Set axis line size
                  font=1
      )
      title(xlab="LD2",    # % variance explained on PC1
            ylab="LD3",    # % variance explained on PC2 
            main=" ",                  # Title
            cex.lab=1,                    # size of label text
            cex.main=1                    # size of title text
      )
      dev.off()
      
      png(filename="lda_biplot_2_4_ellipses_samples.png", units="in",width=2, height=2 , pointsize=7, res=1200)
      par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2)
      ellipseplot(lda_val[,2],                          # data for x-axis
                  lda_val[,4],                          # data for y-axis
                  as.factor(lda_val[,5]),                    # factor with classes
                  elev=0.95,
                  pcol=cols_timeline,    # colors for plotting (must match # of factors)
                  pbgcol=TRUE,                           # point borders black?
                  cexsize=0.1,                            # size of points 
                  ppch=rep(21,5),                          # shape of points (must match # of factors)
                  axissize=0.8,                           # Set axis text size
                  linewidth=1,                          # Set axis line size
                  font=1
      )
      title(xlab="LD2",    # % variance explained on PC1
            ylab="LD4",    # % variance explained on PC2 
            main=" ",                  # Title
            cex.lab=1,                    # size of label text
            cex.main=1                    # size of title text
      )
      dev.off()
      
      png(filename="lda_biplot_3_4_ellipses_samples.png", units="in",width=2, height=2 , pointsize=7, res=1200)
      par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2)
      ellipseplot(lda_val[,3],                          # data for x-axis
                  lda_val[,4],                          # data for y-axis
                  as.factor(lda_val[,5]),                    # factor with classes
                  elev=0.95,
                  pcol=cols_timeline,    # colors for plotting (must match # of factors)
                  pbgcol=TRUE,                           # point borders black?
                  cexsize=0.1,                            # size of points 
                  ppch=rep(21,5),                          # shape of points (must match # of factors)
                  axissize=0.8,                           # Set axis text size
                  linewidth=1,                          # Set axis line size
                  font=1
      )
      title(xlab="LD3",    # % variance explained on PC1
            ylab="LD4",    # % variance explained on PC2 
            main=" ",                  # Title
            cex.lab=1,                    # size of label text
            cex.main=1                    # size of title text
      )
      dev.off()
    }
    #area of intersection
    {
      #1vs2
      {
        edf <- data.frame(LV1 = lda_val[,1], LV2=lda_val[,2], factr = lda_val[,5]) # create data frame with data and factor
        ellipses <- dlply(edf, .(factr), function(x) {
          LV1 <- x[,1]
          LV2 <- x[,2]
          dataEllipse(LV1, LV2, levels=0.95, robust=TRUE, draw=FALSE) # Get confidence ellipse points from dataEllipse() function by factor level
        })
        col1 <- colorRampPalette(c("#7F0000", "#7F0000","#7F0000","#7F0000","#7F0000","#7F0000","#7F0000","#7F0000","#7F0000"
                                   ,"red", "#FF7F00", "yellow", "#7FFF7F",
                                   "cyan", "#007FFF", "blue", "#00007F"))
        
        area_int<-as.data.frame(matrix(NA, nrow=5,ncol=5))
        colnames(area_int)<-rev(attributes(ellipses)$names)
        rownames(area_int)<-(attributes(ellipses)$names)
        for(k in 1:5){
          for(o in 5:1){
            ell1 <- xy2SP(ellipses[[k]])
            ell2 <- xy2SP(ellipses[[6-o]])
            a<-gIntersection(ell1, ell2)
            if(length(a)>0){
              area_int[k,o]<-gArea(gIntersection(ell1, ell2))/(gArea(ell1)+gArea(ell2)-gArea(gIntersection(ell1, ell2)))
            }
            if(length(a)==0){
              area_int[k,o]<-0
            }
            
          }
        }
        png(filename="corrplot_LD1_1_2_overlap_samples.png", units="in",width=2, height=2 , pointsize=5, res=1200)
        par(font.axis = 1,font.lab=1,mar=c(1,1,1,2), font=2)
        corrplot(t(as.matrix((area_int))),tl.col = "black",cl.lim = c(0, 1),method = "color",is.corr=F, col=col1(250),cl.pos = "b",tl.cex = 0.8, cl.ratio = 0.3)
        dev.off()
        
        write.csv(area_int,file="area_intersection_1_2.csv")
        
      }
      #1vs3
      {
        edf <- data.frame(LV1 = lda_val[,1], LV2=lda_val[,3], factr = lda_val[,5]) # create data frame with data and factor
        ellipses <- dlply(edf, .(factr), function(x) {
          LV1 <- x[,1]
          LV2 <- x[,2]
          dataEllipse(LV1, LV2, levels=0.95, robust=TRUE, draw=FALSE) # Get confidence ellipse points from dataEllipse() function by factor level
        })
        col1 <- colorRampPalette(c("#7F0000", "#7F0000","#7F0000","#7F0000","#7F0000","#7F0000","#7F0000","#7F0000","#7F0000"
                                   ,"red", "#FF7F00", "yellow", "#7FFF7F",
                                   "cyan", "#007FFF", "blue", "#00007F"))
        
        area_int<-as.data.frame(matrix(NA, nrow=5,ncol=5))
        colnames(area_int)<-rev(attributes(ellipses)$names)
        rownames(area_int)<-(attributes(ellipses)$names)
        for(k in 1:5){
          for(o in 5:1){
            ell1 <- xy2SP(ellipses[[k]])
            ell2 <- xy2SP(ellipses[[6-o]])
            a<-gIntersection(ell1, ell2)
            if(length(a)>0){
              area_int[k,o]<-gArea(gIntersection(ell1, ell2))/(gArea(ell1)+gArea(ell2)-gArea(gIntersection(ell1, ell2)))
            }
            if(length(a)==0){
              area_int[k,o]<-0
            }
            
          }
        }
        png(filename="corrplot_LD1_1_3_overlap_samples.png", units="in",width=2, height=2 , pointsize=5, res=1200)
        par(font.axis = 1,font.lab=1,mar=c(1,1,1,2), font=2)
        corrplot(t(as.matrix((area_int))),tl.col = "black",cl.lim = c(0, 1),method = "color",is.corr=F, col=col1(250),cl.pos = "b",tl.cex = 0.8, cl.ratio = 0.3)
        dev.off()
        
        write.csv(area_int,file="area_intersection_1_3.csv")
        
      }
      #1vs4
      {
        edf <- data.frame(LV1 = lda_val[,1], LV2=lda_val[,4], factr = lda_val[,5]) # create data frame with data and factor
        ellipses <- dlply(edf, .(factr), function(x) {
          LV1 <- x[,1]
          LV2 <- x[,2]
          dataEllipse(LV1, LV2, levels=0.95, robust=TRUE, draw=FALSE) # Get confidence ellipse points from dataEllipse() function by factor level
        })
        col1 <- colorRampPalette(c("#7F0000", "#7F0000","#7F0000","#7F0000","#7F0000","#7F0000","#7F0000","#7F0000","#7F0000"
                                   ,"red", "#FF7F00", "yellow", "#7FFF7F",
                                   "cyan", "#007FFF", "blue", "#00007F"))
        
        area_int<-as.data.frame(matrix(NA, nrow=5,ncol=5))
        colnames(area_int)<-rev(attributes(ellipses)$names)
        rownames(area_int)<-(attributes(ellipses)$names)
        for(k in 1:5){
          for(o in 5:1){
            ell1 <- xy2SP(ellipses[[k]])
            ell2 <- xy2SP(ellipses[[6-o]])
            a<-gIntersection(ell1, ell2)
            if(length(a)>0){
              area_int[k,o]<-gArea(gIntersection(ell1, ell2))/(gArea(ell1)+gArea(ell2)-gArea(gIntersection(ell1, ell2)))
            }
            if(length(a)==0){
              area_int[k,o]<-0
            }
            
          }
        }
        png(filename="corrplot_LD1_1_4_overlap_samples.png", units="in",width=2, height=2 , pointsize=5, res=1200)
        par(font.axis = 1,font.lab=1,mar=c(1,1,1,2), font=2)
        corrplot(t(as.matrix((area_int))),tl.col = "black",cl.lim = c(0, 1),method = "color",is.corr=F, col=col1(250),cl.pos = "b",tl.cex = 0.8, cl.ratio = 0.3)
        dev.off()
        
        write.csv(area_int,file="area_intersection_1_4.csv")
        
      }
      
      #2vs3
      {
        edf <- data.frame(LV1 = lda_val[,2], LV2=lda_val[,3], factr = lda_val[,5]) # create data frame with data and factor
        ellipses <- dlply(edf, .(factr), function(x) {
          LV1 <- x[,1]
          LV2 <- x[,2]
          dataEllipse(LV1, LV2, levels=0.95, robust=TRUE, draw=FALSE) # Get confidence ellipse points from dataEllipse() function by factor level
        })
        col1 <- colorRampPalette(c("#7F0000", "#7F0000","#7F0000","#7F0000","#7F0000","#7F0000","#7F0000","#7F0000","#7F0000"
                                   ,"red", "#FF7F00", "yellow", "#7FFF7F",
                                   "cyan", "#007FFF", "blue", "#00007F"))
        
        area_int<-as.data.frame(matrix(NA, nrow=5,ncol=5))
        colnames(area_int)<-rev(attributes(ellipses)$names)
        rownames(area_int)<-(attributes(ellipses)$names)
        for(k in 1:5){
          for(o in 5:1){
            ell1 <- xy2SP(ellipses[[k]])
            ell2 <- xy2SP(ellipses[[6-o]])
            a<-gIntersection(ell1, ell2)
            if(length(a)>0){
              area_int[k,o]<-gArea(gIntersection(ell1, ell2))/(gArea(ell1)+gArea(ell2)-gArea(gIntersection(ell1, ell2)))
            }
            if(length(a)==0){
              area_int[k,o]<-0
            }
            
          }
        }
        png(filename="corrplot_LD1_2_3_overlap_samples.png", units="in",width=2, height=2 , pointsize=5, res=1200)
        par(font.axis = 1,font.lab=1,mar=c(1,1,1,2), font=2)
        corrplot(t(as.matrix((area_int))),tl.col = "black",cl.lim = c(0, 1),method = "color",is.corr=F, col=col1(250),cl.pos = "b",tl.cex = 0.8, cl.ratio = 0.3)
        dev.off()
        
        write.csv(area_int,file="area_intersection_2_3.csv")
        
      }
      #2vs4
      {
        edf <- data.frame(LV1 = lda_val[,2], LV2=lda_val[,4], factr = lda_val[,5]) # create data frame with data and factor
        ellipses <- dlply(edf, .(factr), function(x) {
          LV1 <- x[,1]
          LV2 <- x[,2]
          dataEllipse(LV1, LV2, levels=0.95, robust=TRUE, draw=FALSE) # Get confidence ellipse points from dataEllipse() function by factor level
        })
        col1 <- colorRampPalette(c("#7F0000", "#7F0000","#7F0000","#7F0000","#7F0000","#7F0000","#7F0000","#7F0000","#7F0000"
                                   ,"red", "#FF7F00", "yellow", "#7FFF7F",
                                   "cyan", "#007FFF", "blue", "#00007F"))
        
        area_int<-as.data.frame(matrix(NA, nrow=5,ncol=5))
        colnames(area_int)<-rev(attributes(ellipses)$names)
        rownames(area_int)<-(attributes(ellipses)$names)
        for(k in 1:5){
          for(o in 5:1){
            ell1 <- xy2SP(ellipses[[k]])
            ell2 <- xy2SP(ellipses[[6-o]])
            a<-gIntersection(ell1, ell2)
            if(length(a)>0){
              area_int[k,o]<-gArea(gIntersection(ell1, ell2))/(gArea(ell1)+gArea(ell2)-gArea(gIntersection(ell1, ell2)))
            }
            if(length(a)==0){
              area_int[k,o]<-0
            }
            
          }
        }
        png(filename="corrplot_LD1_2_4_overlap_samples.png", units="in",width=2, height=2 , pointsize=5, res=1200)
        par(font.axis = 1,font.lab=1,mar=c(1,1,1,2), font=2)
        corrplot(t(as.matrix((area_int))),tl.col = "black",cl.lim = c(0, 1),method = "color",is.corr=F, col=col1(250),cl.pos = "b",tl.cex = 0.8, cl.ratio = 0.3)
        dev.off()
        
        write.csv(area_int,file="area_intersection_2_4.csv")
        
      }
      
      #3vs4
      {
        edf <- data.frame(LV1 = lda_val[,3], LV2=lda_val[,4], factr = lda_val[,5]) # create data frame with data and factor
        ellipses <- dlply(edf, .(factr), function(x) {
          LV1 <- x[,1]
          LV2 <- x[,2]
          dataEllipse(LV1, LV2, levels=0.95, robust=TRUE, draw=FALSE) # Get confidence ellipse points from dataEllipse() function by factor level
        })
        col1 <- colorRampPalette(c("#7F0000", "#7F0000","#7F0000","#7F0000","#7F0000","#7F0000","#7F0000","#7F0000","#7F0000"
                                   ,"red", "#FF7F00", "yellow", "#7FFF7F",
                                   "cyan", "#007FFF", "blue", "#00007F"))
        
        area_int<-as.data.frame(matrix(NA, nrow=5,ncol=5))
        colnames(area_int)<-rev(attributes(ellipses)$names)
        rownames(area_int)<-(attributes(ellipses)$names)
        for(k in 1:5){
          for(o in 5:1){
            ell1 <- xy2SP(ellipses[[k]])
            ell2 <- xy2SP(ellipses[[6-o]])
            a<-gIntersection(ell1, ell2)
            if(length(a)>0){
              area_int[k,o]<-gArea(gIntersection(ell1, ell2))/(gArea(ell1)+gArea(ell2)-gArea(gIntersection(ell1, ell2)))
            }
            if(length(a)==0){
              area_int[k,o]<-0
            }
            
          }
        }
        png(filename="corrplot_LD1_3_4_overlap_samples.png", units="in",width=2, height=2 , pointsize=5, res=1200)
        par(font.axis = 1,font.lab=1,mar=c(1,1,1,2), font=2)
        corrplot(t(as.matrix((area_int))),tl.col = "black",cl.lim = c(0, 1),method = "color",is.corr=F, col=col1(250),cl.pos = "b",tl.cex = 0.8, cl.ratio = 0.3)
        dev.off()
        
        write.csv(area_int,file="area_intersection_3_4.csv")
        
      }
    }
    
    lda_val$Label<-rownames(re)
    colnames(lda_val)<-c("LD1","LD2","LD3","LD4","sample","Label")
    rm(re,x,edf,k,lengtha,o,ellipses,col1,a,area_int,d1,d2,d3,d4,d5,ell1,ell2,lda_predict,lda_pred_class_train)
  }
  
  #nuclear features that correlate with LD1
  {
    merged_data<-merge(lda_val,nucleus_combined[,which(colnames(nucleus_combined) %in% c(parameters_geo_texture_data,"Label"))], by="Label")
    
    a<-subset(lda_val,lda_val$sample=="D0")
    b<-subset(lda_val,lda_val$sample=="D2")
    c<-subset(lda_val,lda_val$sample=="D4")
    d<-subset(lda_val,lda_val$sample=="D6")
    e<-subset(lda_val,lda_val$sample=="D8")
    
    library(gplots)
    cor_with_ld1=vector()
    for(i in 1:length(parameters_geo_texture_data)){
      cor_with_ld1[i]<-cor(merged_data$LD1,merged_data[,which(colnames(merged_data)==parameters_geo_texture_data[i])], 
                           method="spearman")
      
    }
    names(cor_with_ld1)<-parameters_geo_texture_data
    cor_with_ld2=vector()
    for(i in 1:length(parameters_geo_texture_data)){
      cor_with_ld2[i]<-cor(merged_data$LD2,merged_data[,which(colnames(merged_data)==parameters_geo_texture_data[i])], 
                           method="spearman")
      
    }
    names(cor_with_ld2)<-parameters_geo_texture_data
    cor_with_ld3=vector()
    for(i in 1:length(parameters_geo_texture_data)){
      cor_with_ld3[i]<-cor(merged_data$LD3,merged_data[,which(colnames(merged_data)==parameters_geo_texture_data[i])], 
                           method="spearman")
      
    }
    names(cor_with_ld3)<-parameters_geo_texture_data
    cor_with_ld4=vector()
    for(i in 1:length(parameters_geo_texture_data)){
      cor_with_ld4[i]<-cor(merged_data$LD3,merged_data[,which(colnames(merged_data)==parameters_geo_texture_data[i])], 
                           method="spearman")
      
    }
    names(cor_with_ld4)<-parameters_geo_texture_data
    
    merged_data$Label->rownames(merged_data)
    merged_data_subset<- merged_data[, which(colnames(merged_data) %in% parameters_geo_texture_data)] 
    
    a<-apply(subset(merged_data_subset, rownames(merged_data_subset) %in% D0_lab),2,mean) 
    b<-apply(subset(merged_data_subset, rownames(merged_data_subset) %in% D2_lab),2,mean) 
    c<-apply(subset(merged_data_subset, rownames(merged_data_subset) %in% D4_lab),2,mean) 
    d<-apply(subset(merged_data_subset, rownames(merged_data_subset) %in% D6_lab),2,mean) 
    e<-apply(subset(merged_data_subset, rownames(merged_data_subset) %in% D8_lab),2,mean) 
    
    summary_nuclear_feature<-as.data.frame(t(rbind(a,b,c,d,e)))
    summary_nuclear_feature$Label<-rownames(summary_nuclear_feature)
    cor_with_ld1<-as.data.frame(cor_with_ld1)
    cor_with_ld1$Label<-rownames(cor_with_ld1)
    cor_with_ld2<-as.data.frame(cor_with_ld2)
    cor_with_ld2$Label<-rownames(cor_with_ld2)
    cor_with_ld3<-as.data.frame(cor_with_ld3)
    cor_with_ld3$Label<-rownames(cor_with_ld3)
    cor_with_ld4<-as.data.frame(cor_with_ld4)
    cor_with_ld4$Label<-rownames(cor_with_ld4)
    
    summary_nuclear_feature_merged<-merge(cor_with_ld1,cor_with_ld2, by="Label")
    summary_nuclear_feature_merged<-merge(summary_nuclear_feature_merged,cor_with_ld3, by="Label")
    summary_nuclear_feature_merged<-merge(summary_nuclear_feature_merged,cor_with_ld4, by="Label")
    summary_nuclear_feature_merged<-merge(summary_nuclear_feature_merged,summary_nuclear_feature, by="Label")
    
    colnames(summary_nuclear_feature_merged)<-c("Label","cor_LD1","cor_LD2","cor_LD3","cor_Ld4","D0","D2","D4","D6","D8")
    
    library(plot3D)
    rbPal <- colorRampPalette(c('papayawhip','orchid'))
    summary_nuclear_feature_merged[,6:10]<-(summary_nuclear_feature_merged[,6:10]-apply(summary_nuclear_feature_merged[,6:10],1,mean) )/apply(summary_nuclear_feature_merged[,6:10],1,sd)
    
    cor_cols<-redblue(50)[as.numeric(cut(summary_nuclear_feature_merged$cor_LD1,breaks = 50))]
    cor_cols<-cor_cols[order(summary_nuclear_feature_merged[,10])]
    
    png(filename="nuc_mean_summary_key.png", units="in", width=2, height=2 , pointsize=5, res=1200)
    par(font.axis = 2,font.lab=2,mar=c(4,4,2,2), font=2)
    colkey(col=redblue(50),range(cor_with_ld1$cor_with_ld1),side = 1, add = F,length = 0.6, width = 1,clab = "Correlation with\n LD1",side.clab=1)
    dev.off()
    
    
    summary_nuclear_feature_merged<-summary_nuclear_feature_merged[order(summary_nuclear_feature_merged$cor_LD1),]
    cor_cols<-redblue(100)[as.numeric(cut(summary_nuclear_feature_merged$cor_LD1,breaks = 100))]
    
    png(filename="nuc_mean_summary.png", units="in", width=4, height=8.5 , pointsize=5, res=1200)
    par(font.axis = 2,font.lab=2,mar=c(4,4,2,2), font=2)
    heatmap.2(as.matrix(summary_nuclear_feature_merged[,6:10]), dendrogram ="none",Rowv = NA,Colv = NA,trace="none",density.info="none",srtCol = 0,cexCol = 1.5,
              key.title = "",cexRow = 1,col=rbPal(50),scale = "none",RowSideColors = cor_cols,labRow =summary_nuclear_feature_merged$Label,margins = c(2, 20) )
    dev.off()
    
    summary_nuclear_feature_merged<-summary_nuclear_feature_merged[rev(order(abs(summary_nuclear_feature_merged$cor_LD1))),]
    cor_cols<-redblue(100)[as.numeric(cut(summary_nuclear_feature_merged$cor_LD1[1:20],breaks = 100))]
    
    png(filename="nuc_mean_summary_top20.png", units="in", width=2, height=2 , pointsize=5, res=1200)
    par(font.axis = 2,font.lab=2,mar=c(4,4,2,5), font=2)
    heatmap.2(as.matrix(summary_nuclear_feature_merged[1:20,6:10]), dendrogram ="none",Rowv = NA,Colv = NA,trace="none",density.info="none",srtCol = 0,cexCol = 1.5,
              key.title = "",cexRow = 1,col=rbPal(50),scale = "none",RowSideColors = cor_cols,labRow =summary_nuclear_feature_merged$Label,margins = c(3, 15))
    dev.off()
    
    rm(cor_with_ld1,cor_with_ld2,cor_with_ld3,cor_with_ld4,lda_a,lda_val,merged_data_subset,summary_nuclear_feature,a,b,c,d,e,i)
    
  }
  
}

#oct4deciles
{
  library(plot3D)
  
  #deciles
  {
    oct4_int_selec<-subset(nuc_int_data$oct4, !nuc_int_data$oct4$Label %in% c(label_T1_set2,label_T2_set2,label_T1_set5,label_T1_set1,label_T2_set1))
    decile_oct4<- quantile(oct4_int_selec$Median, prob = seq(0, 1, length = 11), type = 5)
    oct4_labels_0_10<-oct4_int_selec$Label[which(oct4_int_selec$Median >= decile_oct4[1] & oct4_int_selec$Median < decile_oct4[2] )]
    oct4_labels_10_20<-oct4_int_selec$Label[which(oct4_int_selec$Median >= decile_oct4[2] & oct4_int_selec$Median < decile_oct4[3] )]
    oct4_labels_20_30<-oct4_int_selec$Label[which(oct4_int_selec$Median >= decile_oct4[3] & oct4_int_selec$Median < decile_oct4[4] )]
    oct4_labels_30_40<-oct4_int_selec$Label[which(oct4_int_selec$Median >= decile_oct4[4] & oct4_int_selec$Median < decile_oct4[5] )]
    oct4_labels_40_50<-oct4_int_selec$Label[which(oct4_int_selec$Median >= decile_oct4[5] & oct4_int_selec$Median < decile_oct4[6] )]
    oct4_labels_50_60<-oct4_int_selec$Label[which(oct4_int_selec$Median >= decile_oct4[6] & oct4_int_selec$Median < decile_oct4[7] )]
    oct4_labels_60_70<-oct4_int_selec$Label[which(oct4_int_selec$Median >= decile_oct4[7] & oct4_int_selec$Median < decile_oct4[8] )]
    oct4_labels_70_80<-oct4_int_selec$Label[which(oct4_int_selec$Median >= decile_oct4[8] & oct4_int_selec$Median < decile_oct4[9] )]
    oct4_labels_80_90<-oct4_int_selec$Label[which(oct4_int_selec$Median >= decile_oct4[9] & oct4_int_selec$Median < decile_oct4[10] )]
    oct4_labels_90_100<-oct4_int_selec$Label[which(oct4_int_selec$Median >= decile_oct4[10] & oct4_int_selec$Median < decile_oct4[11] )]
    
    colfunc <- colorRampPalette(c("#FF3E96","#FFEC8B","#21840d" ))
    cols<-colfunc(11)
    
    png(filename="intensity_oct4_hist.png", units="in", width=1, height=1 , pointsize=5,res=1200)
    par(font.axis=2,font.lab=2,mar = c(4,4,1,1), font=2)
    d<-density(oct4_int_selec$Mean)
    plot(d, main="",las=1, xlab="Nuclear Oct4 levels",cex.axis=0.7)
    dec<-quantile(oct4_int_selec$Mean, prob = seq(0, 1, length = 11), type = 5)
    # Add the shaded area.
    polygon(c( dec[11],d$x[d$x>=dec[11]], max(d$x) ),  c(0,d$y[d$x>=dec[11]],0), col=cols[11],lwd=0.5,border = "black")
    polygon(c( dec[10],d$x[d$x>=dec[10] & d$x<dec[11]],dec[11] ),  c(0,d$y[d$x>=dec[10] & d$x<dec[11]],0), col=cols[11],lwd=0.5,border = "black")
    polygon(c( dec[9],d$x[d$x>=dec[9] & d$x<dec[10]],dec[10] ),  c(0,d$y[d$x>=dec[9] & d$x<dec[10]],0), col=cols[10],lwd=0.5,border = "black")
    polygon(c( dec[8],d$x[d$x>=dec[8] & d$x<dec[9]],dec[9] ),  c(0,d$y[d$x>=dec[8] & d$x<dec[9]],0), col=cols[9],lwd=0.5,border = "black")
    polygon(c( dec[7],d$x[d$x>=dec[7] & d$x<dec[8]],dec[8] ),  c(0,d$y[d$x>=dec[7] & d$x<dec[8]],0), col=cols[8],lwd=0.5,border = "black")
    polygon(c( dec[6],d$x[d$x>=dec[6] & d$x<dec[7]],dec[7] ),  c(0,d$y[d$x>=dec[6] & d$x<dec[7]],0), col=cols[7],lwd=0.5,border = "black")
    polygon(c( dec[5],d$x[d$x>=dec[5] & d$x<dec[6]],dec[6] ),  c(0,d$y[d$x>=dec[5] & d$x<dec[6]],0), col=cols[6],lwd=0.5,border = "black")
    polygon(c( dec[4],d$x[d$x>=dec[4] & d$x<dec[5]],dec[5] ),  c(0,d$y[d$x>=dec[4] & d$x<dec[5]],0), col=cols[5],lwd=0.5,border = "black")
    polygon(c( dec[3],d$x[d$x>=dec[3] & d$x<dec[4]],dec[4] ),  c(0,d$y[d$x>=dec[3] & d$x<dec[4]],0), col=cols[4],lwd=0.5,border = "black")
    polygon(c( dec[2],d$x[d$x>=dec[2] & d$x<dec[3]],dec[3] ),  c(0,d$y[d$x>=dec[2] & d$x<dec[3]],0), col=cols[3],lwd=0.5,border = "black")
    polygon(c( dec[1],d$x[d$x>=dec[1] & d$x<dec[2]],dec[2] ),  c(0,d$y[d$x>=dec[1] & d$x<dec[2]],0), col=cols[2],lwd=0.5,border = "black")
    polygon(c( min(d$x),d$x[d$x<dec[1]],dec[1]),  c(0,d$y[d$x<dec[1]],0), col=cols[1],lwd=0.5,border = "black")  
    dev.off()
    
    a<-subset(oct4_int_selec,oct4_int_selec$Label %in% D0_lab)
    b<-subset(oct4_int_selec,oct4_int_selec$Label %in% D2_lab)
    c<-subset(oct4_int_selec,oct4_int_selec$Label %in% D4_lab)
    d<-subset(oct4_int_selec,oct4_int_selec$Label %in% D6_lab)
    e<-subset(oct4_int_selec,oct4_int_selec$Label %in% D8_lab)
    
    png(filename="intensity_oct4.png", units="in", width=1.5, height=1.5 , pointsize=5,res=1200)
    par(font.axis=2,font.lab=2,mar = c(2,4,1,1), font=2)
    boxplot(a$Median,b$Median,c$Median,d$Median,e$Median,lty=1,pch=19,las=1,col=cols_timeline,
            names=c("D0","D2","D4","D6","D8"),ylab="Oct4 Nuclear Levels",cex=0.2)
    dev.off()
    rm(a,b,c,d,e)
    
    png(filename="nuc_oct4_summary_key.png", units="in", width=2, height=2 , pointsize=7, res=1200)
    par(font.axis = 2,font.lab=2,mar=c(4,4,2,2), font=2)
    colkey(col=colfunc(11),range(0,1),side = 1, add = F,length = 1, width = 1,clab = "Nuclear Oct4",side.clab=1)
    dev.off()
    
  }
  
  #plots
  {
    merge_ld1_prop_oct4<-merge(merged_data,oct4_int_selec[,which(colnames(oct4_int_selec) %in% c("Median","Label"))],by="Label")
    
    lmMod <- lm(Median ~LD1,data=merge_ld1_prop_oct4)  # build the model
    
    png(filename="oct4_vs_ld1.png", units="in", width=1.5, height=1.5 , pointsize=5,res=1200)
    par(font.axis=2,font.lab=2,mar = c(4,4,1,1), font=2)
    plot(merge_ld1_prop_oct4$LD1,merge_ld1_prop_oct4$Median,las=1,pch=19,ylab="Nuclear Oct4 Levels",xlab="LD1",col=alpha.col(col="gray",alpha=0.5))
    abline(lmMod,lwd=2,col="red")
    dev.off()
    
    d1<-subset(merge_ld1_prop_oct4,merge_ld1_prop_oct4$Label %in% oct4_labels_0_10)
    d2<-subset(merge_ld1_prop_oct4,merge_ld1_prop_oct4$Label %in% oct4_labels_10_20)
    d3<-subset(merge_ld1_prop_oct4,merge_ld1_prop_oct4$Label %in% oct4_labels_20_30)
    d4<-subset(merge_ld1_prop_oct4,merge_ld1_prop_oct4$Label %in% oct4_labels_30_40)
    d5<-subset(merge_ld1_prop_oct4,merge_ld1_prop_oct4$Label %in% oct4_labels_40_50)
    d6<-subset(merge_ld1_prop_oct4,merge_ld1_prop_oct4$Label %in% oct4_labels_50_60)
    d7<-subset(merge_ld1_prop_oct4,merge_ld1_prop_oct4$Label %in% oct4_labels_60_70)
    d8<-subset(merge_ld1_prop_oct4,merge_ld1_prop_oct4$Label %in% oct4_labels_70_80)
    d9<-subset(merge_ld1_prop_oct4,merge_ld1_prop_oct4$Label %in% oct4_labels_80_90)
    d10<-subset(merge_ld1_prop_oct4,merge_ld1_prop_oct4$Label %in% oct4_labels_90_100)
    
    png(filename="oct4_ld1vs_ld2.png", units="in", width=1.5, height=1.5 , pointsize=5,res=1200)
    par(font.axis=2,font.lab=2,mar = c(4,4,1,1), font=2,font.axis=2,font.lab=2)
    plot(merge_ld1_prop_oct4$LD1,merge_ld1_prop_oct4$LD2,col='NA',las=1,xlab="LD1",ylab="LD2")
    points(d1$LD1,d1$LD2,col=cols[1],pch=19)
    points(d2$LD1,d2$LD2,col=cols[2],pch=19)
    points(d3$LD1,d3$LD2,col=cols[3],pch=19)
    points(d4$LD1,d4$LD2,col=cols[4],pch=19)
    points(d5$LD1,d5$LD2,col=cols[5],pch=19)
    points(d6$LD1,d6$LD2,col=cols[6],pch=19)
    points(d7$LD1,d7$LD2,col=cols[7],pch=19)
    points(d8$LD1,d8$LD2,col=cols[8],pch=19)
    points(d9$LD1,d9$LD2,col=cols[9],pch=19)
    points(d10$LD1,d10$LD2,col=cols[10],pch=19)
    dev.off()
    
    d1$samplem="a10th"
    d2$samplem="b20th"
    d3$samplem="c30th"
    d4$samplem="d40th"
    d5$samplem="e50th"
    d6$samplem="f60th"
    d7$samplem="g70th"
    d8$samplem="h80th"
    d9$samplem="i90th"
    d10$samplem="j100th"
    
    merge_ld1_prop_oct4_temp<-rbind(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10)
    
    png(filename="oct4_ld1vs_ld2_ellipse.png", units="in", width=1.5, height=1.5 , pointsize=5,res=1200)
    par(font.axis=2,font.lab=2,mar = c(4,4,1,1), font=2,font.axis=2,font.lab=2)
    ellipseplot(merge_ld1_prop_oct4_temp$LD1,                          # data for x-axis
                merge_ld1_prop_oct4_temp$LD2,                          # data for y-axis
                as.factor(merge_ld1_prop_oct4_temp$samplem),                    # factor with classes
                elev=0.95,
                pcol=cols,    # colors for plotting (must match # of factors)
                pbgcol=TRUE,                           # point borders black?
                cexsize=0.6,                            # size of points 
                ppch=rep(21,10),                          # shape of points (must match # of factors)
                axissize=0.8,                           # Set axis text size
                linewidth=1,                          # Set axis line size
                font=1
    )
    title(xlab="LD1",    # % variance explained on PC1
          ylab="LD2",    # % variance explained on PC2 
          main=" ",                  # Title
          cex.lab=1,                    # size of label text
          cex.main=1 ,                   # size of title text
          font.axis=2
    )
    dev.off()
  }
  
  
}


#h3k9acdeciles
{
  library(plot3D)
  
  #deciles
  {
    h3k9ac_int_selec<-nuc_int_data$h3k9ac
    a<-subset(h3k9ac_int_selec,h3k9ac_int_selec$Label %in% D0_lab)
    b<-subset(h3k9ac_int_selec,h3k9ac_int_selec$Label %in% D2_lab)
    c<-subset(h3k9ac_int_selec,h3k9ac_int_selec$Label %in% D4_lab)
    d<-subset(h3k9ac_int_selec,h3k9ac_int_selec$Label %in% D6_lab)
    e<-subset(h3k9ac_int_selec,h3k9ac_int_selec$Label %in% D8_lab)
    
    png(filename="intensity_h3k9ac.png", units="in", width=1.5, height=1.5 , pointsize=5,res=1200)
    par(font.axis=2,font.lab=2,mar = c(2,4,1,1), font=2)
    boxplot(a$Median,b$Median,c$Median,d$Median,e$Median,lty=1,pch=19,las=1,col=cols_timeline,
            names=c("D0","D2","D4","D6","D8"),ylab="h3k9ac Nuclear Levels",cex=0.2)
    dev.off()
    rm(a,b,c,d,e)
    
    decile_h3k9ac<- quantile(h3k9ac_int_selec$Median, prob = seq(0, 1, length = 11), type = 5)
    h3k9ac_labels_0_10<-h3k9ac_int_selec$Label[which(h3k9ac_int_selec$Median >= decile_h3k9ac[1] & h3k9ac_int_selec$Median < decile_h3k9ac[2] )]
    h3k9ac_labels_10_20<-h3k9ac_int_selec$Label[which(h3k9ac_int_selec$Median >= decile_h3k9ac[2] & h3k9ac_int_selec$Median < decile_h3k9ac[3] )]
    h3k9ac_labels_20_30<-h3k9ac_int_selec$Label[which(h3k9ac_int_selec$Median >= decile_h3k9ac[3] & h3k9ac_int_selec$Median < decile_h3k9ac[4] )]
    h3k9ac_labels_30_40<-h3k9ac_int_selec$Label[which(h3k9ac_int_selec$Median >= decile_h3k9ac[4] & h3k9ac_int_selec$Median < decile_h3k9ac[5] )]
    h3k9ac_labels_40_50<-h3k9ac_int_selec$Label[which(h3k9ac_int_selec$Median >= decile_h3k9ac[5] & h3k9ac_int_selec$Median < decile_h3k9ac[6] )]
    h3k9ac_labels_50_60<-h3k9ac_int_selec$Label[which(h3k9ac_int_selec$Median >= decile_h3k9ac[6] & h3k9ac_int_selec$Median < decile_h3k9ac[7] )]
    h3k9ac_labels_60_70<-h3k9ac_int_selec$Label[which(h3k9ac_int_selec$Median >= decile_h3k9ac[7] & h3k9ac_int_selec$Median < decile_h3k9ac[8] )]
    h3k9ac_labels_70_80<-h3k9ac_int_selec$Label[which(h3k9ac_int_selec$Median >= decile_h3k9ac[8] & h3k9ac_int_selec$Median < decile_h3k9ac[9] )]
    h3k9ac_labels_80_90<-h3k9ac_int_selec$Label[which(h3k9ac_int_selec$Median >= decile_h3k9ac[9] & h3k9ac_int_selec$Median < decile_h3k9ac[10] )]
    h3k9ac_labels_90_100<-h3k9ac_int_selec$Label[which(h3k9ac_int_selec$Median >= decile_h3k9ac[10] & h3k9ac_int_selec$Median < decile_h3k9ac[11] )]
    
    colfunc <- colorRampPalette(c("#FF3E96","#FFEC8B","#21840d" ))
    cols<-colfunc(11)
    
    png(filename="intensity_h3k9ac_hist.png", units="in", width=1, height=1 , pointsize=5,res=1200)
    par(font.axis=2,font.lab=2,mar = c(4,4,1,1), font=2)
    d<-density(h3k9ac_int_selec$Median)
    plot(d, main="",las=1, xlab="Nuclear h3k9ac levels",cex.axis=0.7)
    dec<-quantile(h3k9ac_int_selec$Median, prob = seq(0, 1, length = 11), type = 5)
    # Add the shaded area.
    polygon(c( dec[11],d$x[d$x>=dec[11]], max(d$x) ),  c(0,d$y[d$x>=dec[11]],0), col=cols[11],lwd=0.5,border = "black")
    polygon(c( dec[10],d$x[d$x>=dec[10] & d$x<dec[11]],dec[11] ),  c(0,d$y[d$x>=dec[10] & d$x<dec[11]],0), col=cols[11],lwd=0.5,border = "black")
    polygon(c( dec[9],d$x[d$x>=dec[9] & d$x<dec[10]],dec[10] ),  c(0,d$y[d$x>=dec[9] & d$x<dec[10]],0), col=cols[10],lwd=0.5,border = "black")
    polygon(c( dec[8],d$x[d$x>=dec[8] & d$x<dec[9]],dec[9] ),  c(0,d$y[d$x>=dec[8] & d$x<dec[9]],0), col=cols[9],lwd=0.5,border = "black")
    polygon(c( dec[7],d$x[d$x>=dec[7] & d$x<dec[8]],dec[8] ),  c(0,d$y[d$x>=dec[7] & d$x<dec[8]],0), col=cols[8],lwd=0.5,border = "black")
    polygon(c( dec[6],d$x[d$x>=dec[6] & d$x<dec[7]],dec[7] ),  c(0,d$y[d$x>=dec[6] & d$x<dec[7]],0), col=cols[7],lwd=0.5,border = "black")
    polygon(c( dec[5],d$x[d$x>=dec[5] & d$x<dec[6]],dec[6] ),  c(0,d$y[d$x>=dec[5] & d$x<dec[6]],0), col=cols[6],lwd=0.5,border = "black")
    polygon(c( dec[4],d$x[d$x>=dec[4] & d$x<dec[5]],dec[5] ),  c(0,d$y[d$x>=dec[4] & d$x<dec[5]],0), col=cols[5],lwd=0.5,border = "black")
    polygon(c( dec[3],d$x[d$x>=dec[3] & d$x<dec[4]],dec[4] ),  c(0,d$y[d$x>=dec[3] & d$x<dec[4]],0), col=cols[4],lwd=0.5,border = "black")
    polygon(c( dec[2],d$x[d$x>=dec[2] & d$x<dec[3]],dec[3] ),  c(0,d$y[d$x>=dec[2] & d$x<dec[3]],0), col=cols[3],lwd=0.5,border = "black")
    polygon(c( dec[1],d$x[d$x>=dec[1] & d$x<dec[2]],dec[2] ),  c(0,d$y[d$x>=dec[1] & d$x<dec[2]],0), col=cols[2],lwd=0.5,border = "black")
    polygon(c( min(d$x),d$x[d$x<dec[1]],dec[1]),  c(0,d$y[d$x<dec[1]],0), col=cols[1],lwd=0.5,border = "black")  
    dev.off()
    
   
    
    png(filename="nuc_h3k9ac_summary_key.png", units="in", width=2, height=2 , pointsize=7, res=1200)
    par(font.axis = 2,font.lab=2,mar=c(4,4,2,2), font=2)
    colkey(col=colfunc(11),range(0,1),side = 1, add = F,length = 1, width = 1,clab = "Nuclear H3K9Ac",side.clab=1)
    dev.off()
    
  }
  
  #plots
  {
    merge_ld1_prop_h3k9ac<-merge(merged_data,h3k9ac_int_selec[,which(colnames(h3k9ac_int_selec) %in% c("Median","Label"))],by="Label")
    
    lmMod <- lm(Median ~LD1,data=merge_ld1_prop_h3k9ac)  # build the model
    
    png(filename="h3k9ac_vs_ld1.png", units="in", width=1.5, height=1.5 , pointsize=5,res=1200)
    par(font.axis=2,font.lab=2,mar = c(4,4,1,1), font=2)
    plot(merge_ld1_prop_h3k9ac$LD1,merge_ld1_prop_h3k9ac$Median,las=1,pch=19,ylab="Nuclear h3k9ac Levels",xlab="LD1",col=alpha.col(col="gray",alpha=0.5))
    abline(lmMod,lwd=2,col="red")
    dev.off()
    
    d1<-subset(merge_ld1_prop_h3k9ac,merge_ld1_prop_h3k9ac$Label %in% h3k9ac_labels_0_10)
    d2<-subset(merge_ld1_prop_h3k9ac,merge_ld1_prop_h3k9ac$Label %in% h3k9ac_labels_10_20)
    d3<-subset(merge_ld1_prop_h3k9ac,merge_ld1_prop_h3k9ac$Label %in% h3k9ac_labels_20_30)
    d4<-subset(merge_ld1_prop_h3k9ac,merge_ld1_prop_h3k9ac$Label %in% h3k9ac_labels_30_40)
    d5<-subset(merge_ld1_prop_h3k9ac,merge_ld1_prop_h3k9ac$Label %in% h3k9ac_labels_40_50)
    d6<-subset(merge_ld1_prop_h3k9ac,merge_ld1_prop_h3k9ac$Label %in% h3k9ac_labels_50_60)
    d7<-subset(merge_ld1_prop_h3k9ac,merge_ld1_prop_h3k9ac$Label %in% h3k9ac_labels_60_70)
    d8<-subset(merge_ld1_prop_h3k9ac,merge_ld1_prop_h3k9ac$Label %in% h3k9ac_labels_70_80)
    d9<-subset(merge_ld1_prop_h3k9ac,merge_ld1_prop_h3k9ac$Label %in% h3k9ac_labels_80_90)
    d10<-subset(merge_ld1_prop_h3k9ac,merge_ld1_prop_h3k9ac$Label %in% h3k9ac_labels_90_100)
    
    png(filename="h3k9ac_ld1vs_ld2.png", units="in", width=1.5, height=1.5 , pointsize=5,res=1200)
    par(font.axis=2,font.lab=2,mar = c(4,4,1,1), font=2,font.axis=2,font.lab=2)
    plot(merge_ld1_prop_h3k9ac$LD1,merge_ld1_prop_h3k9ac$LD2,col='NA',las=1,xlab="LD1",ylab="LD2")
    points(d1$LD1,d1$LD2,col=cols[1],pch=19)
    points(d2$LD1,d2$LD2,col=cols[2],pch=19)
    points(d3$LD1,d3$LD2,col=cols[3],pch=19)
    points(d4$LD1,d4$LD2,col=cols[4],pch=19)
    points(d5$LD1,d5$LD2,col=cols[5],pch=19)
    points(d6$LD1,d6$LD2,col=cols[6],pch=19)
    points(d7$LD1,d7$LD2,col=cols[7],pch=19)
    points(d8$LD1,d8$LD2,col=cols[8],pch=19)
    points(d9$LD1,d9$LD2,col=cols[9],pch=19)
    points(d10$LD1,d10$LD2,col=cols[10],pch=19)
    dev.off()
    
    d1$samplem="a10th"
    d2$samplem="b20th"
    d3$samplem="c30th"
    d4$samplem="d40th"
    d5$samplem="e50th"
    d6$samplem="f60th"
    d7$samplem="g70th"
    d8$samplem="h80th"
    d9$samplem="i90th"
    d10$samplem="j100th"
    
    merge_ld1_prop_h3k9ac_temp<-rbind(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10)
    
    png(filename="h3k9ac_ld1vs_ld2_ellipse.png", units="in", width=1.5, height=1.5 , pointsize=5,res=1200)
    par(font.axis=2,font.lab=2,mar = c(4,4,1,1), font=2,font.axis=2,font.lab=2)
    ellipseplot(merge_ld1_prop_h3k9ac_temp$LD1,                          # data for x-axis
                merge_ld1_prop_h3k9ac_temp$LD2,                          # data for y-axis
                as.factor(merge_ld1_prop_h3k9ac_temp$samplem),                    # factor with classes
                elev=0.9,
                pcol=cols,    # colors for plotting (must match # of factors)
                pbgcol=TRUE,                           # point borders black?
                cexsize=0.6,                            # size of points 
                ppch=rep(21,10),                          # shape of points (must match # of factors)
                axissize=0.8,                           # Set axis text size
                linewidth=1,                          # Set axis line size
                font=1
    )
    title(xlab="LD1",    # % variance explained on PC1
          ylab="LD2",    # % variance explained on PC2 
          main=" ",                  # Title
          cex.lab=1,                    # size of label text
          cex.main=1 ,                   # size of title text
          font.axis=2
    )
    dev.off()
  }
  
  
}

#ki67deciles
{
  library(plot3D)
  
  #deciles
  {
    ki67_int_selec<-nuc_int_data$ki67
    a<-subset(ki67_int_selec,ki67_int_selec$Label %in% D0_lab)
    b<-subset(ki67_int_selec,ki67_int_selec$Label %in% D2_lab)
    c<-subset(ki67_int_selec,ki67_int_selec$Label %in% D4_lab)
    d<-subset(ki67_int_selec,ki67_int_selec$Label %in% D6_lab)
    e<-subset(ki67_int_selec,ki67_int_selec$Label %in% D8_lab)
    
    png(filename="intensity_ki67.png", units="in", width=1.5, height=1.5 , pointsize=5,res=1200)
    par(font.axis=2,font.lab=2,mar = c(2,4,1,1), font=2)
    boxplot(a$Mean,b$Mean,c$Mean,d$Mean,e$Mean,lty=1,pch=19,las=1,col=cols_timeline,
            names=c("D0","D2","D4","D6","D8"),ylab="ki67 Nuclear Levels",cex=0.2)
    dev.off()
    rm(a,b,c,d,e)
    
    decile_ki67<- quantile(ki67_int_selec$Mean, prob = seq(0, 1, length = 11), type = 5)
    ki67_labels_0_10<-ki67_int_selec$Label[which(ki67_int_selec$Mean >= decile_ki67[1] & ki67_int_selec$Mean < decile_ki67[2] )]
    ki67_labels_10_20<-ki67_int_selec$Label[which(ki67_int_selec$Mean >= decile_ki67[2] & ki67_int_selec$Mean < decile_ki67[3] )]
    ki67_labels_20_30<-ki67_int_selec$Label[which(ki67_int_selec$Mean >= decile_ki67[3] & ki67_int_selec$Mean < decile_ki67[4] )]
    ki67_labels_30_40<-ki67_int_selec$Label[which(ki67_int_selec$Mean >= decile_ki67[4] & ki67_int_selec$Mean < decile_ki67[5] )]
    ki67_labels_40_50<-ki67_int_selec$Label[which(ki67_int_selec$Mean >= decile_ki67[5] & ki67_int_selec$Mean < decile_ki67[6] )]
    ki67_labels_50_60<-ki67_int_selec$Label[which(ki67_int_selec$Mean >= decile_ki67[6] & ki67_int_selec$Mean < decile_ki67[7] )]
    ki67_labels_60_70<-ki67_int_selec$Label[which(ki67_int_selec$Mean >= decile_ki67[7] & ki67_int_selec$Mean < decile_ki67[8] )]
    ki67_labels_70_80<-ki67_int_selec$Label[which(ki67_int_selec$Mean >= decile_ki67[8] & ki67_int_selec$Mean < decile_ki67[9] )]
    ki67_labels_80_90<-ki67_int_selec$Label[which(ki67_int_selec$Mean >= decile_ki67[9] & ki67_int_selec$Mean < decile_ki67[10] )]
    ki67_labels_90_100<-ki67_int_selec$Label[which(ki67_int_selec$Mean >= decile_ki67[10] & ki67_int_selec$Mean < decile_ki67[11] )]
    
    colfunc <- colorRampPalette(c("#FF3E96","#FFEC8B","#21840d" ))
    cols<-colfunc(11)
    
    png(filename="intensity_ki67_hist.png", units="in", width=1, height=1 , pointsize=5,res=1200)
    par(font.axis=2,font.lab=2,mar = c(4,4,1,1), font=2)
    d<-density(ki67_int_selec$Mean)
    plot(d, main="",las=1, xlab="Nuclear ki67 levels",cex.axis=0.7)
    dec<-quantile(ki67_int_selec$Mean, prob = seq(0, 1, length = 11), type = 5)
    # Add the shaded area.
    polygon(c( dec[11],d$x[d$x>=dec[11]], max(d$x) ),  c(0,d$y[d$x>=dec[11]],0), col=cols[11],lwd=0.5,border = "black")
    polygon(c( dec[10],d$x[d$x>=dec[10] & d$x<dec[11]],dec[11] ),  c(0,d$y[d$x>=dec[10] & d$x<dec[11]],0), col=cols[11],lwd=0.5,border = "black")
    polygon(c( dec[9],d$x[d$x>=dec[9] & d$x<dec[10]],dec[10] ),  c(0,d$y[d$x>=dec[9] & d$x<dec[10]],0), col=cols[10],lwd=0.5,border = "black")
    polygon(c( dec[8],d$x[d$x>=dec[8] & d$x<dec[9]],dec[9] ),  c(0,d$y[d$x>=dec[8] & d$x<dec[9]],0), col=cols[9],lwd=0.5,border = "black")
    polygon(c( dec[7],d$x[d$x>=dec[7] & d$x<dec[8]],dec[8] ),  c(0,d$y[d$x>=dec[7] & d$x<dec[8]],0), col=cols[8],lwd=0.5,border = "black")
    polygon(c( dec[6],d$x[d$x>=dec[6] & d$x<dec[7]],dec[7] ),  c(0,d$y[d$x>=dec[6] & d$x<dec[7]],0), col=cols[7],lwd=0.5,border = "black")
    polygon(c( dec[5],d$x[d$x>=dec[5] & d$x<dec[6]],dec[6] ),  c(0,d$y[d$x>=dec[5] & d$x<dec[6]],0), col=cols[6],lwd=0.5,border = "black")
    polygon(c( dec[4],d$x[d$x>=dec[4] & d$x<dec[5]],dec[5] ),  c(0,d$y[d$x>=dec[4] & d$x<dec[5]],0), col=cols[5],lwd=0.5,border = "black")
    polygon(c( dec[3],d$x[d$x>=dec[3] & d$x<dec[4]],dec[4] ),  c(0,d$y[d$x>=dec[3] & d$x<dec[4]],0), col=cols[4],lwd=0.5,border = "black")
    polygon(c( dec[2],d$x[d$x>=dec[2] & d$x<dec[3]],dec[3] ),  c(0,d$y[d$x>=dec[2] & d$x<dec[3]],0), col=cols[3],lwd=0.5,border = "black")
    polygon(c( dec[1],d$x[d$x>=dec[1] & d$x<dec[2]],dec[2] ),  c(0,d$y[d$x>=dec[1] & d$x<dec[2]],0), col=cols[2],lwd=0.5,border = "black")
    polygon(c( min(d$x),d$x[d$x<dec[1]],dec[1]),  c(0,d$y[d$x<dec[1]],0), col=cols[1],lwd=0.5,border = "black")  
    dev.off()
    
    
    
    png(filename="nuc_ki67_summary_key.png", units="in", width=2, height=2 , pointsize=7, res=1200)
    par(font.axis = 2,font.lab=2,mar=c(4,4,2,2), font=2)
    colkey(col=colfunc(11),range(0,1),side = 1, add = F,length = 1, width = 1,clab = "Nuclear ki67",side.clab=1)
    dev.off()
    
  }
  
  #plots
  {
    merge_ld1_prop_ki67<-merge(merged_data,ki67_int_selec[,which(colnames(ki67_int_selec) %in% c("Mean","Label"))],by="Label")
    
    lmMod <- lm(Mean ~LD1,data=merge_ld1_prop_ki67)  # build the model
    
    png(filename="ki67_vs_ld1.png", units="in", width=1.5, height=1.5 , pointsize=5,res=1200)
    par(font.axis=2,font.lab=2,mar = c(4,4,1,1), font=2)
    plot(merge_ld1_prop_ki67$LD1,merge_ld1_prop_ki67$Mean,las=1,pch=19,ylab="Nuclear ki67 Levels",xlab="LD1",col=alpha.col(col="gray",alpha=0.5))
    abline(lmMod,lwd=2,col="red")
    dev.off()
    
    d1<-subset(merge_ld1_prop_ki67,merge_ld1_prop_ki67$Label %in% ki67_labels_0_10)
    d2<-subset(merge_ld1_prop_ki67,merge_ld1_prop_ki67$Label %in% ki67_labels_10_20)
    d3<-subset(merge_ld1_prop_ki67,merge_ld1_prop_ki67$Label %in% ki67_labels_20_30)
    d4<-subset(merge_ld1_prop_ki67,merge_ld1_prop_ki67$Label %in% ki67_labels_30_40)
    d5<-subset(merge_ld1_prop_ki67,merge_ld1_prop_ki67$Label %in% ki67_labels_40_50)
    d6<-subset(merge_ld1_prop_ki67,merge_ld1_prop_ki67$Label %in% ki67_labels_50_60)
    d7<-subset(merge_ld1_prop_ki67,merge_ld1_prop_ki67$Label %in% ki67_labels_60_70)
    d8<-subset(merge_ld1_prop_ki67,merge_ld1_prop_ki67$Label %in% ki67_labels_70_80)
    d9<-subset(merge_ld1_prop_ki67,merge_ld1_prop_ki67$Label %in% ki67_labels_80_90)
    d10<-subset(merge_ld1_prop_ki67,merge_ld1_prop_ki67$Label %in% ki67_labels_90_100)
    
    png(filename="ki67_ld1vs_ld2.png", units="in", width=1.5, height=1.5 , pointsize=5,res=1200)
    par(font.axis=2,font.lab=2,mar = c(4,4,1,1), font=2,font.axis=2,font.lab=2)
    plot(merge_ld1_prop_ki67$LD1,merge_ld1_prop_ki67$LD2,col='NA',las=1,xlab="LD1",ylab="LD2")
    points(d1$LD1,d1$LD2,col=cols[1],pch=19)
    points(d2$LD1,d2$LD2,col=cols[2],pch=19)
    points(d3$LD1,d3$LD2,col=cols[3],pch=19)
    points(d4$LD1,d4$LD2,col=cols[4],pch=19)
    points(d5$LD1,d5$LD2,col=cols[5],pch=19)
    points(d6$LD1,d6$LD2,col=cols[6],pch=19)
    points(d7$LD1,d7$LD2,col=cols[7],pch=19)
    points(d8$LD1,d8$LD2,col=cols[8],pch=19)
    points(d9$LD1,d9$LD2,col=cols[9],pch=19)
    points(d10$LD1,d10$LD2,col=cols[10],pch=19)
    dev.off()
    
    d1$samplem="a10th"
    d2$samplem="b20th"
    d3$samplem="c30th"
    d4$samplem="d40th"
    d5$samplem="e50th"
    d6$samplem="f60th"
    d7$samplem="g70th"
    d8$samplem="h80th"
    d9$samplem="i90th"
    d10$samplem="j100th"
    
    merge_ld1_prop_ki67_temp<-rbind(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10)
    
    png(filename="ki67_ld1vs_ld2_ellipse.png", units="in", width=1.5, height=1.5 , pointsize=5,res=1200)
    par(font.axis=2,font.lab=2,mar = c(4,4,1,1), font=2,font.axis=2,font.lab=2)
    ellipseplot(merge_ld1_prop_ki67_temp$LD1,                          # data for x-axis
                merge_ld1_prop_ki67_temp$LD2,                          # data for y-axis
                as.factor(merge_ld1_prop_ki67_temp$samplem),                    # factor with classes
                elev=0.9,
                pcol=cols,    # colors for plotting (must match # of factors)
                pbgcol=TRUE,                           # point borders black?
                cexsize=0.6,                            # size of points 
                ppch=rep(21,10),                          # shape of points (must match # of factors)
                axissize=0.8,                           # Set axis text size
                linewidth=1,                          # Set axis line size
                font=1
    )
    title(xlab="LD1",    # % variance explained on PC1
          ylab="LD2",    # % variance explained on PC2 
          main=" ",                  # Title
          cex.lab=1,                    # size of label text
          cex.main=1 ,                   # size of title text
          font.axis=2
    )
    dev.off()
  }
  
  
}
rm(ls())

#required functions
{
  require(plyr)
  library(MASS)
  library(corrplot)
  require(sp)
  require(rgeos)
  require(car)
  plotat <- function(RANGE) {
    if(length(RANGE) != 2) stop("RANGE argument must have a length of 2")
    if(RANGE[1] > RANGE[2]) stop("First element in RANGE must be smaller then second element")
    prettyres <- pretty(sprintf("%.2f",RANGE[1]):sprintf("%.2f",RANGE[2]), 7)
    while((min(prettyres) < RANGE[1]) == FALSE) {
      prdiff <- prettyres[2] - prettyres[1]
      prettyres[length(prettyres) + 1] <- prettyres[1] - prdiff
      prettyres <- sort(prettyres)
    } 
    while((max(prettyres) > RANGE[2]) == FALSE) {
      prdiff <- prettyres[2] - prettyres[1]
      prettyres[length(prettyres) + 1] <- prettyres[length(prettyres)] + prdiff
      prettyres <- sort(prettyres)    
    }   
    plotticks <- as.numeric(sprintf("%.2f",prettyres))
    plotticks
  }
  
  ## ellipseplot function
  ellipseplot <- function(x, y, factr, 
                          elev=0.95, # Ellipse probability level
                          pcol=NULL, # manual addition of colors, must meet length of factors
                          cexsize=1, # point size
                          ppch=21, # Point type, must meet length of factors
                          pbgcol=TRUE,
                          axissize=1, 
                          linewidth=1, 
                          font=1) {
    
    ## Set factor levels
    if(is.factor(factr)) {
      f <- factr
    } 
    else {
      f <- factor(factr, levels=unique(as.character(factr)))
    }
    intfactr <- as.integer(f) # Set integer vector that matches factor levels
    
    # Checking to make sure length of ppch equals number of factor levels
    if((length(ppch) > 1 & length(unique(intfactr)) != length(ppch))) stop("Can only increase point shape if equal to factor levels")
    
    ## Get data for ellipses
    edf <- data.frame(LV1 = x, LV2=y, factr = f) # create data frame with data and factor
    ellipses <- dlply(edf, .(factr), function(x) {
      LV1 <- x[,1]
      LV2 <- x[,2]
      dataEllipse(LV1, LV2, levels=elev, robust=TRUE, draw=FALSE) # Get confidence ellipse points from dataEllipse() function by factor level
    })
    ## Get range of x and y data
    xrange <- plotat(range(c(as.vector(sapply(ellipses, function(x) x[,1])), min(x), max(x))))
    yrange <- plotat(range(c(as.vector(sapply(ellipses, function(x) x[,2])), min(y), max(y))))
    
    
    ## Set colors for plots
    if(is.null(pcol) != TRUE) { # If colors are supplied by user
      ptcol <- pcol
      pgcol <- paste(pcol, "7e", sep="") # adds opaqueness
    } else { # Default
      pgcol <- c("#e41a1c7e","#377eb87e","#4daf4a7e","#984ea37e","#807f7d7e") # Defaults at 5 colors
      ptcol <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#807f7d") # For opaqueness
    }
    # Plotting graphic
    plot(x,y, type="n", xlab="", ylab="", main="", xlim=range(xrange), ylim=range(yrange), axes=FALSE)
    axis(1, at=xrange, labels=xrange, cex.axis=axissize,lwd=linewidth, font=font)
    axis(2, las=2, cex.axis=axissize,lwd=linewidth, font=font)
    box(lwd=linewidth, font=font)
    abline(h=0, v=0, col="gray", lty=2) # Adds lines at 0
    legpch <- c() # vector to collect legend pch data
    legcol <- c() # vector to collect legend col data
    
    ## Adds points, ellipse, and determines color specifications for legend 
    if(pbgcol==TRUE)  {
      for(i in 1:length(unique(intfactr))){
        points(x[intfactr==i], y[intfactr==i], pch=ppch[i], col=ptcol[i], bg=ptcol[i],cex=cexsize)
        polygon(ellipses[[i]], col=pgcol[i], border=ptcol[i])
        legpch[i] <- ppch[i]
        legcol[i] <- ptcol[i]
      }
    } else {
      for(i in 1:length(unique(intfactr))){
        points(x[intfactr==i], y[intfactr==i], pch=ppch[i], col="black", bg=ptcol[i],cex=cexsize)
        polygon(ellipses[[i]], col=pgcol[i], border=ptcol[i])
        legpch[i] <- ppch[i]
        legcol[i] <- ptcol[i]       
      }
    }
    
  }     
  
  ## Function for creating a SpatialPolygons object from data.frame of coords
  xy2SP <- function(xy, ID=NULL) {
    if(is.null(ID)) ID <- sample(1e12, size=1)
    SpatialPolygons(list(Polygons(list(Polygon(xy)), ID=ID)),
                    proj4string=CRS("+proj=merc"))
  }
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  
  normalize_2d_mean_trials <-function(data){
    
    if(length(levels(as.factor(data$trial))==3)){
      T1<-subset(data,data$trial=="T1")
      T1$Mean<-range01(T1$Mean)
      T1$IntDen<-range01(T1$IntDen)
      T1$Mode<-range01(T1$Mode)
      T1$Min<-range01(T1$Min)
      T1$Max<-range01(T1$Max)
      T1$Median<-range01(T1$Median)
      T1$RawIntDen<-range01(T1$RawIntDen)
      
      T2<-subset(data,data$trial=="T2")
      T2$Mean<-range01(T2$Mean)
      T2$IntDen<-range01(T2$IntDen)
      T2$Mode<-range01(T2$Mode)
      T2$Min<-range01(T2$Min)
      T2$Max<-range01(T2$Max)
      T2$Median<-range01(T2$Median)
      T2$RawIntDen<-range01(T2$RawIntDen)
      
      T3<-subset(data,data$trial=="T3")
      T3$Mean<-range01(T3$Mean)
      T3$IntDen<-range01(T3$IntDen)
      T3$Mode<-range01(T3$Mode)
      T3$Min<-range01(T3$Min)
      T3$Max<-range01(T3$Max)
      T3$Median<-range01(T3$Median)
      T3$RawIntDen<-range01(T3$RawIntDen)
      
      return(rbind(T1,T2,T3))
    }
    
    else if(length(levels(as.factor(data$trial))==4)){
      T1<-subset(data,data$trial=="T1")
      T1$Mean<-range01(T1$Mean)
      T1$IntDen<-range01(T1$IntDen)
      T1$Mode<-range01(T1$Mode)
      T1$Min<-range01(T1$Min)
      T1$Max<-range01(T1$Max)
      T1$Median<-range01(T1$Median)
      T1$RawIntDen<-range01(T1$RawIntDen)
      
      T2<-subset(data,data$trial=="T2")
      T2$Mean<-range01(T2$Mean)
      T2$IntDen<-range01(T2$IntDen)
      T2$Mode<-range01(T2$Mode)
      T2$Min<-range01(T2$Min)
      T2$Max<-range01(T2$Max)
      T2$Median<-range01(T2$Median)
      T2$RawIntDen<-range01(T2$RawIntDen)
      
      T3<-subset(data,data$trial=="T3")
      T3$Mean<-range01(T3$Mean)
      T3$IntDen<-range01(T3$IntDen)
      T3$Mode<-range01(T3$Mode)
      T3$Min<-range01(T3$Min)
      T3$Max<-range01(T3$Max)
      T3$Median<-range01(T3$Median)
      T3$RawIntDen<-range01(T3$RawIntDen)
      
      T4<-subset(data,data$trial=="T4")
      T4$Mean<-range01(T4$Mean)
      T4$IntDen<-range01(T4$IntDen)
      T4$Mode<-range01(T4$Mode)
      T4$Min<-range01(T4$Min)
      T4$Max<-range01(T4$Max)
      T4$Median<-range01(T4$Median)
      T4$RawIntDen<-range01(T4$RawIntDen)
      
      
      return(rbind(T1,T2,T3,T4))
    }
  }
  normalize_mean_trials <-function(data){
    
    if(length(levels(as.factor(data$trial))==3)){
      T1<-subset(data,data$trial=="T1")
      T1$Mean<-range01(T1$Mean)
      T1$IntDen<-range01(T1$IntDen)
      T1$Mode<-range01(T1$Mode)
      T1$Min<-range01(T1$Min)
      T1$Max<-range01(T1$Max)
      
      T2<-subset(data,data$trial=="T2")
      T2$Mean<-range01(T2$Mean)
      T2$IntDen<-range01(T2$IntDen)
      T2$Mode<-range01(T2$Mode)
      T2$Min<-range01(T2$Min)
      T2$Max<-range01(T2$Max)
      
      T3<-subset(data,data$trial=="T3")
      T3$Mean<-range01(T3$Mean)
      T3$IntDen<-range01(T3$IntDen)
      T3$Mode<-range01(T3$Mode)
      T3$Min<-range01(T3$Min)
      T3$Max<-range01(T3$Max)
      
      return(rbind(T1,T2,T3))
    }
    
    else if(length(levels(as.factor(data$trial))==4)){
      T1<-subset(data,data$trial=="T1")
      T1$Mean<-range01(T1$Mean)
      T1$IntDen<-range01(T1$IntDen)
      T1$Mode<-range01(T1$Mode)
      T1$Min<-range01(T1$Min)
      T1$Max<-range01(T1$Max)
      
      T2<-subset(data,data$trial=="T2")
      T2$Mean<-range01(T2$Mean)
      T2$IntDen<-range01(T2$IntDen)
      T2$Mode<-range01(T2$Mode)
      T2$Min<-range01(T2$Min)
      T2$Max<-range01(T2$Max)
      
      T3<-subset(data,data$trial=="T3")
      T3$Mean<-range01(T3$Mean)
      T3$IntDen<-range01(T3$IntDen)
      T3$Mode<-range01(T3$Mode)
      T3$Min<-range01(T3$Min)
      T3$Max<-range01(T3$Max)
      
      T4<-subset(data,data$trial=="T4")
      T4$Mean<-range01(T4$Mean)
      T4$IntDen<-range01(T4$IntDen)
      T4$Mode<-range01(T4$Mode)
      T4$Min<-range01(T4$Min)
      T4$Max<-range01(T4$Max)
      
      
      return(rbind(T1,T2,T3,T4))
    }
  }
  calclda <- function(variables,loadings)
  {
    # find the number of samples in the data set
    as.data.frame(variables)
    numsamples <- nrow(variables)
    # make a vector to store the discriminant function
    ld <- numeric(numsamples)
    # find the number of variables
    numvariables <- length(variables)
    # calculate the value of the discriminant function for each sample
    for (i in 1:numsamples)
    {
      valuei <- 0
      for (j in 1:numvariables)
      {
        valueij <- variables[i,j]
        loadingj <- loadings[j]
        valuei <- valuei + (valueij * loadingj)
      }
      ld[i] <- valuei
    }
    # standardise the discriminant function so that its mean value is 0:
    ld <- as.data.frame(scale(ld, center=TRUE, scale=FALSE))
    ld <- ld[[1]]
    return(ld)
  }
  
}
#parameters
{
  
  
  parameters_2d_int_data<-c("Mean","StdDev","Mode","Min","Max","IntDen","Median","Skew","Kurt","RawIntDen")
  labelling_2d_int_data<-paste("Proj.Spheroid Level",parameters_2d_int_data,sep="\n")
  nuc_labelling_2d_int_data<-paste("Proj.Nuclear Level",parameters_2d_int_data,sep="\n")
  
  parameters_int_data<-c("AtCenter","IntDen","Min","Max","Mean","Sigma","Mode","Mode.NonZero")
  labelling_int_data<-paste("Spheroid Level",parameters_int_data,sep="\n")
  nuclabelling_int_data<-paste("Nuclear Level",parameters_int_data,sep="\n")
  celllabelling_int_data<-paste("Cellular Level",parameters_int_data,sep="\n")
  
  parameters_geo_texture_data<-c( "RatioVolbox","Vol..unit.","Surf..unit." ,"Comp..unit.","Spher..unit.",
                                  "Feret..unit.","Ell_MajRad","Ell_Elon","Ell_Flatness","RatioVolEllipsoid",
                                  "Moment1","Moment2","Moment3","Moment4","Moment5","DCMin..unit.",
                                  "DCMax..unit.","DCMean..unit.","DCSD..unit.",
                                  "Area","Perim.","Major","Minor","Circ.","Feret","MinFeret","AR","Round",
                                  
                                  "HCcontent","ECcontent","HCvolume","ECvolume","HC_EC_content","HC_EC_volume",
                                  "EFD1","EFD2","EFD3","EFD4","EFD5","EFD6","EFD7","EFD8","EFD9","EFD10",                                  
                                  
                                  "Deg_0_Angular_Second_Moment._step5","Deg_0_Contrast._step5","Deg_0_Correlation._step5",
                                  "Deg_0_Inverse_Difference_Moment._step5","Deg_0_Entropy._step5",
                                  "Deg_90_Angular_Second_Moment._step5","Deg_90_Contrast._step5","Deg_90_Correlation._step5",
                                  "Deg_90_Inverse_Difference_Moment._step5","Deg_90_Entropy._step5",
                                  "Deg_180_Angular_Second_Moment._step5","Deg_180_Contrast._step5", "Deg_180_Correlation._step5",
                                  "Deg_180_Inverse_Difference_Moment._step5","Deg_180_Entropy._step5",                   
                                  "Deg_270_Angular_Second_Moment._step5","Deg_270_Contrast._step5","Deg_270_Correlation._step5",
                                  "Deg_270_Inverse_Difference_Moment._step5","Deg_270_Entropy._step5",
                                  
                                  "Deg_0_Angular_Second_Moment._step15","Deg_0_Contrast._step15","Deg_0_Correlation._step15",
                                  "Deg_0_Inverse_Difference_Moment._step15" ,"Deg_0_Entropy._step15",
                                  "Deg_90_Angular_Second_Moment._step15","Deg_90_Contrast._step15","Deg_90_Correlation._step15",
                                  "Deg_90_Inverse_Difference_Moment._step15","Deg_90_Entropy._step15",                    
                                  "Deg_180_Angular_Second_Moment._step15","Deg_180_Contrast._step15","Deg_180_Correlation._step15",              
                                  "Deg_180_Inverse_Difference_Moment._step15","Deg_180_Entropy._step15",
                                  "Deg_270_Angular_Second_Moment._step15","Deg_270_Contrast._step15","Deg_270_Correlation._step15",
                                  "Deg_270_Inverse_Difference_Moment._step15","Deg_270_Entropy._step15",
                                  
                                  "Deg_0_Angular_Second_Moment._step25","Deg_0_Contrast._step25","Deg_0_Correlation._step25",
                                  "Deg_0_Inverse_Difference_Moment._step25","Deg_0_Entropy._step25",
                                  "Deg_90_Angular_Second_Moment._step25","Deg_90_Contrast._step25","Deg_90_Correlation._step25",               
                                  "Deg_90_Inverse_Difference_Moment._step25","Deg_90_Entropy._step25",
                                  "Deg_180_Angular_Second_Moment._step25","Deg_180_Contrast._step25","Deg_180_Correlation._step25",
                                  "Deg_180_Inverse_Difference_Moment._step25","Deg_180_Entropy._step25",
                                  "Deg_270_Angular_Second_Moment._step25","Deg_270_Contrast._step25",               
                                  "Deg_270_Correlation._step25","Deg_270_Inverse_Difference_Moment._step25","Deg_270_Entropy._step25",
                                  
                                  "Deg_0_Angular_Second_Moment._step35","Deg_0_Contrast._step35","Deg_0_Correlation._step35",                
                                  "Deg_0_Inverse_Difference_Moment._step35","Deg_0_Entropy._step35",
                                  "Deg_90_Angular_Second_Moment._step35","Deg_90_Contrast._step35","Deg_90_Correlation._step35",  
                                  "Deg_90_Inverse_Difference_Moment._step35","Deg_90_Entropy._step35",
                                  "Deg_180_Angular_Second_Moment._step35","Deg_180_Contrast._step35","Deg_180_Correlation._step35",              
                                  "Deg_180_Inverse_Difference_Moment._step35","Deg_180_Entropy._step35",
                                  "Deg_270_Angular_Second_Moment._step35","Deg_270_Contrast._step35","Deg_270_Correlation._step35",
                                  "Deg_270_Inverse_Difference_Moment._step35","Deg_270_Entropy._step35",
                                  
                                  "Deg_0_Angular_Second_Moment._step45","Deg_0_Contrast._step45","Deg_0_Correlation._step45",
                                  "Deg_0_Inverse_Difference_Moment._step45","Deg_0_Entropy._step45",
                                  "Deg_90_Angular_Second_Moment._step45","Deg_90_Contrast._step45","Deg_90_Correlation._step45",
                                  "Deg_90_Inverse_Difference_Moment._step45","Deg_90_Entropy._step45",                   
                                  "Deg_180_Angular_Second_Moment._step45","Deg_180_Contrast._step45","Deg_180_Correlation._step45",               
                                  "Deg_180_Inverse_Difference_Moment._step45","Deg_180_Entropy._step45",
                                  "Deg_270_Angular_Second_Moment._step45","Deg_270_Contrast._step45","Deg_270_Correlation._step45",
                                  "Deg_270_Inverse_Difference_Moment._step45","Deg_270_Entropy._step45",
                                  
                                  "Deg_0_Angular_Second_Moment._step100","Deg_0_Contrast._step100","Deg_0_Correlation._step100",
                                  "Deg_0_Inverse_Difference_Moment._step100","Deg_0_Entropy._step100",                                     
                                  "Deg_90_Angular_Second_Moment._step100","Deg_90_Contrast._step100","Deg_90_Correlation._step100",
                                  "Deg_90_Inverse_Difference_Moment._step100" ,"Deg_90_Entropy._step100",
                                  "Deg_180_Angular_Second_Moment._step100","Deg_180_Contrast._step100","Deg_180_Correlation._step100",
                                  "Deg_180_Inverse_Difference_Moment._step100","Deg_180_Entropy._step100",               
                                  "Deg_270_Angular_Second_Moment._step100","Deg_270_Contrast._step100","Deg_270_Correlation._step100",            
                                  "Deg_270_Inverse_Difference_Moment._step100","Deg_270_Entropy._step100")
  
  labelling_geo_texture_data<-c("Bounding Box Ratio","Volume", "Surface Area","Compaction","Sphericity",
                                "Feret (3D)", "Ellipsoid R1","Ellipsoid Ellongation","Ellipsoid Flatness", "Ellipsoid Ratio",
                                "Moment1","Moment2","Moment3","Moment4","Moment5", 
                                "DC_Min.","DC_Max","DC_Mean.","DC_SD",
                                "Proj. Area", "Proj. Perimeter","Proj. Major Axis", "Proj. Minor Axis", "Proj. Circularity",
                                "Proj. Feret","Proj. Min Feret","Proj. A.R","Proj. Roundness","Proj. Solidity",
                                
                                "HCcontent","ECcontent","HCvolume","ECvolume","HC_EC_content","HC_EC_volume",
                                "EFD1","EFD2","EFD3","EFD4","EFD5","EFD6","EFD7","EFD8","EFD9","EFD10",                                  
                                
                                "Deg_0_Angular_Second_Moment _step5","Deg_0_Contrast _step5","Deg_0_Correlation _step5",
                                "Deg_0_Inverse_Difference_Moment _step5","Deg_0_Entropy _step5",
                                "Deg_90_Angular_Second_Moment _step5","Deg_90_Contrast _step5","Deg_90_Correlation _step5",
                                "Deg_90_Inverse_Difference_Moment _step5","Deg_90_Entropy _step5",
                                "Deg_180_Angular_Second_Moment _step5","Deg_180_Contrast _step5", "Deg_180_Correlation _step5",
                                "Deg_180_Inverse_Difference_Moment _step5","Deg_180_Entropy _step5",                   
                                "Deg_270_Angular_Second_Moment _step5","Deg_270_Contrast _step5","Deg_270_Correlation _step5",
                                "Deg_270_Inverse_Difference_Moment _step5","Deg_270_Entropy _step5",
                                
                                "Deg_0_Angular_Second_Moment _step15","Deg_0_Contrast _step15","Deg_0_Correlation _step15",
                                "Deg_0_Inverse_Difference_Moment _step15" ,"Deg_0_Entropy _step15",
                                "Deg_90_Angular_Second_Moment _step15","Deg_90_Contrast _step15","Deg_90_Correlation _step15",
                                "Deg_90_Inverse_Difference_Moment _step15","Deg_90_Entropy _step15",                   
                                "Deg_180_Angular_Second_Moment _step15","Deg_180_Contrast _step15","Deg_180_Correlation _step15",              
                                "Deg_180_Inverse_Difference_Moment _step15","Deg_180_Entropy _step15",
                                "Deg_270_Angular_Second_Moment _step15","Deg_270_Contrast _step15","Deg_270_Correlation _step15",
                                "Deg_270_Inverse_Difference_Moment _step15","Deg_270_Entropy _step15",
                                
                                "Deg_0_Angular_Second_Moment _step25","Deg_0_Contrast _step25","Deg_0_Correlation _step25",
                                "Deg_0_Inverse_Difference_Moment _step25","Deg_0_Entropy _step25",
                                "Deg_90_Angular_Second_Moment _step25","Deg_90_Contrast _step25","Deg_90_Correlation _step25",               
                                "Deg_90_Inverse_Difference_Moment _step25","Deg_90_Entropy _step25",
                                "Deg_180_Angular_Second_Moment _step25","Deg_180_Contrast _step25","Deg_180_Correlation _step25",
                                "Deg_180_Inverse_Difference_Moment _step25","Deg_180_Entropy _step25",
                                "Deg_270_Angular_Second_Moment _step25","Deg_270_Contrast _step25",               
                                "Deg_270_Correlation _step25","Deg_270_Inverse_Difference_Moment _step25","Deg_270_Entropy _step25",
                                
                                "Deg_0_Angular_Second_Moment _step35","Deg_0_Contrast _step35","Deg_0_Correlation _step35",                
                                "Deg_0_Inverse_Difference_Moment _step35","Deg_0_Entropy _step35",
                                "Deg_90_Angular_Second_Moment _step35","Deg_90_Contrast _step35","Deg_90_Correlation _step35",  
                                "Deg_90_Inverse_Difference_Moment _step35","Deg_90_Entropy _step35",
                                "Deg_180_Angular_Second_Moment _step35","Deg_180_Contrast _step35","Deg_180_Correlation _step35",              
                                "Deg_180_Inverse_Difference_Moment _step35","Deg_180_Entropy _step35",
                                "Deg_270_Angular_Second_Moment _step35","Deg_270_Contrast _step35","Deg_270_Correlation _step35",
                                "Deg_270_Inverse_Difference_Moment _step35","Deg_270_Entropy _step35",
                                
                                "Deg_0_Angular_Second_Moment _step45","Deg_0_Contrast _step45","Deg_0_Correlation _step45",
                                "Deg_0_Inverse_Difference_Moment _step45","Deg_0_Entropy _step45",
                                "Deg_90_Angular_Second_Moment _step45","Deg_90_Contrast _step45","Deg_90_Correlation _step45",
                                "Deg_90_Inverse_Difference_Moment _step45","Deg_90_Entropy _step45",                   
                                "Deg_180_Angular_Second_Moment _step45","Deg_180_Contrast _step45","Deg_180_Correlation _step45",               
                                "Deg_180_Inverse_Difference_Moment _step45","Deg_180_Entropy _step45",
                                "Deg_270_Angular_Second_Moment _step45","Deg_270_Contrast _step45","Deg_270_Correlation _step45",
                                "Deg_270_Inverse_Difference_Moment _step45","Deg_270_Entropy _step45",
                                
                                "Deg_0_Angular_Second_Moment _step100","Deg_0_Contrast _step100","Deg_0_Correlation _step100",
                                "Deg_0_Inverse_Difference_Moment _step100","Deg_0_Entropy _step100",                                     
                                "Deg_90_Angular_Second_Moment _step100","Deg_90_Contrast _step100","Deg_90_Correlation _step100",
                                "Deg_90_Inverse_Difference_Moment _step100" ,"Deg_90_Entropy _step100",
                                "Deg_180_Angular_Second_Moment _step100","Deg_180_Contrast _step100","Deg_180_Correlation _step100",
                                "Deg_180_Inverse_Difference_Moment _step100","Deg_180_Entropy _step100",               
                                "Deg_270_Angular_Second_Moment _step100","Deg_270_Contrast _step100","Deg_270_Correlation _step100",            
                                "Deg_270_Inverse_Difference_Moment _step100","Deg_270_Entropy _step100")
  
}
#subset labels
{
  label_T1_set1<-subset(nucleus_combined,nucleus_combined$meg_trial=="T1_set1")$Label
  label_T2_set1<-subset(nucleus_combined,nucleus_combined$meg_trial=="T2_set1")$Label
  label_T3_set1<-subset(nucleus_combined,nucleus_combined$meg_trial=="T3_set1")$Label
  label_T1_set2<-subset(nucleus_combined,nucleus_combined$meg_trial=="T1_set2")$Label
  label_T2_set2<-subset(nucleus_combined,nucleus_combined$meg_trial=="T2_set2")$Label
  label_T3_set2<-subset(nucleus_combined,nucleus_combined$meg_trial=="T3_set2")$Label
  label_T1_set3<-subset(nucleus_combined,nucleus_combined$meg_trial=="T1_set3")$Label
  label_T2_set3<-subset(nucleus_combined,nucleus_combined$meg_trial=="T2_set3")$Label
  label_T3_set3<-subset(nucleus_combined,nucleus_combined$meg_trial=="T3_set3")$Label
  label_T1_set4<-subset(nucleus_combined,nucleus_combined$meg_trial=="T1_set4")$Label
  label_T2_set4<-subset(nucleus_combined,nucleus_combined$meg_trial=="T2_set4")$Label
  label_T3_set4<-subset(nucleus_combined,nucleus_combined$meg_trial=="T3_set4")$Label
  label_T1_set5<-subset(nucleus_combined,nucleus_combined$meg_trial=="T1_set5")$Label
  label_T2_set5<-subset(nucleus_combined,nucleus_combined$meg_trial=="T2_set5")$Label
  label_T3_set5<-subset(nucleus_combined,nucleus_combined$meg_trial=="T3_set5")$Label
  label_T4_set5<-subset(nucleus_combined,nucleus_combined$meg_trial=="T4_set5")$Label
  
  
  D0_lab<-subset(nucleus_combined,nucleus_combined$sample=="D0")$Label
  D2_lab<-subset(nucleus_combined,nucleus_combined$sample=="D2")$Label
  D4_lab<-subset(nucleus_combined,nucleus_combined$sample=="D4")$Label
  D6_lab<-subset(nucleus_combined,nucleus_combined$sample=="D6")$Label
  D8_lab<-subset(nucleus_combined,nucleus_combined$sample=="D8")$Label
  
}
cols_timeline<-colorRampPalette(c("red","yellow"))( 5) 
cols_trials<-colorRampPalette(c("gray","black")) (length(levels(as.factor(nucleus_combined$meg_trial))))
#read nucleus combined data from 5 mega sets 
{
  setwd("~/Desktop/Reprogramming/actin_pmlc_oct4/combined_data/")
  
  nucleus_combined_set1 <- read.csv("nucleus_combined.csv",stringsAsFactors = F)
  nucleus_combined_set1<-nucleus_combined_set1[!(duplicated(nucleus_combined_set1$Label)),]
  nucleus_combined_set1$meg_trial<-paste(nucleus_combined_set1$trial,"set1",sep="_")
  
  setwd("~/Desktop/Reprogramming/vimentin_ecad_oct4/combined_data/")
  
  nucleus_combined_set2 <- read.csv("nucleus_combined.csv",stringsAsFactors = F)
  nucleus_combined_set2<-nucleus_combined_set2[!(duplicated(nucleus_combined_set2$Label)),]
  nucleus_combined_set2$meg_trial<-paste(nucleus_combined_set2$trial,"set2",sep="_")
  
  setwd("~/Desktop/Reprogramming/ki67/combined_data/")
  
  nucleus_combined_set3 <- read.csv("nucleus_combined.csv",stringsAsFactors = F)
  nucleus_combined_set3<-nucleus_combined_set3[!(duplicated(nucleus_combined_set3$Label)),]
  nucleus_combined_set3<-nucleus_combined_set3[,which(colnames(nucleus_combined_set3) %in% colnames(nucleus_combined_set2))]
  nucleus_combined_set3$meg_trial<-paste(nucleus_combined_set3$trial,"set3",sep="_")
  
  
  setwd("~/Desktop/Reprogramming/H3k9ac/combined_data/")
  
  nucleus_combined_set4 <- read.csv("nucleus_combined.csv",stringsAsFactors = F)
  nucleus_combined_set4<-nucleus_combined_set4[!(duplicated(nucleus_combined_set4$Label)),]
  nucleus_combined_set4$meg_trial<-paste(nucleus_combined_set4$trial,"set4",sep="_")
  
  setwd("~/Desktop/Reprogramming/actin_laminac_oct4/combined_data/")
  
  nucleus_combined_set5 <- read.csv("nucleus_combined.csv",stringsAsFactors = F)
  nucleus_combined_set5<-nucleus_combined_set5[!(duplicated(nucleus_combined_set5$Label)),]
  nucleus_combined_set5$meg_trial<-paste(nucleus_combined_set5$trial,"set5",sep="_")
  
  nucleus_combined<-rbind(nucleus_combined_set1,nucleus_combined_set2,nucleus_combined_set3,nucleus_combined_set4,nucleus_combined_set5)
  nucleus_combined<-nucleus_combined[!(duplicated(nucleus_combined$Label)),]
  rownames(nucleus_combined)<-nucleus_combined$Label
  
  rm(nucleus_combined_set1,nucleus_combined_set2,nucleus_combined_set3,nucleus_combined_set4,nucleus_combined_set5)
  
}
#nuc_int_data
{
  setwd("~/Desktop/Reprogramming/actin_pmlc_oct4/combined_data/")
  
  actin <- read.csv("cellular_int_actin.csv",stringsAsFactors = F)
  actin$Label<-substring(actin$Label, 5, nchar(actin$Label)-4)
  actin<-actin[!(duplicated(actin$Label)),]
  actin<-normalize_mean_trials(actin)
  
  oct4_set1 <- read.csv("nuclear_2d_int_oct4.csv",stringsAsFactors = F)
  oct4_set1<-oct4_set1[!(duplicated(oct4_set1$Label)),]
  oct4_set1<-normalize_2d_mean_trials(oct4_set1)
  
  pmlc <- read.csv("cellular_int_pmlc.csv",stringsAsFactors = F)
  pmlc<-pmlc[!(duplicated(pmlc$Label)),]
  pmlc<-normalize_mean_trials(pmlc)
  
  setwd("~/Desktop/Reprogramming/ki67/combined_data/")
  ki67 <- read.csv("nuclear_int_KI67.csv",stringsAsFactors = F)
  ki67<-ki67[!(duplicated(ki67$Label)),]
  ki67<-normalize_mean_trials(ki67)
  
  setwd("~/Desktop/Reprogramming/H3k9ac/combined_data/")
  h3k9ac <- read.csv("nuclear_2d_int_H3K9AC.csv",stringsAsFactors = F)
  h3k9ac<-h3k9ac[!(duplicated(h3k9ac$Label)),]
  h3k9ac<-normalize_2d_mean_trials(h3k9ac)
  
  setwd("~/Desktop/Reprogramming/actin_laminac_oct4/combined_data/")
  laminac <- read.csv("cellular_int_LAMINAC.csv",stringsAsFactors = F)
  laminac<-laminac[!(duplicated(laminac$Label)),]
  laminac<-normalize_mean_trials(laminac)
  
  oct4_set2 <- read.csv("nuclear_2d_int_oct4.csv",stringsAsFactors = F)
  oct4_set2<-oct4_set2[!(duplicated(oct4_set2$Label)),]
  oct4_set2<-normalize_2d_mean_trials(oct4_set2)
  
  setwd("~/Desktop/Reprogramming/vimentin_ecad_oct4/combined_data/")
  
  vimentin <- read.csv("cellular_int_VIMENTIN.csv",stringsAsFactors = F)
  vimentin<-vimentin[!(duplicated(vimentin$Label)),]
  vimentin<-normalize_mean_trials(vimentin)
  
  ecad <- read.csv("cellular_int_ECAD.csv",stringsAsFactors = F)
  ecad<-ecad[!(duplicated(ecad$Label)),]
  ecad<-normalize_mean_trials(ecad)
  
  oct4_set3 <- read.csv("nuclear_2d_int_oct4.csv",stringsAsFactors = F)
  oct4_set3<-oct4_set3[!(duplicated(oct4_set3$Label)),]
  oct4_set3<-normalize_2d_mean_trials(oct4_set3)
  
  oct4<-rbind(oct4_set1,oct4_set2,oct4_set3)
  
  nuc_int_data<-list(actin=actin,oct4=oct4,pmlc=pmlc,
                     vimentin=vimentin,ecad=ecad,
                     ki67=ki67,h3k9ac=h3k9ac,
                     laminac=laminac)
  rm(actin,oct4,pmlc,vimentin,ecad,ki67,h3k9ac,laminac,oct4_set1,oct4_set2,oct4_set3)
  
}
dir.create("/Users/saradhavenkatachalapathy/Desktop/reprogramming/timeline_nucleus_analysis_scaled/equal_sampling_oct4/")
setwd("/Users/saradhavenkatachalapathy/Desktop/reprogramming/timeline_nucleus_analysis_scaled/equal_sampling_oct4/")
dird<-"/Users/saradhavenkatachalapathy/Desktop/reprogramming/timeline_nucleus_analysis_scaled/equal_sampling_oct4/"
#LDA oct4 intensity
{
  #deciles
  {
    oct4_int_selec<-subset(nuc_int_data$oct4, !nuc_int_data$oct4$Label %in% c(label_T1_set2,label_T2_set2,label_T1_set5,label_T1_set1,label_T2_set1))
    decile_oct4<- quantile(oct4_int_selec$Median, prob = seq(0, 1, length = 11), type = 5)
    oct4_labels_0_10<-oct4_int_selec$Label[which(oct4_int_selec$Median >= decile_oct4[1] & oct4_int_selec$Median < decile_oct4[2] )]
    oct4_labels_10_20<-oct4_int_selec$Label[which(oct4_int_selec$Median >= decile_oct4[2] & oct4_int_selec$Median < decile_oct4[3] )]
    oct4_labels_20_30<-oct4_int_selec$Label[which(oct4_int_selec$Median >= decile_oct4[3] & oct4_int_selec$Median < decile_oct4[4] )]
    oct4_labels_30_40<-oct4_int_selec$Label[which(oct4_int_selec$Median >= decile_oct4[4] & oct4_int_selec$Median < decile_oct4[5] )]
    oct4_labels_40_50<-oct4_int_selec$Label[which(oct4_int_selec$Median >= decile_oct4[5] & oct4_int_selec$Median < decile_oct4[6] )]
    oct4_labels_50_60<-oct4_int_selec$Label[which(oct4_int_selec$Median >= decile_oct4[6] & oct4_int_selec$Median < decile_oct4[7] )]
    oct4_labels_60_70<-oct4_int_selec$Label[which(oct4_int_selec$Median >= decile_oct4[7] & oct4_int_selec$Median < decile_oct4[8] )]
    oct4_labels_70_80<-oct4_int_selec$Label[which(oct4_int_selec$Median >= decile_oct4[8] & oct4_int_selec$Median < decile_oct4[9] )]
    oct4_labels_80_90<-oct4_int_selec$Label[which(oct4_int_selec$Median >= decile_oct4[9] & oct4_int_selec$Median < decile_oct4[10] )]
    oct4_labels_90_100<-oct4_int_selec$Label[which(oct4_int_selec$Median >= decile_oct4[10] & oct4_int_selec$Median < decile_oct4[11] )]
    
    colfunc <- colorRampPalette(c("#FF3E96","#FFEC8B","#21840d" ))
    cols<-colfunc(11)
    
    png(filename="intensity_oct4.png", units="in", width=1, height=1 , pointsize=5,res=1200)
    par(font.axis=2,font.lab=2,mar = c(4,4,1,1), font=2)
    d<-density(oct4_int_selec$Mean)
    plot(d, main="",las=1, xlab="Nuclear Oct4 levels",cex.axis=0.7)
    dec<-quantile(oct4_int_selec$Mean, prob = seq(0, 1, length = 11), type = 5)
    # Add the shaded area.
    polygon(c( dec[11],d$x[d$x>=dec[11]], max(d$x) ),  c(0,d$y[d$x>=dec[11]],0), col=cols[11],lwd=0.5,border = "black")
    polygon(c( dec[10],d$x[d$x>=dec[10] & d$x<dec[11]],dec[11] ),  c(0,d$y[d$x>=dec[10] & d$x<dec[11]],0), col=cols[11],lwd=0.5,border = "black")
    polygon(c( dec[9],d$x[d$x>=dec[9] & d$x<dec[10]],dec[10] ),  c(0,d$y[d$x>=dec[9] & d$x<dec[10]],0), col=cols[10],lwd=0.5,border = "black")
    polygon(c( dec[8],d$x[d$x>=dec[8] & d$x<dec[9]],dec[9] ),  c(0,d$y[d$x>=dec[8] & d$x<dec[9]],0), col=cols[9],lwd=0.5,border = "black")
    polygon(c( dec[7],d$x[d$x>=dec[7] & d$x<dec[8]],dec[8] ),  c(0,d$y[d$x>=dec[7] & d$x<dec[8]],0), col=cols[8],lwd=0.5,border = "black")
    polygon(c( dec[6],d$x[d$x>=dec[6] & d$x<dec[7]],dec[7] ),  c(0,d$y[d$x>=dec[6] & d$x<dec[7]],0), col=cols[7],lwd=0.5,border = "black")
    polygon(c( dec[5],d$x[d$x>=dec[5] & d$x<dec[6]],dec[6] ),  c(0,d$y[d$x>=dec[5] & d$x<dec[6]],0), col=cols[6],lwd=0.5,border = "black")
    polygon(c( dec[4],d$x[d$x>=dec[4] & d$x<dec[5]],dec[5] ),  c(0,d$y[d$x>=dec[4] & d$x<dec[5]],0), col=cols[5],lwd=0.5,border = "black")
    polygon(c( dec[3],d$x[d$x>=dec[3] & d$x<dec[4]],dec[4] ),  c(0,d$y[d$x>=dec[3] & d$x<dec[4]],0), col=cols[4],lwd=0.5,border = "black")
    polygon(c( dec[2],d$x[d$x>=dec[2] & d$x<dec[3]],dec[3] ),  c(0,d$y[d$x>=dec[2] & d$x<dec[3]],0), col=cols[3],lwd=0.5,border = "black")
    polygon(c( dec[1],d$x[d$x>=dec[1] & d$x<dec[2]],dec[2] ),  c(0,d$y[d$x>=dec[1] & d$x<dec[2]],0), col=cols[2],lwd=0.5,border = "black")
    polygon(c( min(d$x),d$x[d$x<dec[1]],dec[1]),  c(0,d$y[d$x<dec[1]],0), col=cols[1],lwd=0.5,border = "black")  
    dev.off()
    
    a<-subset(oct4_int_selec,oct4_int_selec$Label %in% D0_lab)
    b<-subset(oct4_int_selec,oct4_int_selec$Label %in% D2_lab)
    c<-subset(oct4_int_selec,oct4_int_selec$Label %in% D4_lab)
    d<-subset(oct4_int_selec,oct4_int_selec$Label %in% D6_lab)
    e<-subset(oct4_int_selec,oct4_int_selec$Label %in% D8_lab)
    
    boxplot(a$Median,b$Median,c$Median,d$Median,e$Median)
    rm(a,b,c,d,e)
    
  }
  #PCA/LDA
  {
    nucleus_subset<-nucleus_combined[,which(colnames(nucleus_combined) %in% parameters_geo_texture_data)]
    
    scaled_nucleus_subset<-subset(nucleus_subset,rownames(nucleus_subset) %in% label_T1_set1)
    scaled_nucleus_subset<-scale(scaled_nucleus_subset,center = T, scale = T)
    temp2<-subset(nucleus_subset,rownames(nucleus_subset) %in% label_T2_set1)
    temp2<-scale(temp2,center = T, scale = T)
    scaled_nucleus_subset<-rbind(scaled_nucleus_subset,temp2)
    temp2<-subset(nucleus_subset,rownames(nucleus_subset) %in% label_T3_set1)
    temp2<-scale(temp2,center = T, scale = T)
    scaled_nucleus_subset<-rbind(scaled_nucleus_subset,temp2)
    temp2<-subset(nucleus_subset,rownames(nucleus_subset) %in% label_T1_set2)
    temp2<-scale(temp2,center = T, scale = T)
    scaled_nucleus_subset<-rbind(scaled_nucleus_subset,temp2)
    temp2<-subset(nucleus_subset,rownames(nucleus_subset) %in% label_T2_set2)
    temp2<-scale(temp2,center = T, scale = T)
    scaled_nucleus_subset<-rbind(scaled_nucleus_subset,temp2)
    temp2<-subset(nucleus_subset,rownames(nucleus_subset) %in% label_T3_set2)
    temp2<-scale(temp2,center = T, scale = T)
    scaled_nucleus_subset<-rbind(scaled_nucleus_subset,temp2)
    temp2<-subset(nucleus_subset,rownames(nucleus_subset) %in% label_T1_set3)
    temp2<-scale(temp2,center = T, scale = T)
    scaled_nucleus_subset<-rbind(scaled_nucleus_subset,temp2)
    temp2<-subset(nucleus_subset,rownames(nucleus_subset) %in% label_T2_set3)
    temp2<-scale(temp2,center = T, scale = T)
    scaled_nucleus_subset<-rbind(scaled_nucleus_subset,temp2)
    temp2<-subset(nucleus_subset,rownames(nucleus_subset) %in% label_T3_set3)
    temp2<-scale(temp2,center = T, scale = T)
    scaled_nucleus_subset<-rbind(scaled_nucleus_subset,temp2)
    temp2<-subset(nucleus_subset,rownames(nucleus_subset) %in% label_T1_set4)
    temp2<-scale(temp2,center = T, scale = T)
    scaled_nucleus_subset<-rbind(scaled_nucleus_subset,temp2)
    temp2<-subset(nucleus_subset,rownames(nucleus_subset) %in% label_T2_set4)
    temp2<-scale(temp2,center = T, scale = T)
    scaled_nucleus_subset<-rbind(scaled_nucleus_subset,temp2)
    temp2<-subset(nucleus_subset,rownames(nucleus_subset) %in% label_T3_set4)
    temp2<-scale(temp2,center = T, scale = T)
    scaled_nucleus_subset<-rbind(scaled_nucleus_subset,temp2)
    temp2<-subset(nucleus_subset,rownames(nucleus_subset) %in% label_T1_set5)
    temp2<-scale(temp2,center = T, scale = T)
    scaled_nucleus_subset<-rbind(scaled_nucleus_subset,temp2)
    temp2<-subset(nucleus_subset,rownames(nucleus_subset) %in% label_T2_set5)
    temp2<-scale(temp2,center = T, scale = T)
    scaled_nucleus_subset<-rbind(scaled_nucleus_subset,temp2)
    temp2<-subset(nucleus_subset,rownames(nucleus_subset) %in% label_T3_set5)
    temp2<-scale(temp2,center = T, scale = T)
    scaled_nucleus_subset<-rbind(scaled_nucleus_subset,temp2)
    temp2<-subset(nucleus_subset,rownames(nucleus_subset) %in% label_T4_set5)
    temp2<-scale(temp2,center = T, scale = T)
    scaled_nucleus_subset<-rbind(scaled_nucleus_subset,temp2)
    
    nucleus_subset<-as.matrix(na.omit((scaled_nucleus_subset)))
    
    re<-as.data.frame(nucleus_subset)
    d1<-subset(re, rownames(re) %in% oct4_labels_0_10)
    d2<-subset(re, rownames(re) %in% oct4_labels_10_20)
    d3<-subset(re, rownames(re) %in% oct4_labels_20_30)
    d4<-subset(re, rownames(re) %in% oct4_labels_30_40)
    d5<-subset(re, rownames(re) %in% oct4_labels_40_50)
    d6<-subset(re, rownames(re) %in% oct4_labels_50_60)
    d7<-subset(re, rownames(re) %in% oct4_labels_60_70)
    d8<-subset(re, rownames(re) %in% oct4_labels_70_80)
    d9<-subset(re, rownames(re) %in% oct4_labels_80_90)
    d10<-subset(re, rownames(re) %in% oct4_labels_90_100)
    
    d1$samplem<-"0-10th"
    d2$samplem<-"10-20th"
    d3$samplem<-"20-30th"
    d4$samplem<-"30-40th"
    d5$samplem<-"40-50th"
    d6$samplem<-"50-60th"
    d7$samplem<-"60-70th"
    d8$samplem<-"70-80th"
    d9$samplem<-"80-90th"
    d10$samplem<-"90-100th"
    
    re<-rbind(d1, d2)
    re<-rbind(re, d3)
    re<-rbind(re, d4)
    re<-rbind(re, d5)
    re<-rbind(re, d6)
    re<-rbind(re, d7)
    re<-rbind(re, d8)
    re<-rbind(re, d9)
    re<-rbind(re, d10)
     
    re$Label<-row.names(re)
    lengtha<-min(nrow(d1),nrow(d2),nrow(d3),nrow(d4),nrow(d5),nrow(d6),nrow(d7),nrow(d8),nrow(d9),nrow(d10))
    
    trainingData<-ddply(re,.(samplem),function(x) x[sample(nrow(x),round(lengtha,0)),])
    
    rownames(trainingData)<-trainingData$Label
    nucleus_subset<-trainingData[,which(colnames(trainingData) %in% parameters_geo_texture_data)]
    nucleus_pca <- princomp(nucleus_subset,  cor = TRUE, scores = TRUE)
    pca_scores<-nucleus_pca$scores
    
    re<-as.data.frame(pca_scores)
    d1<-subset(re, rownames(re) %in% oct4_labels_0_10)
    d2<-subset(re, rownames(re) %in% oct4_labels_10_20)
    d3<-subset(re, rownames(re) %in% oct4_labels_20_30)
    d4<-subset(re, rownames(re) %in% oct4_labels_30_40)
    d5<-subset(re, rownames(re) %in% oct4_labels_40_50)
    d6<-subset(re, rownames(re) %in% oct4_labels_50_60)
    d7<-subset(re, rownames(re) %in% oct4_labels_60_70)
    d8<-subset(re, rownames(re) %in% oct4_labels_70_80)
    d9<-subset(re, rownames(re) %in% oct4_labels_80_90)
    d10<-subset(re, rownames(re) %in% oct4_labels_90_100)
    
    d1$samplem<-"0-10th"
    d2$samplem<-"10-20th"
    d3$samplem<-"20-30th"
    d4$samplem<-"30-40th"
    d5$samplem<-"40-50th"
    d6$samplem<-"50-60th"
    d7$samplem<-"60-70th"
    d8$samplem<-"70-80th"
    d9$samplem<-"80-90th"
    d10$samplem<-"90-100th"
    
    re<-rbind(d1, d2)
    re<-rbind(re, d3)
    re<-rbind(re, d4)
    re<-rbind(re, d5)
    re<-rbind(re, d6)
    re<-rbind(re, d7)
    re<-rbind(re, d8)
    re<-rbind(re, d9)
    re<-rbind(re, d10)
    
    png(filename="pca_biplot_1_2_ellipses.png", units="in",width=2, height=2 , pointsize=7, res=1200)
    par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2)
    ellipseplot(re[,1],                          # data for x-axis
                re[,2],                          # data for y-axis
                as.factor(re$samplem),                    # factor with classes
                elev=0.95,
                pcol=cols,    # colors for plotting (must match # of factors)
                pbgcol=TRUE,                           # point borders black?
                cexsize=0.1,                            # size of points 
                ppch=rep(21,10),                          # shape of points (must match # of factors)
                axissize=0.8,                           # Set axis text size
                linewidth=1,                          # Set axis line size
                font=1
    )
    title(xlab="PC1",    # % variance explained on PC1
          ylab="PC2",    # % variance explained on PC2 
          main=" ",                  # Title
          cex.lab=1,                    # size of label text
          cex.main=1                    # size of title text
    )
    dev.off()
    
    png(filename="pca_biplot_1_3_ellipses.png", units="in",width=2, height=2 , pointsize=7, res=1200)
    par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2)
    ellipseplot(re[,1],                          # data for x-axis
                re[,3],                          # data for y-axis
                as.factor(re$samplem),                    # factor with classes
                elev=0.95,
                pcol=cols,    # colors for plotting (must match # of factors)
                pbgcol=TRUE,                           # point borders black?
                cexsize=0.1,                            # size of points 
                ppch=rep(21,10),                          # shape of points (must match # of factors)
                axissize=0.8,                           # Set axis text size
                linewidth=1,                          # Set axis line size
                font=1
    )
    title(xlab="PC1",    # % variance explained on PC1
          ylab="PC3",    # % variance explained on PC2 
          main=" ",                  # Title
          cex.lab=1,                    # size of label text
          cex.main=1                    # size of title text
    )
    dev.off()
    
    png(filename="pca_biplot_2_3_ellipses.png", units="in",width=2, height=2 , pointsize=7, res=1200)
    par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2)
    ellipseplot(re[,2],                          # data for x-axis
                re[,3],                          # data for y-axis
                as.factor(re$samplem),                    # factor with classes
                elev=0.95,
                pcol=cols,    # colors for plotting (must match # of factors)
                pbgcol=TRUE,                           # point borders black?
                cexsize=0.1,                            # size of points 
                ppch=rep(21,10),                          # shape of points (must match # of factors)
                axissize=0.8,                           # Set axis text size
                linewidth=1,                          # Set axis line size
                font=1
    )
    title(xlab="PC2",    # % variance explained on PC1
          ylab="PC3",    # % variance explained on PC2 
          main=" ",                  # Title
          cex.lab=1,                    # size of label text
          cex.main=1                    # size of title text
    )
    dev.off()
    
    {
      {
        edf <- data.frame(LV1 = re[,1], LV2=re[,2], factr = re$samplem) # create data frame with data and factor
        ellipses <- dlply(edf, .(factr), function(x) {
          LV1 <- x[,1]
          LV2 <- x[,2]
          dataEllipse(LV1, LV2, levels=0.95, robust=TRUE, draw=FALSE) # Get confidence ellipse points from dataEllipse() function by factor level
        })
        col1 <- colorRampPalette(c("#7F0000", "#7F0000","#7F0000","#7F0000","#7F0000","#7F0000","#7F0000","#7F0000","#7F0000"
                                   ,"red", "#FF7F00", "yellow", "#7FFF7F",
                                   "cyan", "#007FFF", "blue", "#00007F"))
        
        area_int<-as.data.frame(matrix(NA, nrow=10,ncol=10))
        colnames(area_int)<-rev(attributes(ellipses)$names)
        rownames(area_int)<-(attributes(ellipses)$names)
        for(k in 1:10){
          for(o in 10:1){
            ell1 <- xy2SP(ellipses[[k]])
            ell2 <- xy2SP(ellipses[[11-o]])
            a<-gIntersection(ell1, ell2)
            if(length(a)>0){
              area_int[k,o]<-gArea(gIntersection(ell1, ell2))/(gArea(ell1)+gArea(ell2)-gArea(gIntersection(ell1, ell2)))
            }
            if(length(a)==0){
              area_int[k,o]<-0
            }
            
          }
        }
        png(filename="corrplot_PC_1_2.png", units="in",width=2, height=2 , pointsize=5, res=1200)
        par(font.axis = 1,font.lab=1,mar=c(1,1,1,2), font=2)
        corrplot(t(as.matrix((area_int))),tl.col = "black",cl.lim = c(0, 1),method = "color",is.corr=F, col=col1(250),cl.pos = "b",tl.cex = 0.8, cl.ratio = 0.3)
        dev.off()
        
        write.csv(area_int,file="area_intersection_1_2.csv")
        
      }
      {
        edf <- data.frame(LV1 = re[,1], LV2=re[,3], factr = re$samplem) # create data frame with data and factor
        ellipses <- dlply(edf, .(factr), function(x) {
          LV1 <- x[,1]
          LV2 <- x[,2]
          dataEllipse(LV1, LV2, levels=0.95, robust=TRUE, draw=FALSE) # Get confidence ellipse points from dataEllipse() function by factor level
        })
        col1 <- colorRampPalette(c("#7F0000", "#7F0000","#7F0000","#7F0000","#7F0000","#7F0000","#7F0000","#7F0000","#7F0000"
                                   ,"red", "#FF7F00", "yellow", "#7FFF7F",
                                   "cyan", "#007FFF", "blue", "#00007F"))
        
        area_int<-as.data.frame(matrix(NA, nrow=10,ncol=10))
        colnames(area_int)<-rev(attributes(ellipses)$names)
        rownames(area_int)<-(attributes(ellipses)$names)
        for(k in 1:10){
          for(o in 10:1){
            ell1 <- xy2SP(ellipses[[k]])
            ell2 <- xy2SP(ellipses[[11-o]])
            a<-gIntersection(ell1, ell2)
            if(length(a)>0){
              area_int[k,o]<-gArea(gIntersection(ell1, ell2))/(gArea(ell1)+gArea(ell2)-gArea(gIntersection(ell1, ell2)))
            }
            if(length(a)==0){
              area_int[k,o]<-0
            }
            
          }
        }
        png(filename="corrplot_PC_1_3.png", units="in",width=2, height=2 , pointsize=5, res=1200)
        par(font.axis = 1,font.lab=1,mar=c(1,1,1,2), font=2)
        corrplot(t(as.matrix((area_int))),tl.col = "black",cl.lim = c(0, 1),method = "color",is.corr=F, col=col1(250),cl.pos = "b",tl.cex = 0.8, cl.ratio = 0.3)
        dev.off()
        
        write.csv(area_int,file="area_intersection_1_3.csv")
        
      }
      {
        edf <- data.frame(LV1 = re[,2], LV2=re[,3], factr = re$samplem) # create data frame with data and factor
        ellipses <- dlply(edf, .(factr), function(x) {
          LV1 <- x[,1]
          LV2 <- x[,2]
          dataEllipse(LV1, LV2, levels=0.95, robust=TRUE, draw=FALSE) # Get confidence ellipse points from dataEllipse() function by factor level
        })
        col1 <- colorRampPalette(c("#7F0000", "#7F0000","#7F0000","#7F0000","#7F0000","#7F0000","#7F0000","#7F0000","#7F0000"
                                   ,"red", "#FF7F00", "yellow", "#7FFF7F",
                                   "cyan", "#007FFF", "blue", "#00007F"))
        
        area_int<-as.data.frame(matrix(NA, nrow=10,ncol=10))
        colnames(area_int)<-rev(attributes(ellipses)$names)
        rownames(area_int)<-(attributes(ellipses)$names)
        for(k in 1:10){
          for(o in 10:1){
            ell1 <- xy2SP(ellipses[[k]])
            ell2 <- xy2SP(ellipses[[11-o]])
            a<-gIntersection(ell1, ell2)
            if(length(a)>0){
              area_int[k,o]<-gArea(gIntersection(ell1, ell2))/(gArea(ell1)+gArea(ell2)-gArea(gIntersection(ell1, ell2)))
            }
            if(length(a)==0){
              area_int[k,o]<-0
            }
            
          }
        }
        png(filename="corrplot_PC_2_3.png", units="in",width=2, height=2 , pointsize=5, res=1200)
        par(font.axis = 1,font.lab=1,mar=c(1,1,1,2), font=2)
        corrplot(t(as.matrix((area_int))),tl.col = "black",cl.lim = c(0, 1),method = "color",is.corr=F, col=col1(250),cl.pos = "b",tl.cex = 0.8, cl.ratio = 0.3)
        dev.off()
        
        write.csv(area_int,file="area_intersection_2_3.csv")
        
      }
      
    }
    
    
    lda_a <- lda(samplem ~Comp.1+Comp.2+Comp.3, data = re)
    s <- lda_a
    capture.output(s, file = "lda_traning.txt")
    
    
    lda_val<-as.data.frame(matrix(nrow=nrow(re),ncol=3))
    lda_val[,1]<-calclda(re[,1:3], lda_a$scaling[,1])
    lda_val[,2]<-calclda(re[,1:3], lda_a$scaling[,2])
    lda_val[,3]<-re$samplem
    colorset = cols[as.factor(lda_val$V3)]
    
    lda_val$Label<-rownames(re)
    write.csv(lda_val,file="lda_training.csv")
    
    lda_predict <- predict(lda_a, newdata = re)
    lda_pred_class_train<-as.data.frame(matrix(nrow=nrow(re),ncol=2))
    lda_pred_class_train[,1]<- lda_predict$class
    lda_pred_class_train[,2]<- re$samplem
    colnames(lda_pred_class_train)<-c("predicted","actual")
    x<-table(lda_pred_class_train)  
    col1 <- colorRampPalette(c("#7F0000", "#7F0000","#7F0000","#7F0000","#7F0000","#7F0000","#7F0000","#7F0000","#7F0000"
                               ,"red", "#FF7F00", "yellow", "#7FFF7F",
                               "cyan", "#007FFF", "blue", "#00007F"))
    
    png(filename="LDA_training_pred.png", units="in",width=2, height=2 , pointsize=5, res=1200)
    par(font.axis = 1,font.lab=1,mar=c(1,1,1,2), font=2)
    corrplot(x,tl.col = "black",method = "color",is.corr=F, col=col1(250),cl.pos = "b",tl.cex = 0.8, cl.ratio = 0.3)
    dev.off()
    write.csv(x,file="confusion_matrix_lda_training.csv")
    
    png(filename="lda_biplot_1_2_ellipses.png", units="in",width=2, height=2 , pointsize=7, res=1200)
    par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2)
    ellipseplot(lda_val[,1],                          # data for x-axis
                lda_val[,2],                          # data for y-axis
                as.factor(lda_val[,3]),                    # factor with classes
                elev=0.95,
                pcol=cols,    # colors for plotting (must match # of factors)
                pbgcol=TRUE,                           # point borders black?
                cexsize=0.1,                            # size of points 
                ppch=rep(21,10),                          # shape of points (must match # of factors)
                axissize=0.8,                           # Set axis text size
                linewidth=1,                          # Set axis line size
                font=1
    )
    title(xlab="LD1",    # % variance explained on PC1
          ylab="LD2",    # % variance explained on PC2 
          main=" ",                  # Title
          cex.lab=1,                    # size of label text
          cex.main=1                    # size of title text
    )
    dev.off()
    
    {
      {
        edf <- data.frame(LV1 = lda_val[,1], LV2=lda_val[,2], factr = lda_val[,3]) # create data frame with data and factor
        ellipses <- dlply(edf, .(factr), function(x) {
          LV1 <- x[,1]
          LV2 <- x[,2]
          dataEllipse(LV1, LV2, levels=0.95, robust=TRUE, draw=FALSE) # Get confidence ellipse points from dataEllipse() function by factor level
        })
        col1 <- colorRampPalette(c("#7F0000", "#7F0000","#7F0000","#7F0000","#7F0000","#7F0000","#7F0000","#7F0000","#7F0000"
                                   ,"red", "#FF7F00", "yellow", "#7FFF7F",
                                   "cyan", "#007FFF", "blue", "#00007F"))
        
        area_int<-as.data.frame(matrix(NA, nrow=10,ncol=10))
        colnames(area_int)<-rev(attributes(ellipses)$names)
        rownames(area_int)<-(attributes(ellipses)$names)
        for(k in 1:10){
          for(o in 10:1){
            ell1 <- xy2SP(ellipses[[k]])
            ell2 <- xy2SP(ellipses[[11-o]])
            a<-gIntersection(ell1, ell2)
            if(length(a)>0){
              area_int[k,o]<-gArea(gIntersection(ell1, ell2))/(gArea(ell1)+gArea(ell2)-gArea(gIntersection(ell1, ell2)))
            }
            if(length(a)==0){
              area_int[k,o]<-0
            }
            
          }
        }
        png(filename="corrplot_LD1_1_2.png", units="in",width=2, height=2 , pointsize=5, res=1200)
        par(font.axis = 1,font.lab=1,mar=c(1,1,1,2), font=2)
        corrplot(t(as.matrix((area_int))),tl.col = "black",cl.lim = c(0, 1),method = "color",is.corr=F, col=col1(250),cl.pos = "b",tl.cex = 0.8, cl.ratio = 0.3)
        dev.off()
        
        write.csv(area_int,file="area_intersection_1_2.csv")
        
      }
      
    }
  }
  #Linear regression
  {
    dataset_nn<-merge(lda_val,oct4_int_selec,by="Label")
    lmMod <- lm(Median ~V1,data=dataset_nn)  # build the model
    
    distPred <- predict(lmMod, dataset_nn)  # predict 
    actuals_preds <- data.frame(cbind(actuals=dataset_nn$Mean, predicteds=distPred))  # make actuals_predicteds dataframe.
    correlation_accuracy <- cor(actuals_preds) 
    min_max_accuracy <- mean(apply(actuals_preds, 1, min) / apply(actuals_preds, 1, max))  
    mape <- mean(abs((actuals_preds$predicteds - actuals_preds$actuals))/actuals_preds$actuals)  
    
    s <- summary(lmMod)
    capture.output(s, file = "linear_regression.txt")
    
    png(filename="LD1_vs_mean_intesity.png", units="in",width=2, height=2 , pointsize=7, res=1200)
    par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2)
    plot(dataset_nn$Median~dataset_nn$V1,pch=18,las=1,xlab="LD1",ylab="Oct4 Nuclear Intensity",cex.axis=0.6)
    abline(lmMod,lwd=2,col="red")
    usra <- par( "usr" )
    text(usra[ 2 ]*0.45, usra[ 4 ]*0.95, paste("adj.R.sq=",round(s$adj.r.squared,2)),cex=0.6) 
    dev.off()
    
    png(filename="actual_vs_predicted.png", units="in",width=2, height=2 , pointsize=7, res=1200)
    par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2)
    plot(actuals_preds$actuals~actuals_preds$predicteds,pch=18, las=1, ylab="actual Oct4 levels",xlab="predicted Oct4 levels",cex.axis=0.6)
    abline(0,1,col="red",lwd=2)    
    usra <- par( "usr" )
    text(usra[ 1 ]*2, usra[ 4]*0.95, paste("Cor=",round(correlation_accuracy[1,2],2)),cex=0.6) 
    text(usra[ 1 ]*2, usra[ 4]*0.85, paste("SE=",round(mape,2)),cex=0.6) 
    dev.off()
    write.csv(actuals_preds,file="linear_regression_actual_prediction.csv")
    
    
  }
  
}

#nuc feature changes over oct4
{
  merged_dataset<-merge(lda_val,nucleus_combined,by="Label")
  cor_with_ld1=vector()
  for(i in 1:length(parameters_geo_texture_data)){
    cor_with_ld1[i]<-cor(merged_dataset$V1,merged_dataset[,which(colnames(merged_dataset)==parameters_geo_texture_data[i])], 
                         method="spearman")
    
  }
  names(cor_with_ld1)<-parameters_geo_texture_data
  
  sub_features<-nucleus_combined[,which(colnames(nucleus_combined) %in% parameters_geo_texture_data)]
  rownames(sub_features)<-nucleus_combined$Label
  d1<-apply(sub_features[which(rownames(sub_features) %in% oct4_labels_0_10),],2,mean,na.rm=T)
  d2<-apply(sub_features[which(rownames(sub_features) %in% oct4_labels_10_20),],2,mean,na.rm=T)
  d3<-apply(sub_features[which(rownames(sub_features) %in% oct4_labels_20_30),],2,mean,na.rm=T)
  d4<-apply(sub_features[which(rownames(sub_features) %in% oct4_labels_30_40),],2,mean,na.rm=T)
  d5<-apply(sub_features[which(rownames(sub_features) %in% oct4_labels_40_50),],2,mean,na.rm=T)
  d6<-apply(sub_features[which(rownames(sub_features) %in% oct4_labels_50_60),],2,mean,na.rm=T)
  d7<-apply(sub_features[which(rownames(sub_features) %in% oct4_labels_60_70),],2,mean,na.rm=T)
  d8<-apply(sub_features[which(rownames(sub_features) %in% oct4_labels_70_80),],2,mean,na.rm=T)
  d9<-apply(sub_features[which(rownames(sub_features) %in% oct4_labels_80_90),],2,mean,na.rm=T)
  d10<-apply(sub_features[which(rownames(sub_features) %in% oct4_labels_90_100),],2,mean,na.rm=T)
  
  sub_features<-as.data.frame(cbind(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10))
  sub_features$Label<-rownames(sub_features)
  cor_with_ld1<-as.data.frame(cor_with_ld1)
  cor_with_ld1$Label<-rownames(cor_with_ld1)
  sub_features_merged<-merge(sub_features,cor_with_ld1,by="Label")
  
  library(plot3D)
  rbPal <- colorRampPalette(c('papayawhip','orchid'))
  sub_features_merged[,2:11]<-(sub_features_merged[,2:11]-apply(sub_features_merged[,2:11],1,mean) )/apply(sub_features_merged[,2:11],1,sd)
  
  cor_cols<-redblue(50)[as.numeric(cut(sub_features_merged$cor_with_ld1,breaks = 50))]
  cor_cols<-cor_cols[order(sub_features_merged[,12])]
  
  png(filename="nuc_mean_summary_key.png", units="in", width=2, height=2 , pointsize=5, res=1200)
  par(font.axis = 2,font.lab=2,mar=c(4,4,2,2), font=2)
  colkey(col=redblue(50),range(sub_features_merged$cor_with_ld1),side = 1, add = F,length = 0.6, width = 1,clab = "Correlation with\n LD1",side.clab=1)
  dev.off()
  
  
  sub_features_merged<-sub_features_merged[order(sub_features_merged$cor_with_ld1),]
  cor_cols<-redblue(100)[as.numeric(cut(sub_features_merged$cor_with_ld1,breaks = 100))]
  
  png(filename="nuc_mean_summary.png", units="in", width=4, height=8.5 , pointsize=5, res=1200)
  par(font.axis = 2,font.lab=2,mar=c(4,4,2,2), font=2)
  heatmap.2(as.matrix(sub_features_merged[,2:11]), dendrogram ="none",Rowv = NA,Colv = NA,trace="none",density.info="none",srtCol = 0,cexCol = 1.5,
            key.title = "",cexRow = 1,col=rbPal(50),scale = "none",RowSideColors = cor_cols,labRow =sub_features_merged$Label,margins = c(2, 20) )
  dev.off()
  
  cor_cols<-cor_cols[rev(order(abs(sub_features_merged$cor_with_ld1)))][1:30]
  sub_features_merged<-sub_features_merged[rev(order(abs(sub_features_merged$cor_with_ld1))),]
  
  png(filename="nuc_mean_summary_top20.png", units="in", width=2, height=2 , pointsize=5, res=1200)
  par(font.axis = 2,font.lab=2,mar=c(4,4,2,5), font=2)
  heatmap.2(as.matrix(sub_features_merged[1:30,2:11]), dendrogram ="none",Rowv = NA,Colv = NA,trace="none",density.info="none",srtCol = 0,cexCol = 1.5,
            key.title = "",cexRow = 1,col=rbPal(50),scale = "none",RowSideColors = cor_cols,labRow =sub_features_merged$Label,margins = c(3, 15))
  dev.off()
  
  
  rm(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10)
  
}




