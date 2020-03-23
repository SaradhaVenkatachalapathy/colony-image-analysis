normalize_2d_mean_trials <-function(data){
  
  if(length(levels(as.factor(data$trial))==3)){
    T1<-subset(data,data$trial=="T1")
    T1$Mean<-T1$Mean/max(T1$Mean)
    T1$IntDen<-T1$IntDen/max(T1$IntDen)
    T1$Mode<-T1$Mode/max(T1$Mode)
    T1$Min<-T1$Min/max(T1$Min)
    T1$Max<-T1$Max/max(T1$Max)
    T1$Median<-T1$Max/max(T1$Median)
    T1$RawIntDen<-T1$Max/max(T1$RawIntDen)
    
    T2<-subset(data,data$trial=="T2")
    T2$Mean<-T2$Mean/max(T2$Mean)
    T2$IntDen<-T2$IntDen/max(T2$IntDen)
    T2$Mode<-T2$Mode/max(T2$Mode)
    T2$Min<-T2$Min/max(T2$Min)
    T2$Max<-T2$Max/max(T2$Max)
    T2$Median<-T2$Max/max(T2$Median)
    T2$RawIntDen<-T2$Max/max(T2$RawIntDen)
    
    T3<-subset(data,data$trial=="T3")
    T3$Mean<-T3$Mean/max(T3$Mean)
    T3$IntDen<-T3$IntDen/max(T3$IntDen)
    T3$Mode<-T3$Mode/max(T3$Mode)
    T3$Min<-T3$Min/max(T3$Min)
    T3$Max<-T3$Max/max(T3$Max)
    T3$Median<-T3$Max/max(T3$Median)
    T3$RawIntDen<-T3$Max/max(T3$RawIntDen)
    
    return(rbind(T1,T2,T3))
  }
  
  else if(length(levels(as.factor(data$trial))==4)){
    T1<-subset(data,data$trial=="T1")
    T1$Mean<-T1$Mean/max(T1$Mean)
    T1$IntDen<-T1$IntDen/max(T1$IntDen)
    T1$Mode<-T1$Mode/max(T1$Mode)
    T1$Min<-T1$Min/max(T1$Min)
    T1$Max<-T1$Max/max(T1$Max)
    T1$Median<-T1$Max/max(T1$Median)
    T1$RawIntDen<-T1$Max/max(T1$RawIntDen)
    
    T2<-subset(data,data$trial=="T2")
    T2$Mean<-T2$Mean/max(T2$Mean)
    T2$IntDen<-T2$IntDen/max(T2$IntDen)
    T2$Mode<-T2$Mode/max(T2$Mode)
    T2$Min<-T2$Min/max(T2$Min)
    T2$Max<-T2$Max/max(T2$Max)
    T2$Median<-T2$Max/max(T2$Median)
    T2$RawIntDen<-T2$Max/max(T2$RawIntDen)
    
    T3<-subset(data,data$trial=="T3")
    T3$Mean<-T3$Mean/max(T3$Mean)
    T3$IntDen<-T3$IntDen/max(T3$IntDen)
    T3$Mode<-T3$Mode/max(T3$Mode)
    T3$Min<-T3$Min/max(T3$Min)
    T3$Max<-T3$Max/max(T3$Max)
    T3$Median<-T3$Max/max(T3$Median)
    T3$RawIntDen<-T3$Max/max(T3$RawIntDen)
    
    T4<-subset(data,data$trial=="T4")
    T4$Mean<-T4$Mean/max(T4$Mean)
    T4$IntDen<-T4$IntDen/max(T4$IntDen)
    T4$Mode<-T4$Mode/max(T4$Mode)
    T4$Min<-T4$Min/max(T4$Min)
    T4$Max<-T4$Max/max(T4$Max)
    T4$Median<-T4$Max/max(T4$Median)
    T4$RawIntDen<-T4$Max/max(T4$RawIntDen)
    
    
    return(rbind(T1,T2,T3,T4))
  }
}
normalize_mean_trials <-function(data){
  
  if(length(levels(as.factor(data$trial))==3)){
    T1<-subset(data,data$trial=="T1")
    T1$Mean<-T1$Mean/max(T1$Mean)
    T1$IntDen<-T1$IntDen/max(T1$IntDen)
    T1$Mode<-T1$Mode/max(T1$Mode)
    T1$Min<-T1$Min/max(T1$Min)
    T1$Max<-T1$Max/max(T1$Max)

    T2<-subset(data,data$trial=="T2")
    T2$Mean<-T2$Mean/max(T2$Mean)
    T2$IntDen<-T2$IntDen/max(T2$IntDen)
    T2$Mode<-T2$Mode/max(T2$Mode)
    T2$Min<-T2$Min/max(T2$Min)
    T2$Max<-T2$Max/max(T2$Max)

    T3<-subset(data,data$trial=="T3")
    T3$Mean<-T3$Mean/max(T3$Mean)
    T3$IntDen<-T3$IntDen/max(T3$IntDen)
    T3$Mode<-T3$Mode/max(T3$Mode)
    T3$Min<-T3$Min/max(T3$Min)
    T3$Max<-T3$Max/max(T3$Max)

    return(rbind(T1,T2,T3))
  }
  
  else if(length(levels(as.factor(data$trial))==4)){
    T1<-subset(data,data$trial=="T1")
    T1$Mean<-T1$Mean/max(T1$Mean)
    T1$IntDen<-T1$IntDen/max(T1$IntDen)
    T1$Mode<-T1$Mode/max(T1$Mode)
    T1$Min<-T1$Min/max(T1$Min)
    T1$Max<-T1$Max/max(T1$Max)

    T2<-subset(data,data$trial=="T2")
    T2$Mean<-T2$Mean/max(T2$Mean)
    T2$IntDen<-T2$IntDen/max(T2$IntDen)
    T2$Mode<-T2$Mode/max(T2$Mode)
    T2$Min<-T2$Min/max(T2$Min)
    T2$Max<-T2$Max/max(T2$Max)

    T3<-subset(data,data$trial=="T3")
    T3$Mean<-T3$Mean/max(T3$Mean)
    T3$IntDen<-T3$IntDen/max(T3$IntDen)
    T3$Mode<-T3$Mode/max(T3$Mode)
    T3$Min<-T3$Min/max(T3$Min)
    T3$Max<-T3$Max/max(T3$Max)

    T4<-subset(data,data$trial=="T4")
    T4$Mean<-T4$Mean/max(T4$Mean)
    T4$IntDen<-T4$IntDen/max(T4$IntDen)
    T4$Mode<-T4$Mode/max(T4$Mode)
    T4$Min<-T4$Min/max(T4$Min)
    T4$Max<-T4$Max/max(T4$Max)

    
    return(rbind(T1,T2,T3,T4))
  }
}

# note remove spheroid data for exp without actin
#read spheroid combined data from 5 mega sets 
{
  setwd("~/Desktop/Reprogramming/actin_pmlc_oct4/combined_data/")
  
  spheroid_combined_set1 <- read.csv("spheroid_combined.csv",stringsAsFactors = F)
  spheroid_combined_set1<-spheroid_combined_set1[!(duplicated(spheroid_combined_set1$Label)),]
  spheroid_combined_set1$meg_trial<-paste(spheroid_combined_set1$trial,"set1",sep="_")
  
  setwd("~/Desktop/Reprogramming/vimentin_ecad_oct4/combined_data/")
  
  spheroid_combined_set2 <- read.csv("spheroid_combined.csv",stringsAsFactors = F)
  spheroid_combined_set2<-spheroid_combined_set2[!(duplicated(spheroid_combined_set2$Label)),]
  spheroid_combined_set2$meg_trial<-paste(spheroid_combined_set2$trial,"set2",sep="_")
  
  setwd("~/Desktop/Reprogramming/ki67/combined_data/")
  
  spheroid_combined_set3 <- read.csv("spheroid_combined.csv",stringsAsFactors = F)
  spheroid_combined_set3<-spheroid_combined_set3[!(duplicated(spheroid_combined_set3$Label)),]
  spheroid_combined_set3<-spheroid_combined_set3[,which(colnames(spheroid_combined_set3) %in% colnames(spheroid_combined_set2))]
  spheroid_combined_set3$meg_trial<-paste(spheroid_combined_set3$trial,"set3",sep="_")
  
  
  setwd("~/Desktop/Reprogramming/H3k9ac/combined_data/")
  
  spheroid_combined_set4 <- read.csv("spheroid_combined.csv",stringsAsFactors = F)
  spheroid_combined_set4<-spheroid_combined_set4[!(duplicated(spheroid_combined_set4$Label)),]
  spheroid_combined_set4$meg_trial<-paste(spheroid_combined_set4$trial,"set4",sep="_")
  
  setwd("~/Desktop/Reprogramming/actin_laminac_oct4/combined_data/")
  
  spheroid_combined_set5 <- read.csv("spheroid_combined.csv",stringsAsFactors = F)
  spheroid_combined_set5<-spheroid_combined_set5[!(duplicated(spheroid_combined_set5$Label)),]
  spheroid_combined_set5$meg_trial<-paste(spheroid_combined_set5$trial,"set5",sep="_")
  
  spheroid_combined<-rbind(spheroid_combined_set1,spheroid_combined_set2,spheroid_combined_set3,spheroid_combined_set4,spheroid_combined_set5)
  spheroid_combined<-spheroid_combined[!(duplicated(spheroid_combined$Label)),]
  rownames(spheroid_combined)<-spheroid_combined$Label
  
  rm(spheroid_combined_set1,spheroid_combined_set2,spheroid_combined_set3,spheroid_combined_set4,spheroid_combined_set5)
  
}
#spheroid_2d_data
{
  setwd("~/Desktop/Reprogramming/actin_pmlc_oct4/combined_data/")
  
  actin <- read.csv("spheroid_2d_int_actin.csv",stringsAsFactors = F)
  actin$Label<-substring(actin$Label, 5, nchar(actin$Label)-4)
  actin<-actin[!(duplicated(actin$Label)),]
  actin<-normalize_2d_mean_trials(actin)
  
  oct4_set1 <- read.csv("spheroid_2d_int_oct4.csv",stringsAsFactors = F)
  oct4_set1$Label<-substring(oct4_set1$Label, 5, nchar(oct4_set1$Label)-4)
  oct4_set1<-oct4_set1[!(duplicated(oct4_set1$Label)),]
  oct4_set1<-normalize_2d_mean_trials(oct4_set1)
  
  pmlc <- read.csv("spheroid_2d_int_pmlc.csv",stringsAsFactors = F)
  pmlc$Label<-substring(pmlc$Label, 5, nchar(pmlc$Label)-4)
  pmlc<-pmlc[!(duplicated(pmlc$Label)),]
  pmlc<-normalize_2d_mean_trials(pmlc)
  
  
  
  setwd("~/Desktop/Reprogramming/ki67/combined_data/")
  ki67 <- read.csv("spheroid_2d_int_KI67.csv",stringsAsFactors = F)
  ki67$Label<-substring(ki67$Label, 5, nchar(ki67$Label)-4)
  ki67<-ki67[!(duplicated(ki67$Label)),]
  ki67<-normalize_2d_mean_trials(ki67)
  
  setwd("~/Desktop/Reprogramming/H3k9ac/combined_data/")
  h3k9ac <- read.csv("spheroid_2d_int_H3K9AC.csv",stringsAsFactors = F)
  h3k9ac$Label<-substring(h3k9ac$Label, 5, nchar(h3k9ac$Label)-4)
  h3k9ac<-h3k9ac[!(duplicated(h3k9ac$Label)),]
  h3k9ac<-normalize_2d_mean_trials(h3k9ac)
  
  setwd("~/Desktop/Reprogramming/actin_laminac_oct4/combined_data/")
  laminac <- read.csv("spheroid_2d_int_LAMINAC.csv",stringsAsFactors = F)
  laminac$Label<-substring(laminac$Label, 5, nchar(laminac$Label)-4)
  laminac<-laminac[!(duplicated(laminac$Label)),]
  laminac<-normalize_2d_mean_trials(laminac)
  
  oct4_set2 <- read.csv("spheroid_2d_int_oct4.csv",stringsAsFactors = F)
  oct4_set2$Label<-substring(oct4_set2$Label, 5, nchar(oct4_set2$Label)-4)
  oct4_set2<-oct4_set2[!(duplicated(oct4_set2$Label)),]
  oct4_set2<-normalize_2d_mean_trials(oct4_set2)
  
  setwd("~/Desktop/Reprogramming/vimentin_ecad_oct4/combined_data/")
  
  vimentin <- read.csv("spheroid_2d_int_VIMENTIN.csv",stringsAsFactors = F)
  vimentin$Label<-substring(vimentin$Label, 5, nchar(vimentin$Label)-4)
  vimentin<-vimentin[!(duplicated(vimentin$Label)),]
  vimentin<-normalize_2d_mean_trials(vimentin)
  
  ecad <- read.csv("spheroid_2d_int_ECAD.csv",stringsAsFactors = F)
  ecad$Label<-substring(ecad$Label, 5, nchar(ecad$Label)-4)
  ecad<-ecad[!(duplicated(ecad$Label)),]
  ecad<-normalize_2d_mean_trials(ecad)
  
  oct4_set3 <- read.csv("spheroid_2d_int_oct4.csv",stringsAsFactors = F)
  oct4_set3$Label<-substring(oct4_set3$Label, 5, nchar(oct4_set3$Label)-4)
  oct4_set3<-oct4_set3[!(duplicated(oct4_set3$Label)),]
  oct4_set3<-normalize_2d_mean_trials(oct4_set3)
  
  oct4<-rbind(oct4_set1,oct4_set2)
  
  spheroid_2dint_data<-list(actin=actin,oct4=oct4,pmlc=pmlc,
                            vimentin=vimentin,ecad=ecad,
                            ki67=ki67,h3k9ac=h3k9ac,
                            laminac=laminac)
  rm(actin,oct4,pmlc,vimentin,ecad,ki67,h3k9ac,laminac,oct4_set1,oct4_set2,oct4_set3)
  
}
#spheroid3d
{
  setwd("~/Desktop/Reprogramming/actin_pmlc_oct4/combined_data/")
  
  actin <- read.csv("spheroid_int_actin.csv",stringsAsFactors = F)
  actin$Label<-substring(actin$Label, 1, nchar(actin$Label)-1)
  actin<-actin[!(duplicated(actin$Label)),]
  actin<-normalize_mean_trials(actin)
  
  oct4_set1 <- read.csv("spheroid_int_oct4.csv",stringsAsFactors = F)
  oct4_set1$Label<-substring(oct4_set1$Label, 1, nchar(oct4_set1$Label)-1)
  oct4_set1<-oct4_set1[!(duplicated(oct4_set1$Label)),]
  oct4_set1<-normalize_mean_trials(oct4_set1)
  
  pmlc <- read.csv("spheroid_int_pmlc.csv",stringsAsFactors = F)
  pmlc$Label<-substring(pmlc$Label, 1, nchar(pmlc$Label)-1)
  pmlc<-pmlc[!(duplicated(pmlc$Label)),]
  pmlc<-normalize_mean_trials(pmlc)
  
  
  setwd("~/Desktop/Reprogramming/vimentin_ecad_oct4/combined_data/")
  
  vimentin <- read.csv("spheroid_int_VIMENTIN.csv",stringsAsFactors = F)
  vimentin$Label<-substring(vimentin$Label, 1, nchar(vimentin$Label)-1)
  vimentin<-vimentin[!(duplicated(vimentin$Label)),]
  vimentin<-normalize_mean_trials(vimentin)
  
  ecad <- read.csv("spheroid_int_ECAD.csv",stringsAsFactors = F)
  ecad$Label<-substring(ecad$Label, 1, nchar(ecad$Label)-1)
  ecad<-ecad[!(duplicated(ecad$Label)),]
  ecad<-normalize_mean_trials(ecad)
  
  setwd("~/Desktop/Reprogramming/ki67/combined_data/")
  ki67 <- read.csv("spheroid_int_KI67.csv",stringsAsFactors = F)
  ki67<-ki67[!(duplicated(ki67$Label)),]
  ki67<-normalize_mean_trials(ki67)
  
  setwd("~/Desktop/Reprogramming/H3k9ac/combined_data/")
  h3k9ac <- read.csv("spheroid_int_H3K9AC.csv",stringsAsFactors = F)
  h3k9ac$Label<-substring(h3k9ac$Label, 1, nchar(h3k9ac$Label)-1)
  h3k9ac<-h3k9ac[!(duplicated(h3k9ac$Label)),]
  h3k9ac<-normalize_mean_trials(h3k9ac)
  
  setwd("~/Desktop/Reprogramming/actin_laminac_oct4/combined_data/")
  laminac <- read.csv("spheroid_int_LAMINAC.csv",stringsAsFactors = F)
  laminac$Label<-substring(laminac$Label, 1, nchar(laminac$Label)-1)
  laminac<-laminac[!(duplicated(laminac$Label)),]
  laminac<-normalize_mean_trials(laminac)
  
  oct4_set2 <- read.csv("spheroid_int_oct4.csv",stringsAsFactors = F)
  oct4_set2$Label<-substring(oct4_set2$Label, 1, nchar(oct4_set2$Label)-1)
  oct4_set2<-oct4_set2[!(duplicated(oct4_set2$Label)),]
  oct4_set2<-normalize_mean_trials(oct4_set2)
  
  oct4<-rbind(oct4_set1,oct4_set2)
  
  spheroid_3dint_data<-list(actin=actin,oct4=oct4,pmlc=pmlc,
                            vimentin=vimentin,ecad=ecad,
                            ki67=ki67,h3k9ac=h3k9ac,
                            laminac=laminac)
  rm(actin,oct4,pmlc,vimentin,ecad,ki67,h3k9ac,laminac,oct4_set1,oct4_set2)
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
                                  "Area","Perim.","Major","Minor","Circ.","Feret","MinFeret","AR","Round","Solidity",
                                  
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
  label_T1_set1<-subset(spheroid_combined,spheroid_combined$meg_trial=="T1_set1")$Label
  label_T2_set1<-subset(spheroid_combined,spheroid_combined$meg_trial=="T2_set1")$Label
  label_T3_set1<-subset(spheroid_combined,spheroid_combined$meg_trial=="T3_set1")$Label
  label_T1_set2<-subset(spheroid_combined,spheroid_combined$meg_trial=="T1_set2")$Label
  label_T2_set2<-subset(spheroid_combined,spheroid_combined$meg_trial=="T2_set2")$Label
  label_T3_set2<-subset(spheroid_combined,spheroid_combined$meg_trial=="T3_set2")$Label
  label_T1_set3<-subset(spheroid_combined,spheroid_combined$meg_trial=="T1_set3")$Label
  label_T2_set3<-subset(spheroid_combined,spheroid_combined$meg_trial=="T2_set3")$Label
  label_T3_set3<-subset(spheroid_combined,spheroid_combined$meg_trial=="T3_set3")$Label
  label_T1_set4<-subset(spheroid_combined,spheroid_combined$meg_trial=="T1_set4")$Label
  label_T2_set4<-subset(spheroid_combined,spheroid_combined$meg_trial=="T2_set4")$Label
  label_T3_set4<-subset(spheroid_combined,spheroid_combined$meg_trial=="T3_set4")$Label
  label_T1_set5<-subset(spheroid_combined,spheroid_combined$meg_trial=="T1_set5")$Label
  label_T2_set5<-subset(spheroid_combined,spheroid_combined$meg_trial=="T2_set5")$Label
  label_T3_set5<-subset(spheroid_combined,spheroid_combined$meg_trial=="T3_set5")$Label
  label_T4_set5<-subset(spheroid_combined,spheroid_combined$meg_trial=="T4_set5")$Label
  
  
  D0_lab<-subset(spheroid_combined,spheroid_combined$sample=="D0")$Label
  D2_lab<-subset(spheroid_combined,spheroid_combined$sample=="D2")$Label
  D4_lab<-subset(spheroid_combined,spheroid_combined$sample=="D4")$Label
  D6_lab<-subset(spheroid_combined,spheroid_combined$sample=="D6")$Label
  D8_lab<-subset(spheroid_combined,spheroid_combined$sample=="D8")$Label
  
}
cols_timeline<-colorRampPalette(c("red","yellow"))( 5) 
cols_trials<-colorRampPalette(c("gray","black")) (length(levels(as.factor(spheroid_combined$meg_trial))))

dir.create("/Users/saradhavenkatachalapathy/Desktop/reprogramming/timeline_spheroid_analysis_scaled/")
setwd("/Users/saradhavenkatachalapathy/Desktop/reprogramming/timeline_spheroid_analysis_scaled/")

#pca
{
  spheroid_subset<-spheroid_combined[,which(colnames(spheroid_combined) %in% parameters_geo_texture_data)]
  spheroid_subset <- spheroid_subset[ - as.numeric(which(apply(spheroid_subset, 2, var) == 0))]
  
  scaled_spheroid_subset<-subset(spheroid_subset,rownames(spheroid_subset) %in% label_T1_set1)
  scaled_spheroid_subset<-scale(scaled_spheroid_subset,center = T, scale = T)
  temp2<-subset(spheroid_subset,rownames(spheroid_subset) %in% label_T2_set1)
  temp2<-scale(temp2,center = T, scale = T)
  scaled_spheroid_subset<-rbind(scaled_spheroid_subset,temp2)
  temp2<-subset(spheroid_subset,rownames(spheroid_subset) %in% label_T3_set1)
  temp2<-scale(temp2,center = T, scale = T)
  scaled_spheroid_subset<-rbind(scaled_spheroid_subset,temp2)
  temp2<-subset(spheroid_subset,rownames(spheroid_subset) %in% label_T1_set2)
  temp2<-scale(temp2,center = T, scale = T)
  scaled_spheroid_subset<-rbind(scaled_spheroid_subset,temp2)
  temp2<-subset(spheroid_subset,rownames(spheroid_subset) %in% label_T2_set2)
  temp2<-scale(temp2,center = T, scale = T)
  scaled_spheroid_subset<-rbind(scaled_spheroid_subset,temp2)
  temp2<-subset(spheroid_subset,rownames(spheroid_subset) %in% label_T3_set2)
  temp2<-scale(temp2,center = T, scale = T)
  scaled_spheroid_subset<-rbind(scaled_spheroid_subset,temp2)
  temp2<-subset(spheroid_subset,rownames(spheroid_subset) %in% label_T1_set3)
  temp2<-scale(temp2,center = T, scale = T)
  scaled_spheroid_subset<-rbind(scaled_spheroid_subset,temp2)
  temp2<-subset(spheroid_subset,rownames(spheroid_subset) %in% label_T2_set3)
  temp2<-scale(temp2,center = T, scale = T)
  scaled_spheroid_subset<-rbind(scaled_spheroid_subset,temp2)
  temp2<-subset(spheroid_subset,rownames(spheroid_subset) %in% label_T3_set3)
  temp2<-scale(temp2,center = T, scale = T)
  scaled_spheroid_subset<-rbind(scaled_spheroid_subset,temp2)
  temp2<-subset(spheroid_subset,rownames(spheroid_subset) %in% label_T1_set4)
  temp2<-scale(temp2,center = T, scale = T)
  scaled_spheroid_subset<-rbind(scaled_spheroid_subset,temp2)
  temp2<-subset(spheroid_subset,rownames(spheroid_subset) %in% label_T2_set4)
  temp2<-scale(temp2,center = T, scale = T)
  scaled_spheroid_subset<-rbind(scaled_spheroid_subset,temp2)
  temp2<-subset(spheroid_subset,rownames(spheroid_subset) %in% label_T3_set4)
  temp2<-scale(temp2,center = T, scale = T)
  scaled_spheroid_subset<-rbind(scaled_spheroid_subset,temp2)
  temp2<-subset(spheroid_subset,rownames(spheroid_subset) %in% label_T1_set5)
  temp2<-scale(temp2,center = T, scale = T)
  scaled_spheroid_subset<-rbind(scaled_spheroid_subset,temp2)
  temp2<-subset(spheroid_subset,rownames(spheroid_subset) %in% label_T2_set5)
  temp2<-scale(temp2,center = T, scale = T)
  scaled_spheroid_subset<-rbind(scaled_spheroid_subset,temp2)
  temp2<-subset(spheroid_subset,rownames(spheroid_subset) %in% label_T3_set5)
  temp2<-scale(temp2,center = T, scale = T)
  scaled_spheroid_subset<-rbind(scaled_spheroid_subset,temp2)
  temp2<-subset(spheroid_subset,rownames(spheroid_subset) %in% label_T4_set5)
  temp2<-scale(temp2,center = T, scale = T)
  scaled_spheroid_subset<-rbind(scaled_spheroid_subset,temp2)
  
  spheroid_subset<-as.matrix(na.omit((scaled_spheroid_subset)))
  
  spheroid_pca <- princomp(spheroid_subset,  cor = TRUE, scores = TRUE)
  pca_scores<-spheroid_pca$scores
  
  png(filename="pca_variance_vs_PCs.png", units="in", width=2, height=2 , pointsize=7, res=1200)
  par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2)
  plot(spheroid_pca, type="l", las=1, main="", pch=19, las=1)
  box()
  dev.off()
  
  s<-abs(with(spheroid_pca, unclass(loadings)))
  
  png(filename="pc1_loadings.png", units="in", width=2, height=4 , pointsize=7, res=1200)
  so<-s[order(s[,1]),1]
  par(font.axis = 2,font.lab=2,mar=c(5,9,1,1), font=2)
  barplot(so, horiz=T, las=1, xlab="loading on pC1")
  box()
  dev.off()
  
  png(filename="pc2_loadings.png", units="in", width=2, height=4 , pointsize=7, res=1200)
  so<-s[order(s[,2]),2]
  par(font.axis = 2,font.lab=2,mar=c(5,9,1,1), font=2)
  barplot(so, horiz=T, las=1, xlab="loading on pC2")
  box()
  dev.off()
  
  a<-subset(pca_scores,rownames(pca_scores)%in% D0_lab)
  b<-subset(pca_scores,rownames(pca_scores)%in% D2_lab)
  c<-subset(pca_scores,rownames(pca_scores)%in% D4_lab)
  d<-subset(pca_scores,rownames(pca_scores)%in% D6_lab)
  e<-subset(pca_scores,rownames(pca_scores)%in% D8_lab)
  
  png(filename="pca_samples.png", units="in", width=2, height=2 , pointsize=7, res=1200)
  par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2)
  plot(pca_scores[,1],pca_scores[,2],col=NA,xlab="PC1",ylab="PC2")
  points(a[,1],a[,2],col=cols_timeline[1],pch=19,cex=0.3)
  points(b[,1],b[,2],col=cols_timeline[2],pch=19,cex=0.3)
  points(c[,1],c[,2],col=cols_timeline[3],pch=19,cex=0.3)
  points(d[,1],d[,2],col=cols_timeline[4],pch=19,cex=0.3)
  points(e[,1],e[,2],col=cols_timeline[5],pch=19,cex=0.3)
  dev.off()
  
  png(filename="pca_samplesendpoints.png", units="in", width=2, height=2 , pointsize=7, res=1200)
  par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2)
  plot(pca_scores[,1],pca_scores[,2],col=NA,xlab="PC1",ylab="PC2")
  points(a[,1],a[,2],col=cols_timeline[1],pch=19,cex=0.3)
  points(e[,1],e[,2],col=cols_timeline[5],pch=19,cex=0.3)
  dev.off()
 
  pca_T1_set1<-subset(pca_scores,rownames(pca_scores)%in% label_T1_set1)
  pca_T1_set2<-subset(pca_scores,rownames(pca_scores)%in% label_T1_set2)
  pca_T1_set3<-subset(pca_scores,rownames(pca_scores)%in% label_T1_set3)
  pca_T1_set4<-subset(pca_scores,rownames(pca_scores)%in% label_T1_set4)
  pca_T1_set5<-subset(pca_scores,rownames(pca_scores)%in% label_T1_set5)
  pca_T2_set1<-subset(pca_scores,rownames(pca_scores)%in% label_T2_set1)
  pca_T2_set2<-subset(pca_scores,rownames(pca_scores)%in% label_T2_set2)
  pca_T2_set3<-subset(pca_scores,rownames(pca_scores)%in% label_T2_set3)
  pca_T2_set4<-subset(pca_scores,rownames(pca_scores)%in% label_T2_set4)
  pca_T2_set5<-subset(pca_scores,rownames(pca_scores)%in% label_T2_set5)
  pca_T3_set1<-subset(pca_scores,rownames(pca_scores)%in% label_T3_set1)
  pca_T3_set2<-subset(pca_scores,rownames(pca_scores)%in% label_T3_set2)
  pca_T3_set3<-subset(pca_scores,rownames(pca_scores)%in% label_T3_set3)
  pca_T3_set4<-subset(pca_scores,rownames(pca_scores)%in% label_T3_set4)
  pca_T3_set5<-subset(pca_scores,rownames(pca_scores)%in% label_T3_set5)
  pca_T4_set5<-subset(pca_scores,rownames(pca_scores)%in% label_T4_set5)
  
  png(filename="pca_trials.png", units="in", width=2, height=2 , pointsize=7, res=1200)
  par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2)
  plot(pca_scores[,1],pca_scores[,2],col=NA,xlab="PC1",ylab="PC2")
  points(pca_T1_set1[,1],pca_T1_set1[,2],col=cols_trials[1],pch=19,cex=0.1)
  points(pca_T1_set2[,1],pca_T1_set2[,2],col=cols_trials[2],pch=19,cex=0.1)
  points(pca_T1_set3[,1],pca_T1_set3[,2],col=cols_trials[3],pch=19,cex=0.1)
  points(pca_T1_set4[,1],pca_T1_set4[,2],col=cols_trials[4],pch=19,cex=0.1)
  points(pca_T1_set5[,1],pca_T1_set5[,2],col=cols_trials[5],pch=19,cex=0.1)
  
  points(pca_T2_set1[,1],pca_T2_set1[,2],col=cols_trials[6],pch=19,cex=0.1)
  points(pca_T2_set2[,1],pca_T2_set2[,2],col=cols_trials[7],pch=19,cex=0.1)
  points(pca_T2_set3[,1],pca_T2_set3[,2],col=cols_trials[8],pch=19,cex=0.1)
  points(pca_T2_set4[,1],pca_T2_set4[,2],col=cols_trials[9],pch=19,cex=0.1)
  points(pca_T2_set5[,1],pca_T2_set5[,2],col=cols_trials[10],pch=19,cex=0.1)
  
  points(pca_T3_set1[,1],pca_T3_set1[,2],col=cols_trials[11],pch=19,cex=0.1)
  points(pca_T3_set2[,1],pca_T3_set2[,2],col=cols_trials[12],pch=19,cex=0.1)
  points(pca_T3_set3[,1],pca_T3_set3[,2],col=cols_trials[13],pch=19,cex=0.1)
  points(pca_T3_set4[,1],pca_T3_set4[,2],col=cols_trials[14],pch=19,cex=0.1)
  points(pca_T3_set5[,1],pca_T3_set5[,2],col=cols_trials[15],pch=19,cex=0.1)
  
  points(pca_T4_set5[,1],pca_T4_set5[,2],col=cols_trials[16],pch=19,cex=0.1)
  
  dev.off()
  
  rm(a,b,c,d,e,s,so,temp2,scaled_spheroid_subset,
     pca_T1_set1,pca_T1_set2,pca_T1_set3,pca_T1_set4,pca_T1_set5,
     pca_T2_set1,pca_T2_set2,pca_T2_set3,pca_T2_set4,pca_T2_set5,
     pca_T3_set1,pca_T3_set2,pca_T3_set3,pca_T3_set4,pca_T3_set5,
     pca_T4_set5)
  
  
}

#diffusion map
{
  library(destiny)
  
  png(filename="sigmas.png", units="in", width=2, height=2 , pointsize=7, res=1200)
  sigmas <- find_sigmas(as.matrix(pca_scores[,1:6], verbose = FALSE))
  dev.off()
  
  dm<-DiffusionMap(as.matrix(pca_scores[,1:6]),sigma = 1.0)
  
  #just to select_tips
  dpt <- DPT(dm)
  dmp<-as.data.frame(cbind(dpt$DC1,dpt$DC2))
  rownames(dmp)<-rownames(pca_scores)
  dmp$dpt<-dpt$dpt
  
  plot(dpt)

  a<-dmp[which(rownames(dmp) %in% D0_lab),]
  b<-dmp[which(rownames(dmp) %in% D2_lab),]
  c<-dmp[which(rownames(dmp) %in% D4_lab),]
  d<-dmp[which(rownames(dmp) %in% D6_lab),]
  e<-dmp[which(rownames(dmp) %in% D8_lab),]
  plot(dmp[,1],dmp[,2],col=NA, xlab="Diffusion Component 1", ylab="Diffusion Component 2",las=1)
  points(e[,1],e[,2],col=cols_timeline[5],pch=19,cex=0.5)
  points(d[,1],d[,2],col=cols_timeline[4],pch=19,cex=0.5)
  points(c[,1],c[,2],col=cols_timeline[3],pch=19,cex=0.5)
  points(b[,1],b[,2],col=cols_timeline[2],pch=19,cex=0.5)
  points(a[,1],a[,2],col=cols_timeline[1],pch=19,cex=0.5)
  rm(a,b,c,d,e)
  
  #re do to coorect time
  dpt <- DPT(dm, tips= which(rownames(pca_scores) == "26_spheroid_354_146_0_14"))
  dmp<-as.data.frame(cbind(dpt$DC1,dpt$DC2))
  rownames(dmp)<-rownames(pca_scores)
  dmp$dpt<-dpt$dpt
  
  png(filename="diffusion_map.png", units="in", width=3, height=2 , pointsize=7, res=1200)
  par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2)
  plot(dpt)
  dev.off()
  
  png(filename="diffusion_map_braches.png", units="in", width=3, height=2 , pointsize=6, res=1200)
  par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2)
  plot(dpt, col_by = 'branch')
  dev.off()
  
  a<-dmp[which(rownames(dmp) %in% D0_lab),]
  b<-dmp[which(rownames(dmp) %in% D2_lab),]
  c<-dmp[which(rownames(dmp) %in% D4_lab),]
  d<-dmp[which(rownames(dmp) %in% D6_lab),]
  e<-dmp[which(rownames(dmp) %in% D8_lab),]
  
  png(filename="diffusion_map_samples.png", units="in", width=2, height=2 , pointsize=7, res=1200)
  par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2)
  plot(dmp[,1],dmp[,2],col=NA, xlab="Diffusion Component 1", ylab="Diffusion Component 2",las=1)
  points(e[,1],e[,2],col=cols_timeline[5],pch=19,cex=0.5)
  points(d[,1],d[,2],col=cols_timeline[4],pch=19,cex=0.5)
  points(c[,1],c[,2],col=cols_timeline[3],pch=19,cex=0.5)
  points(b[,1],b[,2],col=cols_timeline[2],pch=19,cex=0.5)
  points(a[,1],a[,2],col=cols_timeline[1],pch=19,cex=0.5)
  dev.off()
  
  diffusion_psuedotime<-dmp
  colnames(diffusion_psuedotime)<-c("DC1","DC2","DPT") 
  spheroid_subset<-as.data.frame(spheroid_subset)
  spheroid_subset$Label=rownames(spheroid_subset)
  diffusion_psuedotime$Label=rownames(diffusion_psuedotime)
  spheroid_diffusion=merge(spheroid_subset,diffusion_psuedotime, by="Label")
  
  rm(a,b,c,d,e,dm,dpt,dmp,pca_scores,sigmas)
  
}

# correlating spheroid features with diffusion analysis
{
  library(gplots)
  parameters_geo_texture_data_spheroid<-parameters_geo_texture_data[which(parameters_geo_texture_data %in% colnames(spheroid_diffusion))]
  
  cor_with_dpt=vector()
  for(i in 1:length(parameters_geo_texture_data_spheroid)){
    cor_with_dpt[i]<-cor(spheroid_diffusion$DPT,spheroid_diffusion[,which(colnames(spheroid_diffusion)==parameters_geo_texture_data_spheroid[i])], 
                         method="spearman")
    
  }
  names(cor_with_dpt)<-parameters_geo_texture_data_spheroid
  
  cor_with_dc2=vector()
  for(i in 1:length(parameters_geo_texture_data_spheroid)){
    cor_with_dc2[i]<-cor(spheroid_diffusion$DC2,spheroid_diffusion[,which(colnames(spheroid_diffusion)==parameters_geo_texture_data_spheroid[i])], 
                         method="spearman")
  }
  names(cor_with_dc2)<-parameters_geo_texture_data_spheroid
  
  cor_with_dc1=vector()
  for(i in 1:length(parameters_geo_texture_data_spheroid)){
    cor_with_dc1[i]<-cor(spheroid_diffusion$DC1,spheroid_diffusion[,which(colnames(spheroid_diffusion)==parameters_geo_texture_data_spheroid[i])], 
                         method="spearman")
  }
  names(cor_with_dc1)<-parameters_geo_texture_data_spheroid
  
  combined_cor=cbind(DC1=cor_with_dc1,DC2=cor_with_dc2,DPT=cor_with_dpt)
  combined_cor<-combined_cor[order(combined_cor[,3]),]
  
  png(filename="correlation_properties_DM.png", units="in", width=2, height=2 , pointsize=5, res=1200)
  par(font.axis = 2,font.lab=2,mar=c(4,4,2,2), font=2)
  heatmap.2(combined_cor, dendrogram ="none",Rowv = NA,Colv = NA,trace="none",density.info="none",srtCol = 0,cexCol = 1,
            key.xlab = "Cor",key.title = "",cexRow = 0.6, col=redblue(50))
  dev.off()
  rm(cor_with_dc1,cor_with_dc2,cor_with_dpt,i,combined_cor)
  
}

#spheroid feature changes over time
{
  rbPal <- colorRampPalette(c('papayawhip','orchid'))
  
  spheroid_diffusion$Label->rownames(spheroid_diffusion)
  spheroid_diffusion_subset<- spheroid_diffusion[, which(colnames(spheroid_diffusion) %in% parameters_geo_texture_data_spheroid)] 
  
  a<-apply(subset(spheroid_diffusion_subset, rownames(spheroid_diffusion_subset) %in% D0_lab),2,mean) 
  b<-apply(subset(spheroid_diffusion_subset, rownames(spheroid_diffusion_subset) %in% D2_lab),2,mean) 
  c<-apply(subset(spheroid_diffusion_subset, rownames(spheroid_diffusion_subset) %in% D4_lab),2,mean) 
  d<-apply(subset(spheroid_diffusion_subset, rownames(spheroid_diffusion_subset) %in% D6_lab),2,mean) 
  e<-apply(subset(spheroid_diffusion_subset, rownames(spheroid_diffusion_subset) %in% D8_lab),2,mean) 
  
  summary_spheroid_feature<-as.data.frame(t(rbind(a,b,c,d,e)))
  colnames(summary_spheroid_feature)<-c("D0","D2","D4","D6","D8")
  
  summary_spheroid_feature<-summary_spheroid_feature[order(summary_spheroid_feature[,5]),]
  png(filename="spheroid_mean_summary.png", units="in", width=2, height=2 , pointsize=5, res=1200)
  par(font.axis = 2,font.lab=2,mar=c(4,4,2,2), font=2)
  heatmap.2(as.matrix(summary_spheroid_feature), dendrogram ="none",Rowv = NA,Colv = NA,trace="none",density.info="none",srtCol = 0,cexCol = 1,
            key.title = "",cexRow = 0.6,col=rbPal(50),scale = "row")
  dev.off()
  rm(spheroid_diffusion_subset,summary_spheroid_feature,a,b,c,d,e)
}

#mapping spheroid features on diffusion maps
{
  rbPal <- colorRampPalette(c('papayawhip','orchid'))
  library(plot3D)
  
 
  png(filename="volume_diffusion.png", units="in", width=2.5, height=2 , pointsize=7, res=1200)
  par(font.axis = 2,font.lab=2,mar=c(4,4,1,5), font=2)
  spheroid_diffusion$Col <- rbPal(50)[as.numeric(cut(spheroid_diffusion$Vol..unit.,breaks = 50))]
  plot(spheroid_diffusion$DC1,spheroid_diffusion$DC2,col=spheroid_diffusion$Col, 
       xlab="Diffusion Component 1", ylab="Diffusion Component 2",las=1,pch=19,cex=0.5)
  colkey(col=rbPal(50),range(spheroid_diffusion$Vol..unit.),side = 4, add = TRUE,length = 0.5, width = 0.5,clab = "Spheroid Volume",side.clab=4)
  dev.off()
  
  png(filename="sphericity_diffusion.png", units="in", width=2.5, height=2 , pointsize=7, res=1200)
  par(font.axis = 2,font.lab=2,mar=c(4,4,1,5), font=2)
  spheroid_diffusion$Col <- rbPal(50)[as.numeric(cut(spheroid_diffusion$Spher..unit.,breaks = 50))]
  plot(spheroid_diffusion$DC1,spheroid_diffusion$DC2,col=spheroid_diffusion$Col, 
       xlab="Diffusion Component 1", ylab="Diffusion Component 2",las=1,pch=19,cex=0.5)
  colkey(col=rbPal(50),range(spheroid_diffusion$Circ.),side = 4, add = TRUE,length = 0.5, width = 0.5,clab = "Spheroid Sphericity",side.clab=4)
  dev.off()
  
  png(filename="Feret_diffusion.png", units="in", width=2.5, height=2 , pointsize=7, res=1200)
  par(font.axis = 2,font.lab=2,mar=c(4,4,1,5), font=2)
  spheroid_diffusion$Col <- rbPal(50)[as.numeric(cut(spheroid_diffusion$Feret..unit.,breaks = 50))]
  plot(spheroid_diffusion$DC1,spheroid_diffusion$DC2,col=spheroid_diffusion$Col, 
       xlab="Diffusion Component 1", ylab="Diffusion Component 2",las=1,pch=19,cex=0.5)
  colkey(col=rbPal(50),range(spheroid_diffusion$Circ.),side = 4, add = TRUE,length = 0.5, width = 0.5,clab = "Spheroid Feret (3D)",side.clab=4)
  dev.off()
  
  png(filename="DC_min_diffusion.png", units="in", width=2.5, height=2 , pointsize=7, res=1200)
  par(font.axis = 2,font.lab=2,mar=c(4,4,1,5), font=2)
  spheroid_diffusion$Col <- rbPal(50)[as.numeric(cut(spheroid_diffusion$DCMin..unit.,breaks = 50))]
  plot(spheroid_diffusion$DC1,spheroid_diffusion$DC2,col=spheroid_diffusion$Col, 
       xlab="Diffusion Component 1", ylab="Diffusion Component 2",las=1,pch=19,cex=0.5)
  colkey(col=rbPal(50),range(spheroid_diffusion$Circ.),side = 4, add = TRUE,length = 0.5, width = 0.5,clab = "Min Dist from center to surface",side.clab=4)
  dev.off()
  
  png(filename="Ellongation_diffusion.png", units="in", width=2.5, height=2 , pointsize=7, res=1200)
  par(font.axis = 2,font.lab=2,mar=c(4,4,1,5), font=2)
  spheroid_diffusion$Col <- rbPal(50)[as.numeric(cut(spheroid_diffusion$Ell_Elon,breaks = 50))]
  plot(spheroid_diffusion$DC1,spheroid_diffusion$DC2,col=spheroid_diffusion$Col, 
       xlab="Diffusion Component 1", ylab="Diffusion Component 2",las=1,pch=19,cex=0.5)
  colkey(col=rbPal(50),range(spheroid_diffusion$Circ.),side = 4, add = TRUE,length = 0.5, width = 0.5,clab = "Spheroid Ellongation",side.clab=4)
  dev.off()
  
  png(filename="DC_SD_diffusion.png", units="in", width=2.5, height=2 , pointsize=7, res=1200)
  par(font.axis = 2,font.lab=2,mar=c(4,4,1,5), font=2)
  spheroid_diffusion$Col <- rbPal(50)[as.numeric(cut(spheroid_diffusion$DCSD..unit.,breaks = 50))]
  plot(spheroid_diffusion$DC1,spheroid_diffusion$DC2,col=spheroid_diffusion$Col, 
       xlab="Diffusion Component 1", ylab="Diffusion Component 2",las=1,pch=19,cex=0.5)
  colkey(col=rbPal(50),range(spheroid_diffusion$Circ.),side = 4, add = TRUE,length = 0.5, width = 0.5,clab = "Stdev Dist from center to surface",side.clab=4)
  dev.off()
  

}

#reprogramming established
{
  rbPal <- colorRampPalette(c('deeppink','papayawhip','green4'))
  library(plot3D)
  
  ymin_max=c(min(spheroid_diffusion$DC2),max(spheroid_diffusion$DC2))
  xmin_max=c(min(spheroid_diffusion$DC1),max(spheroid_diffusion$DC1))
  
  png(filename="oct4_diffusion.png", units="in", width=2.5, height=2 , pointsize=7, res=1200)
  par(font.axis = 2,font.lab=2,mar=c(4,4,1,5), font=2)
  temp<-merge(spheroid_2dint_data$oct4,spheroid_diffusion, by="Label")
  temp$Col <- rbPal(50)[as.numeric(cut(temp$Mean,breaks = 50))]
  plot(temp$DC1,temp$DC2,col=temp$Col, cex=0.5, xlim=xmin_max,ylim=ymin_max,
       xlab="Diffusion Component 1", ylab="Diffusion Component 2",las=1,pch=19)
  colkey(col=rbPal(50),range(temp$Mean),side = 4, add = TRUE,length = 0.5, width = 0.5,clab = "Avg. Oct4 Levels",side.clab=4)
  dev.off()
  
  png(filename="vimentin_diffusion.png", units="in", width=2.5, height=2 , pointsize=7, res=1200)
  par(font.axis = 2,font.lab=2,mar=c(4,4,1,5), font=2)
  temp<-merge(spheroid_3dint_data$vimentin,spheroid_diffusion, by="Label")
  temp$Col <- rbPal(50)[as.numeric(cut(temp$Mean,breaks = 50))]
  plot(temp$DC1,temp$DC2,col=temp$Col, cex=0.5, xlim=xmin_max,ylim=ymin_max,
       xlab="Diffusion Component 1", ylab="Diffusion Component 2",las=1,pch=19)
  colkey(col=rbPal(50),range(temp$Mean),side = 4, add = TRUE,length = 0.5, width = 0.5,clab = "Avg. Vimentin Levels",side.clab=4)
  dev.off()
  rm(temp)
}
#replot_time
{
  rbPal <- colorRampPalette(c('firebrick','tomato2','tomato2','tomato2','yellow','green4','blue3','blue3','blue3'))
  library(plot3D)
  
  png(filename="time.png", units="in", width=2.5, height=2 , pointsize=7, res=1200)
  par(font.axis = 2,font.lab=2,mar=c(4,4,1,5), font=2)
  temp<-spheroid_diffusion
  temp$Col <- rbPal(100)[as.numeric(cut(temp$DPT,breaks = 100))]
  plot(temp$DC1,temp$DC2,col=temp$Col, cex=0.5, 
       xlab="Diffusion Component 1", ylab="Diffusion Component 2",las=1,pch=19)
  colkey(col=rbPal(100),range(temp$DPT),side = 4, add = TRUE,length = 0.5, width = 0.5,clab = "Psuedotime",side.clab=4)
  dev.off()
  
  a<-(subset(spheroid_diffusion,spheroid_diffusion$Label %in% D0_lab))
  b<-(subset(spheroid_diffusion,spheroid_diffusion$Label %in% D2_lab))
  c<-(subset(spheroid_diffusion,spheroid_diffusion$Label %in% D4_lab))
  d<-(subset(spheroid_diffusion,spheroid_diffusion$Label %in% D6_lab))
  e<-(subset(spheroid_diffusion,spheroid_diffusion$Label %in% D8_lab))
  a$sample<-"D0"
  b$sample<-"D2"
  c$sample<-"D4"
  d$sample<-"D6"
  e$sample<-"D8"
  
  temp<-rbind.data.frame(a,b,c,d,e)
  png(filename="time_DPT.png", units="in", width=2, height=2 , pointsize=5, res=1200)
  par(font.axis = 2,font.lab=2,mar=c(4,4,2,2), font=2)
  stripchart(temp$DPT~temp$sample,vertical = F,method = "jitter", add = F, pch = 20, 
             col = cols_timeline,xlab="pseudotime",las=1)
  dev.off()
}

# correlating protein features with diffusion analysis
{
  library(gplots)
  
  {
    temp<-merge(spheroid_2dint_data$oct4,spheroid_diffusion, by="Label")
    cor_oct4<-cor(temp$Mean,temp$DPT,method="spearman")
    temp<-merge(spheroid_3dint_data$actin,spheroid_diffusion, by="Label")
    cor_actin<-cor(temp$Mean,temp$DPT,method="spearman")
    temp<-merge(spheroid_2dint_data$pmlc,spheroid_diffusion, by="Label")
    cor_pmlc<-cor(temp$Min,temp$DPT,method="spearman")
    temp<-merge(spheroid_3dint_data$vimentin,spheroid_diffusion, by="Label")
    cor_vimentin<-cor(temp$Mean,temp$DPT,method="pearson")
    temp<-merge(spheroid_3dint_data$ecad,spheroid_diffusion, by="Label")
    cor_ecad<-cor(temp$Mean,temp$DPT,method="pearson")
    temp<-merge(spheroid_3dint_data$ki67,spheroid_diffusion, by="Label")
    cor_ki67<-cor(temp$Mean,temp$DPT,method="spearman")
    temp<-merge(spheroid_2dint_data$h3k9ac,spheroid_diffusion, by="Label")
    cor_h3k9ac<-cor(temp$IntDen,temp$DPT,method="pearson")
    temp<-merge(spheroid_3dint_data$laminac,spheroid_diffusion, by="Label")
    cor_laminac<-cor(temp$Mean,temp$DPT,method="spearman")
    cor_dot_proteins<-c(cor_oct4,cor_actin,cor_pmlc,cor_vimentin,cor_ecad,cor_ki67,cor_h3k9ac,cor_laminac)
    rm(temp,cor_oct4,cor_actin,cor_pmlc,cor_vimentin,cor_ecad,cor_ki67,cor_h3k9ac,cor_laminac)
    
  }
  
  {
    data<-spheroid_2dint_data$oct4
    a<-subset(data,data$Label %in% D0_lab)$Mean
    b<-subset(data,data$Label %in% D2_lab)$Mean
    c<-subset(data,data$Label %in% D4_lab)$Mean
    d<-subset(data,data$Label %in% D6_lab)$Mean
    e<-subset(data,data$Label %in% D8_lab)$Mean
    oct4<-c(median(a),median(b),median(c),median(d),median(e))
    
    data<-spheroid_3dint_data$actin
    a<-subset(data,data$Label %in% D0_lab)$Mean
    b<-subset(data,data$Label %in% D2_lab)$Mean
    c<-subset(data,data$Label %in% D4_lab)$Mean
    d<-subset(data,data$Label %in% D6_lab)$Mean
    e<-subset(data,data$Label %in% D8_lab)$Mean
    actin<-c(median(a),median(b),median(c),median(d),median(e))
    
    data<-spheroid_2dint_data$pmlc
    a<-subset(data,data$Label %in% D0_lab)$Min
    b<-subset(data,data$Label %in% D2_lab)$Min
    c<-subset(data,data$Label %in% D4_lab)$Min
    d<-subset(data,data$Label %in% D6_lab)$Min
    e<-subset(data,data$Label %in% D8_lab)$Min
    pmlc<-c(median(a),median(b),median(c),median(d),median(e))
    
    data<-spheroid_3dint_data$vimentin
    a<-subset(data,data$Label %in% D0_lab)$Mean
    b<-subset(data,data$Label %in% D2_lab)$Mean
    c<-subset(data,data$Label %in% D4_lab)$Mean
    d<-subset(data,data$Label %in% D6_lab)$Mean
    e<-subset(data,data$Label %in% D8_lab)$Mean
    vimentin<-c(mean(a),mean(b),mean(c),mean(d),mean(e))/mean(a)
    
    data<-spheroid_3dint_data$ecad
    a<-subset(data,data$Label %in% D0_lab)$Mean
    b<-subset(data,data$Label %in% D2_lab)$Mean
    c<-subset(data,data$Label %in% D4_lab)$Mean
    d<-subset(data,data$Label %in% D6_lab)$Mean
    e<-subset(data,data$Label %in% D8_lab)$Mean
    ecad<-c(median(a),median(b),median(c),median(d),median(e))
    
    data<-spheroid_3dint_data$ki67
    a<-subset(data,data$Label %in% D0_lab)$Mean
    b<-subset(data,data$Label %in% D2_lab)$Mean
    c<-subset(data,data$Label %in% D4_lab)$Mean
    d<-subset(data,data$Label %in% D6_lab)$Mean
    e<-subset(data,data$Label %in% D8_lab)$Mean
    ki67<-c(median(a),median(b),median(c),median(d),median(e))
    
    data<-spheroid_2dint_data$h3k9ac
    a<-subset(data,data$Label %in% D0_lab)$IntDen
    b<-subset(data,data$Label %in% D2_lab)$IntDen
    c<-subset(data,data$Label %in% D4_lab)$IntDen
    d<-subset(data,data$Label %in% D6_lab)$IntDen
    e<-subset(data,data$Label %in% D8_lab)$IntDen
    h3k9ac<-c(median(a),median(b),median(c),median(d),median(e))
    
    data<-spheroid_3dint_data$laminac
    a<-subset(data,data$Label %in% D0_lab)$Mean
    b<-subset(data,data$Label %in% D2_lab)$Mean
    c<-subset(data,data$Label %in% D4_lab)$Mean
    d<-subset(data,data$Label %in% D6_lab)$Mean
    e<-subset(data,data$Label %in% D8_lab)$Mean
    laminac<-c(median(a),median(b),median(c),median(d),median(e))
    
    temporal_protein_profile<-rbind(oct4,actin,pmlc,vimentin,ecad,ki67,h3k9ac,laminac)
    colnames(temporal_protein_profile)<-c("D0","D2","D4","D6","D8")
    rownames(temporal_protein_profile)<-c("Oct4","F-Actin","pMLC","Vimentin","ECadherin","KI67","H3K9Ac","LaminAC")
    
    rm(a,b,c,d,e,data,oct4,actin,pmlc,vimentin,ecad,ki67,h3k9ac,laminac)
  }
  
  cor_cols<-redblue(50)[as.numeric(cut(cor_dot_proteins,breaks = 50))]
  png(filename="correlation_DM.png", units="in", width=2, height=2 , pointsize=5, res=1200)
  par(font.axis = 2,font.lab=2,mar=c(4,4,2,2), font=2)
  heatmap.2(temporal_protein_profile, dendrogram ="none",Rowv = NA,Colv = NA,trace="none",density.info="none",srtCol = 0,cexCol = 1,
          key.title = "",cexRow = 0.9, col=redgreen(50),scale='row',RowSideColors = cor_cols)
  dev.off()
  png(filename="correlation_DM_key.png", units="in", width=2, height=2 , pointsize=5, res=1200)
  par(font.axis = 2,font.lab=2,mar=c(4,4,2,2), font=2)
  colkey(col=redblue(50),range(cor_dot_proteins),side = 1, add = F,length = 0.6, width = 1,clab = "Correlation with\n psuedotime",side.clab=1)
  dev.off()
  
  rm(cor_dot_proteins,temporal_protein_profile,cor_cols)
  
  
}


