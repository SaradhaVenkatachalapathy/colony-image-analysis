#functions
{
  open_and_compile <- function(dir, exp, samplen, trial,file_type){
    setwd(dir)
    file_data<-list.files()
    for (j in 1:length(file_data)){
      if (!exists("dataset")){
        if(file_type=="tsv"){
          dataset <- read.table(file_data[j], header=TRUE, sep="\t",stringsAsFactors=FALSE)
        }
        else if(file_type=="csv"){
          dataset <- read.table(file_data[j], header=TRUE, sep=",",stringsAsFactors=FALSE)
        }
      }
      else if (exists("dataset")){
        if(file_type=="tsv"){
          temp_dataset <- read.table(file_data[j], header=TRUE, sep="\t",stringsAsFactors=FALSE)
        }
        else if(file_type=="csv"){
          temp_dataset <- read.table(file_data[j], header=TRUE, sep=",",stringsAsFactors=FALSE)
        }
        if(nrow(temp_dataset)>0){
          dataset<-rbind(dataset, temp_dataset)
        }
        rm(temp_dataset)
      }
    }
    dataset$Exp<-exp
    dataset$sample<-samplen
    dataset$trial<-trial
    return(dataset)
  }
  combine_sample_sets<-function(path_to_ex,samples,data_type,file_type,exp, trial){
    dataset<-open_and_compile(paste(path_to_ex,samples[1],data_type,sep=""), 
                              exp, samples[1],trial,file_type)
    for (j in 2:length(samples)){
      temp<-open_and_compile(paste(path_to_ex,samples[j],data_type,sep=""), 
                             exp, samples[j],trial,file_type)
      dataset<-rbind(dataset,temp)
      rm(temp)
    }
    return(dataset)
    
  }
  combine_sample_sets_2d<-function(path_to_ex,samples,data_type,exp, trial){
    dataset<-read.csv(paste(path_to_ex,samples[1],data_type,sep=""),header=T, stringsAsFactors = F)
    dataset$Exp<-exp
    dataset$sample<-samples[1]
    dataset$trial<-trial
    for (j in 2:length(samples)){
      temp<-read.csv(paste(path_to_ex,samples[j],data_type,sep=""),header=T, stringsAsFactors = F)
      temp$Exp<-exp
      temp$sample<-samples[j]
      temp$trial<-trial
      dataset<-rbind(dataset,temp)
      rm(temp)
    }
    return(dataset)
    
  }
  open_and_compile_shells <- function(dir, exp, samplen, trial,file_type){
    
    setwd(dir)
    file_data<-list.files()
    for (j in 1:length(file_data)){
      if (!exists("dataset")){
        dataset <- read.table(file_data[j], header=TRUE, sep="\t",stringsAsFactors=FALSE)
        dataset$spheroidnum<-(j)
      }
      else if (exists("dataset")){
        temp_dataset <- read.table(file_data[j], header=TRUE, sep="\t",stringsAsFactors=FALSE)
        temp_dataset$spheroidnum<-(j)
        
        if(nrow(temp_dataset)>0){
          dataset<-rbind(dataset, temp_dataset)
        }
        rm(temp_dataset)
      }
    }
    dataset$Exp<-exp
    dataset$sample<-samplen
    dataset$trial<-trial
    return(dataset)
  }
  shell_volfraction<- function(dirg,diri, exp, sample, trial, dirres){
    
    setwd(dirg)
    file_data<-list.files()
    if(length(file_data) >2){
      geometrical_data_spheroid_shells<-open_and_compile_shells(dirg, 
                                                                exp, sample, trial,"tsv")
      oct4_intensity_spheroid_shells<-open_and_compile_shells(diri, 
                                                              exp, sample, trial,"tsv")
      
      intden_nom<-as.data.frame(matrix(nrow=length(file_data),ncol=33))
      colnames(intden_nom)<-c("spheroidnum","Shell0","Shell1","Shell2","Shell3","Shell4","Shell5","Shell6","Shell7","Shell8","Shell9","Shell10",
                              "Shell11","Shell12","Shell13","Shell14","Shell15","Shell16","Shell17","Shell18","Shell19","Shell20",
                              "Shell21","Shell22","Shell23","Shell24","Shell25","Shell26","Shell27","Shell28","Shell29","Shell30", "totint")
      intden_nom$nuc<-file_data
      vol_nom<-as.data.frame(matrix(nrow=length(file_data),ncol=34))
      colnames(vol_nom)<-c("spheroidnum","Shell0","Shell1","Shell2","Shell3","Shell4","Shell5","Shell6","Shell7","Shell8","Shell9","Shell10",
                           "Shell11","Shell12","Shell13","Shell14","Shell15","Shell16","Shell17","Shell18","Shell19","Shell20",
                           "Shell21","Shell22","Shell23","Shell24","Shell25","Shell26","Shell27","Shell28","Shell29","Shell30","totvol", "sphere")
      vol_nom$nuc<-file_data
      volring_nom<-as.data.frame(matrix(nrow=length(file_data),ncol=32))
      colnames(volring_nom)<-c("spheroidnum","Shell0","Shell1","Shell2","Shell3","Shell4","Shell5","Shell6","Shell7","Shell8","Shell9","Shell10",
                               "Shell11","Shell12","Shell13","Shell14","Shell15","Shell16","Shell17","Shell18","Shell19","Shell20",
                               "Shell21","Shell22","Shell23","Shell24","Shell25","Shell26","Shell27","Shell28","Shell29","Shell30")
      for( i in 1:length(file_data)){
        intsub<-subset(oct4_intensity_spheroid_shells,oct4_intensity_spheroid_shells$spheroidnum==i)
        intden_nom[i,33]<-intsub$IntDen[1]
        volsub<-subset(geometrical_data_spheroid_shells,geometrical_data_spheroid_shells$spheroidnum==i)
        vol_nom[i,33]<-volsub$Vol..unit.[1]
        vol_nom[i,34]<-volsub$Spher..unit.[1]
        
      }   
      for( i in 1:length(file_data)){
        intsub<-subset(oct4_intensity_spheroid_shells,oct4_intensity_spheroid_shells$spheroidnum==i)
        if(nrow(intsub)>1){
          for(k in 1:(nrow(intsub)-1)){
            intsub$IntDen[k]<-intsub$IntDen[k]-intsub$IntDen[k+1]
          }
          intden_nom[i,1]<-i
          intden_nom[i,2:(nrow(intsub)+1)]<-intsub$IntDen
          
        }
        
        
        volsub<-subset(geometrical_data_spheroid_shells,geometrical_data_spheroid_shells$spheroidnum==i)
        if(nrow(volsub)>1){
          for(k in 1:(nrow(volsub)-1)){
            volsub$Vol..unit.[k]<-volsub$Vol..unit.[k]-volsub$Vol..unit.[k+1]
          }
          volring_nom[i,1]<-i
          volring_nom[i,2:(nrow(volsub)+1)]<-volsub$Vol..unit.
          
          volsub<-subset(geometrical_data_spheroid_shells,geometrical_data_spheroid_shells$spheroidnum==i)
          volsub$Vol..unit. <-volsub$Vol..unit./volsub$Vol..unit.[1]
          vol_nom[i,1]<-i
          vol_nom[i,2:(nrow(volsub)+1)]<-volsub$Vol..unit.
        }
        
        
      }    
      
      intden_nom[,2]<-intden_nom[,2]/volring_nom[,2]
      intden_nom[,3]<-intden_nom[,3]/volring_nom[,3]
      intden_nom[,4]<-intden_nom[,4]/volring_nom[,4]
      intden_nom[,5]<-intden_nom[,5]/volring_nom[,5]
      intden_nom[,6]<-intden_nom[,6]/volring_nom[,6]
      intden_nom[,7]<-intden_nom[,7]/volring_nom[,7]
      intden_nom[,8]<-intden_nom[,8]/volring_nom[,8]
      intden_nom[,9]<-intden_nom[,9]/volring_nom[,9]
      intden_nom[,10]<-intden_nom[,10]/volring_nom[,10]
      intden_nom[,11]<-intden_nom[,11]/volring_nom[,11]
      intden_nom[,12]<-intden_nom[,12]/volring_nom[,12]
      intden_nom[,13]<-intden_nom[,13]/volring_nom[,13]
      intden_nom[,14]<-intden_nom[,14]/volring_nom[,14]
      intden_nom[,15]<-intden_nom[,15]/volring_nom[,15]
      intden_nom[,16]<-intden_nom[,16]/volring_nom[,16]
      intden_nom[,17]<-intden_nom[,17]/volring_nom[,17]
      intden_nom[,18]<-intden_nom[,18]/volring_nom[,18]
      intden_nom[,19]<-intden_nom[,19]/volring_nom[,19]
      intden_nom[,20]<-intden_nom[,20]/volring_nom[,20]
      intden_nom[,21]<-intden_nom[,21]/volring_nom[,21]
      intden_nom[,22]<-intden_nom[,22]/volring_nom[,22]
      intden_nom[,23]<-intden_nom[,23]/volring_nom[,23]
      intden_nom[,24]<-intden_nom[,24]/volring_nom[,24]
      intden_nom[,25]<-intden_nom[,25]/volring_nom[,25]
      intden_nom[,26]<-intden_nom[,26]/volring_nom[,26]
      intden_nom[,27]<-intden_nom[,27]/volring_nom[,27]
      intden_nom[,28]<-intden_nom[,28]/volring_nom[,28]
      intden_nom[,29]<-intden_nom[,29]/volring_nom[,29]
      intden_nom[,30]<-intden_nom[,30]/volring_nom[,30]
      intden_nom[,31]<-intden_nom[,31]/volring_nom[,31]
      intden_nom[,32]<-intden_nom[,32]/volring_nom[,32]
      
      for( i in 1:nrow(intden_nom)){
        intden_nom[i,2:32]<-intden_nom[i,2:32]/max(intden_nom[i,2:32],na.rm=T)
      }
      
      intden_nom$Label<-file_data
      vol_nom$Label<-file_data
      
      setwd(dirres)
      
      name=paste("Shell_volume_fraction",sample,trial,"lines.png",sep="_")
      png(filename=name, units="in",width=1, height=1 , pointsize=4, res=1200)
      par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2)
      plot(1:0,1:0,col="white",xlim=c(0,1),xlab="Volume fraction",ylab="Intensity fraction",las=1)
      cols<-rgb(runif(length(file_data)),runif(length(file_data)),runif(length(file_data))) 
      for( i in 1:length(file_data)){
        lines(vol_nom[i,2:32],intden_nom[i,2:32],col=cols[i],type="o",pch=18)
      }
      dev.off()
      
      fractions<-as.data.frame(matrix(nrow=length(file_data),ncol=11))
      fractions[,1]<-file_data
      colnames(fractions)<-c("Label","0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%")
      for( i in 1:length(file_data)){
        y<-(na.omit(unlist(intden_nom[i,2:32])))
        x<-(na.omit(unlist(vol_nom[i,2:32])))
        if(length(unique(x))>=4 & length(unique(x))==length(unique(y))){
          lines(y~x, col = cols[i], pch = 19,type="o")
          points(y~x, col = cols[i], pch = 19)
          fit <- smooth.spline(x, y)
          prd<-predict(fit, x)$y
          lines(prd~x,col=cols[i])
          
          fractions[i,2]<-predict(fit, 0.1)$y
          fractions[i,3]<-predict(fit, 0.2)$y
          fractions[i,4]<-predict(fit, 0.3)$y
          fractions[i,5]<-predict(fit, 0.4)$y
          fractions[i,6]<-predict(fit, 0.5)$y
          fractions[i,7]<-predict(fit, 0.6)$y
          fractions[i,8]<-predict(fit, 0.7)$y
          fractions[i,9]<-predict(fit, 0.8)$y
          fractions[i,10]<-predict(fit, 0.9)$y
          fractions[i,11]<-predict(fit, 1)$y
        }
        
      }
      
      mergedfile<-merge(intden_nom[,c(33,35)], fractions, by="Label")
      mergedfile<-merge(vol_nom[,c(33,34,36)], mergedfile, by="Label")
      
      library(gplots)
      library(RColorBrewer)
      Colors=brewer.pal(11,"Spectral")
      
      fractions<-fractions[order(fractions[,2]),]
      row.names(fractions)<-fractions[,1]
      temp<-na.omit(as.matrix(fractions[,2:11]))
      row.names(fractions)<-fractions[,1]
      if(nrow(temp)>2 & ncol(temp)>0){
        name=paste("Shell_volume_fraction",sample,trial,"heatmap.png",sep="_")
        png(filename=name, units="in",width=2.5, height=2 , pointsize=5, res=1200)
        par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2)
        data_heat<-heatmap.2(temp,scale="none",Rowv = TRUE, Colv=FALSE, dendrogram = "row", density.info = "none",trace="none",col=Colors, 
                             labRow = FALSE)
        dev.off()
      }
      
      
      mergedfile$exp<-exp
      mergedfile$sample<-sample
      mergedfile$trial<-trial
      
      return(mergedfile)
    }
    else
      
      return(data.frame(matrix(ncol = 17, nrow = 0)))
    
  }
  shell_analsysis<-function(path_to_exn,samplesn,data_type,expn, trialn, dird){
    temp<-shell_volfraction(paste(path_to_exn,samplesn[1],"/Shells/geometric_measures",sep=""), 
                            paste(path_to_exn,samplesn[1],data_type,sep=""),
                            expn, samplesn[1],trialn,dird)
    sample1<-temp[complete.cases(temp), ]
    for (j in 2:length(samples)){
      temp<-shell_volfraction(paste(path_to_exn,samplesn[j],"/Shells/geometric_measures",sep=""), 
                              paste(path_to_exn,samplesn[j],data_type,sep=""),
                              expn, samplesn[j],trialn,dird)
      tempsample1<-temp[complete.cases(temp), ]
      sample1<-rbind(sample1,tempsample1)
      rm(temp,tempsample1)
    }
    
    sample1$totint<-log(sample1$totint)
    sample1$totvol<-log(sample1$totvol)
    sample1$Label<-sub(".*M_","",sub(".tsv*","",sample1$Label))
    
    df=sample1
    for (i in 1:nrow(df)){
      df$maxint[i]<-which(df[i,5:14] %in% max(df[i,5:14]))
      df$minint[i]<-which(df[i,5:14] %in% min(df[i,5:14]))
      df$peripheral_index_sp[i]<-cor(unlist(df[i,5:14]),c(10,20,30,40,50,60,70,80,90,100),  method="spearman") 
      df$central_index_sp[i]<-cor(unlist(df[i,5:14]),rev(c(10,20,30,40,50,60,70,80,90,100)),  method="spearman") 
      df$peripheral_index_pc[i]<-cor(unlist(df[i,5:14]),c(10,20,30,40,50,60,70,80,90,100),  method="pearson") 
      df$central_index_pc[i]<-cor(unlist(df[i,5:14]),rev(c(10,20,30,40,50,60,70,80,90,100)),  method="pearson") 
    }
    return(df)
    
  }     
  open_and_compile_axial_measures <- function(dir, exp, samplen, trial){
    setwd(dir)
    file_data<-list.files()
    axial_intensity<-(matrix(nrow=45,ncol=length(file_data)))
    for (j in 1:length(file_data)){
      temp<-read.table(file_data[j], header=TRUE, sep=",",stringsAsFactors=FALSE)
      axial_intensity[1:nrow(temp),j]<-temp$Mean
      rm(temp)
    }
    axial_intensity<-as.data.frame(t(axial_intensity))
    colnames(axial_intensity)<-paste("Zposition",seq(1,45,1),sep="_")
    rownames(axial_intensity)<-sub("\\..*", "", file_data)
    
    axial_intensity$Exp<-exp
    axial_intensity$sample<-samplen
    axial_intensity$trial<-trial
    return(axial_intensity)
  }
  axial_analysis<-function(path_to_exn,samplesn,expn, trialn,data_type){
    sample1<-open_and_compile_axial_measures(paste(path_to_exn,samplesn[1],"/stepwise_measures_spheroid/",data_type,sep=""),
                                             expn, samplesn[1],trialn)
    for (j in 2:length(samples)){
      tempsample1<-open_and_compile_axial_measures(paste(path_to_exn,samplesn[j],"/stepwise_measures_spheroid/",data_type,sep=""),
                                                   expn, samplesn[j],trialn)
      sample1<-rbind(sample1,tempsample1)
      rm(tempsample1)
    }
    
    df=sample1
    for (i in 1:nrow(df)){
      x<-na.omit(unlist(df[i,1:45]))
      if(mean(x)!="NaN"){
        df$maxint[i]<-which(x %in% max(x))-floor(length(x)/2)
        df$minint[i]<-which(x %in% min(x))-floor(length(x)/2)
        df$axial_index_biased_sp[i]<-cor(x,seq(1,length(x),1),  method="spearman") 
        df$basal_index_biased_sp[i]<-cor(x,rev(seq(1,length(x),1)),  method="spearman")
        df$axial_index_biased_pc[i]<-cor(x,seq(1,length(x),1),  method="pearson") 
        df$basal_index_biased_pc[i]<-cor(x,rev(seq(1,length(x),1)),  method="pearson")
        x1<-c(x,rep(0, (45-length(x))))
        df$axial_index_sp[i]<-cor(x1,seq(1,length(x1),1),  method="spearman") 
        df$basal_index_sp[i]<-cor(x1,rev(seq(1,length(x1),1)),  method="spearman")
        df$axial_index_pc[i]<-cor(x1,seq(1,length(x1),1),  method="pearson") 
        df$basal_index_pc[i]<-cor(x1,rev(seq(1,length(x1),1)),  method="pearson")
        
      }
    }
    sample1=df
    sample1$Label<-rownames(sample1)
    return(sample1)
    
  }
  open_and_compile_angles <- function(dir, exp, samplen, trial){
    setwd(dir)
    file_data<-list.files()
    angles<-as.data.frame(matrix(ncol=6,nrow=length(file_data)))
    colnames(angles)<-c("Label","max_freq_angle","min_freq_angle","meadian_freq_angle","mean_freq_angle","error_sd_freq")
    for (j in 1:length(file_data)){
      temp<-read.table(file_data[j], header=TRUE, sep=",",stringsAsFactors=FALSE)
      
      max_freq_angle<-temp[which(temp[,2]==max(temp[,2])),1]
      min_freq_angle<-temp[which(temp[,2]==min(temp[,2])),1]
      meadian_freq_angle<-mean(temp[which(abs(temp[,2]-median(temp[,2]))==min(abs(temp[,2]-median(temp[,2])))),1])
      mean_freq_angle<-mean(temp[which(abs(temp[,2]-mean(temp[,2]))==min(abs(temp[,2]-mean(temp[,2])))),1])
      error_sd_freq<-sd(temp[,2])*100
      
      angles[j,1]<-file_data[j]
      angles[j,2]<-max_freq_angle
      angles[j,3]<-min_freq_angle
      angles[j,4]<-meadian_freq_angle
      angles[j,5]<-mean_freq_angle
      angles[j,6]<-error_sd_freq
      
      rm(temp)
    }
    angles$Exp<-exp
    angles$sample<-samplen
    angles$trial<-trial
    return(angles)
  }
  angle_analsysis<-function(path_to_exn,samplesn,data_type,expn, trialn){
    
    sample1<-open_and_compile_angles(paste(path_to_exn,samples[1],"/angles_spheroid/",data_type,sep=""),
                                     expn, samplesn[1],trialn)
    for (j in 2:length(samples)){
      tempsample1<-open_and_compile_angles(paste(path_to_exn,samplesn[j],"/angles_spheroid/",data_type,sep=""),
                                           expn, samplesn[j],trialn)
      sample1<-rbind(sample1,tempsample1)
      rm(tempsample1)
    }
    
    sample1$Label<-substring(sample1$Label,1, nchar(sample1$Label)-4)
    
    return(sample1)
    
  }     
  
  read_glcm<-function(path_to_ex,samples,data_types_2d,exp,trial){
    
    T1_gclm0_5_2D_data<-combine_sample_sets_2d(path_to_ex,samples,data_types_2d[4],exp,trial)
    colnames(T1_gclm0_5_2D_data)[2:6]<-paste(colnames(T1_gclm0_5_2D_data)[2:6],"_step5");
    T1_gclm90_5_2D_data<-combine_sample_sets_2d(path_to_ex,samples,data_types_2d[5],exp,trial)
    colnames(T1_gclm90_5_2D_data)[2:6]<-paste(colnames(T1_gclm90_5_2D_data)[2:6],"_step5");
    T1_gclm180_5_2D_data<-combine_sample_sets_2d(path_to_ex,samples,data_types_2d[6],exp,trial)
    colnames(T1_gclm180_5_2D_data)[2:6]<-paste(colnames(T1_gclm180_5_2D_data)[2:6],"_step5");
    T1_gclm270_5_2D_data<-combine_sample_sets_2d(path_to_ex,samples,data_types_2d[7],exp,trial)
    colnames(T1_gclm270_5_2D_data)[2:6]<-paste(colnames(T1_gclm270_5_2D_data)[2:6],"_step5");
    T1_gclm_5<-merge(T1_gclm0_5_2D_data,T1_gclm90_5_2D_data[,1:6],by="Label")
    T1_gclm_5<-merge(T1_gclm_5,T1_gclm180_5_2D_data[,1:6],by="Label")
    T1_gclm_5<-merge(T1_gclm_5,T1_gclm270_5_2D_data[,1:6],by="Label")
    rm(T1_gclm0_5_2D_data,T1_gclm90_5_2D_data,T1_gclm180_5_2D_data,T1_gclm270_5_2D_data)
    
    T1_gclm0_15_2D_data<-combine_sample_sets_2d(path_to_ex,samples,data_types_2d[8],exp,trial)
    colnames(T1_gclm0_15_2D_data)[2:6]<-paste(colnames(T1_gclm0_15_2D_data)[2:6],"_step15");
    T1_gclm90_15_2D_data<-combine_sample_sets_2d(path_to_ex,samples,data_types_2d[9],exp,trial)
    colnames(T1_gclm90_15_2D_data)[2:6]<-paste(colnames(T1_gclm90_15_2D_data)[2:6],"_step15");
    T1_gclm180_15_2D_data<-combine_sample_sets_2d(path_to_ex,samples,data_types_2d[10],exp,trial)
    colnames(T1_gclm180_15_2D_data)[2:6]<-paste(colnames(T1_gclm180_15_2D_data)[2:6],"_step15");
    T1_gclm270_15_2D_data<-combine_sample_sets_2d(path_to_ex,samples,data_types_2d[11],exp,trial)
    colnames(T1_gclm270_15_2D_data)[2:6]<-paste(colnames(T1_gclm270_15_2D_data)[2:6],"_step15");
    T1_gclm_15<-merge(T1_gclm0_15_2D_data,T1_gclm90_15_2D_data[,1:6],by="Label")
    T1_gclm_15<-merge(T1_gclm_15,T1_gclm180_15_2D_data[,1:6],by="Label")
    T1_gclm_15<-merge(T1_gclm_15,T1_gclm270_15_2D_data[,1:6],by="Label")
    rm(T1_gclm0_15_2D_data,T1_gclm90_15_2D_data,T1_gclm180_15_2D_data,T1_gclm270_15_2D_data)
    
    T1_gclm0_25_2D_data<-combine_sample_sets_2d(path_to_ex,samples,data_types_2d[12],exp,trial)
    colnames(T1_gclm0_25_2D_data)[2:6]<-paste(colnames(T1_gclm0_25_2D_data)[2:6],"_step25");
    T1_gclm90_25_2D_data<-combine_sample_sets_2d(path_to_ex,samples,data_types_2d[13],exp,trial)
    colnames(T1_gclm90_25_2D_data)[2:6]<-paste(colnames(T1_gclm90_25_2D_data)[2:6],"_step25");
    T1_gclm180_25_2D_data<-combine_sample_sets_2d(path_to_ex,samples,data_types_2d[14],exp,trial)
    colnames(T1_gclm180_25_2D_data)[2:6]<-paste(colnames(T1_gclm180_25_2D_data)[2:6],"_step25");
    T1_gclm270_25_2D_data<-combine_sample_sets_2d(path_to_ex,samples,data_types_2d[15],exp,trial)
    colnames(T1_gclm270_25_2D_data)[2:6]<-paste(colnames(T1_gclm270_25_2D_data)[2:6],"_step25");
    T1_gclm_25<-merge(T1_gclm0_25_2D_data,T1_gclm90_25_2D_data[,1:6],by="Label")
    T1_gclm_25<-merge(T1_gclm_25,T1_gclm180_25_2D_data[,1:6],by="Label")
    T1_gclm_25<-merge(T1_gclm_25,T1_gclm270_25_2D_data[,1:6],by="Label")
    rm(T1_gclm0_25_2D_data,T1_gclm90_25_2D_data,T1_gclm180_25_2D_data,T1_gclm270_25_2D_data)
    
    
    T1_gclm0_35_2D_data<-combine_sample_sets_2d(path_to_ex,samples,data_types_2d[16],exp,trial)
    colnames(T1_gclm0_35_2D_data)[2:6]<-paste(colnames(T1_gclm0_35_2D_data)[2:6],"_step35");
    T1_gclm90_35_2D_data<-combine_sample_sets_2d(path_to_ex,samples,data_types_2d[17],exp,trial)
    colnames(T1_gclm90_35_2D_data)[2:6]<-paste(colnames(T1_gclm90_35_2D_data)[2:6],"_step35");
    T1_gclm180_35_2D_data<-combine_sample_sets_2d(path_to_ex,samples,data_types_2d[18],exp,trial)
    colnames(T1_gclm180_35_2D_data)[2:6]<-paste(colnames(T1_gclm180_35_2D_data)[2:6],"_step35");
    T1_gclm270_35_2D_data<-combine_sample_sets_2d(path_to_ex,samples,data_types_2d[19],exp,trial)
    colnames(T1_gclm270_35_2D_data)[2:6]<-paste(colnames(T1_gclm270_35_2D_data)[2:6],"_step35");
    T1_gclm_35<-merge(T1_gclm0_35_2D_data,T1_gclm90_35_2D_data[,1:6],by="Label")
    T1_gclm_35<-merge(T1_gclm_35,T1_gclm180_35_2D_data[,1:6],by="Label")
    T1_gclm_35<-merge(T1_gclm_35,T1_gclm270_35_2D_data[,1:6],by="Label")
    rm(T1_gclm0_35_2D_data,T1_gclm90_35_2D_data,T1_gclm180_35_2D_data,T1_gclm270_35_2D_data)
    
    T1_gclm0_45_2D_data<-combine_sample_sets_2d(path_to_ex,samples,data_types_2d[20],exp,trial)
    colnames(T1_gclm0_45_2D_data)[2:6]<-paste(colnames(T1_gclm0_45_2D_data)[2:6],"_step45");
    T1_gclm90_45_2D_data<-combine_sample_sets_2d(path_to_ex,samples,data_types_2d[21],exp,trial)
    colnames(T1_gclm90_45_2D_data)[2:6]<-paste(colnames(T1_gclm90_45_2D_data)[2:6],"_step45");
    T1_gclm180_45_2D_data<-combine_sample_sets_2d(path_to_ex,samples,data_types_2d[22],exp,trial)
    colnames(T1_gclm180_45_2D_data)[2:6]<-paste(colnames(T1_gclm180_45_2D_data)[2:6],"_step45");
    T1_gclm270_45_2D_data<-combine_sample_sets_2d(path_to_ex,samples,data_types_2d[23],exp,trial)
    colnames(T1_gclm270_45_2D_data)[2:6]<-paste(colnames(T1_gclm270_45_2D_data)[2:6],"_step45");
    T1_gclm_45<-merge(T1_gclm0_45_2D_data,T1_gclm90_45_2D_data[,1:6],by="Label")
    T1_gclm_45<-merge(T1_gclm_45,T1_gclm180_45_2D_data[,1:6],by="Label")
    T1_gclm_45<-merge(T1_gclm_45,T1_gclm270_45_2D_data[,1:6],by="Label")
    rm(T1_gclm0_45_2D_data,T1_gclm90_45_2D_data,T1_gclm180_45_2D_data,T1_gclm270_45_2D_data)
    
    T1_gclm0_100_2D_data<-combine_sample_sets_2d(path_to_ex,samples,data_types_2d[24],exp,trial)
    colnames(T1_gclm0_100_2D_data)[2:6]<-paste(colnames(T1_gclm0_100_2D_data)[2:6],"_step100");
    T1_gclm90_100_2D_data<-combine_sample_sets_2d(path_to_ex,samples,data_types_2d[25],exp,trial)
    colnames(T1_gclm90_100_2D_data)[2:6]<-paste(colnames(T1_gclm90_100_2D_data)[2:6],"_step100");
    T1_gclm180_100_2D_data<-combine_sample_sets_2d(path_to_ex,samples,data_types_2d[26],exp,trial)
    colnames(T1_gclm180_100_2D_data)[2:6]<-paste(colnames(T1_gclm180_100_2D_data)[2:6],"_step100");
    T1_gclm270_100_2D_data<-combine_sample_sets_2d(path_to_ex,samples,data_types_2d[27],exp,trial)
    colnames(T1_gclm270_100_2D_data)[2:6]<-paste(colnames(T1_gclm270_100_2D_data)[2:6],"_step100");
    T1_gclm_100<-merge(T1_gclm0_100_2D_data,T1_gclm90_100_2D_data[,1:6],by="Label")
    T1_gclm_100<-merge(T1_gclm_100,T1_gclm180_100_2D_data[,1:6],by="Label")
    T1_gclm_100<-merge(T1_gclm_100,T1_gclm270_100_2D_data[,1:6],by="Label")
    rm(T1_gclm0_100_2D_data,T1_gclm90_100_2D_data,T1_gclm180_100_2D_data,T1_gclm270_100_2D_data)
    
    T1_gclm<-cbind(T1_gclm_5[,-c(7:10)],T1_gclm_15[,-c(1,7:10),],
                   T1_gclm_25[,-c(1,7:10)],T1_gclm_35[,-c(1,7:10)],
                   T1_gclm_45[,-c(1,7:10)],T1_gclm_100)
    rm(T1_gclm_5,T1_gclm_15,T1_gclm_25,T1_gclm_35,T1_gclm_45,T1_gclm_100)
    return(T1_gclm)
  }
  
  
  plot_black_bg_histogram_timeline<-function(data,colstimeline,parameters,labelling, prefix){
    
    for( i in 1:length(parameters)){
      index<-which(colnames(data) %in% parameters[i])
      filename_png<-paste(prefix,parameters[i],"hist.png",sep="_")
      png(filename=filename_png, units="in",width=2, height=2 , pointsize=6, res=1200)
      par(font.axis = 2,font.lab=2,mar=c(4,4.5,1,1), font=2,bg='black',fg='white',col.lab='white',col.axis='white')
      a<-na.omit(subset(data, data$sample=="D0")[,index])
      b<-na.omit(subset(data, data$sample=="D2")[,index])
      c<-na.omit(subset(data, data$sample=="D4")[,index])
      d<-na.omit(subset(data, data$sample=="D6")[,index])
      e<-na.omit(subset(data, data$sample=="D8")[,index])
      if(length(a)>0) d1<-density(a) else d1<-density(c(0,0,0,0,0,0,0,0,0))
      if(length(b)>0) d2<-density(b) else d2<-density(c(0,0,0,0,0,0,0,0,0))
      if(length(c)>0) d3<-density(c) else d3<-density(c(0,0,0,0,0,0,0,0,0))
      if(length(d)>0) d4<-density(d) else d4<-density(c(0,0,0,0,0,0,0,0,0))
      if(length(e)>0) d5<-density(e) else d5<-density(c(0,0,0,0,0,0,0,0,0))
      
      ymax<-max(c(d1$y,d2$y,d3$y,d4$y,d5$y))*1.2
      xmax<-max(c(d1$x,d2$x,d3$x,d4$x,d5$x))*1.1
      xmin<-min(c(d1$x,d2$x,d3$x,d4$x,d5$x))*0.9
      plot(d3,main="",col=NA,lwd=2,xlab=labelling[i],ylim=c(0,ymax),las=1,cex.axis=0.8, xlim=c(xmin,xmax), ylab="Probability Density")
      lines(d1,col=colstimeline[1],lwd=2)
      lines(d2,col=colstimeline[2],lwd=2)
      lines(d3,col=colstimeline[3],lwd=2)
      lines(d4,col=colstimeline[4],lwd=2)
      lines(d5,col=colstimeline[5],lwd=2)
      legend('topright',legend=c("D0","D2","D4","D6","D8"),fill=colstimeline, horiz=TRUE,cex=0.8)
      dev.off()
    }
    
    
    
  }
  plot_black_bg_boxplot_timeline<-function(data,colstimeline,parameters,labelling, prefix){
    for( i in 1:length(parameters)){
      index<-which(colnames(data) %in% parameters[i])
      filename_png<-paste(prefix,parameters[i],"box.png",sep="_")
      png(filename=filename_png, units="in",width=2, height=2 , pointsize=6, res=1200)
      par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2,bg='black',fg='white',col.lab='white',col.axis='white')
      a<-subset(data, data$sample=="D0")[,index]
      b<-subset(data, data$sample=="D2")[,index]
      c<-subset(data, data$sample=="D4")[,index]
      d<-subset(data, data$sample=="D6")[,index]
      e<-subset(data, data$sample=="D8")[,index]
      xmax<-max(data[,index])*1.2
      xmin<-min(data[,index])*0.9
      boxplot(a,b,c,d,e,main="",col=colstimeline,xlab=labelling[i],las=1,cex.axis=1,horizontal = T, 
              names=c("D0","D2","D4","D6","D8"),lty=1,lwd=0.5, pch=18)
      dev.off()
    }
    
  }
  plot_black_bg_barplot_timeline<-function(data,colstimeline,parameters,labelling, prefix){
    
    for( i in 1:length(parameters)){
      index<-which(colnames(data) %in% parameters[i])
      filename_png<-paste(prefix,parameters[i],"bar.png",sep="_")
      
      T1<-subset(data,data$trial=="T1")
      T1_means<-c(mean(subset(T1, T1$sample=="D0")[,index], na.rm = T),mean(subset(T1, T1$sample=="D2")[,index], na.rm = T),
                  mean(subset(T1, T1$sample=="D4")[,index], na.rm = T),mean(subset(T1, T1$sample=="D6")[,index], na.rm = T),
                  mean(subset(T1, T1$sample=="D8")[,index], na.rm = T))
      
      T1_means[is.nan(T1_means)] <- NA
      T1_means<-as.vector(na.omit(T1_means))
      T1_means<-T1_means/T1_means[1]
      rm(T1)
      
      T1<-subset(data,data$trial=="T2")
      T2_means<-c(mean(subset(T1, T1$sample=="D0")[,index], na.rm = T),mean(subset(T1, T1$sample=="D2")[,index], na.rm = T),
                  mean(subset(T1, T1$sample=="D4")[,index], na.rm = T),mean(subset(T1, T1$sample=="D6")[,index], na.rm = T),
                  mean(subset(T1, T1$sample=="D8")[,index], na.rm = T))
      T2_means[is.nan(T2_means)] <- NA
      T2_means<-as.vector(na.omit(T2_means))
      T2_means<-T2_means/T2_means[1]
      rm(T1)
      T1<-subset(data,data$trial=="T3")
      T3_means<-c(mean(subset(T1, T1$sample=="D0")[,index], na.rm = T),mean(subset(T1, T1$sample=="D2")[,index], na.rm = T),
                  mean(subset(T1, T1$sample=="D4")[,index], na.rm = T),mean(subset(T1, T1$sample=="D6")[,index], na.rm = T),
                  mean(subset(T1, T1$sample=="D8")[,index], na.rm = T))
      T3_means[is.nan(T3_means)] <- NA
      T3_means<-as.vector(na.omit(T3_means))
      T3_means<-T3_means/T3_means[1]
      rm(T1)
      
      if(length(T1_means)==5 & length(T2_means)==5 & length(T3_means)==5){
        pop_means<-c(1,mean(c(T1_means[2],T2_means[2],T3_means[2])),
                     mean(c(T1_means[3],T2_means[3],T3_means[3])),mean(c(T1_means[4],T2_means[4],T3_means[4])),
                     mean(c(T1_means[5],T2_means[5],T3_means[5])))
        SE_means<-c(0,sd(c(T1_means[2],T2_means[2],T3_means[2])),
                    sd(c(T1_means[3],T2_means[3],T3_means[3])),sd(c(T1_means[4],T2_means[4],T3_means[4])),
                    sd(c(T1_means[5],T2_means[5],T3_means[5])))/sqrt(3)
        
        png(filename=filename_png, units="in",width=2, height=2 , pointsize=6, res=1200)
        par(font.axis = 2,font.lab=2,mar=c(4,6,1,1), font=2,bg='black',fg='white',col.lab='white',col.axis='white')
        ymax=max(pop_means+SE_means)*1.2
        ymin=min(pop_means-SE_means)*0.8
        foo<-barplot(pop_means,las=1, col=colstimeline,ylim=c(ymin,ymax),names.arg = c("D0","D2","D4","D6","D8"),ylab=labelling[i],xpd = FALSE)
        arrows(x0=foo,y0=pop_means+SE_means,y1=pop_means-SE_means,angle=90,code=3,length=0.1)
        box()
        dev.off()
      }
      else{
        fill_num<-5-length(T1_means)
        temp<-rep(0,fill_num)
        temp[(fill_num+1):5]<-T1_means
        T1_means<-temp
        fill_num<-5-length(T2_means)
        temp<-rep(0,fill_num)
        temp[(fill_num+1):5]<-T2_means
        T2_means<-temp
        fill_num<-5-length(T3_means)
        temp<-rep(0,fill_num)
        temp[(fill_num+1):5]<-T3_means
        T3_means<-temp
        rm(fill_num,temp)
        
        pop_means<-c(mean(c(T1_means[1],T2_means[1],T3_means[1])),mean(c(T1_means[2],T2_means[2],T3_means[2])),
                     mean(c(T1_means[3],T2_means[3],T3_means[3])),mean(c(T1_means[4],T2_means[4],T3_means[4])),
                     mean(c(T1_means[5],T2_means[5],T3_means[5])))
        SE_means<-c(sd(c(T1_means[1],T2_means[1],T3_means[1])),sd(c(T1_means[2],T2_means[2],T3_means[2])),
                    sd(c(T1_means[3],T2_means[3],T3_means[3])),sd(c(T1_means[4],T2_means[4],T3_means[4])),
                    sd(c(T1_means[5],T2_means[5],T3_means[5])))/sqrt(3)
        
        png(filename=filename_png, units="in",width=2, height=2 , pointsize=6, res=1200)
        par(font.axis = 2,font.lab=2,mar=c(4,6,1,1), font=2,bg='black',fg='white',col.lab='white',col.axis='white')
        ymax=max(pop_means+SE_means)*1.2
        ymin=min(pop_means-SE_means)*0.8
        foo<-barplot(pop_means,las=1, col=colstimeline,ylim=c(ymin,ymax),names.arg = c("D0","D2","D4","D6","D8"),ylab=labelling[i],xpd = FALSE)
        arrows(x0=foo,y0=pop_means+SE_means,y1=pop_means-SE_means,angle=90,code=3,length=0.1)
        box()
        dev.off()
      }
      
      
      
      
    }
  }
  axial_plot_black_bg_boxplot_timeline<-function(data,colstimeline,ylabel, filename_png){
    D0_mean<-colMeans(subset(data,data$sample=="D0")[1:45])
    D0_sd<-apply(subset(data,data$sample=="D0")[1:45],2,sd)
    D2_mean<-colMeans(subset(data,data$sample=="D2")[1:45])
    D2_sd<-apply(subset(data,data$sample=="D2")[1:45],2,sd)
    D4_mean<-colMeans(subset(data,data$sample=="D4")[1:45])
    D4_sd<-apply(subset(data,data$sample=="D4")[1:45],2,sd)
    D6_mean<-colMeans(subset(data,data$sample=="D6")[1:45])
    D6_sd<-apply(subset(data,data$sample=="D6")[1:45],2,sd)
    D8_mean<-colMeans(subset(data,data$sample=="D8")[1:45])
    D8_sd<-apply(subset(data,data$sample=="D8")[1:45],2,sd)
    
    png(filename=filename_png, units="in",width=2, height=2 , pointsize=6, res=1200)
    par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2,bg='black',fg='white',col.lab='white',col.axis='white')
    plot(D0_mean, type="l",col=cols[1],xlim=c(0,30),ylim=c(0,40),las=1,ylab=ylabel,xlab="Zposition")
    lines(D2_mean, type="l",col=cols[2],xlim=c(0,30),ylim=c(0,40),las=1,ylab=ylabel,xlab="Zposition")
    lines(D4_mean, type="l",col=cols[3],xlim=c(0,30),ylim=c(0,40),las=1,ylab=ylabel,xlab="Zposition")
    lines(D6_mean, type="l",col=cols[4],xlim=c(0,30),ylim=c(0,40),las=1,ylab=ylabel,xlab="Zposition")
    lines(D8_mean, type="l",col=cols[5],xlim=c(0,30),ylim=c(0,40),las=1,ylab=ylabel,xlab="Zposition")
    segments(1:45, D0_mean - D0_sd, 1:45, D0_mean + D0_sd,col=cols[1]);
    segments(1:45, D2_mean - D2_sd, 1:45, D2_mean + D2_sd,col=cols[2]);
    segments(1:45, D4_mean - D4_sd, 1:45, D4_mean + D4_sd,col=cols[3]);
    segments(1:45, D6_mean - D6_sd, 1:45, D6_mean + D6_sd,col=cols[4]);
    segments(1:45, D8_mean - D8_sd, 1:45, D8_mean + D8_sd,col=cols[5]);
    dev.off()
    
    
  }
  radial_plot_black_bg_boxplot_timeline<-function(data,colstimeline,ylabel, filename_png){
    D0_mean<-colMeans(subset(data,data$sample=="D0")[5:14])
    D0_sd<-apply(subset(data,data$sample=="D0")[5:14],2,sd)
    D2_mean<-colMeans(subset(data,data$sample=="D2")[5:14])
    D2_sd<-apply(subset(data,data$sample=="D2")[5:14],2,sd)
    D4_mean<-colMeans(subset(data,data$sample=="D4")[5:14])
    D4_sd<-apply(subset(data,data$sample=="D4")[5:14],2,sd)
    D6_mean<-colMeans(subset(data,data$sample=="D6")[5:14])
    D6_sd<-apply(subset(data,data$sample=="D6")[5:14],2,sd)
    D8_mean<-colMeans(subset(data,data$sample=="D8")[5:14])
    D8_sd<-apply(subset(data,data$sample=="D8")[5:14],2,sd)
    
    png(filename=filename_png, units="in",width=2, height=2 , pointsize=6, res=1200)
    par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2,bg='black',fg='white',col.lab='white',col.axis='white')
    plot(D0_mean, type="l",col=cols[1],xlim=c(0,10),ylim=c(0,1),las=1,ylab=ylabel,xlab="Radial position")
    lines(D2_mean, type="l",col=cols[2],xlim=c(0,10),ylim=c(0,1),las=1,ylab=ylabel,xlab="Radial position")
    lines(D4_mean, type="l",col=cols[3],xlim=c(0,10),ylim=c(0,1),las=1,ylab=ylabel,xlab="Radial position")
    lines(D6_mean, type="l",col=cols[4],xlim=c(0,10),ylim=c(0,1),las=1,ylab=ylabel,xlab="Radial position")
    lines(D8_mean, type="l",col=cols[5],xlim=c(0,10),ylim=c(0,1),las=1,ylab=ylabel,xlab="Radial position")
    segments(1:10, D0_mean - D0_sd, 1:10, D0_mean + D0_sd,col=cols[1]);
    segments(1:10, D2_mean - D2_sd, 1:10, D2_mean + D2_sd,col=cols[2]);
    segments(1:10, D4_mean - D4_sd, 1:10, D4_mean + D4_sd,col=cols[3]);
    segments(1:10, D6_mean - D6_sd, 1:10, D6_mean + D6_sd,col=cols[4]);
    segments(1:10, D8_mean - D8_sd, 1:10, D8_mean + D8_sd,col=cols[5]);
    dev.off()
    
    
  }
  angles_plot_black_bg_timeline <- function(path_to_exn, samplen,data_type,filename_png,dire){
    dir=paste(path_to_exn,samplen[1],"/angles_spheroid/",data_type,sep="")
    setwd(dir)
    file_data<-list.files()
    D0_mean_angles<-as.data.frame(matrix(ncol=length(file_data),nrow=91))
    for (j in 1:length(file_data)){
      temp<-read.table(file_data[j], header=TRUE, sep=",",stringsAsFactors=FALSE)
      D0_mean_angles[,j]<-temp[,2]/max(temp[,2])
      rm(temp)
    }
    dir=paste(path_to_exn,samplen[2],"/angles_spheroid/",data_type,sep="")
    setwd(dir)
    file_data<-list.files()
    D2_mean_angles<-as.data.frame(matrix(ncol=length(file_data),nrow=91))
    for (j in 1:length(file_data)){
      temp<-read.table(file_data[j], header=TRUE, sep=",",stringsAsFactors=FALSE)
      D2_mean_angles[,j]<-temp[,2]/max(temp[,2])
      rm(temp)
    }
    dir=paste(path_to_exn,samplen[3],"/angles_spheroid/",data_type,sep="")
    setwd(dir)
    file_data<-list.files()
    D4_mean_angles<-as.data.frame(matrix(ncol=length(file_data),nrow=91))
    for (j in 1:length(file_data)){
      temp<-read.table(file_data[j], header=TRUE, sep=",",stringsAsFactors=FALSE)
      D4_mean_angles[,j]<-temp[,2]/max(temp[,2])
      rm(temp)
    }
    dir=paste(path_to_exn,samplen[4],"/angles_spheroid/",data_type,sep="")
    setwd(dir)
    file_data<-list.files()
    D6_mean_angles<-as.data.frame(matrix(ncol=length(file_data),nrow=91))
    for (j in 1:length(file_data)){
      temp<-read.table(file_data[j], header=TRUE, sep=",",stringsAsFactors=FALSE)
      D6_mean_angles[,j]<-temp[,2]/max(temp[,2])
      rm(temp)
    }
    dir=paste(path_to_exn,samplen[5],"/angles_spheroid/",data_type,sep="")
    setwd(dir)
    file_data<-list.files()
    D8_mean_angles<-as.data.frame(matrix(ncol=length(file_data),nrow=91))
    for (j in 1:length(file_data)){
      temp<-read.table(file_data[j], header=TRUE, sep=",",stringsAsFactors=FALSE)
      D8_mean_angles[,j]<-temp[,2]/max(temp[,2])
      rm(temp)
    }
    
    
    D0_mean<-rowMeans(D0_mean_angles)
    D0_sd<-apply(as.matrix(D0_mean_angles),1,sd)
    D2_mean<-rowMeans(D2_mean_angles)
    D2_sd<-apply(as.matrix(D2_mean_angles),1,sd)
    D4_mean<-rowMeans(D4_mean_angles)
    D4_sd<-apply(as.matrix(D4_mean_angles),1,sd)
    D6_mean<-rowMeans(D6_mean_angles)
    D6_sd<-apply(as.matrix(D6_mean_angles),1,sd)
    D8_mean<-rowMeans(D8_mean_angles)
    D8_sd<-apply(as.matrix(D8_mean_angles),1,sd)
    
    cols<-colorRampPalette(c("red","yellow"))( 5) 
    setwd(dire)
    
    png(filename=filename_png, units="in",width=2, height=2 , pointsize=6, res=1200)
    par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2,bg='black',fg='white',col.lab='white',col.axis='white')
    plot(NA,xlim=c(0,90),ylim=c(0,1),ylab="Norm. Frequency",xlab="Angle",las=1)
    lines(D0_mean,col=cols[1])
    lines(D2_mean,col=cols[1])
    lines(D4_mean,col=cols[1])
    lines(D6_mean,col=cols[1])
    lines(D8_mean,col=cols[1])
    segments(1:91,D0_mean+D0_sd,1:91,D0_mean-D0_sd, col=cols[1])
    segments(1:91,D2_mean+D2_sd,1:91,D2_mean-D2_sd, col=cols[2])
    segments(1:91,D4_mean+D4_sd,1:91,D4_mean-D4_sd, col=cols[3])
    segments(1:91,D6_mean+D6_sd,1:91,D6_mean-D6_sd, col=cols[4])
    segments(1:91,D8_mean+D8_sd,1:91,D8_mean-D8_sd, col=cols[5])
    dev.off()
    
  }
  
  plot_white_bg_histogram_timeline<-function(data,colstimeline,parameters,labelling, prefix){
    
    for( i in 1:length(parameters)){
      index<-which(colnames(data) %in% parameters[i])
      filename_png<-paste(prefix,parameters[i],"hist.png",sep="_")
      png(filename=filename_png, units="in",width=2, height=2 , pointsize=6, res=1200)
      par(font.axis = 2,font.lab=2,mar=c(4,4.5,1,1), font=2,bg='white',fg='black',col.lab='black',col.axis='black')
      a<-subset(data, data$sample=="D0")[,index]
      b<-subset(data, data$sample=="D2")[,index]
      c<-subset(data, data$sample=="D4")[,index]
      d<-subset(data, data$sample=="D6")[,index]
      e<-subset(data, data$sample=="D8")[,index]
      if(length(a)>0) d1<-density(a) else d1<-density(c(0,0,0,0,0,0,0,0,0))
      if(length(b)>0) d2<-density(b) else d2<-density(c(0,0,0,0,0,0,0,0,0))
      if(length(c)>0) d3<-density(c) else d3<-density(c(0,0,0,0,0,0,0,0,0))
      if(length(d)>0) d4<-density(d) else d4<-density(c(0,0,0,0,0,0,0,0,0))
      if(length(e)>0) d5<-density(e) else d5<-density(c(0,0,0,0,0,0,0,0,0))
      
      ymax<-max(c(d1$y,d2$y,d3$y,d4$y,d5$y))*1.2
      xmax<-max(c(d1$x,d2$x,d3$x,d4$x,d5$x))*1.1
      xmin<-min(c(d1$x,d2$x,d3$x,d4$x,d5$x))*0.9
      plot(d3,main="",col=NA,lwd=2,xlab=labelling[i],ylim=c(0,ymax),las=1,cex.axis=0.8, xlim=c(xmin,xmax), ylab="Probability Density")
      lines(d1,col=colstimeline[1],lwd=2)
      lines(d2,col=colstimeline[2],lwd=2)
      lines(d3,col=colstimeline[3],lwd=2)
      lines(d4,col=colstimeline[4],lwd=2)
      lines(d5,col=colstimeline[5],lwd=2)
      legend('topright',legend=c("D0","D2","D4","D6","D8"),fill=colstimeline, horiz=TRUE,cex=0.8)
      dev.off()
    }
    
    
    
  }
  plot_white_bg_boxplot_timeline<-function(data,colstimeline,parameters,labelling, prefix){
    for( i in 1:length(parameters)){
      index<-which(colnames(data) %in% parameters[i])
      filename_png<-paste(prefix,parameters[i],"box.png",sep="_")
      png(filename=filename_png, units="in",width=2, height=2 , pointsize=6, res=1200)
      par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2,bg='white',fg='black',col.lab='black',col.axis='black')
      a<-subset(data, data$sample=="D0")[,index]
      b<-subset(data, data$sample=="D2")[,index]
      c<-subset(data, data$sample=="D4")[,index]
      d<-subset(data, data$sample=="D6")[,index]
      e<-subset(data, data$sample=="D8")[,index]
      xmax<-max(data[,index])*1.2
      xmin<-min(data[,index])*0.9
      boxplot(a,b,c,d,e,main="",col=colstimeline,xlab=labelling[i],las=1,cex.axis=1,horizontal = T, 
              names=c("D0","D2","D4","D6","D8"),lty=1,lwd=0.5, pch=18)
      dev.off()
    }
    
  }
  plot_white_bg_barplot_timeline<-function(data,colstimeline,parameters,labelling, prefix){
    
    for( i in 1:length(parameters)){
      index<-which(colnames(data) %in% parameters[i])
      filename_png<-paste(prefix,parameters[i],"bar.png",sep="_")
      
      T1<-subset(data,data$trial=="T1")
      T1_means<-c(mean(subset(T1, T1$sample=="D0")[,index], na.rm = T),mean(subset(T1, T1$sample=="D2")[,index], na.rm = T),
                  mean(subset(T1, T1$sample=="D4")[,index], na.rm = T),mean(subset(T1, T1$sample=="D6")[,index], na.rm = T),
                  mean(subset(T1, T1$sample=="D8")[,index], na.rm = T))
      
      T1_means[is.nan(T1_means)] <- NA
      T1_means<-as.vector(na.omit(T1_means))
      T1_means<-T1_means/T1_means[1]
      rm(T1)
      
      T1<-subset(data,data$trial=="T2")
      T2_means<-c(mean(subset(T1, T1$sample=="D0")[,index], na.rm = T),mean(subset(T1, T1$sample=="D2")[,index], na.rm = T),
                  mean(subset(T1, T1$sample=="D4")[,index], na.rm = T),mean(subset(T1, T1$sample=="D6")[,index], na.rm = T),
                  mean(subset(T1, T1$sample=="D8")[,index], na.rm = T))
      T2_means[is.nan(T2_means)] <- NA
      T2_means<-as.vector(na.omit(T2_means))
      T2_means<-T2_means/T2_means[1]
      rm(T1)
      T1<-subset(data,data$trial=="T3")
      T3_means<-c(mean(subset(T1, T1$sample=="D0")[,index], na.rm = T),mean(subset(T1, T1$sample=="D2")[,index], na.rm = T),
                  mean(subset(T1, T1$sample=="D4")[,index], na.rm = T),mean(subset(T1, T1$sample=="D6")[,index], na.rm = T),
                  mean(subset(T1, T1$sample=="D8")[,index], na.rm = T))
      T3_means[is.nan(T3_means)] <- NA
      T3_means<-as.vector(na.omit(T3_means))
      T3_means<-T3_means/T3_means[1]
      rm(T1)
      
      if(length(T1_means)==5 & length(T2_means)==5 & length(T3_means)==5){
        pop_means<-c(1,mean(c(T1_means[2],T2_means[2],T3_means[2])),
                     mean(c(T1_means[3],T2_means[3],T3_means[3])),mean(c(T1_means[4],T2_means[4],T3_means[4])),
                     mean(c(T1_means[5],T2_means[5],T3_means[5])))
        SE_means<-c(0,sd(c(T1_means[2],T2_means[2],T3_means[2])),
                    sd(c(T1_means[3],T2_means[3],T3_means[3])),sd(c(T1_means[4],T2_means[4],T3_means[4])),
                    sd(c(T1_means[5],T2_means[5],T3_means[5])))/sqrt(3)
        
        png(filename=filename_png, units="in",width=2, height=2 , pointsize=6, res=1200)
        par(font.axis = 2,font.lab=2,mar=c(4,6,1,1), font=2,bg='white',fg='black',col.lab='black',col.axis='black')
        ymax=max(pop_means+SE_means)*1.2
        ymin=min(pop_means-SE_means)*0.8
        foo<-barplot(pop_means,las=1, col=colstimeline,ylim=c(ymin,ymax),names.arg = c("D0","D2","D4","D6","D8"),ylab=labelling[i],xpd = FALSE)
        arrows(x0=foo,y0=pop_means+SE_means,y1=pop_means-SE_means,angle=90,code=3,length=0.1)
        box()
        dev.off()
      }
      else{
        fill_num<-5-length(T1_means)
        temp<-rep(0,fill_num)
        temp[(fill_num+1):5]<-T1_means
        T1_means<-temp
        fill_num<-5-length(T2_means)
        temp<-rep(0,fill_num)
        temp[(fill_num+1):5]<-T2_means
        T2_means<-temp
        fill_num<-5-length(T3_means)
        temp<-rep(0,fill_num)
        temp[(fill_num+1):5]<-T3_means
        T3_means<-temp
        rm(fill_num,temp)
        
        pop_means<-c(mean(c(T1_means[1],T2_means[1],T3_means[1])),mean(c(T1_means[2],T2_means[2],T3_means[2])),
                     mean(c(T1_means[3],T2_means[3],T3_means[3])),mean(c(T1_means[4],T2_means[4],T3_means[4])),
                     mean(c(T1_means[5],T2_means[5],T3_means[5])))
        SE_means<-c(sd(c(T1_means[1],T2_means[1],T3_means[1])),sd(c(T1_means[2],T2_means[2],T3_means[2])),
                    sd(c(T1_means[3],T2_means[3],T3_means[3])),sd(c(T1_means[4],T2_means[4],T3_means[4])),
                    sd(c(T1_means[5],T2_means[5],T3_means[5])))/sqrt(3)
        
        png(filename=filename_png, units="in",width=2, height=2 , pointsize=6, res=1200)
        par(font.axis = 2,font.lab=2,mar=c(4,6,1,1), font=2,bg='white',fg='black',col.lab='black',col.axis='black')
        ymax=max(pop_means+SE_means)*1.2
        ymin=min(pop_means-SE_means)*0.8
        foo<-barplot(pop_means,las=1, col=colstimeline,ylim=c(ymin,ymax),names.arg = c("D0","D2","D4","D6","D8"),ylab=labelling[i],xpd = FALSE)
        arrows(x0=foo,y0=pop_means+SE_means,y1=pop_means-SE_means,angle=90,code=3,length=0.1)
        box()
        dev.off()
      }
      
      
      
      
    }
  }
  axial_plot_white_bg_boxplot_timeline<-function(data,colstimeline,ylabel, filename_png){
    D0_mean<-colMeans(subset(data,data$sample=="D0")[1:45])
    D0_sd<-apply(subset(data,data$sample=="D0")[1:45],2,sd)
    D2_mean<-colMeans(subset(data,data$sample=="D2")[1:45])
    D2_sd<-apply(subset(data,data$sample=="D2")[1:45],2,sd)
    D4_mean<-colMeans(subset(data,data$sample=="D4")[1:45])
    D4_sd<-apply(subset(data,data$sample=="D4")[1:45],2,sd)
    D6_mean<-colMeans(subset(data,data$sample=="D6")[1:45])
    D6_sd<-apply(subset(data,data$sample=="D6")[1:45],2,sd)
    D8_mean<-colMeans(subset(data,data$sample=="D8")[1:45])
    D8_sd<-apply(subset(data,data$sample=="D8")[1:45],2,sd)
    
    png(filename=filename_png, units="in",width=2, height=2 , pointsize=6, res=1200)
    par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2,bg='white',fg='black',col.lab='black',col.axis='black')
    plot(D0_mean, type="l",col=cols[1],xlim=c(0,30),ylim=c(0,40),las=1,ylab=ylabel,xlab="Zposition")
    lines(D2_mean, type="l",col=cols[2],xlim=c(0,30),ylim=c(0,40),las=1,ylab=ylabel,xlab="Zposition")
    lines(D4_mean, type="l",col=cols[3],xlim=c(0,30),ylim=c(0,40),las=1,ylab=ylabel,xlab="Zposition")
    lines(D6_mean, type="l",col=cols[4],xlim=c(0,30),ylim=c(0,40),las=1,ylab=ylabel,xlab="Zposition")
    lines(D8_mean, type="l",col=cols[5],xlim=c(0,30),ylim=c(0,40),las=1,ylab=ylabel,xlab="Zposition")
    segments(1:45, D0_mean - D0_sd, 1:45, D0_mean + D0_sd,col=cols[1]);
    segments(1:45, D2_mean - D2_sd, 1:45, D2_mean + D2_sd,col=cols[2]);
    segments(1:45, D4_mean - D4_sd, 1:45, D4_mean + D4_sd,col=cols[3]);
    segments(1:45, D6_mean - D6_sd, 1:45, D6_mean + D6_sd,col=cols[4]);
    segments(1:45, D8_mean - D8_sd, 1:45, D8_mean + D8_sd,col=cols[5]);
    dev.off()
    
    
  }
  radial_plot_white_bg_boxplot_timeline<-function(data,colstimeline,ylabel, filename_png){
    D0_mean<-colMeans(subset(data,data$sample=="D0")[5:14])
    D0_sd<-apply(subset(data,data$sample=="D0")[5:14],2,sd)
    D2_mean<-colMeans(subset(data,data$sample=="D2")[5:14])
    D2_sd<-apply(subset(data,data$sample=="D2")[5:14],2,sd)
    D4_mean<-colMeans(subset(data,data$sample=="D4")[5:14])
    D4_sd<-apply(subset(data,data$sample=="D4")[5:14],2,sd)
    D6_mean<-colMeans(subset(data,data$sample=="D6")[5:14])
    D6_sd<-apply(subset(data,data$sample=="D6")[5:14],2,sd)
    D8_mean<-colMeans(subset(data,data$sample=="D8")[5:14])
    D8_sd<-apply(subset(data,data$sample=="D8")[5:14],2,sd)
    
    png(filename=filename_png, units="in",width=2, height=2 , pointsize=6, res=1200)
    par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2,bg='white',fg='black',col.lab='black',col.axis='black')
    plot(D0_mean, type="l",col=cols[1],xlim=c(0,10),ylim=c(0,1),las=1,ylab=ylabel,xlab="Radial position")
    lines(D2_mean, type="l",col=cols[2],xlim=c(0,10),ylim=c(0,1),las=1,ylab=ylabel,xlab="Radial position")
    lines(D4_mean, type="l",col=cols[3],xlim=c(0,10),ylim=c(0,1),las=1,ylab=ylabel,xlab="Radial position")
    lines(D6_mean, type="l",col=cols[4],xlim=c(0,10),ylim=c(0,1),las=1,ylab=ylabel,xlab="Radial position")
    lines(D8_mean, type="l",col=cols[5],xlim=c(0,10),ylim=c(0,1),las=1,ylab=ylabel,xlab="Radial position")
    segments(1:10, D0_mean - D0_sd, 1:10, D0_mean + D0_sd,col=cols[1]);
    segments(1:10, D2_mean - D2_sd, 1:10, D2_mean + D2_sd,col=cols[2]);
    segments(1:10, D4_mean - D4_sd, 1:10, D4_mean + D4_sd,col=cols[3]);
    segments(1:10, D6_mean - D6_sd, 1:10, D6_mean + D6_sd,col=cols[4]);
    segments(1:10, D8_mean - D8_sd, 1:10, D8_mean + D8_sd,col=cols[5]);
    dev.off()
    
    
  }
  angles_plot_white_bg_timeline <- function(path_to_exn, samplen,data_type,filename_png,dire){
    dir=paste(path_to_exn,samplen[1],"/angles_spheroid/",data_type,sep="")
    setwd(dir)
    file_data<-list.files()
    D0_mean_angles<-as.data.frame(matrix(ncol=length(file_data),nrow=91))
    for (j in 1:length(file_data)){
      temp<-read.table(file_data[j], header=TRUE, sep=",",stringsAsFactors=FALSE)
      D0_mean_angles[,j]<-temp[,2]/max(temp[,2])
      rm(temp)
    }
    dir=paste(path_to_exn,samplen[2],"/angles_spheroid/",data_type,sep="")
    setwd(dir)
    file_data<-list.files()
    D2_mean_angles<-as.data.frame(matrix(ncol=length(file_data),nrow=91))
    for (j in 1:length(file_data)){
      temp<-read.table(file_data[j], header=TRUE, sep=",",stringsAsFactors=FALSE)
      D2_mean_angles[,j]<-temp[,2]/max(temp[,2])
      rm(temp)
    }
    dir=paste(path_to_exn,samplen[3],"/angles_spheroid/",data_type,sep="")
    setwd(dir)
    file_data<-list.files()
    D4_mean_angles<-as.data.frame(matrix(ncol=length(file_data),nrow=91))
    for (j in 1:length(file_data)){
      temp<-read.table(file_data[j], header=TRUE, sep=",",stringsAsFactors=FALSE)
      D4_mean_angles[,j]<-temp[,2]/max(temp[,2])
      rm(temp)
    }
    dir=paste(path_to_exn,samplen[4],"/angles_spheroid/",data_type,sep="")
    setwd(dir)
    file_data<-list.files()
    D6_mean_angles<-as.data.frame(matrix(ncol=length(file_data),nrow=91))
    for (j in 1:length(file_data)){
      temp<-read.table(file_data[j], header=TRUE, sep=",",stringsAsFactors=FALSE)
      D6_mean_angles[,j]<-temp[,2]/max(temp[,2])
      rm(temp)
    }
    dir=paste(path_to_exn,samplen[5],"/angles_spheroid/",data_type,sep="")
    setwd(dir)
    file_data<-list.files()
    D8_mean_angles<-as.data.frame(matrix(ncol=length(file_data),nrow=91))
    for (j in 1:length(file_data)){
      temp<-read.table(file_data[j], header=TRUE, sep=",",stringsAsFactors=FALSE)
      D8_mean_angles[,j]<-temp[,2]/max(temp[,2])
      rm(temp)
    }
    
    
    D0_mean<-rowMeans(D0_mean_angles)
    D0_sd<-apply(as.matrix(D0_mean_angles),1,sd)
    D2_mean<-rowMeans(D2_mean_angles)
    D2_sd<-apply(as.matrix(D2_mean_angles),1,sd)
    D4_mean<-rowMeans(D4_mean_angles)
    D4_sd<-apply(as.matrix(D4_mean_angles),1,sd)
    D6_mean<-rowMeans(D6_mean_angles)
    D6_sd<-apply(as.matrix(D6_mean_angles),1,sd)
    D8_mean<-rowMeans(D8_mean_angles)
    D8_sd<-apply(as.matrix(D8_mean_angles),1,sd)
    
    cols<-colorRampPalette(c("red","yellow"))( 5) 
    setwd(dire)
    
    png(filename=filename_png, units="in",width=2, height=2 , pointsize=6, res=1200)
    par(font.axis = 2,font.lab=2,mar=c(4,4,1,1), font=2,bg='white',fg='black',col.lab='black',col.axis='black')
    plot(NA,xlim=c(0,90),ylim=c(0,1),ylab="Norm. Frequency",xlab="Angle",las=1)
    lines(D0_mean,col=cols[1])
    lines(D2_mean,col=cols[1])
    lines(D4_mean,col=cols[1])
    lines(D6_mean,col=cols[1])
    lines(D8_mean,col=cols[1])
    segments(1:91,D0_mean+D0_sd,1:91,D0_mean-D0_sd, col=cols[1])
    segments(1:91,D2_mean+D2_sd,1:91,D2_mean-D2_sd, col=cols[2])
    segments(1:91,D4_mean+D4_sd,1:91,D4_mean-D4_sd, col=cols[3])
    segments(1:91,D6_mean+D6_sd,1:91,D6_mean-D6_sd, col=cols[4])
    segments(1:91,D8_mean+D8_sd,1:91,D8_mean-D8_sd, col=cols[5])
    dev.off()
    
  }
  
}
# initialise the experiment details and the subdirectories required
{
  path_to_ex<-c("E:/HMF3A_reprogramming/Timeline_H3K4ME3/20200110_hmf3a_reprogramming_H3k4me3_488/","none")
  
  exp<-c("20200110_hmf3a_reprogramming_H3k4me3_488","none")
  
  data_types<-c("/3D geometrical data spheroid/",
                "/3D int_data spheroid/DNA/","/3D int_data spheroid/H3K4ME3/")
  data_types_2d<-c("/2D_measures_spheroid/2D_spheroid_DNA.csv","/2D_measures_spheroid/2D_spheroid_H3K4ME3.csv")
  nuc_data_types<-c("/3D geometrical data/","/3D ellipsoid/","/3D geometerical_simple/","/3D shape measure/",
                    "/3D int_data/DNA/","/3D int_data/H3K4ME3/")
  nuc_data_types_2d<-c("/2D_measures_nuclei/2D_nucleus_DNA.csv","/2D_measures_nuclei/compaction.csv","/2D_measures_nuclei/edf.CSV",
                       "/2D_measures_nuclei/gclm_5_deg0.csv","/2D_measures_nuclei/gclm_5_deg90.csv",
                       "/2D_measures_nuclei/gclm_5_deg180.csv","/2D_measures_nuclei/gclm_5_deg270.csv",
                       "/2D_measures_nuclei/gclm_15_deg0.csv","/2D_measures_nuclei/gclm_15_deg90.csv",
                       "/2D_measures_nuclei/gclm_15_deg180.csv","/2D_measures_nuclei/gclm_15_deg270.csv",
                       "/2D_measures_nuclei/gclm_25_deg0.csv","/2D_measures_nuclei/gclm_25_deg90.csv",
                       "/2D_measures_nuclei/gclm_25_deg180.csv","/2D_measures_nuclei/gclm_25_deg270.csv",
                       "/2D_measures_nuclei/gclm_35_deg0.csv","/2D_measures_nuclei/gclm_35_deg90.csv",
                       "/2D_measures_nuclei/gclm_35_deg180.csv","/2D_measures_nuclei/gclm_35_deg270.csv",
                       "/2D_measures_nuclei/gclm_45_deg0.csv","/2D_measures_nuclei/gclm_45_deg90.csv",
                       "/2D_measures_nuclei/gclm_45_deg180.csv","/2D_measures_nuclei/gclm_45_deg270.csv",
                       "/2D_measures_nuclei/gclm_100_deg0.csv","/2D_measures_nuclei/gclm_100_deg90.csv",
                       "/2D_measures_nuclei/gclm_100_deg180.csv","/2D_measures_nuclei/gclm_100_deg270.csv",
                       "/2D_measures_nuclei/2D_nucleusH3K4ME3.csv")
  samples<-c("D0","D2","D4","D6","D8")
  trial<-c("T1","T2","T3")
  
  dird_apo<-"E:/HMF3A_reprogramming/R_analysis/H3K4ME3/"
  dir.create(dird_apo)
  dird_H3K4ME3<-paste(dird_apo,"spheroid_shells_H3K4ME3", sep="")
  dird_dna<-paste(dird_apo,"spheroid_shells_dna", sep="")
  dir.create(dird_H3K4ME3)
  dir.create(dird_dna)
}

# read in all the required features for spheroids
{
  plot(0:1,0:1)
  #T1
  {
    geometrical_data<-combine_sample_sets(path_to_ex[1],samples,data_types[1],"tsv",exp[1],trial[1])
    T1DNA_int_data<-combine_sample_sets(path_to_ex[1],samples,data_types[2],"tsv",exp[1],trial[1])
    T1H3K4ME3_int_data<-combine_sample_sets(path_to_ex[1],samples,data_types[3],"tsv",exp[1],trial[1])
    geo_int_2D_data<-combine_sample_sets_2d(path_to_ex[1],samples,data_types_2d[1],exp[1],trial[1])
    
    T1_2d_data_H3K4ME3<-combine_sample_sets_2d(path_to_ex[1],samples,data_types_2d[2],exp[1],trial[1])
    
    shell_H3K4ME3_T1<-shell_analsysis(path_to_ex[1],samples,"/Shells/H3K4ME3",exp[1],trial[1],dird_H3K4ME3)
    shell_dna_T1<-shell_analsysis(path_to_ex[1],samples,"/Shells/DNA",exp[1],trial[1],dird_dna)
    
    T1_dna_axial<-axial_analysis(path_to_ex[1],samples,exp[1],trial[1],"/DNA/")
    T1_H3K4ME3_axial<-axial_analysis(path_to_ex[1],samples,exp[1],trial[1],"/H3K4ME3/")
    
    T1_H3K4ME3_angle<-angle_analsysis(path_to_ex[1],samples,"/H3K4ME3",exp[1],trial[1])
    T1_dna_angle<-angle_analsysis(path_to_ex[1],samples,"/DNA",exp[1],trial[1])
    
    geo_int_2D_data$Label<-substring(geo_int_2D_data$Label,5, nchar(geo_int_2D_data$Label)-4)
    geometrical_data$Label<-substring(geometrical_data$Label,0, nchar(geometrical_data$Label)-1)
    
    T1combined<-merge(geometrical_data,geo_int_2D_data, by="Label")
    
    rm(geometrical_data,geo_int_2D_data)
    
  }
  
  # #T2
  # {
  #   geometrical_data<-combine_sample_sets(path_to_ex[2],samples,data_types[1],"tsv",exp[2],trial[2])
  #   T2DNA_int_data<-combine_sample_sets(path_to_ex[2],samples,data_types[2],"tsv",exp[2],trial[2])
  #   T2H3K4ME3_int_data<-combine_sample_sets(path_to_ex[2],samples,data_types[3],"tsv",exp[2],trial[2])
  #   geo_int_2D_data<-combine_sample_sets_2d(path_to_ex[2],samples,data_types_2d[1],exp[2],trial[2])
  #   
  #   T2_2d_data_H3K4ME3<-combine_sample_sets_2d(path_to_ex[2],samples,data_types_2d[2],exp[2],trial[2])
  #   
  #   shell_H3K4ME3_T2<-shell_analsysis(path_to_ex[2],samples,"/Shells/H3K4ME3",exp[2],trial[2],dird_H3K4ME3)
  #   shell_dna_T2<-shell_analsysis(path_to_ex[2],samples,"/Shells/DNA",exp[2],trial[2],dird_dna)
  #   
  #   T2_dna_axial<-axial_analysis(path_to_ex[2],samples,exp[2],trial[2],"/DNA/")
  #   T2_H3K4ME3_axial<-axial_analysis(path_to_ex[2],samples,exp[2],trial[2],"/H3K4ME3/")
  #   
  #   T2_H3K4ME3_angle<-angle_analsysis(path_to_ex[2],samples,"/H3K4ME3",exp[2],trial[2])
  #   T2_dna_angle<-angle_analsysis(path_to_ex[2],samples,"/DNA",exp[2],trial[2])
  #   
  #   geo_int_2D_data$Label<-substring(geo_int_2D_data$Label,5, nchar(geo_int_2D_data$Label)-4)
  #   geometrical_data$Label<-substring(geometrical_data$Label,0, nchar(geometrical_data$Label)-1)
  #   T2combined<-merge(geometrical_data,geo_int_2D_data, by="Label")
  #   
  #   rm(geometrical_data,geo_int_2D_data)
  #   
  # }
  # 
  # #T3
  # {
  #   geometrical_data<-combine_sample_sets(path_to_ex[3],samples,data_types[1],"tsv",exp[3],trial[3])
  #   T3DNA_int_data<-combine_sample_sets(path_to_ex[3],samples,data_types[2],"tsv",exp[3],trial[3])
  #   T3H3K4ME3_int_data<-combine_sample_sets(path_to_ex[3],samples,data_types[3],"tsv",exp[3],trial[3])
  #   geo_int_2D_data<-combine_sample_sets_2d(path_to_ex[3],samples,data_types_2d[1],exp[3],trial[3])
  #   
  #   T3_2d_data_H3K4ME3<-combine_sample_sets_2d(path_to_ex[3],samples,data_types_2d[2],exp[3],trial[3])
  #   
  #   shell_H3K4ME3_T3<-shell_analsysis(path_to_ex[3],samples,"/Shells/H3K4ME3",exp[3],trial[3],dird_H3K4ME3)
  #   shell_dna_T3<-shell_analsysis(path_to_ex[3],samples,"/Shells/DNA",exp[3],trial[3],dird_dna)
  #   
  #   T3_dna_axial<-axial_analysis(path_to_ex[3],samples,exp[3],trial[3],"/DNA/")
  #   T3_H3K4ME3_axial<-axial_analysis(path_to_ex[3],samples,exp[3],trial[3],"/H3K4ME3/")
  #   
  #   T3_H3K4ME3_angle<-angle_analsysis(path_to_ex[3],samples,"/H3K4ME3",exp[3],trial[3])
  #   T3_dna_angle<-angle_analsysis(path_to_ex[3],samples,"/DNA",exp[3],trial[3])
  #   
  #   geo_int_2D_data$Label<-substring(geo_int_2D_data$Label,5, nchar(geo_int_2D_data$Label)-4)
  #   geometrical_data$Label<-substring(geometrical_data$Label,0, nchar(geometrical_data$Label)-1)
  #   T3combined<-merge(geometrical_data,geo_int_2D_data, by="Label")
  #   
  #   rm(geometrical_data,geo_int_2D_data)
  #   
  # }
  # 
  # spheroid_2dint_data_H3K4ME3<-rbind(T1_2d_data_H3K4ME3,T2_2d_data_H3K4ME3,T3_2d_data_H3K4ME3)
  # 
  # spheroid_int_data_H3K4ME3<-rbind(T1H3K4ME3_int_data,T2H3K4ME3_int_data,T3H3K4ME3_int_data)
  # spheroid_int_data_dna<-rbind(T1DNA_int_data,T2DNA_int_data,T3DNA_int_data)
  # 
  # spheroid_radial_H3K4ME3<-rbind(shell_H3K4ME3_T1,shell_H3K4ME3_T2,shell_H3K4ME3_T3)
  # spheroid_radial_dna<-rbind(shell_dna_T1,shell_dna_T2,shell_dna_T3)
  # 
  # spheroid_axial_dna<-rbind(T1_dna_axial,T2_dna_axial,T3_dna_axial)
  # spheroid_axial_H3K4ME3<-rbind(T1_H3K4ME3_axial,T2_H3K4ME3_axial,T3_H3K4ME3_axial)
  # 
  # spheroid_angle_H3K4ME3<-rbind(T1_H3K4ME3_angle,T2_H3K4ME3_angle,T3_H3K4ME3_angle)
  # spheroid_angle_dna<-rbind(T1_dna_angle,T2_dna_angle,T3_dna_angle)
  # 
  # spheroid_int_data<-list( H3K4ME3=spheroid_int_data_H3K4ME3,
  #                          dna=spheroid_int_data_dna)
  # spheroid_radial<-list(H3K4ME3=spheroid_radial_H3K4ME3,
  #                       dna=spheroid_radial_dna)
  # spheroid_axial<-list( H3K4ME3=spheroid_axial_H3K4ME3,
  #                       dna=spheroid_axial_dna)
  # spheroid_angle<-list(H3K4ME3=spheroid_angle_H3K4ME3,
  #                      dna=spheroid_angle_dna)
  # spheroid_2dint_data<-list( H3K4ME3=spheroid_2dint_data_H3K4ME3)
  # 
  # 
  # spheroid_combined<-rbind(T1combined,T2combined,T3combined)
  # colnames(spheroid_combined)[86:87]<-c("sample","trial")
  # 
  # rm(T1combined,T2combined,T3combined,T1VIMENTIN_int_data,T2VIMENTIN_int_data,T3VIMENTIN_int_data,
  #    T1H3K4ME3_int_data,T2H3K4ME3_int_data,T3H3K4ME3_int_data,T1ECAD_int_data,T2ECAD_int_data,T3ECAD_int_data,
  #    shell_H3K4ME3_T1,shell_H3K4ME3_T2,shell_H3K4ME3_T3,shell_vimentin_T1,shell_vimentin_T2,shell_vimentin_T3,
  #    shell_ecad_T1,shell_ecad_T2,shell_ecad_T3,shell_dna_T1,shell_dna_T2,shell_dna_T3,
  #    T1_vimentin_axial,T2_vimentin_axial,T3_vimentin_axial,T1_dna_axial,T2_dna_axial,T3_dna_axial,
  #    T1_H3K4ME3_axial,T2_H3K4ME3_axial,T3_H3K4ME3_axial,T1_ecad_axial,T2_ecad_axial,T3_ecad_axial,
  #    T1_vimentin_angle,T2_vimentin_angle,T3_vimentin_angle,T1_H3K4ME3_angle,T2_H3K4ME3_angle,T3_H3K4ME3_angle,
  #    T1_ecad_angle,T2_ecad_angle,T3_ecad_angle,T1_dna_angle,T2_dna_angle,T3_dna_angle,
  #    T1DNA_int_data,T2DNA_int_data,T3DNA_int_data,
  #    spheroid_angle_vimentin, spheroid_angle_H3K4ME3,spheroid_angle_ecad,spheroid_angle_dna,
  #    spheroid_axial_vimentin, spheroid_axial_H3K4ME3,spheroid_axial_ecad,spheroid_axial_dna,
  #    spheroid_radial_vimentin, spheroid_radial_H3K4ME3,spheroid_radial_ecad,spheroid_radial_dna,
  #    spheroid_int_data_vimentin, spheroid_int_data_H3K4ME3,spheroid_int_data_ecad,spheroid_int_data_dna,
  #    T1_2d_data_VIMENTIN,T2_2d_data_VIMENTIN,T3_2d_data_VIMENTIN,T1_2d_data_H3K4ME3,T2_2d_data_H3K4ME3,T3_2d_data_H3K4ME3,
  #    T1_2d_data_ECAD,T2_2d_data_ECAD,T3_2d_data_ECAD,spheroid_2dint_data_vimentin, spheroid_2dint_data_H3K4ME3,
  #    spheroid_2dint_data_ecad)
  
}

# read in all the required features for nuclei
{
  #T1
  {
    T1_geometrical_data<-combine_sample_sets(path_to_ex[1],samples,nuc_data_types[1],"tsv",exp[1],trial[1])
    T1_DNA_int_data<-combine_sample_sets(path_to_ex[1],samples,nuc_data_types[5],"tsv",exp[1],trial[1])
    T1_geo_int_2D_data<-combine_sample_sets_2d(path_to_ex[1],samples,nuc_data_types_2d[1],exp[1],trial[1])
    T1_compact_2D_data<-combine_sample_sets_2d(path_to_ex[1],samples,nuc_data_types_2d[2],exp[1],trial[1])
    T1_edf_2D_data<-combine_sample_sets_2d(path_to_ex[1],samples,nuc_data_types_2d[3],exp[1],trial[1])
    T1_gclm<-read_glcm(path_to_ex[1],samples,nuc_data_types_2d,exp[1],trial[1])
    
    T1_H3K4ME3_int_data<-combine_sample_sets(path_to_ex[1],samples,nuc_data_types[6],"tsv",exp[1],trial[1])
    
    T1_2d_H3K4ME3_int_data<-combine_sample_sets_2d(path_to_ex[1],samples,nuc_data_types_2d[28],exp[1],trial[1])
    
    
    colnames(T1_compact_2D_data)[2]<-"Label"
    T1_compact_2D_data$Label<-substring(T1_compact_2D_data$Label,0, nchar(T1_compact_2D_data$Label)-4)
    T1_geo_int_2D_data$Label<-substring(T1_geo_int_2D_data$Label,5, nchar(T1_geo_int_2D_data$Label)-4)
    T1_geometrical_data$Label<-substring(T1_geometrical_data$Label,0, nchar(T1_geometrical_data$Label)-1)
    T1_DNA_int_data$Label<-substring(T1_DNA_int_data$Label,0, nchar(T1_DNA_int_data$Label)-1)
    T1_2d_H3K4ME3_int_data$Label<-substring(T1_2d_H3K4ME3_int_data$Label,5, nchar(T1_2d_H3K4ME3_int_data$Label)-4)
    T1_H3K4ME3_int_data$Label<-substring(T1_H3K4ME3_int_data$Label,0, nchar(T1_H3K4ME3_int_data$Label)-1)
    
    T1_combined<-merge(T1_geometrical_data,T1_geo_int_2D_data[,2:36], by="Label")
    T1_combined<-merge(T1_combined,T1_compact_2D_data[,2:8], by="Label")
    T1_combined<-merge(T1_combined,T1_DNA_int_data[,c(3:5,12:18)], by="Label")
    T1_combined<-merge(T1_combined,T1_edf_2D_data[,1:11], by="Label")
    T1_combined<-cbind(T1_combined,T1_gclm[,-c(1,128:131)])
    
    rm(T1_geometrical_data,T1_geo_int_2D_data,T1_compact_2D_data,T1_edf_2D_data,T1_gclm)
  }
  # #T2
  # {
  #   T2_geometrical_data<-combine_sample_sets(path_to_ex[2],samples,nuc_data_types[1],"tsv",exp[2],trial[2])
  #   T2_DNA_int_data<-combine_sample_sets(path_to_ex[2],samples,nuc_data_types[5],"tsv",exp[2],trial[2])
  #   T2_geo_int_2D_data<-combine_sample_sets_2d(path_to_ex[2],samples,nuc_data_types_2d[1],exp[2],trial[2])
  #   T2_compact_2D_data<-combine_sample_sets_2d(path_to_ex[2],samples,nuc_data_types_2d[2],exp[2],trial[2])
  #   T2_edf_2D_data<-combine_sample_sets_2d(path_to_ex[2],samples,nuc_data_types_2d[3],exp[2],trial[2])
  #   T2_gclm<-read_glcm(path_to_ex[2],samples,nuc_data_types_2d,exp[2],trial[2])
  #   
  #   T2_H3K4ME3_int_data<-combine_sample_sets(path_to_ex[2],samples,nuc_data_types[6],"tsv",exp[2],trial[2])
  #   
  #   T2_2d_H3K4ME3_int_data<-combine_sample_sets_2d(path_to_ex[2],samples,nuc_data_types_2d[28],exp[2],trial[2])
  #   
  #   
  #   colnames(T2_compact_2D_data)[2]<-"Label"
  #   T2_compact_2D_data$Label<-substring(T2_compact_2D_data$Label,0, nchar(T2_compact_2D_data$Label)-4)
  #   T2_geo_int_2D_data$Label<-substring(T2_geo_int_2D_data$Label,5, nchar(T2_geo_int_2D_data$Label)-4)
  #   T2_geometrical_data$Label<-substring(T2_geometrical_data$Label,0, nchar(T2_geometrical_data$Label)-1)
  #   T2_2d_H3K4ME3_int_data$Label<-substring(T2_2d_H3K4ME3_int_data$Label,5, nchar(T2_2d_H3K4ME3_int_data$Label)-4)
  #   T2_DNA_int_data$Label<-substring(T2_DNA_int_data$Label,0, nchar(T2_DNA_int_data$Label)-1)
  #   T2_H3K4ME3_int_data$Label<-substring(T2_H3K4ME3_int_data$Label,0, nchar(T2_H3K4ME3_int_data$Label)-1)
  #   
  #   T2_combined<-merge(T2_geometrical_data,T2_geo_int_2D_data[,2:36], by="Label")
  #   T2_combined<-merge(T2_combined,T2_compact_2D_data[,2:8], by="Label")
  #   T2_combined<-merge(T2_combined,T2_DNA_int_data[,c(3:5,12:18)], by="Label")
  #   T2_combined<-merge(T2_combined,T2_edf_2D_data[,1:11], by="Label")
  #   T2_combined<-cbind(T2_combined,T2_gclm[,-c(1,128:131)])
  #   
  #   rm(T2_geometrical_data,T2_geo_int_2D_data,T2_compact_2D_data,T2_edf_2D_data,T2_gclm)
  # }
  # #T3
  # {
  #   T3_geometrical_data<-combine_sample_sets(path_to_ex[3],samples,nuc_data_types[1],"tsv",exp[3],trial[3])
  #   T3_DNA_int_data<-combine_sample_sets(path_to_ex[3],samples,nuc_data_types[5],"tsv",exp[3],trial[3])
  #   T3_geo_int_2D_data<-combine_sample_sets_2d(path_to_ex[3],samples,nuc_data_types_2d[1],exp[3],trial[3])
  #   T3_compact_2D_data<-combine_sample_sets_2d(path_to_ex[3],samples,nuc_data_types_2d[2],exp[3],trial[3])
  #   T3_edf_2D_data<-combine_sample_sets_2d(path_to_ex[3],samples,nuc_data_types_2d[3],exp[3],trial[3])
  #   T3_gclm<-read_glcm(path_to_ex[3],samples,nuc_data_types_2d,exp[3],trial[3])
  #   
  #   T3_H3K4ME3_int_data<-combine_sample_sets(path_to_ex[3],samples,nuc_data_types[6],"tsv",exp[3],trial[3])
  #   
  #   T3_2d_H3K4ME3_int_data<-combine_sample_sets_2d(path_to_ex[3],samples,nuc_data_types_2d[28],exp[3],trial[3])
  #   
  #   
  #   colnames(T3_compact_2D_data)[2]<-"Label"
  #   T3_compact_2D_data$Label<-substring(T3_compact_2D_data$Label,0, nchar(T3_compact_2D_data$Label)-4)
  #   T3_geo_int_2D_data$Label<-substring(T3_geo_int_2D_data$Label,5, nchar(T3_geo_int_2D_data$Label)-4)
  #   T3_geometrical_data$Label<-substring(T3_geometrical_data$Label,0, nchar(T3_geometrical_data$Label)-1)
  #   T3_2d_H3K4ME3_int_data$Label<-substring(T3_2d_H3K4ME3_int_data$Label,5, nchar(T3_2d_H3K4ME3_int_data$Label)-4)
  #   T3_DNA_int_data$Label<-substring(T3_DNA_int_data$Label,0, nchar(T3_DNA_int_data$Label)-1)
  #   T3_H3K4ME3_int_data$Label<-substring(T3_H3K4ME3_int_data$Label,0, nchar(T3_H3K4ME3_int_data$Label)-1)
  #   
  #   T3_combined<-merge(T3_geometrical_data,T3_geo_int_2D_data[,2:36], by="Label")
  #   T3_combined<-merge(T3_combined,T3_compact_2D_data[,2:8], by="Label")
  #   T3_combined<-merge(T3_combined,T3_DNA_int_data[,c(3:5,12:18)], by="Label")
  #   T3_combined<-merge(T3_combined,T3_edf_2D_data[,1:11], by="Label")
  #   T3_combined<-cbind(T3_combined,T3_gclm[,-c(1,128:131)])
  #   
  #   rm(T3_geometrical_data,T3_geo_int_2D_data,T3_compact_2D_data,T3_edf_2D_data,T3_gclm)
  # }
  # 
  # nuclei_2dint_data_H3K4ME3<-rbind(T1_2d_H3K4ME3_int_data,T2_2d_H3K4ME3_int_data,T3_2d_H3K4ME3_int_data)
  # 
  # nuclei_int_data_H3K4ME3<-rbind(T1_H3K4ME3_int_data,T2_H3K4ME3_int_data,T3_H3K4ME3_int_data)
  # nuclei_int_data_dna<-rbind(T1_DNA_int_data,T2_DNA_int_data,T3_DNA_int_data)
  # 
  # 
  # nuclei_2dint_data<-list(H3K4ME3=nuclei_2dint_data_H3K4ME3)
  # nuclei_3dint_data<-list(H3K4ME3=nuclei_int_data_H3K4ME3,dna=nuclei_int_data_dna)
  # 
  # nucleus_combined<-rbind(T1_combined,T2_combined,T3_combined)
  # 
  # rm(T1_2d_H3K4ME3_int_data,T2_2d_H3K4ME3_int_data,T3_2d_H3K4ME3_int_data,T1_2d_OCT4_int_data,T2_2d_OCT4_int_data,T3_2d_OCT4_int_data,
  #    T1_2d_PMLC_int_data,T2_2d_PMLC_int_data,T3_2d_PMLC_int_data,T1_H3K4ME3_int_data,T2_H3K4ME3_int_data,T3_H3K4ME3_int_data,
  #    T1_OCT4_int_data,T2_OCT4_int_data,T3_OCT4_int_data,T1_PMLC_int_data,T2_PMLC_int_data,T3_PMLC_int_data,
  #    T1_DNA_int_data,T2_DNA_int_data,T3_DNA_int_data,T1_cell_H3K4ME3_int_data,T2_cell_H3K4ME3_int_data,T3_cell_H3K4ME3_int_data,
  #    T1_cell_OCT4_int_data,T2_cell_OCT4_int_data,T3_cell_OCT4_int_data,T1_cell_PMLC_int_data,T2_cell_PMLC_int_data,T3_cell_PMLC_int_data,
  #    T1_cell_DNA_int_data,T2_cell_DNA_int_data,T3_cell_DNA_int_data,T1_combined,T2_combined,T3_combined,
  #    nuclei_2dint_data_H3K4ME3, nuclei_2dint_data_oct4,nuclei_2dint_data_pmlc,nuclei_int_data_H3K4ME3, nuclei_int_data_oct4,
  #    nuclei_int_data_pmlc,nuclei_int_data_dna,cell_int_data_H3K4ME3, cell_int_data_oct4,cell_int_data_pmlc,cell_int_data_dna)
  # 
}
#assigning plotting parameters
{
  
  
  parameters_2d_int_data<-c("Mean","StdDev","Mode","Min","Max","IntDen","Median","Skew","Kurt","RawIntDen")
  labelling_2d_int_data<-paste("Proj.Spheroid Level",parameters_2d_int_data,sep="\n")
  nuc_labelling_2d_int_data<-paste("Proj.Nuclear Level",parameters_2d_int_data,sep="\n")
  
  parameters_int_data<-c("AtCenter","IntDen","Min","Max","Mean","Sigma","Mode","Mode.NonZero")
  labelling_int_data<-paste("Spheroid Level",parameters_int_data,sep="\n")
  nuclabelling_int_data<-paste("Nuclear Level",parameters_int_data,sep="\n")
  celllabelling_int_data<-paste("Cellular Level",parameters_int_data,sep="\n")
  
  parameters_axial_data<-c("axial_index_biased_sp","basal_index_biased_sp","axial_index_biased_pc","basal_index_biased_pc",
                           "axial_index_sp","basal_index_sp","axial_index_pc","basal_index_pc")
  labelling_axial_data<-paste("Distribution",c("axial index bisp","basal index bisp","axial index bipc","basal index bipc",
                                               "axial index sp","basal index sp","axial index pc","basal index pc"),sep="\n")
  
  parameters_radial_data<-c("peripheral_index_sp","central_index_sp","peripheral_index_pc","central_index_pc")
  labelling_radial_data<-paste("Distribution",c("peripheral index sp","central index sp","peripheral index pc","central index pc"),sep="\n")
  
  parameters_angle_data<-c("max_freq_angle","min_freq_angle","meadian_freq_angle","mean_freq_angle","error_sd_freq")
  labelling_angle_data<-paste(c("Most Frequent","Least Frequent","Median Frequency","Average Frequency","Std.Dev. Angles"),"angle",sep="\n")
  
  parameters_geometrical_data<-c( "RatioVolbox","Vol..unit.","Surf..unit." ,"Comp..unit.","Spher..unit.","Feret..unit.","Ell_MajRad","Ell_Elon","Ell_Flatness",
                                  "RatioVolEllipsoid","Moment1","Moment2","Moment3","Moment4","Moment5","DCMin..unit.","DCMax..unit.","DCMean..unit.","DCSD..unit.",
                                  "Area","Perim.","Major","Minor","Circ.","Feret","MinFeret","AR","Round")
  
  labelling_geometrical_data<-c("Bounding Box Ratio","Volume", "Surface Area","Compaction","Sphericity","Feret (3D)", "Ellipsoid R1","Ellipsoid Ellongation",
                                "Ellipsoid Flatness", "Ellipsoid Ratio","Moment1","Moment2","Moment3","Moment4","Moment5", "DC_Min.","DC_Max","DC_Mean.","DC_SD",
                                "Proj. Area", "Proj. Perimeter","Proj. Major Axis", "Proj. Minor Axis", "Proj. Circularity","Proj. Feret","Proj. Min Feret", 
                                "Proj. A.R","Proj. Roundness")
  
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

#plotting
{
  dird<-"E:/HMF3A_reprogramming/R_analysis/H3K4ME3/black_bg_plots/"
  dir.create(dird)
  setwd(dird)
  cols<-colorRampPalette(c("red","yellow"))( 5) 
  
  #Spheoid Proj_Int
  {
    
    dire<-paste(dird,"spheroid_proj_int/",sep="")
    dir.create(dire)
    setwd(dire)
    
    plot_black_bg_histogram_timeline(T1_2d_data_H3K4ME3,cols,parameters_2d_int_data,paste("H3K4ME3 ",labelling_2d_int_data,sep=""),"H3K4ME3_T1_")
    plot_black_bg_boxplot_timeline(T1_2d_data_H3K4ME3,cols,parameters_2d_int_data,paste("H3K4ME3 ",labelling_2d_int_data,sep=""),"H3K4ME3_T1_")
    
    plot_black_bg_histogram_timeline(subset(spheroid_2dint_data$H3K4ME3,spheroid_2dint_data$H3K4ME3$trial=="T1"), cols,parameters_2d_int_data,paste("H3K4ME3 ",labelling_2d_int_data,sep=""),"H3K4ME3_T1_")
    plot_black_bg_histogram_timeline(subset(spheroid_2dint_data$H3K4ME3,spheroid_2dint_data$H3K4ME3$trial=="T2"), cols,parameters_2d_int_data,paste("H3K4ME3 ",labelling_2d_int_data,sep=""),"H3K4ME3_T2_")
    plot_black_bg_histogram_timeline(subset(spheroid_2dint_data$H3K4ME3,spheroid_2dint_data$H3K4ME3$trial=="T3"), cols,parameters_2d_int_data,paste("H3K4ME3 ",labelling_2d_int_data,sep=""),"H3K4ME3_T3_")
    
    plot_black_bg_boxplot_timeline(subset(spheroid_2dint_data$H3K4ME3,spheroid_2dint_data$H3K4ME3$trial=="T1"), cols,parameters_2d_int_data,paste("H3K4ME3 ",labelling_2d_int_data,sep=""),"H3K4ME3_T1_")
    plot_black_bg_boxplot_timeline(subset(spheroid_2dint_data$H3K4ME3,spheroid_2dint_data$H3K4ME3$trial=="T2"), cols,parameters_2d_int_data,paste("H3K4ME3 ",labelling_2d_int_data,sep=""),"H3K4ME3_T2_")
    plot_black_bg_boxplot_timeline(subset(spheroid_2dint_data$H3K4ME3,spheroid_2dint_data$H3K4ME3$trial=="T3"), cols,parameters_2d_int_data,paste("H3K4ME3 ",labelling_2d_int_data,sep=""),"H3K4ME3_T3_")
    
    plot_black_bg_barplot_timeline(spheroid_2dint_data$H3K4ME3, cols,parameters_2d_int_data,paste("H3K4ME3 ",labelling_2d_int_data,sep=""),"H3K4ME3_")
    
    
  }
  #Spheoid 3D_Int
  {
    dire<-paste(dird,"spheroid_3d_int/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_black_bg_histogram_timeline(T1_H3K4ME3_int_data,cols,parameters_2d_int_data,paste("H3K4ME3 ",labelling_2d_int_data,sep=""),"H3K4ME3_T1_")
    plot_black_bg_boxplot_timeline(T1_H3K4ME3_int_data,cols,parameters_2d_int_data,paste("H3K4ME3 ",labelling_2d_int_data,sep=""),"H3K4ME3_T1_")
    
    plot_black_bg_histogram_timeline(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T1_")
    plot_black_bg_histogram_timeline(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T2_")
    plot_black_bg_histogram_timeline(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T3_")
    
    plot_black_bg_boxplot_timeline(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T1_")
    plot_black_bg_boxplot_timeline(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T2_")
    plot_black_bg_boxplot_timeline(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T3_")
    
    plot_black_bg_barplot_timeline(spheroid_int_data$dna, cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_")
    
    plot_black_bg_histogram_timeline(subset(spheroid_int_data$H3K4ME3,spheroid_int_data$H3K4ME3$trial=="T1"), cols,parameters_int_data,paste("H3K4ME3 ",labelling_int_data,sep=""),"H3K4ME3_T1_")
    plot_black_bg_histogram_timeline(subset(spheroid_int_data$H3K4ME3,spheroid_int_data$H3K4ME3$trial=="T2"), cols,parameters_int_data,paste("H3K4ME3 ",labelling_int_data,sep=""),"H3K4ME3_T2_")
    plot_black_bg_histogram_timeline(subset(spheroid_int_data$H3K4ME3,spheroid_int_data$H3K4ME3$trial=="T3"), cols,parameters_int_data,paste("H3K4ME3 ",labelling_int_data,sep=""),"H3K4ME3_T3_")
    
    plot_black_bg_boxplot_timeline(subset(spheroid_int_data$H3K4ME3,spheroid_int_data$H3K4ME3$trial=="T1"), cols,parameters_int_data,paste("H3K4ME3 ",labelling_int_data,sep=""),"H3K4ME3_T1_")
    plot_black_bg_boxplot_timeline(subset(spheroid_int_data$H3K4ME3,spheroid_int_data$H3K4ME3$trial=="T2"), cols,parameters_int_data,paste("H3K4ME3 ",labelling_int_data,sep=""),"H3K4ME3_T2_")
    plot_black_bg_boxplot_timeline(subset(spheroid_int_data$H3K4ME3,spheroid_int_data$H3K4ME3$trial=="T3"), cols,parameters_int_data,paste("H3K4ME3 ",labelling_int_data,sep=""),"H3K4ME3_T3_")
    
    plot_black_bg_barplot_timeline(spheroid_int_data$H3K4ME3, cols,parameters_int_data,paste("H3K4ME3 ",labelling_int_data,sep=""),"H3K4ME3_")
    
    
  }
  #Spheoid angles
  {
    dire<-paste(dird,"spheroid_angles/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_black_bg_histogram_timeline(subset(spheroid_angle$dna,spheroid_angle$dna$trial=="T1"), cols,parameters_angle_data,paste("DNA ",labelling_angle_data,sep=""),"DNA_T1_")
    plot_black_bg_histogram_timeline(subset(spheroid_angle$dna,spheroid_angle$dna$trial=="T2"), cols,parameters_angle_data,paste("DNA ",labelling_angle_data,sep=""),"DNA_T2_")
    plot_black_bg_histogram_timeline(subset(spheroid_angle$dna,spheroid_angle$dna$trial=="T3"), cols,parameters_angle_data,paste("DNA ",labelling_angle_data,sep=""),"DNA_T3_")
    
    plot_black_bg_boxplot_timeline(subset(spheroid_angle$dna,spheroid_angle$dna$trial=="T1"), cols,parameters_angle_data,paste("DNA ",labelling_angle_data,sep=""),"DNA_T1_")
    plot_black_bg_boxplot_timeline(subset(spheroid_angle$dna,spheroid_angle$dna$trial=="T2"), cols,parameters_angle_data,paste("DNA ",labelling_angle_data,sep=""),"DNA_T2_")
    plot_black_bg_boxplot_timeline(subset(spheroid_angle$dna,spheroid_angle$dna$trial=="T3"), cols,parameters_angle_data,paste("DNA ",labelling_angle_data,sep=""),"DNA_T3_")
    
    
    plot_black_bg_histogram_timeline(subset(spheroid_angle$H3K4ME3,spheroid_angle$H3K4ME3$trial=="T1"), cols,parameters_angle_data,paste("H3K4ME3 ",labelling_angle_data,sep=""),"H3K4ME3_T1_")
    plot_black_bg_histogram_timeline(subset(spheroid_angle$H3K4ME3,spheroid_angle$H3K4ME3$trial=="T2"), cols,parameters_angle_data,paste("H3K4ME3 ",labelling_angle_data,sep=""),"H3K4ME3_T2_")
    plot_black_bg_histogram_timeline(subset(spheroid_angle$H3K4ME3,spheroid_angle$H3K4ME3$trial=="T3"), cols,parameters_angle_data,paste("H3K4ME3 ",labelling_angle_data,sep=""),"H3K4ME3_T3_")
    
    plot_black_bg_boxplot_timeline(subset(spheroid_angle$H3K4ME3,spheroid_angle$H3K4ME3$trial=="T1"), cols,parameters_angle_data,paste("H3K4ME3 ",labelling_angle_data,sep=""),"H3K4ME3_T1_")
    plot_black_bg_boxplot_timeline(subset(spheroid_angle$H3K4ME3,spheroid_angle$H3K4ME3$trial=="T2"), cols,parameters_angle_data,paste("H3K4ME3 ",labelling_angle_data,sep=""),"H3K4ME3_T2_")
    plot_black_bg_boxplot_timeline(subset(spheroid_angle$H3K4ME3,spheroid_angle$H3K4ME3$trial=="T3"), cols,parameters_angle_data,paste("H3K4ME3 ",labelling_angle_data,sep=""),"H3K4ME3_T3_")
    
    plot_black_bg_barplot_timeline(spheroid_angle$H3K4ME3, cols,parameters_angle_data,paste("H3K4ME3 ",labelling_angle_data,sep=""),"H3K4ME3_")
    
    
  }
  #Spheoid axial distribution
  {
    dire<-paste(dird,"spheroid_axial_distribution/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_black_bg_histogram_timeline(subset(spheroid_axial$dna,spheroid_axial$dna$trial=="T1"), cols,parameters_axial_data,paste("DNA ",labelling_axial_data,sep=""),"DNA_T1_")
    plot_black_bg_histogram_timeline(subset(spheroid_axial$dna,spheroid_axial$dna$trial=="T2"), cols,parameters_axial_data,paste("DNA ",labelling_axial_data,sep=""),"DNA_T2_")
    plot_black_bg_histogram_timeline(subset(spheroid_axial$dna,spheroid_axial$dna$trial=="T3"), cols,parameters_axial_data,paste("DNA ",labelling_axial_data,sep=""),"DNA_T3_")
    
    plot_black_bg_boxplot_timeline(subset(spheroid_axial$dna,spheroid_axial$dna$trial=="T1"), cols,parameters_axial_data,paste("DNA ",labelling_axial_data,sep=""),"DNA_T1_")
    plot_black_bg_boxplot_timeline(subset(spheroid_axial$dna,spheroid_axial$dna$trial=="T2"), cols,parameters_axial_data,paste("DNA ",labelling_axial_data,sep=""),"DNA_T2_")
    plot_black_bg_boxplot_timeline(subset(spheroid_axial$dna,spheroid_axial$dna$trial=="T3"), cols,parameters_axial_data,paste("DNA ",labelling_axial_data,sep=""),"DNA_T3_")
    
    plot_black_bg_barplot_timeline(spheroid_axial$dna, cols,parameters_axial_data,paste("DNA ",labelling_axial_data,sep=""),"DNA_")
    
    plot_black_bg_histogram_timeline(subset(spheroid_axial$H3K4ME3,spheroid_axial$H3K4ME3$trial=="T1"), cols,parameters_axial_data,paste("H3K4ME3 ",labelling_axial_data,sep=""),"H3K4ME3_T1_")
    plot_black_bg_histogram_timeline(subset(spheroid_axial$H3K4ME3,spheroid_axial$H3K4ME3$trial=="T2"), cols,parameters_axial_data,paste("H3K4ME3 ",labelling_axial_data,sep=""),"H3K4ME3_T2_")
    plot_black_bg_histogram_timeline(subset(spheroid_axial$H3K4ME3,spheroid_axial$H3K4ME3$trial=="T3"), cols,parameters_axial_data,paste("H3K4ME3 ",labelling_axial_data,sep=""),"H3K4ME3_T3_")
    
    plot_black_bg_boxplot_timeline(subset(spheroid_axial$H3K4ME3,spheroid_axial$H3K4ME3$trial=="T1"), cols,parameters_axial_data,paste("H3K4ME3 ",labelling_axial_data,sep=""),"H3K4ME3_T1_")
    plot_black_bg_boxplot_timeline(subset(spheroid_axial$H3K4ME3,spheroid_axial$H3K4ME3$trial=="T2"), cols,parameters_axial_data,paste("H3K4ME3 ",labelling_axial_data,sep=""),"H3K4ME3_T2_")
    plot_black_bg_boxplot_timeline(subset(spheroid_axial$H3K4ME3,spheroid_axial$H3K4ME3$trial=="T3"), cols,parameters_axial_data,paste("H3K4ME3 ",labelling_axial_data,sep=""),"H3K4ME3_T3_")
    
    plot_black_bg_barplot_timeline(spheroid_axial$H3K4ME3, cols,parameters_axial_data,paste("H3K4ME3 ",labelling_axial_data,sep=""),"H3K4ME3_")
    
    
  }
  #Spheoid radial distribution
  {
    dire<-paste(dird,"spheroid_radial_distribution/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_black_bg_histogram_timeline(subset(spheroid_radial$dna,spheroid_radial$dna$trial=="T1"), cols,parameters_radial_data,paste("DNA ",labelling_radial_data,sep=""),"DNA_T1_")
    plot_black_bg_histogram_timeline(subset(spheroid_radial$dna,spheroid_radial$dna$trial=="T2"), cols,parameters_radial_data,paste("DNA ",labelling_radial_data,sep=""),"DNA_T2_")
    plot_black_bg_histogram_timeline(subset(spheroid_radial$dna,spheroid_radial$dna$trial=="T3"), cols,parameters_radial_data,paste("DNA ",labelling_radial_data,sep=""),"DNA_T3_")
    
    plot_black_bg_boxplot_timeline(subset(spheroid_radial$dna,spheroid_radial$dna$trial=="T1"), cols,parameters_radial_data,paste("DNA ",labelling_radial_data,sep=""),"DNA_T1_")
    plot_black_bg_boxplot_timeline(subset(spheroid_radial$dna,spheroid_radial$dna$trial=="T2"), cols,parameters_radial_data,paste("DNA ",labelling_radial_data,sep=""),"DNA_T2_")
    plot_black_bg_boxplot_timeline(subset(spheroid_radial$dna,spheroid_radial$dna$trial=="T3"), cols,parameters_radial_data,paste("DNA ",labelling_radial_data,sep=""),"DNA_T3_")
    
    plot_black_bg_barplot_timeline(spheroid_radial$dna, cols,parameters_radial_data,paste("DNA ",labelling_radial_data,sep=""),"DNA_")
    
    plot_black_bg_histogram_timeline(subset(spheroid_radial$H3K4ME3,spheroid_radial$H3K4ME3$trial=="T1"), cols,parameters_radial_data,paste("H3K4ME3 ",labelling_radial_data,sep=""),"H3K4ME3_T1_")
    plot_black_bg_histogram_timeline(subset(spheroid_radial$H3K4ME3,spheroid_radial$H3K4ME3$trial=="T2"), cols,parameters_radial_data,paste("H3K4ME3 ",labelling_radial_data,sep=""),"H3K4ME3_T2_")
    plot_black_bg_histogram_timeline(subset(spheroid_radial$H3K4ME3,spheroid_radial$H3K4ME3$trial=="T3"), cols,parameters_radial_data,paste("H3K4ME3 ",labelling_radial_data,sep=""),"H3K4ME3_T3_")
    
    plot_black_bg_boxplot_timeline(subset(spheroid_radial$H3K4ME3,spheroid_radial$H3K4ME3$trial=="T1"), cols,parameters_radial_data,paste("H3K4ME3 ",labelling_radial_data,sep=""),"H3K4ME3_T1_")
    plot_black_bg_boxplot_timeline(subset(spheroid_radial$H3K4ME3,spheroid_radial$H3K4ME3$trial=="T2"), cols,parameters_radial_data,paste("H3K4ME3 ",labelling_radial_data,sep=""),"H3K4ME3_T2_")
    plot_black_bg_boxplot_timeline(subset(spheroid_radial$H3K4ME3,spheroid_radial$H3K4ME3$trial=="T3"), cols,parameters_radial_data,paste("H3K4ME3 ",labelling_radial_data,sep=""),"H3K4ME3_T3_")
    
    plot_black_bg_barplot_timeline(spheroid_radial$H3K4ME3, cols,parameters_radial_data,paste("H3K4ME3 ",labelling_radial_data,sep=""),"H3K4ME3_")
    
    
  }
  
  #Spheoid geometeric properties
  {
    dire<-paste(dird,"spheroid_properties/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_black_bg_histogram_timeline(subset(spheroid_combined,spheroid_combined$trial=="T1"), cols,parameters_geometrical_data,paste("Spheroid ",labelling_geometrical_data,sep=""),"geo_T1_")
    plot_black_bg_histogram_timeline(subset(spheroid_combined,spheroid_combined$trial=="T2"), cols,parameters_geometrical_data,paste("Spheroid ",labelling_geometrical_data,sep=""),"geo_T2_")
    plot_black_bg_histogram_timeline(subset(spheroid_combined,spheroid_combined$trial=="T3"), cols,parameters_geometrical_data,paste("Spheroid ",labelling_geometrical_data,sep=""),"geo_T3_")
    
    plot_black_bg_boxplot_timeline(subset(spheroid_combined,spheroid_combined$trial=="T1"), cols,parameters_geometrical_data,labelling_geometrical_data,"geo_T1_")
    plot_black_bg_boxplot_timeline(subset(spheroid_combined,spheroid_combined$trial=="T2"), cols,parameters_geometrical_data,labelling_geometrical_data,"geo_T2_")
    plot_black_bg_boxplot_timeline(subset(spheroid_combined,spheroid_combined$trial=="T3"), cols,parameters_geometrical_data,labelling_geometrical_data,"geo_T3_")
    
    plot_black_bg_barplot_timeline(spheroid_combined, cols,parameters_geometrical_data,paste("Spheroid ",labelling_geometrical_data,sep=""),"geo_")
    
  }
  #Nuclear geometeric properties
  {
    dire<-paste(dird,"nuclear_properties/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_black_bg_histogram_timeline(subset(nucleus_combined,nucleus_combined$trial=="T1"), cols,parameters_geo_texture_data,paste("Nuclear ",labelling_geo_texture_data,sep=""),"geo_T1_")
    plot_black_bg_histogram_timeline(subset(nucleus_combined,nucleus_combined$trial=="T2"), cols,parameters_geo_texture_data,paste("Nuclear ",labelling_geo_texture_data,sep=""),"geo_T2_")
    plot_black_bg_histogram_timeline(subset(nucleus_combined,nucleus_combined$trial=="T3"), cols,parameters_geo_texture_data,paste("Nuclear ",labelling_geo_texture_data,sep=""),"geo_T3_")
    
    plot_black_bg_boxplot_timeline(subset(nucleus_combined,nucleus_combined$trial=="T1"), cols,parameters_geo_texture_data,labelling_geo_texture_data,"geo_T1_")
    plot_black_bg_boxplot_timeline(subset(nucleus_combined,nucleus_combined$trial=="T2"), cols,parameters_geo_texture_data,labelling_geo_texture_data,"geo_T2_")
    plot_black_bg_boxplot_timeline(subset(nucleus_combined,nucleus_combined$trial=="T3"), cols,parameters_geo_texture_data,labelling_geo_texture_data,"geo_T3_")
    
    plot_black_bg_barplot_timeline(nucleus_combined, cols,parameters_geo_texture_data,paste("Nuclear ",labelling_geo_texture_data,sep=""),"geo_")
    
  }
  
  #Nuclear Proj_Int
  {
    
    dire<-paste(dird,"nuclear_proj_int/",sep="")
    dir.create(dire)
    setwd(dire)
    
    plot_black_bg_histogram_timeline(subset(nuclei_2dint_data$H3K4ME3,nuclei_2dint_data$H3K4ME3$trial=="T1"), cols,parameters_2d_int_data,paste("H3K4ME3 ",nuc_labelling_2d_int_data,sep=""),"H3K4ME3_T1_")
    plot_black_bg_histogram_timeline(subset(nuclei_2dint_data$H3K4ME3,nuclei_2dint_data$H3K4ME3$trial=="T2"), cols,parameters_2d_int_data,paste("H3K4ME3 ",nuc_labelling_2d_int_data,sep=""),"H3K4ME3_T2_")
    plot_black_bg_histogram_timeline(subset(nuclei_2dint_data$H3K4ME3,nuclei_2dint_data$H3K4ME3$trial=="T3"), cols,parameters_2d_int_data,paste("H3K4ME3 ",nuc_labelling_2d_int_data,sep=""),"H3K4ME3_T3_")
    
    plot_black_bg_boxplot_timeline(subset(nuclei_2dint_data$H3K4ME3,nuclei_2dint_data$H3K4ME3$trial=="T1"), cols,parameters_2d_int_data,paste("H3K4ME3 ",nuc_labelling_2d_int_data,sep=""),"H3K4ME3_T1_")
    plot_black_bg_boxplot_timeline(subset(nuclei_2dint_data$H3K4ME3,nuclei_2dint_data$H3K4ME3$trial=="T2"), cols,parameters_2d_int_data,paste("H3K4ME3 ",nuc_labelling_2d_int_data,sep=""),"H3K4ME3_T2_")
    plot_black_bg_boxplot_timeline(subset(nuclei_2dint_data$H3K4ME3,nuclei_2dint_data$H3K4ME3$trial=="T3"), cols,parameters_2d_int_data,paste("H3K4ME3 ",nuc_labelling_2d_int_data,sep=""),"H3K4ME3_T3_")
    
    plot_black_bg_barplot_timeline(nuclei_2dint_data$H3K4ME3, cols,parameters_2d_int_data,paste("H3K4ME3 ",nuc_labelling_2d_int_data,sep=""),"H3K4ME3_")
    
    
    
  }
  #Nuclear 3D_Int
  {
    dire<-paste(dird,"nuclear_3d_int/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_black_bg_histogram_timeline(subset(nuclei_3dint_data$dna,nuclei_3dint_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T1_")
    plot_black_bg_histogram_timeline(subset(nuclei_3dint_data$dna,nuclei_3dint_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T2_")
    plot_black_bg_histogram_timeline(subset(nuclei_3dint_data$dna,nuclei_3dint_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T3_")
    
    plot_black_bg_boxplot_timeline(subset(nuclei_3dint_data$dna,nuclei_3dint_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T1_")
    plot_black_bg_boxplot_timeline(subset(nuclei_3dint_data$dna,nuclei_3dint_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T2_")
    plot_black_bg_boxplot_timeline(subset(nuclei_3dint_data$dna,nuclei_3dint_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T3_")
    
    plot_black_bg_barplot_timeline(nuclei_3dint_data$dna, cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_")
    
    plot_black_bg_histogram_timeline(subset(nuclei_3dint_data$H3K4ME3,nuclei_3dint_data$H3K4ME3$trial=="T1"), cols,parameters_int_data,paste("H3K4ME3 ",nuclabelling_int_data,sep=""),"H3K4ME3_T1_")
    plot_black_bg_histogram_timeline(subset(nuclei_3dint_data$H3K4ME3,nuclei_3dint_data$H3K4ME3$trial=="T2"), cols,parameters_int_data,paste("H3K4ME3 ",nuclabelling_int_data,sep=""),"H3K4ME3_T2_")
    plot_black_bg_histogram_timeline(subset(nuclei_3dint_data$H3K4ME3,nuclei_3dint_data$H3K4ME3$trial=="T3"), cols,parameters_int_data,paste("H3K4ME3 ",nuclabelling_int_data,sep=""),"H3K4ME3_T3_")
    
    plot_black_bg_boxplot_timeline(subset(nuclei_3dint_data$H3K4ME3,nuclei_3dint_data$H3K4ME3$trial=="T1"), cols,parameters_int_data,paste("H3K4ME3 ",nuclabelling_int_data,sep=""),"H3K4ME3_T1_")
    plot_black_bg_boxplot_timeline(subset(nuclei_3dint_data$H3K4ME3,nuclei_3dint_data$H3K4ME3$trial=="T2"), cols,parameters_int_data,paste("H3K4ME3 ",nuclabelling_int_data,sep=""),"H3K4ME3_T2_")
    plot_black_bg_boxplot_timeline(subset(nuclei_3dint_data$H3K4ME3,nuclei_3dint_data$H3K4ME3$trial=="T3"), cols,parameters_int_data,paste("H3K4ME3 ",nuclabelling_int_data,sep=""),"H3K4ME3_T3_")
    
    plot_black_bg_barplot_timeline(nuclei_3dint_data$H3K4ME3, cols,parameters_int_data,paste("H3K4ME3 ",nuclabelling_int_data,sep=""),"H3K4ME3_")
    
    
    
  }
}

#write the collected data
{
  dir.create("E:/HMF3A_reprogramming/R_analysis/H3K4ME3/combined_data/")
  setwd("E:/HMF3A_reprogramming/R_analysis/H3K4ME3/combined_data/")
  
  write.csv(spheroid_combined, file="spheroid_combined.csv")
  write.csv(nucleus_combined, file="nucleus_combined.csv")
  
  write.csv(spheroid_radial$H3K4ME3, file="spheroid_radial_dist_H3K4ME3.csv")
  
  write.csv(spheroid_2dint_data$H3K4ME3, file="spheroid_2d_int_H3K4ME3.csv")
  
  write.csv(spheroid_axial$H3K4ME3, file="spheroid_axial_dist_H3K4ME3.csv")
  write.csv(spheroid_axial$dna, file="spheroid_axial_dist_dna.csv")
  
  write.csv(spheroid_angle$H3K4ME3, file="spheroid_angle_dist_H3K4ME3.csv")
  write.csv(spheroid_angle$dna, file="spheroid_angle_dist_dna.csv")
  
  write.csv(spheroid_int_data$H3K4ME3, file="spheroid_int_H3K4ME3.csv")
  write.csv(spheroid_int_data$dna, file="spheroid_int_dna.csv")
  
  write.csv(nuclei_2dint_data$H3K4ME3, file="nuclear_2d_int_H3K4ME3.csv")
  
  write.csv(nuclei_3dint_data$H3K4ME3, file="nuclear_int_H3K4ME3.csv")
  write.csv(nuclei_3dint_data$dna, file="nuclear_int_dna.csv")
  
  
}