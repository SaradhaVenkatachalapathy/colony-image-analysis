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
  
  

  
}
source("E:/HMF3A_reprogramming/R_analysis/functions.R")

# initialise the experiment details and the subdirectories required
{
  path_to_ex<-c("E:/HMF3A_reprogramming/Rock_inhibition/20191215_hmf3a_ycompund_oct4_488_phalloidin_568_pmlc_647_T2/",
                "E:/HMF3A_reprogramming/Rock_inhibition/20191216_hmf3a_ycompund_oct4_488_phalloidin_568_pmlc_647_T3/",
                "E:/HMF3A_reprogramming/Rock_inhibition/20191220_hmf_s15_ycomp_pmlc_488_actin_568_oct4_647/")
  
  exp<-c("20191215_hmf3a_ycompund_oct4_488_phalloidin_568_pmlc_647_T2","20191216_hmf3a_ycompund_oct4_488_phalloidin_568_pmlc_647_T3","20191220_hmf_s15_ycomp_pmlc_488_actin_568_oct4_647")
  
  data_types<-c("/3D geometrical data spheroid/",
                "/3D int_data spheroid/DNA/","/3D int_data spheroid/ACTIN/","/3D int_data spheroid/PMLC/",
                "/3D int_data spheroid/OCT4")
  data_types_2d<-c("/2D_measures_spheroid/2D_spheroid_DNA.csv","/2D_measures_spheroid/2D_spheroid_ACTIN.csv",
                   "/2D_measures_spheroid/2D_spheroid_PMLC.csv","/2D_measures_spheroid/2D_spheroid_OCT4.csv")
  nuc_data_types<-c("/3D geometrical data/","/3D ellipsoid/","/3D geometerical_simple/","/3D shape measure/",
                    "/3D int_data/DNA/","/3D int_data/ACTIN/","/3D int_data/PMLC/","/3D int_data/OCT4/",
                    "/cell_2microns_measure/ACTIN/","/cell_2microns_measure/PMLC/","/cell_2microns_measure/OCT4/","/cell_2microns_measure/DNA/")
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
                       "/2D_measures_nuclei/2D_nucleusOCT4.csv","/2D_measures_nuclei/2D_nucleus_ACTIN.csv","/2D_measures_nuclei/2D_nucleus_PMLC.csv")
  
  
  samples<-c("cont","rocki_D6","rocki_D2")
  trial<-c("T1","T2","T3")
  
  dird_apo<-"E:/HMF3A_reprogramming/R_analysis/Ycompund_treatment/"
  dir.create(dird_apo)
  dird_oct4<-paste(dird_apo,"spheroid_shells_oct4", sep="")
  dird_actin<-paste(dird_apo,"spheroid_shells_actin", sep="")
  dird_pmlc<-paste(dird_apo,"spheroid_shells_pmlc", sep="")
  dird_dna<-paste(dird_apo,"spheroid_shells_dna", sep="")
  
  dir.create(dird_oct4)
  dir.create(dird_actin)
  dir.create(dird_pmlc)
  dir.create(dird_dna)
}

# read in all the required features for spheroids
{
  plot(0:1,0:1)
  #T1
  {
    geometrical_data<-combine_sample_sets(path_to_ex[1],samples,data_types[1],"tsv",exp[1],trial[1])
    T1DNA_int_data<-combine_sample_sets(path_to_ex[1],samples,data_types[2],"tsv",exp[1],trial[1])
    T1ACTIN_int_data<-combine_sample_sets(path_to_ex[1],samples,data_types[3],"tsv",exp[1],trial[1])
    T1PMLC_int_data<-combine_sample_sets(path_to_ex[1],samples,data_types[4],"tsv",exp[1],trial[1])
    T1OCT4_int_data<-combine_sample_sets(path_to_ex[1],samples,data_types[5],"tsv",exp[1],trial[1])
    geo_int_2D_data<-combine_sample_sets_2d(path_to_ex[1],samples,data_types_2d[1],exp[1],trial[1])
    
    T1_2d_data_ACTIN<-combine_sample_sets_2d(path_to_ex[1],samples,data_types_2d[2],exp[1],trial[1])
    T1_2d_data_PMLC<-combine_sample_sets_2d(path_to_ex[1],samples,data_types_2d[3],exp[1],trial[1])
    T1_2d_data_OCT4<-combine_sample_sets_2d(path_to_ex[1],samples,data_types_2d[4],exp[1],trial[1])
    
    shell_oct4_T1<-shell_analsysis(path_to_ex[1],samples,"/Shells/OCT4",exp[1],trial[1],dird_oct4)
    shell_actin_T1<-shell_analsysis(path_to_ex[1],samples,"/Shells/ACTIN",exp[1],trial[1],dird_actin)
    shell_pmlc_T1<-shell_analsysis(path_to_ex[1],samples,"/Shells/PMLC",exp[1],trial[1],dird_pmlc)
    shell_dna_T1<-shell_analsysis(path_to_ex[1],samples,"/Shells/DNA",exp[1],trial[1],dird_dna)
    
    T1_actin_axial<-axial_analysis(path_to_ex[1],samples,exp[1],trial[1],"/ACTIN/")
    T1_dna_axial<-axial_analysis(path_to_ex[1],samples,exp[1],trial[1],"/DNA/")
    T1_oct4_axial<-axial_analysis(path_to_ex[1],samples,exp[1],trial[1],"/OCT4/")
    T1_pmlc_axial<-axial_analysis(path_to_ex[1],samples,exp[1],trial[1],"/PMLC/")
    
    T1_actin_angle<-angle_analsysis(path_to_ex[1],samples,"/ACTIN",exp[1],trial[1])
    T1_oct4_angle<-angle_analsysis(path_to_ex[1],samples,"/OCT4",exp[1],trial[1])
    T1_pmlc_angle<-angle_analsysis(path_to_ex[1],samples,"/PMLC",exp[1],trial[1])
    T1_dna_angle<-angle_analsysis(path_to_ex[1],samples,"/DNA",exp[1],trial[1])
    
    geo_int_2D_data$Label<-substring(geo_int_2D_data$Label,5, nchar(geo_int_2D_data$Label)-4)
    geometrical_data$Label<-substring(geometrical_data$Label,0, nchar(geometrical_data$Label)-1)
    
    T1combined<-merge(geometrical_data,geo_int_2D_data, by="Label")
    
    rm(geometrical_data,geo_int_2D_data)
    
  }
  
  #T2
  {
    geometrical_data<-combine_sample_sets(path_to_ex[2],samples,data_types[1],"tsv",exp[2],trial[2])
    T2DNA_int_data<-combine_sample_sets(path_to_ex[2],samples,data_types[2],"tsv",exp[2],trial[2])
    T2ACTIN_int_data<-combine_sample_sets(path_to_ex[2],samples,data_types[3],"tsv",exp[2],trial[2])
    T2PMLC_int_data<-combine_sample_sets(path_to_ex[2],samples,data_types[4],"tsv",exp[2],trial[2])
    T2OCT4_int_data<-combine_sample_sets(path_to_ex[2],samples,data_types[5],"tsv",exp[2],trial[2])
    geo_int_2D_data<-combine_sample_sets_2d(path_to_ex[2],samples,data_types_2d[1],exp[2],trial[2])
    
    T2_2d_data_ACTIN<-combine_sample_sets_2d(path_to_ex[2],samples,data_types_2d[2],exp[2],trial[2])
    T2_2d_data_PMLC<-combine_sample_sets_2d(path_to_ex[2],samples,data_types_2d[3],exp[2],trial[2])
    T2_2d_data_OCT4<-combine_sample_sets_2d(path_to_ex[2],samples,data_types_2d[4],exp[2],trial[2])
    
    shell_oct4_T2<-shell_analsysis(path_to_ex[2],samples,"/Shells/OCT4",exp[2],trial[2],dird_oct4)
    shell_actin_T2<-shell_analsysis(path_to_ex[2],samples,"/Shells/ACTIN",exp[2],trial[2],dird_actin)
    shell_pmlc_T2<-shell_analsysis(path_to_ex[2],samples,"/Shells/PMLC",exp[2],trial[2],dird_pmlc)
    shell_dna_T2<-shell_analsysis(path_to_ex[2],samples,"/Shells/DNA",exp[2],trial[2],dird_dna)
    
    T2_actin_axial<-axial_analysis(path_to_ex[2],samples,exp[2],trial[2],"/ACTIN/")
    T2_dna_axial<-axial_analysis(path_to_ex[2],samples,exp[2],trial[2],"/DNA/")
    T2_oct4_axial<-axial_analysis(path_to_ex[2],samples,exp[2],trial[2],"/OCT4/")
    T2_pmlc_axial<-axial_analysis(path_to_ex[2],samples,exp[2],trial[2],"/PMLC/")
    
    T2_actin_angle<-angle_analsysis(path_to_ex[2],samples,"/ACTIN",exp[2],trial[2])
    T2_oct4_angle<-angle_analsysis(path_to_ex[2],samples,"/OCT4",exp[2],trial[2])
    T2_pmlc_angle<-angle_analsysis(path_to_ex[2],samples,"/PMLC",exp[2],trial[2])
    T2_dna_angle<-angle_analsysis(path_to_ex[2],samples,"/DNA",exp[2],trial[2])
    
    geo_int_2D_data$Label<-substring(geo_int_2D_data$Label,5, nchar(geo_int_2D_data$Label)-4)
    geometrical_data$Label<-substring(geometrical_data$Label,0, nchar(geometrical_data$Label)-1)
    T2combined<-merge(geometrical_data,geo_int_2D_data, by="Label")
    
    rm(geometrical_data,geo_int_2D_data)
    
  }
  
  #T3
  {
    geometrical_data<-combine_sample_sets(path_to_ex[3],samples,data_types[1],"tsv",exp[3],trial[3])
    T3DNA_int_data<-combine_sample_sets(path_to_ex[3],samples,data_types[2],"tsv",exp[3],trial[3])
    T3ACTIN_int_data<-combine_sample_sets(path_to_ex[3],samples,data_types[3],"tsv",exp[3],trial[3])
    T3PMLC_int_data<-combine_sample_sets(path_to_ex[3],samples,data_types[4],"tsv",exp[3],trial[3])
    T3OCT4_int_data<-combine_sample_sets(path_to_ex[3],samples,data_types[5],"tsv",exp[3],trial[3])
    geo_int_2D_data<-combine_sample_sets_2d(path_to_ex[3],samples,data_types_2d[1],exp[3],trial[3])
    
    T3_2d_data_ACTIN<-combine_sample_sets_2d(path_to_ex[3],samples,data_types_2d[2],exp[3],trial[3])
    T3_2d_data_PMLC<-combine_sample_sets_2d(path_to_ex[3],samples,data_types_2d[3],exp[3],trial[3])
    T3_2d_data_OCT4<-combine_sample_sets_2d(path_to_ex[3],samples,data_types_2d[4],exp[3],trial[3])
    
    shell_oct4_T3<-shell_analsysis(path_to_ex[3],samples,"/Shells/OCT4",exp[3],trial[3],dird_oct4)
    shell_actin_T3<-shell_analsysis(path_to_ex[3],samples,"/Shells/ACTIN",exp[3],trial[3],dird_actin)
    shell_pmlc_T3<-shell_analsysis(path_to_ex[3],samples,"/Shells/PMLC",exp[3],trial[3],dird_pmlc)
    shell_dna_T3<-shell_analsysis(path_to_ex[3],samples,"/Shells/DNA",exp[3],trial[3],dird_dna)
    
    T3_actin_axial<-axial_analysis(path_to_ex[3],samples,exp[3],trial[3],"/ACTIN/")
    T3_dna_axial<-axial_analysis(path_to_ex[3],samples,exp[3],trial[3],"/DNA/")
    T3_oct4_axial<-axial_analysis(path_to_ex[3],samples,exp[3],trial[3],"/OCT4/")
    T3_pmlc_axial<-axial_analysis(path_to_ex[3],samples,exp[3],trial[3],"/PMLC/")
    
    T3_actin_angle<-angle_analsysis(path_to_ex[3],samples,"/ACTIN",exp[3],trial[3])
    T3_oct4_angle<-angle_analsysis(path_to_ex[3],samples,"/OCT4",exp[3],trial[3])
    T3_pmlc_angle<-angle_analsysis(path_to_ex[3],samples,"/PMLC",exp[3],trial[3])
    T3_dna_angle<-angle_analsysis(path_to_ex[3],samples,"/DNA",exp[3],trial[3])
    
    geo_int_2D_data$Label<-substring(geo_int_2D_data$Label,5, nchar(geo_int_2D_data$Label)-4)
    geometrical_data$Label<-substring(geometrical_data$Label,0, nchar(geometrical_data$Label)-1)
    T3combined<-merge(geometrical_data,geo_int_2D_data, by="Label")
    
    rm(geometrical_data,geo_int_2D_data)
    
  }
  
  spheroid_2dint_data_actin<-rbind(T1_2d_data_ACTIN,T2_2d_data_ACTIN,T3_2d_data_ACTIN)
  spheroid_2dint_data_oct4<-rbind(T1_2d_data_OCT4,T2_2d_data_OCT4,T3_2d_data_OCT4)
  spheroid_2dint_data_pmlc<-rbind(T1_2d_data_PMLC,T2_2d_data_PMLC,T3_2d_data_PMLC)
  
  spheroid_int_data_actin<-rbind(T1ACTIN_int_data,T2ACTIN_int_data,T3ACTIN_int_data)
  spheroid_int_data_oct4<-rbind(T1OCT4_int_data,T2OCT4_int_data,T3OCT4_int_data)
  spheroid_int_data_pmlc<-rbind(T1PMLC_int_data,T2PMLC_int_data,T3PMLC_int_data)
  spheroid_int_data_dna<-rbind(T1DNA_int_data,T2DNA_int_data,T3DNA_int_data)
  
  spheroid_radial_oct4<-rbind(shell_oct4_T1,shell_oct4_T2,shell_oct4_T3)
  spheroid_radial_actin<-rbind(shell_actin_T1,shell_actin_T2,shell_actin_T3)
  spheroid_radial_pmlc<-rbind(shell_pmlc_T1,shell_pmlc_T2,shell_pmlc_T3)
  spheroid_radial_dna<-rbind(shell_dna_T1,shell_dna_T2,shell_dna_T3)
  
  spheroid_axial_actin<-rbind(T1_actin_axial,T2_actin_axial,T3_actin_axial)
  spheroid_axial_dna<-rbind(T1_dna_axial,T2_dna_axial,T3_dna_axial)
  spheroid_axial_oct4<-rbind(T1_oct4_axial,T2_oct4_axial,T3_oct4_axial)
  spheroid_axial_pmlc<-rbind(T1_pmlc_axial,T2_pmlc_axial,T3_pmlc_axial)
  
  spheroid_angle_actin<-rbind(T1_actin_angle,T2_actin_angle,T3_actin_angle)
  spheroid_angle_oct4<-rbind(T1_oct4_angle,T2_oct4_angle,T3_oct4_angle)
  spheroid_angle_pmlc<-rbind(T1_pmlc_angle,T2_pmlc_angle,T3_pmlc_angle)
  spheroid_angle_dna<-rbind(T1_dna_angle,T2_dna_angle,T3_dna_angle)
  
  spheroid_int_data<-list(actin=spheroid_int_data_actin, oct4=spheroid_int_data_oct4,
                          pmlc=spheroid_int_data_pmlc,dna=spheroid_int_data_dna)
  spheroid_radial<-list(actin=spheroid_radial_actin, oct4=spheroid_radial_oct4,
                        pmlc=spheroid_radial_pmlc,dna=spheroid_radial_dna)
  spheroid_axial<-list(actin=spheroid_axial_actin, oct4=spheroid_axial_oct4,
                       pmlc=spheroid_axial_pmlc,dna=spheroid_axial_dna)
  spheroid_angle<-list(actin=spheroid_angle_actin, oct4=spheroid_angle_oct4,
                       pmlc=spheroid_angle_pmlc,dna=spheroid_angle_dna)
  spheroid_2dint_data<-list(actin=spheroid_2dint_data_actin, oct4=spheroid_2dint_data_oct4,
                            pmlc=spheroid_2dint_data_pmlc)
  
  
  spheroid_combined<-rbind(T1combined,T2combined,T3combined)
  colnames(spheroid_combined)[86:87]<-c("sample","trial")
  rm(T1combined,T2combined,T3combined,T1ACTIN_int_data,T2ACTIN_int_data,T3ACTIN_int_data,
     T1OCT4_int_data,T2OCT4_int_data,T3OCT4_int_data,T1PMLC_int_data,T2PMLC_int_data,T3PMLC_int_data,
     shell_oct4_T1,shell_oct4_T2,shell_oct4_T3,shell_actin_T1,shell_actin_T2,shell_actin_T3,
     shell_pmlc_T1,shell_pmlc_T2,shell_pmlc_T3,shell_dna_T1,shell_dna_T2,shell_dna_T3,
     T1_actin_axial,T2_actin_axial,T3_actin_axial,T1_dna_axial,T2_dna_axial,T3_dna_axial,
     T1_oct4_axial,T2_oct4_axial,T3_oct4_axial,T1_pmlc_axial,T2_pmlc_axial,T3_pmlc_axial,
     T1_actin_angle,T2_actin_angle,T3_actin_angle,T1_oct4_angle,T2_oct4_angle,T3_oct4_angle,
     T1_pmlc_angle,T2_pmlc_angle,T3_pmlc_angle,T1_dna_angle,T2_dna_angle,T3_dna_angle,
     T1DNA_int_data,T2DNA_int_data,T3DNA_int_data,
     spheroid_angle_actin, spheroid_angle_oct4,spheroid_angle_pmlc,spheroid_angle_dna,
     spheroid_axial_actin, spheroid_axial_oct4,spheroid_axial_pmlc,spheroid_axial_dna,
     spheroid_radial_actin, spheroid_radial_oct4,spheroid_radial_pmlc,spheroid_radial_dna,
     spheroid_int_data_actin, spheroid_int_data_oct4,spheroid_int_data_pmlc,spheroid_int_data_dna,
     T1_2d_data_ACTIN,T2_2d_data_ACTIN,T3_2d_data_ACTIN,T1_2d_data_OCT4,T2_2d_data_OCT4,T3_2d_data_OCT4,
     T1_2d_data_PMLC,T2_2d_data_PMLC,T3_2d_data_PMLC,spheroid_2dint_data_actin, spheroid_2dint_data_oct4,
     spheroid_2dint_data_pmlc)
  
}

# read in all the required features for nucleus
{
  #T1
  {
    T1_geometrical_data<-combine_sample_sets(path_to_ex[1],samples,nuc_data_types[1],"tsv",exp[1],trial[1])
    T1_DNA_int_data<-combine_sample_sets(path_to_ex[1],samples,nuc_data_types[5],"tsv",exp[1],trial[1])
    T1_geo_int_2D_data<-combine_sample_sets_2d(path_to_ex[1],samples,nuc_data_types_2d[1],exp[1],trial[1])
    T1_compact_2D_data<-combine_sample_sets_2d(path_to_ex[1],samples,nuc_data_types_2d[2],exp[1],trial[1])
    T1_edf_2D_data<-combine_sample_sets_2d(path_to_ex[1],samples,nuc_data_types_2d[3],exp[1],trial[1])
    T1_gclm<-read_glcm(path_to_ex[1],samples,nuc_data_types_2d,exp[1],trial[1])
    
    T1_ACTIN_int_data<-combine_sample_sets(path_to_ex[1],samples,nuc_data_types[6],"tsv",exp[1],trial[1])
    T1_PMLC_int_data<-combine_sample_sets(path_to_ex[1],samples,nuc_data_types[7],"tsv",exp[1],trial[1])
    T1_OCT4_int_data<-combine_sample_sets(path_to_ex[1],samples,nuc_data_types[8],"tsv",exp[1],trial[1])
    
    T1_2d_ACTIN_int_data<-combine_sample_sets_2d(path_to_ex[1],samples,nuc_data_types_2d[29],exp[1],trial[1])
    T1_2d_PMLC_int_data<-combine_sample_sets_2d(path_to_ex[1],samples,nuc_data_types_2d[30],exp[1],trial[1])
    T1_2d_OCT4_int_data<-combine_sample_sets_2d(path_to_ex[1],samples,nuc_data_types_2d[28],exp[1],trial[1])
    
    T1_cell_ACTIN_int_data<-combine_sample_sets(path_to_ex[1],samples,nuc_data_types[9],"tsv",exp[1],trial[1])
    T1_cell_PMLC_int_data<-combine_sample_sets(path_to_ex[1],samples,nuc_data_types[10],"tsv",exp[1],trial[1])
    T1_cell_OCT4_int_data<-combine_sample_sets(path_to_ex[1],samples,nuc_data_types[11],"tsv",exp[1],trial[1])
    T1_cell_DNA_int_data<-combine_sample_sets(path_to_ex[1],samples,nuc_data_types[12],"tsv",exp[1],trial[1])
    
    colnames(T1_compact_2D_data)[2]<-"Label"
    T1_compact_2D_data$Label<-substring(T1_compact_2D_data$Label,0, nchar(T1_compact_2D_data$Label)-4)
    T1_geo_int_2D_data$Label<-substring(T1_geo_int_2D_data$Label,5, nchar(T1_geo_int_2D_data$Label)-4)
    T1_geometrical_data$Label<-substring(T1_geometrical_data$Label,0, nchar(T1_geometrical_data$Label)-1)
    T1_DNA_int_data$Label<-substring(T1_DNA_int_data$Label,0, nchar(T1_DNA_int_data$Label)-1)
    T1_2d_PMLC_int_data$Label<-substring(T1_2d_PMLC_int_data$Label,5, nchar(T1_2d_PMLC_int_data$Label)-4)
    T1_2d_ACTIN_int_data$Label<-substring(T1_2d_ACTIN_int_data$Label,5, nchar(T1_2d_ACTIN_int_data$Label)-4)
    T1_2d_OCT4_int_data$Label<-substring(T1_2d_OCT4_int_data$Label,5, nchar(T1_2d_OCT4_int_data$Label)-4)
    T1_ACTIN_int_data$Label<-substring(T1_ACTIN_int_data$Label,0, nchar(T1_ACTIN_int_data$Label)-1)
    T1_PMLC_int_data$Label<-substring(T1_PMLC_int_data$Label,0, nchar(T1_PMLC_int_data$Label)-1)
    T1_OCT4_int_data$Label<-substring(T1_OCT4_int_data$Label,0, nchar(T1_OCT4_int_data$Label)-1)
    T1_cell_ACTIN_int_data$Label<-substring(T1_cell_ACTIN_int_data$Label,0, nchar(T1_cell_ACTIN_int_data$Label)-1)
    T1_cell_PMLC_int_data$Label<-substring(T1_cell_PMLC_int_data$Label,0, nchar(T1_cell_PMLC_int_data$Label)-1)
    T1_cell_OCT4_int_data$Label<-substring(T1_cell_OCT4_int_data$Label,0, nchar(T1_cell_OCT4_int_data$Label)-1)
    T1_cell_DNA_int_data$Label<-substring(T1_cell_DNA_int_data$Label,0, nchar(T1_cell_DNA_int_data$Label)-1)
    
    T1_combined<-merge(T1_geometrical_data,T1_geo_int_2D_data[,2:36], by="Label")
    T1_combined<-merge(T1_combined,T1_compact_2D_data[,2:8], by="Label")
    T1_combined<-merge(T1_combined,T1_DNA_int_data[,c(3:5,12:18)], by="Label")
    T1_combined<-merge(T1_combined,T1_edf_2D_data[,1:11], by="Label")
    T1_combined<-cbind(T1_combined,T1_gclm[,-c(1,128:131)])
    
    rm(T1_geometrical_data,T1_geo_int_2D_data,T1_compact_2D_data,T1_edf_2D_data,T1_gclm)
  }
  #T2
  {
    T2_geometrical_data<-combine_sample_sets(path_to_ex[2],samples,nuc_data_types[1],"tsv",exp[2],trial[2])
    T2_DNA_int_data<-combine_sample_sets(path_to_ex[2],samples,nuc_data_types[5],"tsv",exp[2],trial[2])
    T2_geo_int_2D_data<-combine_sample_sets_2d(path_to_ex[2],samples,nuc_data_types_2d[1],exp[2],trial[2])
    T2_compact_2D_data<-combine_sample_sets_2d(path_to_ex[2],samples,nuc_data_types_2d[2],exp[2],trial[2])
    T2_edf_2D_data<-combine_sample_sets_2d(path_to_ex[2],samples,nuc_data_types_2d[3],exp[2],trial[2])
    T2_gclm<-read_glcm(path_to_ex[2],samples,nuc_data_types_2d,exp[2],trial[2])
    
    T2_ACTIN_int_data<-combine_sample_sets(path_to_ex[2],samples,nuc_data_types[6],"tsv",exp[2],trial[2])
    T2_PMLC_int_data<-combine_sample_sets(path_to_ex[2],samples,nuc_data_types[7],"tsv",exp[2],trial[2])
    T2_OCT4_int_data<-combine_sample_sets(path_to_ex[2],samples,nuc_data_types[8],"tsv",exp[2],trial[2])
    
    T2_2d_ACTIN_int_data<-combine_sample_sets_2d(path_to_ex[2],samples,nuc_data_types_2d[29],exp[2],trial[2])
    T2_2d_PMLC_int_data<-combine_sample_sets_2d(path_to_ex[2],samples,nuc_data_types_2d[30],exp[2],trial[2])
    T2_2d_OCT4_int_data<-combine_sample_sets_2d(path_to_ex[2],samples,nuc_data_types_2d[28],exp[2],trial[2])
    
    T2_cell_ACTIN_int_data<-combine_sample_sets(path_to_ex[2],samples,nuc_data_types[9],"tsv",exp[2],trial[2])
    T2_cell_PMLC_int_data<-combine_sample_sets(path_to_ex[2],samples,nuc_data_types[10],"tsv",exp[2],trial[2])
    T2_cell_OCT4_int_data<-combine_sample_sets(path_to_ex[2],samples,nuc_data_types[11],"tsv",exp[2],trial[2])
    T2_cell_DNA_int_data<-combine_sample_sets(path_to_ex[2],samples,nuc_data_types[12],"tsv",exp[2],trial[2])
    
    colnames(T2_compact_2D_data)[2]<-"Label"
    T2_compact_2D_data$Label<-substring(T2_compact_2D_data$Label,0, nchar(T2_compact_2D_data$Label)-4)
    T2_geo_int_2D_data$Label<-substring(T2_geo_int_2D_data$Label,5, nchar(T2_geo_int_2D_data$Label)-4)
    T2_geometrical_data$Label<-substring(T2_geometrical_data$Label,0, nchar(T2_geometrical_data$Label)-1)
    T2_2d_PMLC_int_data$Label<-substring(T2_2d_PMLC_int_data$Label,5, nchar(T2_2d_PMLC_int_data$Label)-4)
    T2_2d_ACTIN_int_data$Label<-substring(T2_2d_ACTIN_int_data$Label,5, nchar(T2_2d_ACTIN_int_data$Label)-4)
    T2_2d_OCT4_int_data$Label<-substring(T2_2d_OCT4_int_data$Label,5, nchar(T2_2d_OCT4_int_data$Label)-4)
    T2_DNA_int_data$Label<-substring(T2_DNA_int_data$Label,0, nchar(T2_DNA_int_data$Label)-1)
    T2_ACTIN_int_data$Label<-substring(T2_ACTIN_int_data$Label,0, nchar(T2_ACTIN_int_data$Label)-1)
    T2_PMLC_int_data$Label<-substring(T2_PMLC_int_data$Label,0, nchar(T2_PMLC_int_data$Label)-1)
    T2_OCT4_int_data$Label<-substring(T2_OCT4_int_data$Label,0, nchar(T2_OCT4_int_data$Label)-1)
    T2_cell_ACTIN_int_data$Label<-substring(T2_cell_ACTIN_int_data$Label,0, nchar(T2_cell_ACTIN_int_data$Label)-1)
    T2_cell_PMLC_int_data$Label<-substring(T2_cell_PMLC_int_data$Label,0, nchar(T2_cell_PMLC_int_data$Label)-1)
    T2_cell_OCT4_int_data$Label<-substring(T2_cell_OCT4_int_data$Label,0, nchar(T2_cell_OCT4_int_data$Label)-1)
    T2_cell_DNA_int_data$Label<-substring(T2_cell_DNA_int_data$Label,0, nchar(T2_cell_DNA_int_data$Label)-1)
    
    T2_combined<-merge(T2_geometrical_data,T2_geo_int_2D_data[,2:36], by="Label")
    T2_combined<-merge(T2_combined,T2_compact_2D_data[,2:8], by="Label")
    T2_combined<-merge(T2_combined,T2_DNA_int_data[,c(3:5,12:18)], by="Label")
    T2_combined<-merge(T2_combined,T2_edf_2D_data[,1:11], by="Label")
    T2_combined<-cbind(T2_combined,T2_gclm[,-c(1,128:131)])
    
    rm(T2_geometrical_data,T2_geo_int_2D_data,T2_compact_2D_data,T2_edf_2D_data,T2_gclm)
  }
  #T3
  {
    T3_geometrical_data<-combine_sample_sets(path_to_ex[3],samples,nuc_data_types[1],"tsv",exp[3],trial[3])
    T3_DNA_int_data<-combine_sample_sets(path_to_ex[3],samples,nuc_data_types[5],"tsv",exp[3],trial[3])
    T3_geo_int_2D_data<-combine_sample_sets_2d(path_to_ex[3],samples,nuc_data_types_2d[1],exp[3],trial[3])
    T3_compact_2D_data<-combine_sample_sets_2d(path_to_ex[3],samples,nuc_data_types_2d[2],exp[3],trial[3])
    T3_edf_2D_data<-combine_sample_sets_2d(path_to_ex[3],samples,nuc_data_types_2d[3],exp[3],trial[3])
    T3_gclm<-read_glcm(path_to_ex[3],samples,nuc_data_types_2d,exp[3],trial[3])
    
    T3_ACTIN_int_data<-combine_sample_sets(path_to_ex[3],samples,nuc_data_types[6],"tsv",exp[3],trial[3])
    T3_PMLC_int_data<-combine_sample_sets(path_to_ex[3],samples,nuc_data_types[7],"tsv",exp[3],trial[3])
    T3_OCT4_int_data<-combine_sample_sets(path_to_ex[3],samples,nuc_data_types[8],"tsv",exp[3],trial[3])
    
    T3_2d_ACTIN_int_data<-combine_sample_sets_2d(path_to_ex[3],samples,nuc_data_types_2d[29],exp[3],trial[3])
    T3_2d_PMLC_int_data<-combine_sample_sets_2d(path_to_ex[3],samples,nuc_data_types_2d[30],exp[3],trial[3])
    T3_2d_OCT4_int_data<-combine_sample_sets_2d(path_to_ex[3],samples,nuc_data_types_2d[28],exp[3],trial[3])
    
    T3_cell_ACTIN_int_data<-combine_sample_sets(path_to_ex[3],samples,nuc_data_types[9],"tsv",exp[3],trial[3])
    T3_cell_PMLC_int_data<-combine_sample_sets(path_to_ex[3],samples,nuc_data_types[10],"tsv",exp[3],trial[3])
    T3_cell_OCT4_int_data<-combine_sample_sets(path_to_ex[3],samples,nuc_data_types[11],"tsv",exp[3],trial[3])
    T3_cell_DNA_int_data<-combine_sample_sets(path_to_ex[3],samples,nuc_data_types[12],"tsv",exp[3],trial[3])
    
    colnames(T3_compact_2D_data)[2]<-"Label"
    T3_compact_2D_data$Label<-substring(T3_compact_2D_data$Label,0, nchar(T3_compact_2D_data$Label)-4)
    T3_geo_int_2D_data$Label<-substring(T3_geo_int_2D_data$Label,5, nchar(T3_geo_int_2D_data$Label)-4)
    T3_geometrical_data$Label<-substring(T3_geometrical_data$Label,0, nchar(T3_geometrical_data$Label)-1)
    T3_2d_PMLC_int_data$Label<-substring(T3_2d_PMLC_int_data$Label,5, nchar(T3_2d_PMLC_int_data$Label)-4)
    T3_2d_ACTIN_int_data$Label<-substring(T3_2d_ACTIN_int_data$Label,5, nchar(T3_2d_ACTIN_int_data$Label)-4)
    T3_2d_OCT4_int_data$Label<-substring(T3_2d_OCT4_int_data$Label,5, nchar(T3_2d_OCT4_int_data$Label)-4)
    T3_DNA_int_data$Label<-substring(T3_DNA_int_data$Label,0, nchar(T3_DNA_int_data$Label)-1)
    T3_ACTIN_int_data$Label<-substring(T3_ACTIN_int_data$Label,0, nchar(T3_ACTIN_int_data$Label)-1)
    T3_PMLC_int_data$Label<-substring(T3_PMLC_int_data$Label,0, nchar(T3_PMLC_int_data$Label)-1)
    T3_OCT4_int_data$Label<-substring(T3_OCT4_int_data$Label,0, nchar(T3_OCT4_int_data$Label)-1)
    T3_cell_ACTIN_int_data$Label<-substring(T3_cell_ACTIN_int_data$Label,0, nchar(T3_cell_ACTIN_int_data$Label)-1)
    T3_cell_PMLC_int_data$Label<-substring(T3_cell_PMLC_int_data$Label,0, nchar(T3_cell_PMLC_int_data$Label)-1)
    T3_cell_OCT4_int_data$Label<-substring(T3_cell_OCT4_int_data$Label,0, nchar(T3_cell_OCT4_int_data$Label)-1)
    T3_cell_DNA_int_data$Label<-substring(T3_cell_DNA_int_data$Label,0, nchar(T3_cell_DNA_int_data$Label)-1)
    
    T3_combined<-merge(T3_geometrical_data,T3_geo_int_2D_data[,2:36], by="Label")
    T3_combined<-merge(T3_combined,T3_compact_2D_data[,2:8], by="Label")
    T3_combined<-merge(T3_combined,T3_DNA_int_data[,c(3:5,12:18)], by="Label")
    T3_combined<-merge(T3_combined,T3_edf_2D_data[,1:11], by="Label")
    T3_combined<-cbind(T3_combined,T3_gclm[,-c(1,128:131)])
    
    rm(T3_geometrical_data,T3_geo_int_2D_data,T3_compact_2D_data,T3_edf_2D_data,T3_gclm)
  }
  
  nuclei_2dint_data_actin<-rbind(T1_2d_ACTIN_int_data,T2_2d_ACTIN_int_data,T3_2d_ACTIN_int_data)
  nuclei_2dint_data_oct4<-rbind(T1_2d_OCT4_int_data,T2_2d_OCT4_int_data,T3_2d_OCT4_int_data)
  nuclei_2dint_data_pmlc<-rbind(T1_2d_PMLC_int_data,T2_2d_PMLC_int_data,T3_2d_PMLC_int_data)
  
  nuclei_int_data_actin<-rbind(T1_ACTIN_int_data,T2_ACTIN_int_data,T3_ACTIN_int_data)
  nuclei_int_data_oct4<-rbind(T1_OCT4_int_data,T2_OCT4_int_data,T3_OCT4_int_data)
  nuclei_int_data_pmlc<-rbind(T1_PMLC_int_data,T2_PMLC_int_data,T3_PMLC_int_data)
  nuclei_int_data_dna<-rbind(T1_DNA_int_data,T2_DNA_int_data,T3_DNA_int_data)
  
  cell_int_data_actin<-rbind(T1_cell_ACTIN_int_data,T2_cell_ACTIN_int_data,T3_cell_ACTIN_int_data)
  cell_int_data_oct4<-rbind(T1_cell_OCT4_int_data,T2_cell_OCT4_int_data,T3_cell_OCT4_int_data)
  cell_int_data_pmlc<-rbind(T1_cell_PMLC_int_data,T2_cell_PMLC_int_data,T3_cell_PMLC_int_data)
  cell_int_data_dna<-rbind(T1_cell_DNA_int_data,T2_cell_DNA_int_data,T3_cell_DNA_int_data)
  
  nuclei_2dint_data<-list(actin=nuclei_2dint_data_actin, oct4=nuclei_2dint_data_oct4,
                          pmlc=nuclei_2dint_data_pmlc)
  nuclei_3dint_data<-list(actin=nuclei_int_data_actin, oct4=nuclei_int_data_oct4,
                          pmlc=nuclei_int_data_pmlc,dna=nuclei_int_data_dna)
  cell_3dint_data<-list(actin=cell_int_data_actin, oct4=cell_int_data_oct4,
                        pmlc=cell_int_data_pmlc,dna=cell_int_data_dna)
  
  nucleus_combined<-rbind(T1_combined,T2_combined,T3_combined)
  
  rm(T1_2d_ACTIN_int_data,T2_2d_ACTIN_int_data,T3_2d_ACTIN_int_data,T1_2d_OCT4_int_data,T2_2d_OCT4_int_data,T3_2d_OCT4_int_data,
     T1_2d_PMLC_int_data,T2_2d_PMLC_int_data,T3_2d_PMLC_int_data,T1_ACTIN_int_data,T2_ACTIN_int_data,T3_ACTIN_int_data,
     T1_OCT4_int_data,T2_OCT4_int_data,T3_OCT4_int_data,T1_PMLC_int_data,T2_PMLC_int_data,T3_PMLC_int_data,
     T1_DNA_int_data,T2_DNA_int_data,T3_DNA_int_data,T1_cell_ACTIN_int_data,T2_cell_ACTIN_int_data,T3_cell_ACTIN_int_data,
     T1_cell_OCT4_int_data,T2_cell_OCT4_int_data,T3_cell_OCT4_int_data,T1_cell_PMLC_int_data,T2_cell_PMLC_int_data,T3_cell_PMLC_int_data,
     T1_cell_DNA_int_data,T2_cell_DNA_int_data,T3_cell_DNA_int_data,T1_combined,T2_combined,T3_combined,
     nuclei_2dint_data_actin, nuclei_2dint_data_oct4,nuclei_2dint_data_pmlc,nuclei_int_data_actin, nuclei_int_data_oct4,
     nuclei_int_data_pmlc,nuclei_int_data_dna,cell_int_data_actin, cell_int_data_oct4,cell_int_data_pmlc,cell_int_data_dna)
  
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
  dird<-"E:/HMF3A_reprogramming/R_analysis/Ycompund_treatment/black_bg_plots/"
  dir.create(dird)
  setwd(dird)
  cols<-c("yellow","blue","violet")
  
  #Spheoid Proj_Int
  {
    
    dire<-paste(dird,"spheroid_proj_int/",sep="")
    dir.create(dire)
    setwd(dire)
    
    plot_black_bg_histogram_rocki(subset(spheroid_2dint_data$actin,spheroid_2dint_data$actin$trial=="T1"), cols,parameters_2d_int_data,paste("F-Actin ",labelling_2d_int_data,sep=""),"ACTIN_T1_")
    plot_black_bg_histogram_rocki(subset(spheroid_2dint_data$actin,spheroid_2dint_data$actin$trial=="T2"), cols,parameters_2d_int_data,paste("F-Actin ",labelling_2d_int_data,sep=""),"ACTIN_T2_")
    plot_black_bg_histogram_rocki(subset(spheroid_2dint_data$actin,spheroid_2dint_data$actin$trial=="T3"), cols,parameters_2d_int_data,paste("F-Actin ",labelling_2d_int_data,sep=""),"ACTIN_T3_")
    
    plot_black_bg_boxplot_rocki(subset(spheroid_2dint_data$actin,spheroid_2dint_data$actin$trial=="T1"), cols,parameters_2d_int_data,paste("F-Actin ",labelling_2d_int_data,sep=""),"ACTIN_T1_")
    plot_black_bg_boxplot_rocki(subset(spheroid_2dint_data$actin,spheroid_2dint_data$actin$trial=="T2"), cols,parameters_2d_int_data,paste("F-Actin ",labelling_2d_int_data,sep=""),"ACTIN_T2_")
    plot_black_bg_boxplot_rocki(subset(spheroid_2dint_data$actin,spheroid_2dint_data$actin$trial=="T3"), cols,parameters_2d_int_data,paste("F-Actin ",labelling_2d_int_data,sep=""),"ACTIN_T3_")
    
    plot_black_bg_barplot_rocki(spheroid_2dint_data$actin, cols,parameters_2d_int_data,paste("F-Actin ",labelling_2d_int_data,sep=""),"ACTIN_")
    
    
    plot_black_bg_histogram_rocki(subset(spheroid_2dint_data$pmlc,spheroid_2dint_data$pmlc$trial=="T1"), cols,parameters_2d_int_data,paste("pMLC ",labelling_2d_int_data,sep=""),"pMLC_T1_")
    plot_black_bg_histogram_rocki(subset(spheroid_2dint_data$pmlc,spheroid_2dint_data$pmlc$trial=="T2"), cols,parameters_2d_int_data,paste("pMLC ",labelling_2d_int_data,sep=""),"pMLC_T2_")
    plot_black_bg_histogram_rocki(subset(spheroid_2dint_data$pmlc,spheroid_2dint_data$pmlc$trial=="T3"), cols,parameters_2d_int_data,paste("pMLC ",labelling_2d_int_data,sep=""),"pMLC_T3_")
    
    plot_black_bg_boxplot_rocki(subset(spheroid_2dint_data$pmlc,spheroid_2dint_data$pmlc$trial=="T1"), cols,parameters_2d_int_data,paste("pMLC ",labelling_2d_int_data,sep=""),"pMLC_T1_")
    plot_black_bg_boxplot_rocki(subset(spheroid_2dint_data$pmlc,spheroid_2dint_data$pmlc$trial=="T2"), cols,parameters_2d_int_data,paste("pMLC ",labelling_2d_int_data,sep=""),"pMLC_T2_")
    plot_black_bg_boxplot_rocki(subset(spheroid_2dint_data$pmlc,spheroid_2dint_data$pmlc$trial=="T3"), cols,parameters_2d_int_data,paste("pMLC ",labelling_2d_int_data,sep=""),"pMLC_T3_")
    
    plot_black_bg_barplot_rocki(spheroid_2dint_data$pmlc, cols,parameters_2d_int_data,paste("pMLC ",labelling_2d_int_data,sep=""),"pMLC_")
    
    plot_black_bg_histogram_rocki(subset(spheroid_2dint_data$oct4,spheroid_2dint_data$oct4$trial=="T1"), cols,parameters_2d_int_data,paste("Oct4 ",labelling_2d_int_data,sep=""),"Oct4_T1_")
    plot_black_bg_histogram_rocki(subset(spheroid_2dint_data$oct4,spheroid_2dint_data$oct4$trial=="T2"), cols,parameters_2d_int_data,paste("Oct4 ",labelling_2d_int_data,sep=""),"Oct4_T2_")
    plot_black_bg_histogram_rocki(subset(spheroid_2dint_data$oct4,spheroid_2dint_data$oct4$trial=="T3"), cols,parameters_2d_int_data,paste("Oct4 ",labelling_2d_int_data,sep=""),"Oct4_T3_")
    
    plot_black_bg_boxplot_rocki(subset(spheroid_2dint_data$oct4,spheroid_2dint_data$oct4$trial=="T1"), cols,parameters_2d_int_data,paste("Oct4 ",labelling_2d_int_data,sep=""),"Oct4_T1_")
    plot_black_bg_boxplot_rocki(subset(spheroid_2dint_data$oct4,spheroid_2dint_data$oct4$trial=="T2"), cols,parameters_2d_int_data,paste("Oct4 ",labelling_2d_int_data,sep=""),"Oct4_T2_")
    plot_black_bg_boxplot_rocki(subset(spheroid_2dint_data$oct4,spheroid_2dint_data$oct4$trial=="T3"), cols,parameters_2d_int_data,paste("Oct4 ",labelling_2d_int_data,sep=""),"Oct4_T3_")
    
    plot_black_bg_barplot_rocki(spheroid_2dint_data$oct4, cols,parameters_2d_int_data,paste("Oct4 ",labelling_2d_int_data,sep=""),"Oct4_")
    
    
  }
  #Spheoid 3D_Int
  {
    dire<-paste(dird,"spheroid_3d_int/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_black_bg_histogram_rocki(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T1_")
    plot_black_bg_histogram_rocki(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T2_")
    plot_black_bg_histogram_rocki(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T3_")
    
    plot_black_bg_boxplot_rocki(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T1_")
    plot_black_bg_boxplot_rocki(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T2_")
    plot_black_bg_boxplot_rocki(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T3_")
    
    plot_black_bg_barplot_rocki(spheroid_int_data$dna, cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_")
    
    plot_black_bg_histogram_rocki(subset(spheroid_int_data$actin,spheroid_int_data$actin$trial=="T1"), cols,parameters_int_data,paste("F-Actin ",labelling_int_data,sep=""),"ACTIN_T1_")
    plot_black_bg_histogram_rocki(subset(spheroid_int_data$actin,spheroid_int_data$actin$trial=="T2"), cols,parameters_int_data,paste("F-Actin ",labelling_int_data,sep=""),"ACTIN_T2_")
    plot_black_bg_histogram_rocki(subset(spheroid_int_data$actin,spheroid_int_data$actin$trial=="T3"), cols,parameters_int_data,paste("F-Actin ",labelling_int_data,sep=""),"ACTIN_T3_")
    
    plot_black_bg_boxplot_rocki(subset(spheroid_int_data$actin,spheroid_int_data$actin$trial=="T1"), cols,parameters_int_data,paste("F-Actin ",labelling_int_data,sep=""),"ACTIN_T1_")
    plot_black_bg_boxplot_rocki(subset(spheroid_int_data$actin,spheroid_int_data$actin$trial=="T2"), cols,parameters_int_data,paste("F-Actin ",labelling_int_data,sep=""),"ACTIN_T2_")
    plot_black_bg_boxplot_rocki(subset(spheroid_int_data$actin,spheroid_int_data$actin$trial=="T3"), cols,parameters_int_data,paste("F-Actin ",labelling_int_data,sep=""),"ACTIN_T3_")
    
    plot_black_bg_barplot_rocki(spheroid_int_data$actin, cols,parameters_int_data,paste("F-Actin ",labelling_int_data,sep=""),"ACTIN_")
    
    
    plot_black_bg_histogram_rocki(subset(spheroid_int_data$pmlc,spheroid_int_data$pmlc$trial=="T1"), cols,parameters_int_data,paste("pMLC ",labelling_int_data,sep=""),"pMLC_T1_")
    plot_black_bg_histogram_rocki(subset(spheroid_int_data$pmlc,spheroid_int_data$pmlc$trial=="T2"), cols,parameters_int_data,paste("pMLC ",labelling_int_data,sep=""),"pMLC_T2_")
    plot_black_bg_histogram_rocki(subset(spheroid_int_data$pmlc,spheroid_int_data$pmlc$trial=="T3"), cols,parameters_int_data,paste("pMLC ",labelling_int_data,sep=""),"pMLC_T3_")
    
    plot_black_bg_boxplot_rocki(subset(spheroid_int_data$pmlc,spheroid_int_data$pmlc$trial=="T1"), cols,parameters_int_data,paste("pMLC ",labelling_int_data,sep=""),"pMLC_T1_")
    plot_black_bg_boxplot_rocki(subset(spheroid_int_data$pmlc,spheroid_int_data$pmlc$trial=="T2"), cols,parameters_int_data,paste("pMLC ",labelling_int_data,sep=""),"pMLC_T2_")
    plot_black_bg_boxplot_rocki(subset(spheroid_int_data$pmlc,spheroid_int_data$pmlc$trial=="T3"), cols,parameters_int_data,paste("pMLC ",labelling_int_data,sep=""),"pMLC_T3_")
    
    plot_black_bg_barplot_rocki(spheroid_int_data$pmlc, cols,parameters_int_data,paste("pMLC ",labelling_int_data,sep=""),"pMLC_")
    
    plot_black_bg_histogram_rocki(subset(spheroid_int_data$oct4,spheroid_int_data$oct4$trial=="T1"), cols,parameters_int_data,paste("Oct4 ",labelling_int_data,sep=""),"Oct4_T1_")
    plot_black_bg_histogram_rocki(subset(spheroid_int_data$oct4,spheroid_int_data$oct4$trial=="T2"), cols,parameters_int_data,paste("Oct4 ",labelling_int_data,sep=""),"Oct4_T2_")
    plot_black_bg_histogram_rocki(subset(spheroid_int_data$oct4,spheroid_int_data$oct4$trial=="T3"), cols,parameters_int_data,paste("Oct4 ",labelling_int_data,sep=""),"Oct4_T3_")
    
    plot_black_bg_boxplot_rocki(subset(spheroid_int_data$oct4,spheroid_int_data$oct4$trial=="T1"), cols,parameters_int_data,paste("Oct4 ",labelling_int_data,sep=""),"Oct4_T1_")
    plot_black_bg_boxplot_rocki(subset(spheroid_int_data$oct4,spheroid_int_data$oct4$trial=="T2"), cols,parameters_int_data,paste("Oct4 ",labelling_int_data,sep=""),"Oct4_T2_")
    plot_black_bg_boxplot_rocki(subset(spheroid_int_data$oct4,spheroid_int_data$oct4$trial=="T3"), cols,parameters_int_data,paste("Oct4 ",labelling_int_data,sep=""),"Oct4_T3_")
    
    plot_black_bg_barplot_rocki(spheroid_int_data$oct4, cols,parameters_int_data,paste("Oct4 ",labelling_int_data,sep=""),"Oct4_")
    
    
  }
  #Spheoid angles
  {
    dire<-paste(dird,"spheroid_angles/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_black_bg_histogram_rocki(subset(spheroid_angle$dna,spheroid_angle$dna$trial=="T1"), cols,parameters_angle_data,paste("DNA ",labelling_angle_data,sep=""),"DNA_T1_")
    plot_black_bg_histogram_rocki(subset(spheroid_angle$dna,spheroid_angle$dna$trial=="T2"), cols,parameters_angle_data,paste("DNA ",labelling_angle_data,sep=""),"DNA_T2_")
    plot_black_bg_histogram_rocki(subset(spheroid_angle$dna,spheroid_angle$dna$trial=="T3"), cols,parameters_angle_data,paste("DNA ",labelling_angle_data,sep=""),"DNA_T3_")
    
    plot_black_bg_boxplot_rocki(subset(spheroid_angle$dna,spheroid_angle$dna$trial=="T1"), cols,parameters_angle_data,paste("DNA ",labelling_angle_data,sep=""),"DNA_T1_")
    plot_black_bg_boxplot_rocki(subset(spheroid_angle$dna,spheroid_angle$dna$trial=="T2"), cols,parameters_angle_data,paste("DNA ",labelling_angle_data,sep=""),"DNA_T2_")
    plot_black_bg_boxplot_rocki(subset(spheroid_angle$dna,spheroid_angle$dna$trial=="T3"), cols,parameters_angle_data,paste("DNA ",labelling_angle_data,sep=""),"DNA_T3_")
    
    
    plot_black_bg_histogram_rocki(subset(spheroid_angle$actin,spheroid_angle$actin$trial=="T1"), cols,parameters_angle_data,paste("F-Actin ",labelling_angle_data,sep=""),"ACTIN_T1_")
    plot_black_bg_histogram_rocki(subset(spheroid_angle$actin,spheroid_angle$actin$trial=="T2"), cols,parameters_angle_data,paste("F-Actin ",labelling_angle_data,sep=""),"ACTIN_T2_")
    plot_black_bg_histogram_rocki(subset(spheroid_angle$actin,spheroid_angle$actin$trial=="T3"), cols,parameters_angle_data,paste("F-Actin ",labelling_angle_data,sep=""),"ACTIN_T3_")
    
    plot_black_bg_boxplot_rocki(subset(spheroid_angle$actin,spheroid_angle$actin$trial=="T1"), cols,parameters_angle_data,paste("F-Actin ",labelling_angle_data,sep=""),"ACTIN_T1_")
    plot_black_bg_boxplot_rocki(subset(spheroid_angle$actin,spheroid_angle$actin$trial=="T2"), cols,parameters_angle_data,paste("F-Actin ",labelling_angle_data,sep=""),"ACTIN_T2_")
    plot_black_bg_boxplot_rocki(subset(spheroid_angle$actin,spheroid_angle$actin$trial=="T3"), cols,parameters_angle_data,paste("F-Actin ",labelling_angle_data,sep=""),"ACTIN_T3_")
    
    plot_black_bg_barplot_rocki(spheroid_angle$actin, cols,parameters_angle_data,paste("F-Actin ",labelling_angle_data,sep=""),"ACTIN_")
    
    
    plot_black_bg_histogram_rocki(subset(spheroid_angle$pmlc,spheroid_angle$pmlc$trial=="T1"), cols,parameters_angle_data,paste("pMLC ",labelling_angle_data,sep=""),"pMLC_T1_")
    plot_black_bg_histogram_rocki(subset(spheroid_angle$pmlc,spheroid_angle$pmlc$trial=="T2"), cols,parameters_angle_data,paste("pMLC ",labelling_angle_data,sep=""),"pMLC_T2_")
    plot_black_bg_histogram_rocki(subset(spheroid_angle$pmlc,spheroid_angle$pmlc$trial=="T3"), cols,parameters_angle_data,paste("pMLC ",labelling_angle_data,sep=""),"pMLC_T3_")
    
    plot_black_bg_boxplot_rocki(subset(spheroid_angle$pmlc,spheroid_angle$pmlc$trial=="T1"), cols,parameters_angle_data,paste("pMLC ",labelling_angle_data,sep=""),"pMLC_T1_")
    plot_black_bg_boxplot_rocki(subset(spheroid_angle$pmlc,spheroid_angle$pmlc$trial=="T2"), cols,parameters_angle_data,paste("pMLC ",labelling_angle_data,sep=""),"pMLC_T2_")
    plot_black_bg_boxplot_rocki(subset(spheroid_angle$pmlc,spheroid_angle$pmlc$trial=="T3"), cols,parameters_angle_data,paste("pMLC ",labelling_angle_data,sep=""),"pMLC_T3_")
    
    plot_black_bg_barplot_rocki(spheroid_angle$pmlc, cols,parameters_angle_data,paste("pMLC ",labelling_angle_data,sep=""),"pMLC_")
    
    plot_black_bg_histogram_rocki(subset(spheroid_angle$oct4,spheroid_angle$oct4$trial=="T1"), cols,parameters_angle_data,paste("Oct4 ",labelling_angle_data,sep=""),"Oct4_T1_")
    plot_black_bg_histogram_rocki(subset(spheroid_angle$oct4,spheroid_angle$oct4$trial=="T2"), cols,parameters_angle_data,paste("Oct4 ",labelling_angle_data,sep=""),"Oct4_T2_")
    plot_black_bg_histogram_rocki(subset(spheroid_angle$oct4,spheroid_angle$oct4$trial=="T3"), cols,parameters_angle_data,paste("Oct4 ",labelling_angle_data,sep=""),"Oct4_T3_")
    
    plot_black_bg_boxplot_rocki(subset(spheroid_angle$oct4,spheroid_angle$oct4$trial=="T1"), cols,parameters_angle_data,paste("Oct4 ",labelling_angle_data,sep=""),"Oct4_T1_")
    plot_black_bg_boxplot_rocki(subset(spheroid_angle$oct4,spheroid_angle$oct4$trial=="T2"), cols,parameters_angle_data,paste("Oct4 ",labelling_angle_data,sep=""),"Oct4_T2_")
    plot_black_bg_boxplot_rocki(subset(spheroid_angle$oct4,spheroid_angle$oct4$trial=="T3"), cols,parameters_angle_data,paste("Oct4 ",labelling_angle_data,sep=""),"Oct4_T3_")
    
    plot_black_bg_barplot_rocki(spheroid_angle$oct4, cols,parameters_angle_data,paste("Oct4 ",labelling_angle_data,sep=""),"Oct4_")
    
    angles_plot_black_bg_rocki(path_to_ex[1],samples,"/ACTIN/","T1_actin.png",dire,cols)
    angles_plot_black_bg_rocki(path_to_ex[2],samples,"ACTIN/","T2_actin.png",dire,cols)
    angles_plot_black_bg_rocki(path_to_ex[3],samples,"ACTIN/","T3_actin.png",dire,cols)
    
    angles_plot_black_bg_rocki(path_to_ex[1],samples,"/DNA/","T1_dna.png",dire,cols)
    angles_plot_black_bg_rocki(path_to_ex[2],samples,"/DNA/","T2_dna.png",dire,cols)
    angles_plot_black_bg_rocki(path_to_ex[3],samples,"/DNA/","T3_dna.png",dire,cols)
    
    angles_plot_black_bg_rocki(path_to_ex[1],samples,"/PMLC/","T1_pmlc.png",dire,cols)
    angles_plot_black_bg_rocki(path_to_ex[2],samples,"/PMLC/","T2_pmlc.png",dire,cols)
    angles_plot_black_bg_rocki(path_to_ex[3],samples,"/PMLC/","T3_pmlc.png",dire,cols)
    
    angles_plot_black_bg_rocki(path_to_ex[1],samples,"/OCT4/","T1_oct4.png",dire,cols)
    angles_plot_black_bg_rocki(path_to_ex[2],samples,"/OCT4/","T2_oct4.png",dire,cols)
    angles_plot_black_bg_rocki(path_to_ex[3],samples,"/OCT4/","T3_oct4.png",dire,cols)
    
  }
  #Spheoid axial distribution
  {
    dire<-paste(dird,"spheroid_axial_distribution/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_black_bg_histogram_rocki(subset(spheroid_axial$dna,spheroid_axial$dna$trial=="T1"), cols,parameters_axial_data,paste("DNA ",labelling_axial_data,sep=""),"DNA_T1_")
    plot_black_bg_histogram_rocki(subset(spheroid_axial$dna,spheroid_axial$dna$trial=="T2"), cols,parameters_axial_data,paste("DNA ",labelling_axial_data,sep=""),"DNA_T2_")
    plot_black_bg_histogram_rocki(subset(spheroid_axial$dna,spheroid_axial$dna$trial=="T3"), cols,parameters_axial_data,paste("DNA ",labelling_axial_data,sep=""),"DNA_T3_")
    
    plot_black_bg_boxplot_rocki(subset(spheroid_axial$dna,spheroid_axial$dna$trial=="T1"), cols,parameters_axial_data,paste("DNA ",labelling_axial_data,sep=""),"DNA_T1_")
    plot_black_bg_boxplot_rocki(subset(spheroid_axial$dna,spheroid_axial$dna$trial=="T2"), cols,parameters_axial_data,paste("DNA ",labelling_axial_data,sep=""),"DNA_T2_")
    plot_black_bg_boxplot_rocki(subset(spheroid_axial$dna,spheroid_axial$dna$trial=="T3"), cols,parameters_axial_data,paste("DNA ",labelling_axial_data,sep=""),"DNA_T3_")
    
    plot_black_bg_barplot_rocki(spheroid_axial$dna, cols,parameters_axial_data,paste("DNA ",labelling_axial_data,sep=""),"DNA_")
    
    plot_black_bg_histogram_rocki(subset(spheroid_axial$actin,spheroid_axial$actin$trial=="T1"), cols,parameters_axial_data,paste("F-Actin ",labelling_axial_data,sep=""),"ACTIN_T1_")
    plot_black_bg_histogram_rocki(subset(spheroid_axial$actin,spheroid_axial$actin$trial=="T2"), cols,parameters_axial_data,paste("F-Actin ",labelling_axial_data,sep=""),"ACTIN_T2_")
    plot_black_bg_histogram_rocki(subset(spheroid_axial$actin,spheroid_axial$actin$trial=="T3"), cols,parameters_axial_data,paste("F-Actin ",labelling_axial_data,sep=""),"ACTIN_T3_")
    
    plot_black_bg_boxplot_rocki(subset(spheroid_axial$actin,spheroid_axial$actin$trial=="T1"), cols,parameters_axial_data,paste("F-Actin ",labelling_axial_data,sep=""),"ACTIN_T1_")
    plot_black_bg_boxplot_rocki(subset(spheroid_axial$actin,spheroid_axial$actin$trial=="T2"), cols,parameters_axial_data,paste("F-Actin ",labelling_axial_data,sep=""),"ACTIN_T2_")
    plot_black_bg_boxplot_rocki(subset(spheroid_axial$actin,spheroid_axial$actin$trial=="T3"), cols,parameters_axial_data,paste("F-Actin ",labelling_axial_data,sep=""),"ACTIN_T3_")
    
    plot_black_bg_barplot_rocki(spheroid_axial$actin, cols,parameters_axial_data,paste("F-Actin ",labelling_axial_data,sep=""),"ACTIN_")
    
    
    plot_black_bg_histogram_rocki(subset(spheroid_axial$pmlc,spheroid_axial$pmlc$trial=="T1"), cols,parameters_axial_data,paste("pMLC ",labelling_axial_data,sep=""),"pMLC_T1_")
    plot_black_bg_histogram_rocki(subset(spheroid_axial$pmlc,spheroid_axial$pmlc$trial=="T2"), cols,parameters_axial_data,paste("pMLC ",labelling_axial_data,sep=""),"pMLC_T2_")
    plot_black_bg_histogram_rocki(subset(spheroid_axial$pmlc,spheroid_axial$pmlc$trial=="T3"), cols,parameters_axial_data,paste("pMLC ",labelling_axial_data,sep=""),"pMLC_T3_")
    
    plot_black_bg_boxplot_rocki(subset(spheroid_axial$pmlc,spheroid_axial$pmlc$trial=="T1"), cols,parameters_axial_data,paste("pMLC ",labelling_axial_data,sep=""),"pMLC_T1_")
    plot_black_bg_boxplot_rocki(subset(spheroid_axial$pmlc,spheroid_axial$pmlc$trial=="T2"), cols,parameters_axial_data,paste("pMLC ",labelling_axial_data,sep=""),"pMLC_T2_")
    plot_black_bg_boxplot_rocki(subset(spheroid_axial$pmlc,spheroid_axial$pmlc$trial=="T3"), cols,parameters_axial_data,paste("pMLC ",labelling_axial_data,sep=""),"pMLC_T3_")
    
    plot_black_bg_barplot_rocki(spheroid_axial$pmlc, cols,parameters_axial_data,paste("pMLC ",labelling_axial_data,sep=""),"pMLC_")
    
    plot_black_bg_histogram_rocki(subset(spheroid_axial$oct4,spheroid_axial$oct4$trial=="T1"), cols,parameters_axial_data,paste("Oct4 ",labelling_axial_data,sep=""),"Oct4_T1_")
    plot_black_bg_histogram_rocki(subset(spheroid_axial$oct4,spheroid_axial$oct4$trial=="T2"), cols,parameters_axial_data,paste("Oct4 ",labelling_axial_data,sep=""),"Oct4_T2_")
    plot_black_bg_histogram_rocki(subset(spheroid_axial$oct4,spheroid_axial$oct4$trial=="T3"), cols,parameters_axial_data,paste("Oct4 ",labelling_axial_data,sep=""),"Oct4_T3_")
    
    plot_black_bg_boxplot_rocki(subset(spheroid_axial$oct4,spheroid_axial$oct4$trial=="T1"), cols,parameters_axial_data,paste("Oct4 ",labelling_axial_data,sep=""),"Oct4_T1_")
    plot_black_bg_boxplot_rocki(subset(spheroid_axial$oct4,spheroid_axial$oct4$trial=="T2"), cols,parameters_axial_data,paste("Oct4 ",labelling_axial_data,sep=""),"Oct4_T2_")
    plot_black_bg_boxplot_rocki(subset(spheroid_axial$oct4,spheroid_axial$oct4$trial=="T3"), cols,parameters_axial_data,paste("Oct4 ",labelling_axial_data,sep=""),"Oct4_T3_")
    
    plot_black_bg_barplot_rocki(spheroid_axial$oct4, cols,parameters_axial_data,paste("Oct4 ",labelling_axial_data,sep=""),"Oct4_")
    
    axial_plot_black_bg_boxplot_rocki(subset(spheroid_axial$dna,spheroid_axial$dna$trial=="T1"),cols,"Mean DNA","T1_dna.png")
    axial_plot_black_bg_boxplot_rocki(subset(spheroid_axial$dna,spheroid_axial$dna$trial=="T2"),cols,"Mean DNA","T2_dna.png")
    axial_plot_black_bg_boxplot_rocki(subset(spheroid_axial$dna,spheroid_axial$dna$trial=="T3"),cols,"Mean DNA","T3_dna.png")
    
    axial_plot_black_bg_boxplot_rocki(subset(spheroid_axial$actin,spheroid_axial$actin$trial=="T1"),cols,"Mean Actin","T1_Actin.png")
    axial_plot_black_bg_boxplot_rocki(subset(spheroid_axial$actin,spheroid_axial$actin$trial=="T2"),cols,"Mean Actin","T2_Actin.png")
    axial_plot_black_bg_boxplot_rocki(subset(spheroid_axial$actin,spheroid_axial$actin$trial=="T3"),cols,"Mean Actin","T3_Actin.png")
    
    axial_plot_black_bg_boxplot_rocki(subset(spheroid_axial$oct4,spheroid_axial$oct4$trial=="T1"),cols,"Mean Oct4","T1_oct4.png")
    axial_plot_black_bg_boxplot_rocki(subset(spheroid_axial$oct4,spheroid_axial$oct4$trial=="T2"),cols,"Mean Oct4","T2_oct4.png")
    axial_plot_black_bg_boxplot_rocki(subset(spheroid_axial$oct4,spheroid_axial$oct4$trial=="T3"),cols,"Mean Oct4","T3_oct4.png")
    
    axial_plot_black_bg_boxplot_rocki(subset(spheroid_axial$oct4,spheroid_axial$oct4$trial=="T1"),cols,"Mean pMLC","T1_PMLC.png")
    axial_plot_black_bg_boxplot_rocki(subset(spheroid_axial$oct4,spheroid_axial$oct4$trial=="T2"),cols,"Mean pMLC","T2_PMLC.png")
    axial_plot_black_bg_boxplot_rocki(subset(spheroid_axial$oct4,spheroid_axial$oct4$trial=="T3"),cols,"Mean pMLC","T3_PMLC.png")
    
    
  }
  #Spheoid radial distribution
  {
    dire<-paste(dird,"spheroid_radial_distribution/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_black_bg_histogram_rocki(subset(spheroid_radial$dna,spheroid_radial$dna$trial=="T1"), cols,parameters_radial_data,paste("DNA ",labelling_radial_data,sep=""),"DNA_T1_")
    plot_black_bg_histogram_rocki(subset(spheroid_radial$dna,spheroid_radial$dna$trial=="T2"), cols,parameters_radial_data,paste("DNA ",labelling_radial_data,sep=""),"DNA_T2_")
    plot_black_bg_histogram_rocki(subset(spheroid_radial$dna,spheroid_radial$dna$trial=="T3"), cols,parameters_radial_data,paste("DNA ",labelling_radial_data,sep=""),"DNA_T3_")
    
    plot_black_bg_boxplot_rocki(subset(spheroid_radial$dna,spheroid_radial$dna$trial=="T1"), cols,parameters_radial_data,paste("DNA ",labelling_radial_data,sep=""),"DNA_T1_")
    plot_black_bg_boxplot_rocki(subset(spheroid_radial$dna,spheroid_radial$dna$trial=="T2"), cols,parameters_radial_data,paste("DNA ",labelling_radial_data,sep=""),"DNA_T2_")
    plot_black_bg_boxplot_rocki(subset(spheroid_radial$dna,spheroid_radial$dna$trial=="T3"), cols,parameters_radial_data,paste("DNA ",labelling_radial_data,sep=""),"DNA_T3_")
    
    plot_black_bg_barplot_rocki(spheroid_radial$dna, cols,parameters_radial_data,paste("DNA ",labelling_radial_data,sep=""),"DNA_")
    
    plot_black_bg_histogram_rocki(subset(spheroid_radial$actin,spheroid_radial$actin$trial=="T1"), cols,parameters_radial_data,paste("F-Actin ",labelling_radial_data,sep=""),"ACTIN_T1_")
    plot_black_bg_histogram_rocki(subset(spheroid_radial$actin,spheroid_radial$actin$trial=="T2"), cols,parameters_radial_data,paste("F-Actin ",labelling_radial_data,sep=""),"ACTIN_T2_")
    plot_black_bg_histogram_rocki(subset(spheroid_radial$actin,spheroid_radial$actin$trial=="T3"), cols,parameters_radial_data,paste("F-Actin ",labelling_radial_data,sep=""),"ACTIN_T3_")
    
    plot_black_bg_boxplot_rocki(subset(spheroid_radial$actin,spheroid_radial$actin$trial=="T1"), cols,parameters_radial_data,paste("F-Actin ",labelling_radial_data,sep=""),"ACTIN_T1_")
    plot_black_bg_boxplot_rocki(subset(spheroid_radial$actin,spheroid_radial$actin$trial=="T2"), cols,parameters_radial_data,paste("F-Actin ",labelling_radial_data,sep=""),"ACTIN_T2_")
    plot_black_bg_boxplot_rocki(subset(spheroid_radial$actin,spheroid_radial$actin$trial=="T3"), cols,parameters_radial_data,paste("F-Actin ",labelling_radial_data,sep=""),"ACTIN_T3_")
    
    plot_black_bg_barplot_rocki(spheroid_radial$actin, cols,parameters_radial_data,paste("F-Actin ",labelling_radial_data,sep=""),"ACTIN_")
    
    
    plot_black_bg_histogram_rocki(subset(spheroid_radial$pmlc,spheroid_radial$pmlc$trial=="T1"), cols,parameters_radial_data,paste("pMLC ",labelling_radial_data,sep=""),"pMLC_T1_")
    plot_black_bg_histogram_rocki(subset(spheroid_radial$pmlc,spheroid_radial$pmlc$trial=="T2"), cols,parameters_radial_data,paste("pMLC ",labelling_radial_data,sep=""),"pMLC_T2_")
    plot_black_bg_histogram_rocki(subset(spheroid_radial$pmlc,spheroid_radial$pmlc$trial=="T3"), cols,parameters_radial_data,paste("pMLC ",labelling_radial_data,sep=""),"pMLC_T3_")
    
    plot_black_bg_boxplot_rocki(subset(spheroid_radial$pmlc,spheroid_radial$pmlc$trial=="T1"), cols,parameters_radial_data,paste("pMLC ",labelling_radial_data,sep=""),"pMLC_T1_")
    plot_black_bg_boxplot_rocki(subset(spheroid_radial$pmlc,spheroid_radial$pmlc$trial=="T2"), cols,parameters_radial_data,paste("pMLC ",labelling_radial_data,sep=""),"pMLC_T2_")
    plot_black_bg_boxplot_rocki(subset(spheroid_radial$pmlc,spheroid_radial$pmlc$trial=="T3"), cols,parameters_radial_data,paste("pMLC ",labelling_radial_data,sep=""),"pMLC_T3_")
    
    plot_black_bg_barplot_rocki(spheroid_radial$pmlc, cols,parameters_radial_data,paste("pMLC ",labelling_radial_data,sep=""),"pMLC_")
    
    plot_black_bg_histogram_rocki(subset(spheroid_radial$oct4,spheroid_radial$oct4$trial=="T1"), cols,parameters_radial_data,paste("Oct4 ",labelling_radial_data,sep=""),"Oct4_T1_")
    plot_black_bg_histogram_rocki(subset(spheroid_radial$oct4,spheroid_radial$oct4$trial=="T2"), cols,parameters_radial_data,paste("Oct4 ",labelling_radial_data,sep=""),"Oct4_T2_")
    plot_black_bg_histogram_rocki(subset(spheroid_radial$oct4,spheroid_radial$oct4$trial=="T3"), cols,parameters_radial_data,paste("Oct4 ",labelling_radial_data,sep=""),"Oct4_T3_")
    
    plot_black_bg_boxplot_rocki(subset(spheroid_radial$oct4,spheroid_radial$oct4$trial=="T1"), cols,parameters_radial_data,paste("Oct4 ",labelling_radial_data,sep=""),"Oct4_T1_")
    plot_black_bg_boxplot_rocki(subset(spheroid_radial$oct4,spheroid_radial$oct4$trial=="T2"), cols,parameters_radial_data,paste("Oct4 ",labelling_radial_data,sep=""),"Oct4_T2_")
    plot_black_bg_boxplot_rocki(subset(spheroid_radial$oct4,spheroid_radial$oct4$trial=="T3"), cols,parameters_radial_data,paste("Oct4 ",labelling_radial_data,sep=""),"Oct4_T3_")
    
    plot_black_bg_barplot_rocki(spheroid_radial$oct4, cols,parameters_radial_data,paste("Oct4 ",labelling_radial_data,sep=""),"Oct4_")
    
    radial_plot_black_bg_boxplot_rocki(subset(spheroid_radial$dna,spheroid_radial$dna$trial=="T1"),cols,"Mean DNA","T1_dna.png")
    radial_plot_black_bg_boxplot_rocki(subset(spheroid_radial$dna,spheroid_radial$dna$trial=="T2"),cols,"Mean DNA","T2_dna.png")
    radial_plot_black_bg_boxplot_rocki(subset(spheroid_radial$dna,spheroid_radial$dna$trial=="T3"),cols,"Mean DNA","T3_dna.png")
    
    radial_plot_black_bg_boxplot_rocki(subset(spheroid_radial$actin,spheroid_radial$actin$trial=="T1"),cols,"Mean Actin","T1_Actin.png")
    radial_plot_black_bg_boxplot_rocki(subset(spheroid_radial$actin,spheroid_radial$actin$trial=="T2"),cols,"Mean Actin","T2_Actin.png")
    radial_plot_black_bg_boxplot_rocki(subset(spheroid_radial$actin,spheroid_radial$actin$trial=="T3"),cols,"Mean Actin","T3_Actin.png")
    
    radial_plot_black_bg_boxplot_rocki(subset(spheroid_radial$oct4,spheroid_radial$oct4$trial=="T1"),cols,"Mean Oct4","T1_oct4.png")
    radial_plot_black_bg_boxplot_rocki(subset(spheroid_radial$oct4,spheroid_radial$oct4$trial=="T2"),cols,"Mean Oct4","T2_oct4.png")
    radial_plot_black_bg_boxplot_rocki(subset(spheroid_radial$oct4,spheroid_radial$oct4$trial=="T3"),cols,"Mean Oct4","T3_oct4.png")
    
    radial_plot_black_bg_boxplot_rocki(subset(spheroid_radial$oct4,spheroid_radial$oct4$trial=="T1"),cols,"Mean pMLC","T1_PMLC.png")
    radial_plot_black_bg_boxplot_rocki(subset(spheroid_radial$oct4,spheroid_radial$oct4$trial=="T2"),cols,"Mean pMLC","T2_PMLC.png")
    radial_plot_black_bg_boxplot_rocki(subset(spheroid_radial$oct4,spheroid_radial$oct4$trial=="T3"),cols,"Mean pMLC","T3_PMLC.png")
    
  }
  #Spheoid geometeric properties
  {
    dire<-paste(dird,"spheroid_properties/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_black_bg_histogram_rocki(subset(spheroid_combined,spheroid_combined$trial=="T1"), cols,parameters_geometrical_data,paste("Spheroid ",labelling_geometrical_data,sep=""),"geo_T1_")
    plot_black_bg_histogram_rocki(subset(spheroid_combined,spheroid_combined$trial=="T2"), cols,parameters_geometrical_data,paste("Spheroid ",labelling_geometrical_data,sep=""),"geo_T2_")
    plot_black_bg_histogram_rocki(subset(spheroid_combined,spheroid_combined$trial=="T3"), cols,parameters_geometrical_data,paste("Spheroid ",labelling_geometrical_data,sep=""),"geo_T3_")
    
    plot_black_bg_boxplot_rocki(subset(spheroid_combined,spheroid_combined$trial=="T1"), cols,parameters_geometrical_data,labelling_geometrical_data,"geo_T1_")
    plot_black_bg_boxplot_rocki(subset(spheroid_combined,spheroid_combined$trial=="T2"), cols,parameters_geometrical_data,labelling_geometrical_data,"geo_T2_")
    plot_black_bg_boxplot_rocki(subset(spheroid_combined,spheroid_combined$trial=="T3"), cols,parameters_geometrical_data,labelling_geometrical_data,"geo_T3_")
    
    plot_black_bg_barplot_rocki(spheroid_combined, cols,parameters_geometrical_data,paste("Spheroid ",labelling_geometrical_data,sep=""),"geo_")
    
  }
  #Nuclear geometeric properties
  {
    dire<-paste(dird,"nuclear_properties/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_black_bg_histogram_rocki(subset(nucleus_combined,nucleus_combined$trial=="T1"), cols,parameters_geo_texture_data,paste("Nuclear ",labelling_geo_texture_data,sep=""),"geo_T1_")
    plot_black_bg_histogram_rocki(subset(nucleus_combined,nucleus_combined$trial=="T2"), cols,parameters_geo_texture_data,paste("Nuclear ",labelling_geo_texture_data,sep=""),"geo_T2_")
    plot_black_bg_histogram_rocki(subset(nucleus_combined,nucleus_combined$trial=="T3"), cols,parameters_geo_texture_data,paste("Nuclear ",labelling_geo_texture_data,sep=""),"geo_T3_")
    
    plot_black_bg_boxplot_rocki(subset(nucleus_combined,nucleus_combined$trial=="T1"), cols,parameters_geo_texture_data,labelling_geo_texture_data,"geo_T1_")
    plot_black_bg_boxplot_rocki(subset(nucleus_combined,nucleus_combined$trial=="T2"), cols,parameters_geo_texture_data,labelling_geo_texture_data,"geo_T2_")
    plot_black_bg_boxplot_rocki(subset(nucleus_combined,nucleus_combined$trial=="T3"), cols,parameters_geo_texture_data,labelling_geo_texture_data,"geo_T3_")
    
    plot_black_bg_barplot_rocki(nucleus_combined, cols,parameters_geo_texture_data,paste("Nuclear ",labelling_geo_texture_data,sep=""),"geo_")
    
  }
  #Nuclear Proj_Int
  {
    
    dire<-paste(dird,"nuclear_proj_int/",sep="")
    dir.create(dire)
    setwd(dire)
    
    plot_black_bg_histogram_rocki(subset(nuclei_2dint_data$actin,nuclei_2dint_data$actin$trial=="T1"), cols,parameters_2d_int_data,paste("F-Actin ",nuc_labelling_2d_int_data,sep=""),"ACTIN_T1_")
    plot_black_bg_histogram_rocki(subset(nuclei_2dint_data$actin,nuclei_2dint_data$actin$trial=="T2"), cols,parameters_2d_int_data,paste("F-Actin ",nuc_labelling_2d_int_data,sep=""),"ACTIN_T2_")
    plot_black_bg_histogram_rocki(subset(nuclei_2dint_data$actin,nuclei_2dint_data$actin$trial=="T3"), cols,parameters_2d_int_data,paste("F-Actin ",nuc_labelling_2d_int_data,sep=""),"ACTIN_T3_")
    
    plot_black_bg_boxplot_rocki(subset(nuclei_2dint_data$actin,nuclei_2dint_data$actin$trial=="T1"), cols,parameters_2d_int_data,paste("F-Actin ",nuc_labelling_2d_int_data,sep=""),"ACTIN_T1_")
    plot_black_bg_boxplot_rocki(subset(nuclei_2dint_data$actin,nuclei_2dint_data$actin$trial=="T2"), cols,parameters_2d_int_data,paste("F-Actin ",nuc_labelling_2d_int_data,sep=""),"ACTIN_T2_")
    plot_black_bg_boxplot_rocki(subset(nuclei_2dint_data$actin,nuclei_2dint_data$actin$trial=="T3"), cols,parameters_2d_int_data,paste("F-Actin ",nuc_labelling_2d_int_data,sep=""),"ACTIN_T3_")
    
    plot_black_bg_barplot_rocki(nuclei_2dint_data$actin, cols,parameters_2d_int_data,paste("F-Actin ",nuc_labelling_2d_int_data,sep=""),"ACTIN_")
    
    
    plot_black_bg_histogram_rocki(subset(nuclei_2dint_data$pmlc,nuclei_2dint_data$pmlc$trial=="T1"), cols,parameters_2d_int_data,paste("pMLC ",nuc_labelling_2d_int_data,sep=""),"pMLC_T1_")
    plot_black_bg_histogram_rocki(subset(nuclei_2dint_data$pmlc,nuclei_2dint_data$pmlc$trial=="T2"), cols,parameters_2d_int_data,paste("pMLC ",nuc_labelling_2d_int_data,sep=""),"pMLC_T2_")
    plot_black_bg_histogram_rocki(subset(nuclei_2dint_data$pmlc,nuclei_2dint_data$pmlc$trial=="T3"), cols,parameters_2d_int_data,paste("pMLC ",nuc_labelling_2d_int_data,sep=""),"pMLC_T3_")
    
    plot_black_bg_boxplot_rocki(subset(nuclei_2dint_data$pmlc,nuclei_2dint_data$pmlc$trial=="T1"), cols,parameters_2d_int_data,paste("pMLC ",nuc_labelling_2d_int_data,sep=""),"pMLC_T1_")
    plot_black_bg_boxplot_rocki(subset(nuclei_2dint_data$pmlc,nuclei_2dint_data$pmlc$trial=="T2"), cols,parameters_2d_int_data,paste("pMLC ",nuc_labelling_2d_int_data,sep=""),"pMLC_T2_")
    plot_black_bg_boxplot_rocki(subset(nuclei_2dint_data$pmlc,nuclei_2dint_data$pmlc$trial=="T3"), cols,parameters_2d_int_data,paste("pMLC ",nuc_labelling_2d_int_data,sep=""),"pMLC_T3_")
    
    plot_black_bg_barplot_rocki(nuclei_2dint_data$pmlc, cols,parameters_2d_int_data,paste("pMLC ",nuc_labelling_2d_int_data,sep=""),"pMLC_")
    
    plot_black_bg_histogram_rocki(subset(nuclei_2dint_data$oct4,nuclei_2dint_data$oct4$trial=="T1"), cols,parameters_2d_int_data,paste("Oct4 ",nuc_labelling_2d_int_data,sep=""),"Oct4_T1_")
    plot_black_bg_histogram_rocki(subset(nuclei_2dint_data$oct4,nuclei_2dint_data$oct4$trial=="T2"), cols,parameters_2d_int_data,paste("Oct4 ",nuc_labelling_2d_int_data,sep=""),"Oct4_T2_")
    plot_black_bg_histogram_rocki(subset(nuclei_2dint_data$oct4,nuclei_2dint_data$oct4$trial=="T3"), cols,parameters_2d_int_data,paste("Oct4 ",nuc_labelling_2d_int_data,sep=""),"Oct4_T3_")
    
    plot_black_bg_boxplot_rocki(subset(nuclei_2dint_data$oct4,nuclei_2dint_data$oct4$trial=="T1"), cols,parameters_2d_int_data,paste("Oct4 ",nuc_labelling_2d_int_data,sep=""),"Oct4_T1_")
    plot_black_bg_boxplot_rocki(subset(nuclei_2dint_data$oct4,nuclei_2dint_data$oct4$trial=="T2"), cols,parameters_2d_int_data,paste("Oct4 ",nuc_labelling_2d_int_data,sep=""),"Oct4_T2_")
    plot_black_bg_boxplot_rocki(subset(nuclei_2dint_data$oct4,nuclei_2dint_data$oct4$trial=="T3"), cols,parameters_2d_int_data,paste("Oct4 ",nuc_labelling_2d_int_data,sep=""),"Oct4_T3_")
    
    plot_black_bg_barplot_rocki(nuclei_2dint_data$oct4, cols,parameters_2d_int_data,paste("Oct4 ",nuc_labelling_2d_int_data,sep=""),"Oct4_")
    
    
  }
  #Nuclear 3D_Int
  {
    dire<-paste(dird,"nuclear_3d_int/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_black_bg_histogram_rocki(subset(nuclei_3dint_data$dna,nuclei_3dint_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T1_")
    plot_black_bg_histogram_rocki(subset(nuclei_3dint_data$dna,nuclei_3dint_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T2_")
    plot_black_bg_histogram_rocki(subset(nuclei_3dint_data$dna,nuclei_3dint_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T3_")
    
    plot_black_bg_boxplot_rocki(subset(nuclei_3dint_data$dna,nuclei_3dint_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T1_")
    plot_black_bg_boxplot_rocki(subset(nuclei_3dint_data$dna,nuclei_3dint_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T2_")
    plot_black_bg_boxplot_rocki(subset(nuclei_3dint_data$dna,nuclei_3dint_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T3_")
    
    plot_black_bg_barplot_rocki(nuclei_3dint_data$dna, cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_")
    
    plot_black_bg_histogram_rocki(subset(nuclei_3dint_data$actin,nuclei_3dint_data$actin$trial=="T1"), cols,parameters_int_data,paste("F-Actin ",nuclabelling_int_data,sep=""),"ACTIN_T1_")
    plot_black_bg_histogram_rocki(subset(nuclei_3dint_data$actin,nuclei_3dint_data$actin$trial=="T2"), cols,parameters_int_data,paste("F-Actin ",nuclabelling_int_data,sep=""),"ACTIN_T2_")
    plot_black_bg_histogram_rocki(subset(nuclei_3dint_data$actin,nuclei_3dint_data$actin$trial=="T3"), cols,parameters_int_data,paste("F-Actin ",nuclabelling_int_data,sep=""),"ACTIN_T3_")
    
    plot_black_bg_boxplot_rocki(subset(nuclei_3dint_data$actin,nuclei_3dint_data$actin$trial=="T1"), cols,parameters_int_data,paste("F-Actin ",nuclabelling_int_data,sep=""),"ACTIN_T1_")
    plot_black_bg_boxplot_rocki(subset(nuclei_3dint_data$actin,nuclei_3dint_data$actin$trial=="T2"), cols,parameters_int_data,paste("F-Actin ",nuclabelling_int_data,sep=""),"ACTIN_T2_")
    plot_black_bg_boxplot_rocki(subset(nuclei_3dint_data$actin,nuclei_3dint_data$actin$trial=="T3"), cols,parameters_int_data,paste("F-Actin ",nuclabelling_int_data,sep=""),"ACTIN_T3_")
    
    plot_black_bg_barplot_rocki(nuclei_3dint_data$actin, cols,parameters_int_data,paste("F-Actin ",nuclabelling_int_data,sep=""),"ACTIN_")
    
    
    plot_black_bg_histogram_rocki(subset(nuclei_3dint_data$pmlc,nuclei_3dint_data$pmlc$trial=="T1"), cols,parameters_int_data,paste("pMLC ",nuclabelling_int_data,sep=""),"pMLC_T1_")
    plot_black_bg_histogram_rocki(subset(nuclei_3dint_data$pmlc,nuclei_3dint_data$pmlc$trial=="T2"), cols,parameters_int_data,paste("pMLC ",nuclabelling_int_data,sep=""),"pMLC_T2_")
    plot_black_bg_histogram_rocki(subset(nuclei_3dint_data$pmlc,nuclei_3dint_data$pmlc$trial=="T3"), cols,parameters_int_data,paste("pMLC ",nuclabelling_int_data,sep=""),"pMLC_T3_")
    
    plot_black_bg_boxplot_rocki(subset(nuclei_3dint_data$pmlc,nuclei_3dint_data$pmlc$trial=="T1"), cols,parameters_int_data,paste("pMLC ",nuclabelling_int_data,sep=""),"pMLC_T1_")
    plot_black_bg_boxplot_rocki(subset(nuclei_3dint_data$pmlc,nuclei_3dint_data$pmlc$trial=="T2"), cols,parameters_int_data,paste("pMLC ",nuclabelling_int_data,sep=""),"pMLC_T2_")
    plot_black_bg_boxplot_rocki(subset(nuclei_3dint_data$pmlc,nuclei_3dint_data$pmlc$trial=="T3"), cols,parameters_int_data,paste("pMLC ",nuclabelling_int_data,sep=""),"pMLC_T3_")
    
    plot_black_bg_barplot_rocki(nuclei_3dint_data$pmlc, cols,parameters_int_data,paste("pMLC ",nuclabelling_int_data,sep=""),"pMLC_")
    
    plot_black_bg_histogram_rocki(subset(nuclei_3dint_data$oct4,nuclei_3dint_data$oct4$trial=="T1"), cols,parameters_int_data,paste("Oct4 ",nuclabelling_int_data,sep=""),"Oct4_T1_")
    plot_black_bg_histogram_rocki(subset(nuclei_3dint_data$oct4,nuclei_3dint_data$oct4$trial=="T2"), cols,parameters_int_data,paste("Oct4 ",nuclabelling_int_data,sep=""),"Oct4_T2_")
    plot_black_bg_histogram_rocki(subset(nuclei_3dint_data$oct4,nuclei_3dint_data$oct4$trial=="T3"), cols,parameters_int_data,paste("Oct4 ",nuclabelling_int_data,sep=""),"Oct4_T3_")
    
    plot_black_bg_boxplot_rocki(subset(nuclei_3dint_data$oct4,nuclei_3dint_data$oct4$trial=="T1"), cols,parameters_int_data,paste("Oct4 ",nuclabelling_int_data,sep=""),"Oct4_T1_")
    plot_black_bg_boxplot_rocki(subset(nuclei_3dint_data$oct4,nuclei_3dint_data$oct4$trial=="T2"), cols,parameters_int_data,paste("Oct4 ",nuclabelling_int_data,sep=""),"Oct4_T2_")
    plot_black_bg_boxplot_rocki(subset(nuclei_3dint_data$oct4,nuclei_3dint_data$oct4$trial=="T3"), cols,parameters_int_data,paste("Oct4 ",nuclabelling_int_data,sep=""),"Oct4_T3_")
    
    plot_black_bg_barplot_rocki(nuclei_3dint_data$oct4, cols,parameters_int_data,paste("Oct4 ",nuclabelling_int_data,sep=""),"Oct4_")
    
    
  }
  #Cellular 3D_Int
  {
    dire<-paste(dird,"cellular_3d_int/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_black_bg_histogram_rocki(subset(cell_3dint_data$dna,cell_3dint_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T1_")
    plot_black_bg_histogram_rocki(subset(cell_3dint_data$dna,cell_3dint_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T2_")
    plot_black_bg_histogram_rocki(subset(cell_3dint_data$dna,cell_3dint_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T3_")
    
    plot_black_bg_boxplot_rocki(subset(cell_3dint_data$dna,cell_3dint_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T1_")
    plot_black_bg_boxplot_rocki(subset(cell_3dint_data$dna,cell_3dint_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T2_")
    plot_black_bg_boxplot_rocki(subset(cell_3dint_data$dna,cell_3dint_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T3_")
    
    plot_black_bg_barplot_rocki(cell_3dint_data$dna, cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_")
    
    plot_black_bg_histogram_rocki(subset(cell_3dint_data$actin,cell_3dint_data$actin$trial=="T1"), cols,parameters_int_data,paste("F-Actin ",celllabelling_int_data,sep=""),"ACTIN_T1_")
    plot_black_bg_histogram_rocki(subset(cell_3dint_data$actin,cell_3dint_data$actin$trial=="T2"), cols,parameters_int_data,paste("F-Actin ",celllabelling_int_data,sep=""),"ACTIN_T2_")
    plot_black_bg_histogram_rocki(subset(cell_3dint_data$actin,cell_3dint_data$actin$trial=="T3"), cols,parameters_int_data,paste("F-Actin ",celllabelling_int_data,sep=""),"ACTIN_T3_")
    
    plot_black_bg_boxplot_rocki(subset(cell_3dint_data$actin,cell_3dint_data$actin$trial=="T1"), cols,parameters_int_data,paste("F-Actin ",celllabelling_int_data,sep=""),"ACTIN_T1_")
    plot_black_bg_boxplot_rocki(subset(cell_3dint_data$actin,cell_3dint_data$actin$trial=="T2"), cols,parameters_int_data,paste("F-Actin ",celllabelling_int_data,sep=""),"ACTIN_T2_")
    plot_black_bg_boxplot_rocki(subset(cell_3dint_data$actin,cell_3dint_data$actin$trial=="T3"), cols,parameters_int_data,paste("F-Actin ",celllabelling_int_data,sep=""),"ACTIN_T3_")
    
    plot_black_bg_barplot_rocki(cell_3dint_data$actin, cols,parameters_int_data,paste("F-Actin ",celllabelling_int_data,sep=""),"ACTIN_")
    
    
    plot_black_bg_histogram_rocki(subset(cell_3dint_data$pmlc,cell_3dint_data$pmlc$trial=="T1"), cols,parameters_int_data,paste("pMLC ",celllabelling_int_data,sep=""),"pMLC_T1_")
    plot_black_bg_histogram_rocki(subset(cell_3dint_data$pmlc,cell_3dint_data$pmlc$trial=="T2"), cols,parameters_int_data,paste("pMLC ",celllabelling_int_data,sep=""),"pMLC_T2_")
    plot_black_bg_histogram_rocki(subset(cell_3dint_data$pmlc,cell_3dint_data$pmlc$trial=="T3"), cols,parameters_int_data,paste("pMLC ",celllabelling_int_data,sep=""),"pMLC_T3_")
    
    plot_black_bg_boxplot_rocki(subset(cell_3dint_data$pmlc,cell_3dint_data$pmlc$trial=="T1"), cols,parameters_int_data,paste("pMLC ",celllabelling_int_data,sep=""),"pMLC_T1_")
    plot_black_bg_boxplot_rocki(subset(cell_3dint_data$pmlc,cell_3dint_data$pmlc$trial=="T2"), cols,parameters_int_data,paste("pMLC ",celllabelling_int_data,sep=""),"pMLC_T2_")
    plot_black_bg_boxplot_rocki(subset(cell_3dint_data$pmlc,cell_3dint_data$pmlc$trial=="T3"), cols,parameters_int_data,paste("pMLC ",celllabelling_int_data,sep=""),"pMLC_T3_")
    
    plot_black_bg_barplot_rocki(cell_3dint_data$pmlc, cols,parameters_int_data,paste("pMLC ",celllabelling_int_data,sep=""),"pMLC_")
    
    plot_black_bg_histogram_rocki(subset(cell_3dint_data$oct4,cell_3dint_data$oct4$trial=="T1"), cols,parameters_int_data,paste("Oct4 ",celllabelling_int_data,sep=""),"Oct4_T1_")
    plot_black_bg_histogram_rocki(subset(cell_3dint_data$oct4,cell_3dint_data$oct4$trial=="T2"), cols,parameters_int_data,paste("Oct4 ",celllabelling_int_data,sep=""),"Oct4_T2_")
    plot_black_bg_histogram_rocki(subset(cell_3dint_data$oct4,cell_3dint_data$oct4$trial=="T3"), cols,parameters_int_data,paste("Oct4 ",celllabelling_int_data,sep=""),"Oct4_T3_")
    
    plot_black_bg_boxplot_rocki(subset(cell_3dint_data$oct4,cell_3dint_data$oct4$trial=="T1"), cols,parameters_int_data,paste("Oct4 ",celllabelling_int_data,sep=""),"Oct4_T1_")
    plot_black_bg_boxplot_rocki(subset(cell_3dint_data$oct4,cell_3dint_data$oct4$trial=="T2"), cols,parameters_int_data,paste("Oct4 ",celllabelling_int_data,sep=""),"Oct4_T2_")
    plot_black_bg_boxplot_rocki(subset(cell_3dint_data$oct4,cell_3dint_data$oct4$trial=="T3"), cols,parameters_int_data,paste("Oct4 ",celllabelling_int_data,sep=""),"Oct4_T3_")
    
    plot_black_bg_barplot_rocki(cell_3dint_data$oct4, cols,parameters_int_data,paste("Oct4 ",celllabelling_int_data,sep=""),"Oct4_")
    
    
  }
  
  dird<-"E:/HMF3A_reprogramming/R_analysis/Ycompund_treatment/white_bg_plots/"
  dir.create(dird)
  setwd(dird)
  cols<-c("yellow","gold","blue","violet")
  
  #Spheoid Proj_Int
  {
    
    dire<-paste(dird,"spheroid_proj_int/",sep="")
    dir.create(dire)
    setwd(dire)
    
    plot_white_bg_histogram_rocki(subset(spheroid_2dint_data$actin,spheroid_2dint_data$actin$trial=="T1"), cols,parameters_2d_int_data,paste("F-Actin ",labelling_2d_int_data,sep=""),"ACTIN_T1_")
    plot_white_bg_histogram_rocki(subset(spheroid_2dint_data$actin,spheroid_2dint_data$actin$trial=="T2"), cols,parameters_2d_int_data,paste("F-Actin ",labelling_2d_int_data,sep=""),"ACTIN_T2_")
    plot_white_bg_histogram_rocki(subset(spheroid_2dint_data$actin,spheroid_2dint_data$actin$trial=="T3"), cols,parameters_2d_int_data,paste("F-Actin ",labelling_2d_int_data,sep=""),"ACTIN_T3_")
    
    plot_white_bg_boxplot_rocki(subset(spheroid_2dint_data$actin,spheroid_2dint_data$actin$trial=="T1"), cols,parameters_2d_int_data,paste("F-Actin ",labelling_2d_int_data,sep=""),"ACTIN_T1_")
    plot_white_bg_boxplot_rocki(subset(spheroid_2dint_data$actin,spheroid_2dint_data$actin$trial=="T2"), cols,parameters_2d_int_data,paste("F-Actin ",labelling_2d_int_data,sep=""),"ACTIN_T2_")
    plot_white_bg_boxplot_rocki(subset(spheroid_2dint_data$actin,spheroid_2dint_data$actin$trial=="T3"), cols,parameters_2d_int_data,paste("F-Actin ",labelling_2d_int_data,sep=""),"ACTIN_T3_")
    
    plot_white_bg_barplot_rocki(spheroid_2dint_data$actin, cols,parameters_2d_int_data,paste("F-Actin ",labelling_2d_int_data,sep=""),"ACTIN_")
    
    
    plot_white_bg_histogram_rocki(subset(spheroid_2dint_data$pmlc,spheroid_2dint_data$pmlc$trial=="T1"), cols,parameters_2d_int_data,paste("pMLC ",labelling_2d_int_data,sep=""),"pMLC_T1_")
    plot_white_bg_histogram_rocki(subset(spheroid_2dint_data$pmlc,spheroid_2dint_data$pmlc$trial=="T2"), cols,parameters_2d_int_data,paste("pMLC ",labelling_2d_int_data,sep=""),"pMLC_T2_")
    plot_white_bg_histogram_rocki(subset(spheroid_2dint_data$pmlc,spheroid_2dint_data$pmlc$trial=="T3"), cols,parameters_2d_int_data,paste("pMLC ",labelling_2d_int_data,sep=""),"pMLC_T3_")
    
    plot_white_bg_boxplot_rocki(subset(spheroid_2dint_data$pmlc,spheroid_2dint_data$pmlc$trial=="T1"), cols,parameters_2d_int_data,paste("pMLC ",labelling_2d_int_data,sep=""),"pMLC_T1_")
    plot_white_bg_boxplot_rocki(subset(spheroid_2dint_data$pmlc,spheroid_2dint_data$pmlc$trial=="T2"), cols,parameters_2d_int_data,paste("pMLC ",labelling_2d_int_data,sep=""),"pMLC_T2_")
    plot_white_bg_boxplot_rocki(subset(spheroid_2dint_data$pmlc,spheroid_2dint_data$pmlc$trial=="T3"), cols,parameters_2d_int_data,paste("pMLC ",labelling_2d_int_data,sep=""),"pMLC_T3_")
    
    plot_white_bg_barplot_rocki(spheroid_2dint_data$pmlc, cols,parameters_2d_int_data,paste("pMLC ",labelling_2d_int_data,sep=""),"pMLC_")
    
    plot_white_bg_histogram_rocki(subset(spheroid_2dint_data$oct4,spheroid_2dint_data$oct4$trial=="T1"), cols,parameters_2d_int_data,paste("Oct4 ",labelling_2d_int_data,sep=""),"Oct4_T1_")
    plot_white_bg_histogram_rocki(subset(spheroid_2dint_data$oct4,spheroid_2dint_data$oct4$trial=="T2"), cols,parameters_2d_int_data,paste("Oct4 ",labelling_2d_int_data,sep=""),"Oct4_T2_")
    plot_white_bg_histogram_rocki(subset(spheroid_2dint_data$oct4,spheroid_2dint_data$oct4$trial=="T3"), cols,parameters_2d_int_data,paste("Oct4 ",labelling_2d_int_data,sep=""),"Oct4_T3_")
    
    plot_white_bg_boxplot_rocki(subset(spheroid_2dint_data$oct4,spheroid_2dint_data$oct4$trial=="T1"), cols,parameters_2d_int_data,paste("Oct4 ",labelling_2d_int_data,sep=""),"Oct4_T1_")
    plot_white_bg_boxplot_rocki(subset(spheroid_2dint_data$oct4,spheroid_2dint_data$oct4$trial=="T2"), cols,parameters_2d_int_data,paste("Oct4 ",labelling_2d_int_data,sep=""),"Oct4_T2_")
    plot_white_bg_boxplot_rocki(subset(spheroid_2dint_data$oct4,spheroid_2dint_data$oct4$trial=="T3"), cols,parameters_2d_int_data,paste("Oct4 ",labelling_2d_int_data,sep=""),"Oct4_T3_")
    
    plot_white_bg_barplot_rocki(spheroid_2dint_data$oct4, cols,parameters_2d_int_data,paste("Oct4 ",labelling_2d_int_data,sep=""),"Oct4_")
    
    
  }
  #Spheoid 3D_Int
  {
    dire<-paste(dird,"spheroid_3d_int/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_white_bg_histogram_rocki(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T1_")
    plot_white_bg_histogram_rocki(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T2_")
    plot_white_bg_histogram_rocki(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T3_")
    
    plot_white_bg_boxplot_rocki(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T1_")
    plot_white_bg_boxplot_rocki(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T2_")
    plot_white_bg_boxplot_rocki(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T3_")
    
    plot_white_bg_barplot_rocki(spheroid_int_data$dna, cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_")
    
    plot_white_bg_histogram_rocki(subset(spheroid_int_data$actin,spheroid_int_data$actin$trial=="T1"), cols,parameters_int_data,paste("F-Actin ",labelling_int_data,sep=""),"ACTIN_T1_")
    plot_white_bg_histogram_rocki(subset(spheroid_int_data$actin,spheroid_int_data$actin$trial=="T2"), cols,parameters_int_data,paste("F-Actin ",labelling_int_data,sep=""),"ACTIN_T2_")
    plot_white_bg_histogram_rocki(subset(spheroid_int_data$actin,spheroid_int_data$actin$trial=="T3"), cols,parameters_int_data,paste("F-Actin ",labelling_int_data,sep=""),"ACTIN_T3_")
    
    plot_white_bg_boxplot_rocki(subset(spheroid_int_data$actin,spheroid_int_data$actin$trial=="T1"), cols,parameters_int_data,paste("F-Actin ",labelling_int_data,sep=""),"ACTIN_T1_")
    plot_white_bg_boxplot_rocki(subset(spheroid_int_data$actin,spheroid_int_data$actin$trial=="T2"), cols,parameters_int_data,paste("F-Actin ",labelling_int_data,sep=""),"ACTIN_T2_")
    plot_white_bg_boxplot_rocki(subset(spheroid_int_data$actin,spheroid_int_data$actin$trial=="T3"), cols,parameters_int_data,paste("F-Actin ",labelling_int_data,sep=""),"ACTIN_T3_")
    
    plot_white_bg_barplot_rocki(spheroid_int_data$actin, cols,parameters_int_data,paste("F-Actin ",labelling_int_data,sep=""),"ACTIN_")
    
    
    plot_white_bg_histogram_rocki(subset(spheroid_int_data$pmlc,spheroid_int_data$pmlc$trial=="T1"), cols,parameters_int_data,paste("pMLC ",labelling_int_data,sep=""),"pMLC_T1_")
    plot_white_bg_histogram_rocki(subset(spheroid_int_data$pmlc,spheroid_int_data$pmlc$trial=="T2"), cols,parameters_int_data,paste("pMLC ",labelling_int_data,sep=""),"pMLC_T2_")
    plot_white_bg_histogram_rocki(subset(spheroid_int_data$pmlc,spheroid_int_data$pmlc$trial=="T3"), cols,parameters_int_data,paste("pMLC ",labelling_int_data,sep=""),"pMLC_T3_")
    
    plot_white_bg_boxplot_rocki(subset(spheroid_int_data$pmlc,spheroid_int_data$pmlc$trial=="T1"), cols,parameters_int_data,paste("pMLC ",labelling_int_data,sep=""),"pMLC_T1_")
    plot_white_bg_boxplot_rocki(subset(spheroid_int_data$pmlc,spheroid_int_data$pmlc$trial=="T2"), cols,parameters_int_data,paste("pMLC ",labelling_int_data,sep=""),"pMLC_T2_")
    plot_white_bg_boxplot_rocki(subset(spheroid_int_data$pmlc,spheroid_int_data$pmlc$trial=="T3"), cols,parameters_int_data,paste("pMLC ",labelling_int_data,sep=""),"pMLC_T3_")
    
    plot_white_bg_barplot_rocki(spheroid_int_data$pmlc, cols,parameters_int_data,paste("pMLC ",labelling_int_data,sep=""),"pMLC_")
    
    plot_white_bg_histogram_rocki(subset(spheroid_int_data$oct4,spheroid_int_data$oct4$trial=="T1"), cols,parameters_int_data,paste("Oct4 ",labelling_int_data,sep=""),"Oct4_T1_")
    plot_white_bg_histogram_rocki(subset(spheroid_int_data$oct4,spheroid_int_data$oct4$trial=="T2"), cols,parameters_int_data,paste("Oct4 ",labelling_int_data,sep=""),"Oct4_T2_")
    plot_white_bg_histogram_rocki(subset(spheroid_int_data$oct4,spheroid_int_data$oct4$trial=="T3"), cols,parameters_int_data,paste("Oct4 ",labelling_int_data,sep=""),"Oct4_T3_")
    
    plot_white_bg_boxplot_rocki(subset(spheroid_int_data$oct4,spheroid_int_data$oct4$trial=="T1"), cols,parameters_int_data,paste("Oct4 ",labelling_int_data,sep=""),"Oct4_T1_")
    plot_white_bg_boxplot_rocki(subset(spheroid_int_data$oct4,spheroid_int_data$oct4$trial=="T2"), cols,parameters_int_data,paste("Oct4 ",labelling_int_data,sep=""),"Oct4_T2_")
    plot_white_bg_boxplot_rocki(subset(spheroid_int_data$oct4,spheroid_int_data$oct4$trial=="T3"), cols,parameters_int_data,paste("Oct4 ",labelling_int_data,sep=""),"Oct4_T3_")
    
    plot_white_bg_barplot_rocki(spheroid_int_data$oct4, cols,parameters_int_data,paste("Oct4 ",labelling_int_data,sep=""),"Oct4_")
    
    
  }
  #Spheoid angles
  {
    dire<-paste(dird,"spheroid_angles/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_white_bg_histogram_rocki(subset(spheroid_angle$dna,spheroid_angle$dna$trial=="T1"), cols,parameters_angle_data,paste("DNA ",labelling_angle_data,sep=""),"DNA_T1_")
    plot_white_bg_histogram_rocki(subset(spheroid_angle$dna,spheroid_angle$dna$trial=="T2"), cols,parameters_angle_data,paste("DNA ",labelling_angle_data,sep=""),"DNA_T2_")
    plot_white_bg_histogram_rocki(subset(spheroid_angle$dna,spheroid_angle$dna$trial=="T3"), cols,parameters_angle_data,paste("DNA ",labelling_angle_data,sep=""),"DNA_T3_")
    
    plot_white_bg_boxplot_rocki(subset(spheroid_angle$dna,spheroid_angle$dna$trial=="T1"), cols,parameters_angle_data,paste("DNA ",labelling_angle_data,sep=""),"DNA_T1_")
    plot_white_bg_boxplot_rocki(subset(spheroid_angle$dna,spheroid_angle$dna$trial=="T2"), cols,parameters_angle_data,paste("DNA ",labelling_angle_data,sep=""),"DNA_T2_")
    plot_white_bg_boxplot_rocki(subset(spheroid_angle$dna,spheroid_angle$dna$trial=="T3"), cols,parameters_angle_data,paste("DNA ",labelling_angle_data,sep=""),"DNA_T3_")
    
    
    plot_white_bg_histogram_rocki(subset(spheroid_angle$actin,spheroid_angle$actin$trial=="T1"), cols,parameters_angle_data,paste("F-Actin ",labelling_angle_data,sep=""),"ACTIN_T1_")
    plot_white_bg_histogram_rocki(subset(spheroid_angle$actin,spheroid_angle$actin$trial=="T2"), cols,parameters_angle_data,paste("F-Actin ",labelling_angle_data,sep=""),"ACTIN_T2_")
    plot_white_bg_histogram_rocki(subset(spheroid_angle$actin,spheroid_angle$actin$trial=="T3"), cols,parameters_angle_data,paste("F-Actin ",labelling_angle_data,sep=""),"ACTIN_T3_")
    
    plot_white_bg_boxplot_rocki(subset(spheroid_angle$actin,spheroid_angle$actin$trial=="T1"), cols,parameters_angle_data,paste("F-Actin ",labelling_angle_data,sep=""),"ACTIN_T1_")
    plot_white_bg_boxplot_rocki(subset(spheroid_angle$actin,spheroid_angle$actin$trial=="T2"), cols,parameters_angle_data,paste("F-Actin ",labelling_angle_data,sep=""),"ACTIN_T2_")
    plot_white_bg_boxplot_rocki(subset(spheroid_angle$actin,spheroid_angle$actin$trial=="T3"), cols,parameters_angle_data,paste("F-Actin ",labelling_angle_data,sep=""),"ACTIN_T3_")
    
    plot_white_bg_barplot_rocki(spheroid_angle$actin, cols,parameters_angle_data,paste("F-Actin ",labelling_angle_data,sep=""),"ACTIN_")
    
    
    plot_white_bg_histogram_rocki(subset(spheroid_angle$pmlc,spheroid_angle$pmlc$trial=="T1"), cols,parameters_angle_data,paste("pMLC ",labelling_angle_data,sep=""),"pMLC_T1_")
    plot_white_bg_histogram_rocki(subset(spheroid_angle$pmlc,spheroid_angle$pmlc$trial=="T2"), cols,parameters_angle_data,paste("pMLC ",labelling_angle_data,sep=""),"pMLC_T2_")
    plot_white_bg_histogram_rocki(subset(spheroid_angle$pmlc,spheroid_angle$pmlc$trial=="T3"), cols,parameters_angle_data,paste("pMLC ",labelling_angle_data,sep=""),"pMLC_T3_")
    
    plot_white_bg_boxplot_rocki(subset(spheroid_angle$pmlc,spheroid_angle$pmlc$trial=="T1"), cols,parameters_angle_data,paste("pMLC ",labelling_angle_data,sep=""),"pMLC_T1_")
    plot_white_bg_boxplot_rocki(subset(spheroid_angle$pmlc,spheroid_angle$pmlc$trial=="T2"), cols,parameters_angle_data,paste("pMLC ",labelling_angle_data,sep=""),"pMLC_T2_")
    plot_white_bg_boxplot_rocki(subset(spheroid_angle$pmlc,spheroid_angle$pmlc$trial=="T3"), cols,parameters_angle_data,paste("pMLC ",labelling_angle_data,sep=""),"pMLC_T3_")
    
    plot_white_bg_barplot_rocki(spheroid_angle$pmlc, cols,parameters_angle_data,paste("pMLC ",labelling_angle_data,sep=""),"pMLC_")
    
    plot_white_bg_histogram_rocki(subset(spheroid_angle$oct4,spheroid_angle$oct4$trial=="T1"), cols,parameters_angle_data,paste("Oct4 ",labelling_angle_data,sep=""),"Oct4_T1_")
    plot_white_bg_histogram_rocki(subset(spheroid_angle$oct4,spheroid_angle$oct4$trial=="T2"), cols,parameters_angle_data,paste("Oct4 ",labelling_angle_data,sep=""),"Oct4_T2_")
    plot_white_bg_histogram_rocki(subset(spheroid_angle$oct4,spheroid_angle$oct4$trial=="T3"), cols,parameters_angle_data,paste("Oct4 ",labelling_angle_data,sep=""),"Oct4_T3_")
    
    plot_white_bg_boxplot_rocki(subset(spheroid_angle$oct4,spheroid_angle$oct4$trial=="T1"), cols,parameters_angle_data,paste("Oct4 ",labelling_angle_data,sep=""),"Oct4_T1_")
    plot_white_bg_boxplot_rocki(subset(spheroid_angle$oct4,spheroid_angle$oct4$trial=="T2"), cols,parameters_angle_data,paste("Oct4 ",labelling_angle_data,sep=""),"Oct4_T2_")
    plot_white_bg_boxplot_rocki(subset(spheroid_angle$oct4,spheroid_angle$oct4$trial=="T3"), cols,parameters_angle_data,paste("Oct4 ",labelling_angle_data,sep=""),"Oct4_T3_")
    
    plot_white_bg_barplot_rocki(spheroid_angle$oct4, cols,parameters_angle_data,paste("Oct4 ",labelling_angle_data,sep=""),"Oct4_")
    
    
    
    }
  #Spheoid axial distribution
  {
    dire<-paste(dird,"spheroid_axial_distribution/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_white_bg_histogram_rocki(subset(spheroid_axial$dna,spheroid_axial$dna$trial=="T1"), cols,parameters_axial_data,paste("DNA ",labelling_axial_data,sep=""),"DNA_T1_")
    plot_white_bg_histogram_rocki(subset(spheroid_axial$dna,spheroid_axial$dna$trial=="T2"), cols,parameters_axial_data,paste("DNA ",labelling_axial_data,sep=""),"DNA_T2_")
    plot_white_bg_histogram_rocki(subset(spheroid_axial$dna,spheroid_axial$dna$trial=="T3"), cols,parameters_axial_data,paste("DNA ",labelling_axial_data,sep=""),"DNA_T3_")
    
    plot_white_bg_boxplot_rocki(subset(spheroid_axial$dna,spheroid_axial$dna$trial=="T1"), cols,parameters_axial_data,paste("DNA ",labelling_axial_data,sep=""),"DNA_T1_")
    plot_white_bg_boxplot_rocki(subset(spheroid_axial$dna,spheroid_axial$dna$trial=="T2"), cols,parameters_axial_data,paste("DNA ",labelling_axial_data,sep=""),"DNA_T2_")
    plot_white_bg_boxplot_rocki(subset(spheroid_axial$dna,spheroid_axial$dna$trial=="T3"), cols,parameters_axial_data,paste("DNA ",labelling_axial_data,sep=""),"DNA_T3_")
    
    plot_white_bg_barplot_rocki(spheroid_axial$dna, cols,parameters_axial_data,paste("DNA ",labelling_axial_data,sep=""),"DNA_")
    
    plot_white_bg_histogram_rocki(subset(spheroid_axial$actin,spheroid_axial$actin$trial=="T1"), cols,parameters_axial_data,paste("F-Actin ",labelling_axial_data,sep=""),"ACTIN_T1_")
    plot_white_bg_histogram_rocki(subset(spheroid_axial$actin,spheroid_axial$actin$trial=="T2"), cols,parameters_axial_data,paste("F-Actin ",labelling_axial_data,sep=""),"ACTIN_T2_")
    plot_white_bg_histogram_rocki(subset(spheroid_axial$actin,spheroid_axial$actin$trial=="T3"), cols,parameters_axial_data,paste("F-Actin ",labelling_axial_data,sep=""),"ACTIN_T3_")
    
    plot_white_bg_boxplot_rocki(subset(spheroid_axial$actin,spheroid_axial$actin$trial=="T1"), cols,parameters_axial_data,paste("F-Actin ",labelling_axial_data,sep=""),"ACTIN_T1_")
    plot_white_bg_boxplot_rocki(subset(spheroid_axial$actin,spheroid_axial$actin$trial=="T2"), cols,parameters_axial_data,paste("F-Actin ",labelling_axial_data,sep=""),"ACTIN_T2_")
    plot_white_bg_boxplot_rocki(subset(spheroid_axial$actin,spheroid_axial$actin$trial=="T3"), cols,parameters_axial_data,paste("F-Actin ",labelling_axial_data,sep=""),"ACTIN_T3_")
    
    plot_white_bg_barplot_rocki(spheroid_axial$actin, cols,parameters_axial_data,paste("F-Actin ",labelling_axial_data,sep=""),"ACTIN_")
    
    
    plot_white_bg_histogram_rocki(subset(spheroid_axial$pmlc,spheroid_axial$pmlc$trial=="T1"), cols,parameters_axial_data,paste("pMLC ",labelling_axial_data,sep=""),"pMLC_T1_")
    plot_white_bg_histogram_rocki(subset(spheroid_axial$pmlc,spheroid_axial$pmlc$trial=="T2"), cols,parameters_axial_data,paste("pMLC ",labelling_axial_data,sep=""),"pMLC_T2_")
    plot_white_bg_histogram_rocki(subset(spheroid_axial$pmlc,spheroid_axial$pmlc$trial=="T3"), cols,parameters_axial_data,paste("pMLC ",labelling_axial_data,sep=""),"pMLC_T3_")
    
    plot_white_bg_boxplot_rocki(subset(spheroid_axial$pmlc,spheroid_axial$pmlc$trial=="T1"), cols,parameters_axial_data,paste("pMLC ",labelling_axial_data,sep=""),"pMLC_T1_")
    plot_white_bg_boxplot_rocki(subset(spheroid_axial$pmlc,spheroid_axial$pmlc$trial=="T2"), cols,parameters_axial_data,paste("pMLC ",labelling_axial_data,sep=""),"pMLC_T2_")
    plot_white_bg_boxplot_rocki(subset(spheroid_axial$pmlc,spheroid_axial$pmlc$trial=="T3"), cols,parameters_axial_data,paste("pMLC ",labelling_axial_data,sep=""),"pMLC_T3_")
    
    plot_white_bg_barplot_rocki(spheroid_axial$pmlc, cols,parameters_axial_data,paste("pMLC ",labelling_axial_data,sep=""),"pMLC_")
    
    plot_white_bg_histogram_rocki(subset(spheroid_axial$oct4,spheroid_axial$oct4$trial=="T1"), cols,parameters_axial_data,paste("Oct4 ",labelling_axial_data,sep=""),"Oct4_T1_")
    plot_white_bg_histogram_rocki(subset(spheroid_axial$oct4,spheroid_axial$oct4$trial=="T2"), cols,parameters_axial_data,paste("Oct4 ",labelling_axial_data,sep=""),"Oct4_T2_")
    plot_white_bg_histogram_rocki(subset(spheroid_axial$oct4,spheroid_axial$oct4$trial=="T3"), cols,parameters_axial_data,paste("Oct4 ",labelling_axial_data,sep=""),"Oct4_T3_")
    
    plot_white_bg_boxplot_rocki(subset(spheroid_axial$oct4,spheroid_axial$oct4$trial=="T1"), cols,parameters_axial_data,paste("Oct4 ",labelling_axial_data,sep=""),"Oct4_T1_")
    plot_white_bg_boxplot_rocki(subset(spheroid_axial$oct4,spheroid_axial$oct4$trial=="T2"), cols,parameters_axial_data,paste("Oct4 ",labelling_axial_data,sep=""),"Oct4_T2_")
    plot_white_bg_boxplot_rocki(subset(spheroid_axial$oct4,spheroid_axial$oct4$trial=="T3"), cols,parameters_axial_data,paste("Oct4 ",labelling_axial_data,sep=""),"Oct4_T3_")
    
    plot_white_bg_barplot_rocki(spheroid_axial$oct4, cols,parameters_axial_data,paste("Oct4 ",labelling_axial_data,sep=""),"Oct4_")
    
       
  }
  #Spheoid radial distribution
  {
    dire<-paste(dird,"spheroid_radial_distribution/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_white_bg_histogram_rocki(subset(spheroid_radial$dna,spheroid_radial$dna$trial=="T1"), cols,parameters_radial_data,paste("DNA ",labelling_radial_data,sep=""),"DNA_T1_")
    plot_white_bg_histogram_rocki(subset(spheroid_radial$dna,spheroid_radial$dna$trial=="T2"), cols,parameters_radial_data,paste("DNA ",labelling_radial_data,sep=""),"DNA_T2_")
    plot_white_bg_histogram_rocki(subset(spheroid_radial$dna,spheroid_radial$dna$trial=="T3"), cols,parameters_radial_data,paste("DNA ",labelling_radial_data,sep=""),"DNA_T3_")
    
    plot_white_bg_boxplot_rocki(subset(spheroid_radial$dna,spheroid_radial$dna$trial=="T1"), cols,parameters_radial_data,paste("DNA ",labelling_radial_data,sep=""),"DNA_T1_")
    plot_white_bg_boxplot_rocki(subset(spheroid_radial$dna,spheroid_radial$dna$trial=="T2"), cols,parameters_radial_data,paste("DNA ",labelling_radial_data,sep=""),"DNA_T2_")
    plot_white_bg_boxplot_rocki(subset(spheroid_radial$dna,spheroid_radial$dna$trial=="T3"), cols,parameters_radial_data,paste("DNA ",labelling_radial_data,sep=""),"DNA_T3_")
    
    plot_white_bg_barplot_rocki(spheroid_radial$dna, cols,parameters_radial_data,paste("DNA ",labelling_radial_data,sep=""),"DNA_")
    
    plot_white_bg_histogram_rocki(subset(spheroid_radial$actin,spheroid_radial$actin$trial=="T1"), cols,parameters_radial_data,paste("F-Actin ",labelling_radial_data,sep=""),"ACTIN_T1_")
    plot_white_bg_histogram_rocki(subset(spheroid_radial$actin,spheroid_radial$actin$trial=="T2"), cols,parameters_radial_data,paste("F-Actin ",labelling_radial_data,sep=""),"ACTIN_T2_")
    plot_white_bg_histogram_rocki(subset(spheroid_radial$actin,spheroid_radial$actin$trial=="T3"), cols,parameters_radial_data,paste("F-Actin ",labelling_radial_data,sep=""),"ACTIN_T3_")
    
    plot_white_bg_boxplot_rocki(subset(spheroid_radial$actin,spheroid_radial$actin$trial=="T1"), cols,parameters_radial_data,paste("F-Actin ",labelling_radial_data,sep=""),"ACTIN_T1_")
    plot_white_bg_boxplot_rocki(subset(spheroid_radial$actin,spheroid_radial$actin$trial=="T2"), cols,parameters_radial_data,paste("F-Actin ",labelling_radial_data,sep=""),"ACTIN_T2_")
    plot_white_bg_boxplot_rocki(subset(spheroid_radial$actin,spheroid_radial$actin$trial=="T3"), cols,parameters_radial_data,paste("F-Actin ",labelling_radial_data,sep=""),"ACTIN_T3_")
    
    plot_white_bg_barplot_rocki(spheroid_radial$actin, cols,parameters_radial_data,paste("F-Actin ",labelling_radial_data,sep=""),"ACTIN_")
    
    
    plot_white_bg_histogram_rocki(subset(spheroid_radial$pmlc,spheroid_radial$pmlc$trial=="T1"), cols,parameters_radial_data,paste("pMLC ",labelling_radial_data,sep=""),"pMLC_T1_")
    plot_white_bg_histogram_rocki(subset(spheroid_radial$pmlc,spheroid_radial$pmlc$trial=="T2"), cols,parameters_radial_data,paste("pMLC ",labelling_radial_data,sep=""),"pMLC_T2_")
    plot_white_bg_histogram_rocki(subset(spheroid_radial$pmlc,spheroid_radial$pmlc$trial=="T3"), cols,parameters_radial_data,paste("pMLC ",labelling_radial_data,sep=""),"pMLC_T3_")
    
    plot_white_bg_boxplot_rocki(subset(spheroid_radial$pmlc,spheroid_radial$pmlc$trial=="T1"), cols,parameters_radial_data,paste("pMLC ",labelling_radial_data,sep=""),"pMLC_T1_")
    plot_white_bg_boxplot_rocki(subset(spheroid_radial$pmlc,spheroid_radial$pmlc$trial=="T2"), cols,parameters_radial_data,paste("pMLC ",labelling_radial_data,sep=""),"pMLC_T2_")
    plot_white_bg_boxplot_rocki(subset(spheroid_radial$pmlc,spheroid_radial$pmlc$trial=="T3"), cols,parameters_radial_data,paste("pMLC ",labelling_radial_data,sep=""),"pMLC_T3_")
    
    plot_white_bg_barplot_rocki(spheroid_radial$pmlc, cols,parameters_radial_data,paste("pMLC ",labelling_radial_data,sep=""),"pMLC_")
    
    plot_white_bg_histogram_rocki(subset(spheroid_radial$oct4,spheroid_radial$oct4$trial=="T1"), cols,parameters_radial_data,paste("Oct4 ",labelling_radial_data,sep=""),"Oct4_T1_")
    plot_white_bg_histogram_rocki(subset(spheroid_radial$oct4,spheroid_radial$oct4$trial=="T2"), cols,parameters_radial_data,paste("Oct4 ",labelling_radial_data,sep=""),"Oct4_T2_")
    plot_white_bg_histogram_rocki(subset(spheroid_radial$oct4,spheroid_radial$oct4$trial=="T3"), cols,parameters_radial_data,paste("Oct4 ",labelling_radial_data,sep=""),"Oct4_T3_")
    
    plot_white_bg_boxplot_rocki(subset(spheroid_radial$oct4,spheroid_radial$oct4$trial=="T1"), cols,parameters_radial_data,paste("Oct4 ",labelling_radial_data,sep=""),"Oct4_T1_")
    plot_white_bg_boxplot_rocki(subset(spheroid_radial$oct4,spheroid_radial$oct4$trial=="T2"), cols,parameters_radial_data,paste("Oct4 ",labelling_radial_data,sep=""),"Oct4_T2_")
    plot_white_bg_boxplot_rocki(subset(spheroid_radial$oct4,spheroid_radial$oct4$trial=="T3"), cols,parameters_radial_data,paste("Oct4 ",labelling_radial_data,sep=""),"Oct4_T3_")
    
    plot_white_bg_barplot_rocki(spheroid_radial$oct4, cols,parameters_radial_data,paste("Oct4 ",labelling_radial_data,sep=""),"Oct4_")
    
    
  }
  #Spheoid geometeric properties
  {
    dire<-paste(dird,"spheroid_properties/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_white_bg_histogram_rocki(subset(spheroid_combined,spheroid_combined$trial.x=="T1"), cols,parameters_geometrical_data,paste("geo ",labelling_geometrical_data,sep=""),"geo_T1_")
    plot_white_bg_histogram_rocki(subset(spheroid_combined,spheroid_combined$trial.x=="T2"), cols,parameters_geometrical_data,paste("geo ",labelling_geometrical_data,sep=""),"geo_T2_")
    plot_white_bg_histogram_rocki(subset(spheroid_combined,spheroid_combined$trial.x=="T3"), cols,parameters_geometrical_data,paste("geo ",labelling_geometrical_data,sep=""),"geo_T3_")
    
    plot_white_bg_boxplot_rocki(subset(spheroid_combined,spheroid_combined$trial.x=="T1"), cols,parameters_geometrical_data,paste("geo ",labelling_geometrical_data,sep=""),"geo_T1_")
    plot_white_bg_boxplot_rocki(subset(spheroid_combined,spheroid_combined$trial.x=="T2"), cols,parameters_geometrical_data,paste("geo ",labelling_geometrical_data,sep=""),"geo_T1_")
    plot_white_bg_boxplot_rocki(subset(spheroid_combined,spheroid_combined$trial.x=="T3"), cols,parameters_geometrical_data,paste("geo ",labelling_geometrical_data,sep=""),"geo_T1_")
    
    plot_white_bg_barplot_rocki(spheroid_combined, cols,parameters_geometrical_data,paste("geo ",labelling_geometrical_data,sep=""),"geo_")
    
  }
  #Nuclear geometeric properties
  {
    dire<-paste(dird,"nuclear_properties/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_white_bg_histogram_rocki(subset(nucleus_combined,nucleus_combined$trial=="T1"), cols,parameters_geo_texture_data,paste("Nuclear ",labelling_geo_texture_data,sep=""),"geo_T1_")
    plot_white_bg_histogram_rocki(subset(nucleus_combined,nucleus_combined$trial=="T2"), cols,parameters_geo_texture_data,paste("Nuclear ",labelling_geo_texture_data,sep=""),"geo_T2_")
    plot_white_bg_histogram_rocki(subset(nucleus_combined,nucleus_combined$trial=="T3"), cols,parameters_geo_texture_data,paste("Nuclear ",labelling_geo_texture_data,sep=""),"geo_T3_")
    
    plot_white_bg_boxplot_rocki(subset(nucleus_combined,nucleus_combined$trial=="T1"), cols,parameters_geo_texture_data,labelling_geo_texture_data,"geo_T1_")
    plot_white_bg_boxplot_rocki(subset(nucleus_combined,nucleus_combined$trial=="T2"), cols,parameters_geo_texture_data,labelling_geo_texture_data,"geo_T2_")
    plot_white_bg_boxplot_rocki(subset(nucleus_combined,nucleus_combined$trial=="T3"), cols,parameters_geo_texture_data,labelling_geo_texture_data,"geo_T3_")
    
    plot_white_bg_barplot_rocki(nucleus_combined, cols,parameters_geo_texture_data,paste("Nuclear ",labelling_geo_texture_data,sep=""),"geo_")
    
  }
  #Nuclear Proj_Int
  {
    
    dire<-paste(dird,"nuclear_proj_int/",sep="")
    dir.create(dire)
    setwd(dire)
    
    plot_white_bg_histogram_rocki(subset(nuclei_2dint_data$actin,nuclei_2dint_data$actin$trial=="T1"), cols,parameters_2d_int_data,paste("F-Actin ",nuc_labelling_2d_int_data,sep=""),"ACTIN_T1_")
    plot_white_bg_histogram_rocki(subset(nuclei_2dint_data$actin,nuclei_2dint_data$actin$trial=="T2"), cols,parameters_2d_int_data,paste("F-Actin ",nuc_labelling_2d_int_data,sep=""),"ACTIN_T2_")
    plot_white_bg_histogram_rocki(subset(nuclei_2dint_data$actin,nuclei_2dint_data$actin$trial=="T3"), cols,parameters_2d_int_data,paste("F-Actin ",nuc_labelling_2d_int_data,sep=""),"ACTIN_T3_")
    
    plot_white_bg_boxplot_rocki(subset(nuclei_2dint_data$actin,nuclei_2dint_data$actin$trial=="T1"), cols,parameters_2d_int_data,paste("F-Actin ",nuc_labelling_2d_int_data,sep=""),"ACTIN_T1_")
    plot_white_bg_boxplot_rocki(subset(nuclei_2dint_data$actin,nuclei_2dint_data$actin$trial=="T2"), cols,parameters_2d_int_data,paste("F-Actin ",nuc_labelling_2d_int_data,sep=""),"ACTIN_T2_")
    plot_white_bg_boxplot_rocki(subset(nuclei_2dint_data$actin,nuclei_2dint_data$actin$trial=="T3"), cols,parameters_2d_int_data,paste("F-Actin ",nuc_labelling_2d_int_data,sep=""),"ACTIN_T3_")
    
    plot_white_bg_barplot_rocki(nuclei_2dint_data$actin, cols,parameters_2d_int_data,paste("F-Actin ",nuc_labelling_2d_int_data,sep=""),"ACTIN_")
    
    
    plot_white_bg_histogram_rocki(subset(nuclei_2dint_data$pmlc,nuclei_2dint_data$pmlc$trial=="T1"), cols,parameters_2d_int_data,paste("pMLC ",nuc_labelling_2d_int_data,sep=""),"pMLC_T1_")
    plot_white_bg_histogram_rocki(subset(nuclei_2dint_data$pmlc,nuclei_2dint_data$pmlc$trial=="T2"), cols,parameters_2d_int_data,paste("pMLC ",nuc_labelling_2d_int_data,sep=""),"pMLC_T2_")
    plot_white_bg_histogram_rocki(subset(nuclei_2dint_data$pmlc,nuclei_2dint_data$pmlc$trial=="T3"), cols,parameters_2d_int_data,paste("pMLC ",nuc_labelling_2d_int_data,sep=""),"pMLC_T3_")
    
    plot_white_bg_boxplot_rocki(subset(nuclei_2dint_data$pmlc,nuclei_2dint_data$pmlc$trial=="T1"), cols,parameters_2d_int_data,paste("pMLC ",nuc_labelling_2d_int_data,sep=""),"pMLC_T1_")
    plot_white_bg_boxplot_rocki(subset(nuclei_2dint_data$pmlc,nuclei_2dint_data$pmlc$trial=="T2"), cols,parameters_2d_int_data,paste("pMLC ",nuc_labelling_2d_int_data,sep=""),"pMLC_T2_")
    plot_white_bg_boxplot_rocki(subset(nuclei_2dint_data$pmlc,nuclei_2dint_data$pmlc$trial=="T3"), cols,parameters_2d_int_data,paste("pMLC ",nuc_labelling_2d_int_data,sep=""),"pMLC_T3_")
    
    plot_white_bg_barplot_rocki(nuclei_2dint_data$pmlc, cols,parameters_2d_int_data,paste("pMLC ",nuc_labelling_2d_int_data,sep=""),"pMLC_")
    
    plot_white_bg_histogram_rocki(subset(nuclei_2dint_data$oct4,nuclei_2dint_data$oct4$trial=="T1"), cols,parameters_2d_int_data,paste("Oct4 ",nuc_labelling_2d_int_data,sep=""),"Oct4_T1_")
    plot_white_bg_histogram_rocki(subset(nuclei_2dint_data$oct4,nuclei_2dint_data$oct4$trial=="T2"), cols,parameters_2d_int_data,paste("Oct4 ",nuc_labelling_2d_int_data,sep=""),"Oct4_T2_")
    plot_white_bg_histogram_rocki(subset(nuclei_2dint_data$oct4,nuclei_2dint_data$oct4$trial=="T3"), cols,parameters_2d_int_data,paste("Oct4 ",nuc_labelling_2d_int_data,sep=""),"Oct4_T3_")
    
    plot_white_bg_boxplot_rocki(subset(nuclei_2dint_data$oct4,nuclei_2dint_data$oct4$trial=="T1"), cols,parameters_2d_int_data,paste("Oct4 ",nuc_labelling_2d_int_data,sep=""),"Oct4_T1_")
    plot_white_bg_boxplot_rocki(subset(nuclei_2dint_data$oct4,nuclei_2dint_data$oct4$trial=="T2"), cols,parameters_2d_int_data,paste("Oct4 ",nuc_labelling_2d_int_data,sep=""),"Oct4_T2_")
    plot_white_bg_boxplot_rocki(subset(nuclei_2dint_data$oct4,nuclei_2dint_data$oct4$trial=="T3"), cols,parameters_2d_int_data,paste("Oct4 ",nuc_labelling_2d_int_data,sep=""),"Oct4_T3_")
    
    plot_white_bg_barplot_rocki(nuclei_2dint_data$oct4, cols,parameters_2d_int_data,paste("Oct4 ",nuc_labelling_2d_int_data,sep=""),"Oct4_")
    
    
  }
  #Nuclear 3D_Int
  {
    dire<-paste(dird,"nuclear_3d_int/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_white_bg_histogram_rocki(subset(nuclei_3dint_data$dna,nuclei_3dint_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T1_")
    plot_white_bg_histogram_rocki(subset(nuclei_3dint_data$dna,nuclei_3dint_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T2_")
    plot_white_bg_histogram_rocki(subset(nuclei_3dint_data$dna,nuclei_3dint_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T3_")
    
    plot_white_bg_boxplot_rocki(subset(nuclei_3dint_data$dna,nuclei_3dint_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T1_")
    plot_white_bg_boxplot_rocki(subset(nuclei_3dint_data$dna,nuclei_3dint_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T2_")
    plot_white_bg_boxplot_rocki(subset(nuclei_3dint_data$dna,nuclei_3dint_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T3_")
    
    plot_white_bg_barplot_rocki(nuclei_3dint_data$dna, cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_")
    
    plot_white_bg_histogram_rocki(subset(nuclei_3dint_data$actin,nuclei_3dint_data$actin$trial=="T1"), cols,parameters_int_data,paste("F-Actin ",nuclabelling_int_data,sep=""),"ACTIN_T1_")
    plot_white_bg_histogram_rocki(subset(nuclei_3dint_data$actin,nuclei_3dint_data$actin$trial=="T2"), cols,parameters_int_data,paste("F-Actin ",nuclabelling_int_data,sep=""),"ACTIN_T2_")
    plot_white_bg_histogram_rocki(subset(nuclei_3dint_data$actin,nuclei_3dint_data$actin$trial=="T3"), cols,parameters_int_data,paste("F-Actin ",nuclabelling_int_data,sep=""),"ACTIN_T3_")
    
    plot_white_bg_boxplot_rocki(subset(nuclei_3dint_data$actin,nuclei_3dint_data$actin$trial=="T1"), cols,parameters_int_data,paste("F-Actin ",nuclabelling_int_data,sep=""),"ACTIN_T1_")
    plot_white_bg_boxplot_rocki(subset(nuclei_3dint_data$actin,nuclei_3dint_data$actin$trial=="T2"), cols,parameters_int_data,paste("F-Actin ",nuclabelling_int_data,sep=""),"ACTIN_T2_")
    plot_white_bg_boxplot_rocki(subset(nuclei_3dint_data$actin,nuclei_3dint_data$actin$trial=="T3"), cols,parameters_int_data,paste("F-Actin ",nuclabelling_int_data,sep=""),"ACTIN_T3_")
    
    plot_white_bg_barplot_rocki(nuclei_3dint_data$actin, cols,parameters_int_data,paste("F-Actin ",nuclabelling_int_data,sep=""),"ACTIN_")
    
    
    plot_white_bg_histogram_rocki(subset(nuclei_3dint_data$pmlc,nuclei_3dint_data$pmlc$trial=="T1"), cols,parameters_int_data,paste("pMLC ",nuclabelling_int_data,sep=""),"pMLC_T1_")
    plot_white_bg_histogram_rocki(subset(nuclei_3dint_data$pmlc,nuclei_3dint_data$pmlc$trial=="T2"), cols,parameters_int_data,paste("pMLC ",nuclabelling_int_data,sep=""),"pMLC_T2_")
    plot_white_bg_histogram_rocki(subset(nuclei_3dint_data$pmlc,nuclei_3dint_data$pmlc$trial=="T3"), cols,parameters_int_data,paste("pMLC ",nuclabelling_int_data,sep=""),"pMLC_T3_")
    
    plot_white_bg_boxplot_rocki(subset(nuclei_3dint_data$pmlc,nuclei_3dint_data$pmlc$trial=="T1"), cols,parameters_int_data,paste("pMLC ",nuclabelling_int_data,sep=""),"pMLC_T1_")
    plot_white_bg_boxplot_rocki(subset(nuclei_3dint_data$pmlc,nuclei_3dint_data$pmlc$trial=="T2"), cols,parameters_int_data,paste("pMLC ",nuclabelling_int_data,sep=""),"pMLC_T2_")
    plot_white_bg_boxplot_rocki(subset(nuclei_3dint_data$pmlc,nuclei_3dint_data$pmlc$trial=="T3"), cols,parameters_int_data,paste("pMLC ",nuclabelling_int_data,sep=""),"pMLC_T3_")
    
    plot_white_bg_barplot_rocki(nuclei_3dint_data$pmlc, cols,parameters_int_data,paste("pMLC ",nuclabelling_int_data,sep=""),"pMLC_")
    
    plot_white_bg_histogram_rocki(subset(nuclei_3dint_data$oct4,nuclei_3dint_data$oct4$trial=="T1"), cols,parameters_int_data,paste("Oct4 ",nuclabelling_int_data,sep=""),"Oct4_T1_")
    plot_white_bg_histogram_rocki(subset(nuclei_3dint_data$oct4,nuclei_3dint_data$oct4$trial=="T2"), cols,parameters_int_data,paste("Oct4 ",nuclabelling_int_data,sep=""),"Oct4_T2_")
    plot_white_bg_histogram_rocki(subset(nuclei_3dint_data$oct4,nuclei_3dint_data$oct4$trial=="T3"), cols,parameters_int_data,paste("Oct4 ",nuclabelling_int_data,sep=""),"Oct4_T3_")
    
    plot_white_bg_boxplot_rocki(subset(nuclei_3dint_data$oct4,nuclei_3dint_data$oct4$trial=="T1"), cols,parameters_int_data,paste("Oct4 ",nuclabelling_int_data,sep=""),"Oct4_T1_")
    plot_white_bg_boxplot_rocki(subset(nuclei_3dint_data$oct4,nuclei_3dint_data$oct4$trial=="T2"), cols,parameters_int_data,paste("Oct4 ",nuclabelling_int_data,sep=""),"Oct4_T2_")
    plot_white_bg_boxplot_rocki(subset(nuclei_3dint_data$oct4,nuclei_3dint_data$oct4$trial=="T3"), cols,parameters_int_data,paste("Oct4 ",nuclabelling_int_data,sep=""),"Oct4_T3_")
    
    plot_white_bg_barplot_rocki(nuclei_3dint_data$oct4, cols,parameters_int_data,paste("Oct4 ",nuclabelling_int_data,sep=""),"Oct4_")
    
    
  }
  #Cellular 3D_Int
  {
    dire<-paste(dird,"cellular_3d_int/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_white_bg_histogram_rocki(subset(cell_3dint_data$dna,cell_3dint_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T1_")
    plot_white_bg_histogram_rocki(subset(cell_3dint_data$dna,cell_3dint_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T2_")
    plot_white_bg_histogram_rocki(subset(cell_3dint_data$dna,cell_3dint_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T3_")
    
    plot_white_bg_boxplot_rocki(subset(cell_3dint_data$dna,cell_3dint_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T1_")
    plot_white_bg_boxplot_rocki(subset(cell_3dint_data$dna,cell_3dint_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T2_")
    plot_white_bg_boxplot_rocki(subset(cell_3dint_data$dna,cell_3dint_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T3_")
    
    plot_white_bg_barplot_rocki(cell_3dint_data$dna, cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_")
    
    plot_white_bg_histogram_rocki(subset(cell_3dint_data$actin,cell_3dint_data$actin$trial=="T1"), cols,parameters_int_data,paste("F-Actin ",celllabelling_int_data,sep=""),"ACTIN_T1_")
    plot_white_bg_histogram_rocki(subset(cell_3dint_data$actin,cell_3dint_data$actin$trial=="T2"), cols,parameters_int_data,paste("F-Actin ",celllabelling_int_data,sep=""),"ACTIN_T2_")
    plot_white_bg_histogram_rocki(subset(cell_3dint_data$actin,cell_3dint_data$actin$trial=="T3"), cols,parameters_int_data,paste("F-Actin ",celllabelling_int_data,sep=""),"ACTIN_T3_")
    
    plot_white_bg_boxplot_rocki(subset(cell_3dint_data$actin,cell_3dint_data$actin$trial=="T1"), cols,parameters_int_data,paste("F-Actin ",celllabelling_int_data,sep=""),"ACTIN_T1_")
    plot_white_bg_boxplot_rocki(subset(cell_3dint_data$actin,cell_3dint_data$actin$trial=="T2"), cols,parameters_int_data,paste("F-Actin ",celllabelling_int_data,sep=""),"ACTIN_T2_")
    plot_white_bg_boxplot_rocki(subset(cell_3dint_data$actin,cell_3dint_data$actin$trial=="T3"), cols,parameters_int_data,paste("F-Actin ",celllabelling_int_data,sep=""),"ACTIN_T3_")
    
    plot_white_bg_barplot_rocki(cell_3dint_data$actin, cols,parameters_int_data,paste("F-Actin ",celllabelling_int_data,sep=""),"ACTIN_")
    
    
    plot_white_bg_histogram_rocki(subset(cell_3dint_data$pmlc,cell_3dint_data$pmlc$trial=="T1"), cols,parameters_int_data,paste("pMLC ",celllabelling_int_data,sep=""),"pMLC_T1_")
    plot_white_bg_histogram_rocki(subset(cell_3dint_data$pmlc,cell_3dint_data$pmlc$trial=="T2"), cols,parameters_int_data,paste("pMLC ",celllabelling_int_data,sep=""),"pMLC_T2_")
    plot_white_bg_histogram_rocki(subset(cell_3dint_data$pmlc,cell_3dint_data$pmlc$trial=="T3"), cols,parameters_int_data,paste("pMLC ",celllabelling_int_data,sep=""),"pMLC_T3_")
    
    plot_white_bg_boxplot_rocki(subset(cell_3dint_data$pmlc,cell_3dint_data$pmlc$trial=="T1"), cols,parameters_int_data,paste("pMLC ",celllabelling_int_data,sep=""),"pMLC_T1_")
    plot_white_bg_boxplot_rocki(subset(cell_3dint_data$pmlc,cell_3dint_data$pmlc$trial=="T2"), cols,parameters_int_data,paste("pMLC ",celllabelling_int_data,sep=""),"pMLC_T2_")
    plot_white_bg_boxplot_rocki(subset(cell_3dint_data$pmlc,cell_3dint_data$pmlc$trial=="T3"), cols,parameters_int_data,paste("pMLC ",celllabelling_int_data,sep=""),"pMLC_T3_")
    
    plot_white_bg_barplot_rocki(cell_3dint_data$pmlc, cols,parameters_int_data,paste("pMLC ",celllabelling_int_data,sep=""),"pMLC_")
    
    plot_white_bg_histogram_rocki(subset(cell_3dint_data$oct4,cell_3dint_data$oct4$trial=="T1"), cols,parameters_int_data,paste("Oct4 ",celllabelling_int_data,sep=""),"Oct4_T1_")
    plot_white_bg_histogram_rocki(subset(cell_3dint_data$oct4,cell_3dint_data$oct4$trial=="T2"), cols,parameters_int_data,paste("Oct4 ",celllabelling_int_data,sep=""),"Oct4_T2_")
    plot_white_bg_histogram_rocki(subset(cell_3dint_data$oct4,cell_3dint_data$oct4$trial=="T3"), cols,parameters_int_data,paste("Oct4 ",celllabelling_int_data,sep=""),"Oct4_T3_")
    
    plot_white_bg_boxplot_rocki(subset(cell_3dint_data$oct4,cell_3dint_data$oct4$trial=="T1"), cols,parameters_int_data,paste("Oct4 ",celllabelling_int_data,sep=""),"Oct4_T1_")
    plot_white_bg_boxplot_rocki(subset(cell_3dint_data$oct4,cell_3dint_data$oct4$trial=="T2"), cols,parameters_int_data,paste("Oct4 ",celllabelling_int_data,sep=""),"Oct4_T2_")
    plot_white_bg_boxplot_rocki(subset(cell_3dint_data$oct4,cell_3dint_data$oct4$trial=="T3"), cols,parameters_int_data,paste("Oct4 ",celllabelling_int_data,sep=""),"Oct4_T3_")
    
    plot_white_bg_barplot_rocki(cell_3dint_data$oct4, cols,parameters_int_data,paste("Oct4 ",celllabelling_int_data,sep=""),"Oct4_")
    
    
  }
}

#write the collected data
{
  dir.create("E:/HMF3A_reprogramming/R_analysis/Ycompund_treatment/combined_data/")
  setwd("E:/HMF3A_reprogramming/R_analysis/Ycompund_treatment/combined_data/")
  
  write.csv(spheroid_combined, file="spheroid_combined.csv")
  write.csv(nucleus_combined, file="nucleus_combined.csv")
  
  write.csv(spheroid_radial$actin, file="spheroid_radial_dist_actin.csv")
  write.csv(spheroid_radial$oct4, file="spheroid_radial_dist_oct4.csv")
  write.csv(spheroid_radial$pmlc, file="spheroid_radial_dist_pmlc.csv")
  write.csv(spheroid_radial$dna, file="spheroid_radial_dist_dna.csv")
  
  write.csv(spheroid_2dint_data$actin, file="spheroid_2d_int_actin.csv")
  write.csv(spheroid_2dint_data$oct4, file="spheroid_2d_int_oct4.csv")
  write.csv(spheroid_2dint_data$pmlc, file="spheroid_2d_int_pmlc.csv")
  
  write.csv(spheroid_axial$actin, file="spheroid_axial_dist_actin.csv")
  write.csv(spheroid_axial$oct4, file="spheroid_axial_dist_oct4.csv")
  write.csv(spheroid_axial$pmlc, file="spheroid_axial_dist_pmlc.csv")
  write.csv(spheroid_axial$dna, file="spheroid_axial_dist_dna.csv")
  
  write.csv(spheroid_angle$actin, file="spheroid_angle_dist_actin.csv")
  write.csv(spheroid_angle$oct4, file="spheroid_angle_dist_oct4.csv")
  write.csv(spheroid_angle$pmlc, file="spheroid_angle_dist_pmlc.csv")
  write.csv(spheroid_angle$dna, file="spheroid_angle_dist_dna.csv")
  
  write.csv(spheroid_int_data$actin, file="spheroid_int_actin.csv")
  write.csv(spheroid_int_data$oct4, file="spheroid_int_oct4.csv")
  write.csv(spheroid_int_data$pmlc, file="spheroid_int_pmlc.csv")
  write.csv(spheroid_int_data$dna, file="spheroid_int_dna.csv")
  
  write.csv(nuclei_2dint_data$actin, file="nuclear_2d_int_actin.csv")
  write.csv(nuclei_2dint_data$oct4, file="nuclear_2d_int_oct4.csv")
  write.csv(nuclei_2dint_data$pmlc, file="nuclear_2d_int_pmlc.csv")
  
  write.csv(nuclei_3dint_data$actin, file="nuclear_int_actin.csv")
  write.csv(nuclei_3dint_data$oct4, file="nuclear_int_oct4.csv")
  write.csv(nuclei_3dint_data$pmlc, file="nuclear_int_pmlc.csv")
  write.csv(nuclei_3dint_data$dna, file="nuclear_int_dna.csv")
  
  write.csv(cell_3dint_data$actin, file="cellular_int_actin.csv")
  write.csv(cell_3dint_data$oct4, file="cellular_int_oct4.csv")
  write.csv(cell_3dint_data$pmlc, file="cellular_int_pmlc.csv")
  write.csv(cell_3dint_data$dna, file="cellular_int_dna.csv")
  
}