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
                                  "DCMax..unit.","DCMean..unit.","DCSD..unit.","FeretAngle",
                                  "CX..unit.","CY..unit.","CZ..unit.",
                                  "Area","Perim.","Major","Minor","Circ.","Feret","MinFeret","AR","Round","Solidity")
  
  labelling_geo_texture_data<-c("Bounding Box Ratio","Volume", "Surface Area","Compaction","Sphericity",
                                "Feret (3D)", "Ellipsoid R1","Ellipsoid Ellongation","Ellipsoid Flatness", "Ellipsoid Ratio",
                                "Moment1","Moment2","Moment3","Moment4","Moment5", 
                                "DC_Min.","DC_Max","DC_Mean.","DC_SD","FeretAngle",
                                "CX","CY","CZ",
                                "Proj. Area", "Proj. Perimeter","Proj. Major Axis", "Proj. Minor Axis", "Proj. Circularity",
                                "Proj. Feret","Proj. Min Feret","Proj. A.R","Proj. Roundness","Proj. Solidity")
  
}

source("E:/HMF3A_reprogramming/R_analysis/functions.R")
# initialise the experiment details and the subdirectories required
{
  path_to_ex<-c("E:/HMF3A_reprogramming/Timeline_Edu/20200123_hmf_edu_T1_high_res/",
                "E:/HMF3A_reprogramming/Timeline_Edu/20200123_hmf_edu_T1_high_res/",
                "E:/HMF3A_reprogramming/Timeline_Edu/20200130_hmf_edu_T2_high_res/")
  
  exp<-c("20200123_hmf_edu_T1_high_res","20200123_hmf_edu_T1_high_res","20200130_hmf_edu_T2_high_res")
  
  data_types<-c("/3D geometrical data spheroid/",
                "/3D int_data spheroid/DNA/","/3D int_data spheroid/EDU/")
  data_types_2d<-c("/2D_measures_spheroid/2D_spheroid_DNA.csv","/2D_measures_spheroid/2D_spheroid_EDU.csv")
  nuc_data_types<-c("/3D geometrical data/","/3D ellipsoid/","/3D geometerical_simple/","/3D shape measure/",
                    "/3D int_data/DNA/","/3D int_data/EDU/",
                    "/cell_2microns_measure/EDU/","/cell_2microns_measure/DNA/")
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
                       "/2D_measures_nuclei/2D_nucleusEDU.csv")
  
  samples<-c("D2","D4","D6","D8")
  trial<-c("T1","T2","T3")
  
  dird_apo<-"E:/HMF3A_reprogramming/R_analysis/Edu/"
  dir.create(dird_apo)
  dird_Edu<-paste(dird_apo,"spheroid_shells_Edu", sep="")
  dird_dna<-paste(dird_apo,"spheroid_shells_dna", sep="")
  dir.create(dird_Edu)
  dir.create(dird_dna)
}

# read in all the required features for spheroids
{
  plot(0:1,0:1)
  #T1
  {
    geometrical_data<-combine_sample_sets(path_to_ex[1],samples,data_types[1],"tsv",exp[1],trial[1])
    T1DNA_int_data<-combine_sample_sets(path_to_ex[1],samples,data_types[2],"tsv",exp[1],trial[1])
    T1EDU_int_data<-combine_sample_sets(path_to_ex[1],samples,data_types[3],"tsv",exp[1],trial[1])
    geo_int_2D_data<-combine_sample_sets_2d(path_to_ex[1],samples,data_types_2d[1],exp[1],trial[1])
    
    T1_2d_data_EDU<-combine_sample_sets_2d(path_to_ex[1],samples,data_types_2d[2],exp[1],trial[1])
    
    shell_EDU_T1<-shell_analsysis(path_to_ex[1],samples,"/Shells/EDU",exp[1],trial[1],dird_Edu)
    shell_dna_T1<-shell_analsysis(path_to_ex[1],samples,"/Shells/DNA",exp[1],trial[1],dird_dna)
    
    T1_dna_axial<-axial_analysis(path_to_ex[1],samples,exp[1],trial[1],"/DNA/")
    T1_EDU_axial<-axial_analysis(path_to_ex[1],samples,exp[1],trial[1],"/EDU/")
    
    T1_EDU_angle<-angle_analsysis(path_to_ex[1],samples,"/EDU",exp[1],trial[1])
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
    T2EDU_int_data<-combine_sample_sets(path_to_ex[2],samples,data_types[3],"tsv",exp[2],trial[2])
    geo_int_2D_data<-combine_sample_sets_2d(path_to_ex[2],samples,data_types_2d[1],exp[2],trial[2])
    
    T2_2d_data_EDU<-combine_sample_sets_2d(path_to_ex[2],samples,data_types_2d[2],exp[2],trial[2])
    
    shell_EDU_T2<-shell_analsysis(path_to_ex[2],samples,"/Shells/EDU",exp[2],trial[2],dird_Edu)
    shell_dna_T2<-shell_analsysis(path_to_ex[2],samples,"/Shells/DNA",exp[2],trial[2],dird_dna)
    
    T2_dna_axial<-axial_analysis(path_to_ex[2],samples,exp[2],trial[2],"/DNA/")
    T2_EDU_axial<-axial_analysis(path_to_ex[2],samples,exp[2],trial[2],"/EDU/")
    
    T2_EDU_angle<-angle_analsysis(path_to_ex[2],samples,"/EDU",exp[2],trial[2])
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
    T3EDU_int_data<-combine_sample_sets(path_to_ex[3],samples,data_types[3],"tsv",exp[3],trial[3])
    geo_int_2D_data<-combine_sample_sets_2d(path_to_ex[3],samples,data_types_2d[1],exp[3],trial[3])
    
    T3_2d_data_EDU<-combine_sample_sets_2d(path_to_ex[3],samples,data_types_2d[2],exp[3],trial[3])
    
    shell_EDU_T3<-shell_analsysis(path_to_ex[3],samples,"/Shells/EDU",exp[3],trial[3],dird_Edu)
    shell_dna_T3<-shell_analsysis(path_to_ex[3],samples,"/Shells/DNA",exp[3],trial[3],dird_dna)
    
    T3_dna_axial<-axial_analysis(path_to_ex[3],samples,exp[3],trial[3],"/DNA/")
    T3_EDU_axial<-axial_analysis(path_to_ex[3],samples,exp[3],trial[3],"/EDU/")
    
    T3_EDU_angle<-angle_analsysis(path_to_ex[3],samples,"/EDU",exp[3],trial[3])
    T3_dna_angle<-angle_analsysis(path_to_ex[3],samples,"/DNA",exp[3],trial[3])
    
    geo_int_2D_data$Label<-substring(geo_int_2D_data$Label,5, nchar(geo_int_2D_data$Label)-4)
    geometrical_data$Label<-substring(geometrical_data$Label,0, nchar(geometrical_data$Label)-1)
    T3combined<-merge(geometrical_data,geo_int_2D_data, by="Label")
    
    rm(geometrical_data,geo_int_2D_data)
    
  }
  
  spheroid_2dint_data_EDU<-rbind(T1_2d_data_EDU,T2_2d_data_EDU,T3_2d_data_EDU)
  spheroid_2dint_data_EDU$Label<-substring(spheroid_2dint_data_EDU$Label,5,(nchar(spheroid_2dint_data_EDU$Label)-4))
  spheroid_2dint_data_EDU$Spheroid_id<-paste(spheroid_2dint_data_EDU$trial,spheroid_2dint_data_EDU$sample,
                                             substring(spheroid_2dint_data_EDU$Label,1,2),sep="_")
  
  spheroid_int_data_EDU<-rbind(T1EDU_int_data,T2EDU_int_data,T3EDU_int_data)
  spheroid_int_data_EDU$Label<-substring(spheroid_int_data_EDU$Label,1,(nchar(spheroid_int_data_EDU$Label)-1))
  spheroid_int_data_EDU$Spheroid_id<-paste(spheroid_int_data_EDU$trial,spheroid_int_data_EDU$sample,
                                           substring(spheroid_int_data_EDU$Label,1,2),sep="_")
  
  spheroid_int_data_dna<-rbind(T1DNA_int_data,T2DNA_int_data,T3DNA_int_data)
  spheroid_int_data_dna$Label<-substring(spheroid_int_data_dna$Label,1,(nchar(spheroid_int_data_dna$Label)-1))
  spheroid_int_data_dna$Spheroid_id<-paste(spheroid_int_data_dna$trial,spheroid_int_data_dna$sample,
                                           substring(spheroid_int_data_dna$Label,1,2),sep="_")
  
  spheroid_radial_EDU<-rbind(shell_EDU_T1,shell_EDU_T2,shell_EDU_T3)
  spheroid_radial_EDU$Spheroid_id<-paste(spheroid_radial_EDU$trial,spheroid_radial_EDU$sample,
                                         substring(spheroid_radial_EDU$Label,1,2),sep="_")
  
  spheroid_radial_dna<-rbind(shell_dna_T1,shell_dna_T2,shell_dna_T3)
  spheroid_radial_dna$Spheroid_id<-paste(spheroid_radial_dna$trial,spheroid_radial_dna$sample,
                                         substring(spheroid_radial_dna$Label,1,2),sep="_")
  
  spheroid_axial_dna<-rbind(T1_dna_axial,T2_dna_axial,T3_dna_axial)
  spheroid_axial_dna$Spheroid_id<-paste(spheroid_axial_dna$trial,spheroid_axial_dna$sample,
                                        substring(spheroid_axial_dna$Label,1,2),sep="_")
  
  spheroid_axial_EDU<-rbind(T1_EDU_axial,T2_EDU_axial,T3_EDU_axial)
  spheroid_axial_EDU$Spheroid_id<-paste(spheroid_axial_EDU$trial,spheroid_axial_EDU$sample,
                                        substring(spheroid_axial_EDU$Label,1,2),sep="_")
  
  spheroid_angle_EDU<-rbind(T1_EDU_angle,T2_EDU_angle,T3_EDU_angle)
  spheroid_angle_EDU$Spheroid_id<-paste(spheroid_angle_EDU$trial,spheroid_angle_EDU$sample,
                                        substring(spheroid_angle_EDU$Label,1,2),sep="_")
  
  spheroid_angle_dna<-rbind(T1_dna_angle,T2_dna_angle,T3_dna_angle)
  spheroid_angle_dna$Spheroid_id<-paste(spheroid_angle_dna$trial,spheroid_angle_dna$sample,
                                        substring(spheroid_angle_dna$Label,1,2),sep="_")
  
  spheroid_combined<-rbind(T1combined,T2combined,T3combined)
  colnames(spheroid_combined)[86:87]<-c("sample","trial")
  spheroid_combined$Spheroid_id<-paste(spheroid_combined$trial,spheroid_combined$sample,
                                       substring(spheroid_combined$Label,1,2),sep="_")
  
  spheroid_int_data<-list( EDU=spheroid_int_data_EDU,
                           dna=spheroid_int_data_dna)
  spheroid_radial<-list(EDU=spheroid_radial_EDU,
                        dna=spheroid_radial_dna)
  spheroid_axial<-list( EDU=spheroid_axial_EDU,
                        dna=spheroid_axial_dna)
  spheroid_angle<-list(EDU=spheroid_angle_EDU,
                       dna=spheroid_angle_dna)
  spheroid_2dint_data<-list( EDU=spheroid_2dint_data_EDU)
  
  
  rm(T1combined,T2combined,T3combined,
     T1EDU_int_data,T2EDU_int_data,T3EDU_int_data,
     shell_EDU_T1,shell_EDU_T2,shell_EDU_T3,
     shell_dna_T1,shell_dna_T2,shell_dna_T3,
     T1_dna_axial,T2_dna_axial,T3_dna_axial,
     T1_EDU_axial,T2_EDU_axial,T3_EDU_axial,
     T1_EDU_angle,T2_EDU_angle,T3_EDU_angle,
     T1_dna_angle,T2_dna_angle,T3_dna_angle,
     T1DNA_int_data,T2DNA_int_data,T3DNA_int_data,
     spheroid_angle_EDU,spheroid_angle_dna,
     spheroid_axial_EDU,spheroid_axial_dna,
     spheroid_radial_EDU,spheroid_radial_dna,
     spheroid_int_data_EDU,spheroid_int_data_dna,
     T1_2d_data_EDU,T2_2d_data_EDU,T3_2d_data_EDU,
     spheroid_2dint_data_EDU)
  
}

# read in all the required features for nuclei
{
  #T1
  {
    T1_geometrical_data<-combine_sample_sets(path_to_ex[1],samples,nuc_data_types[1],"tsv",exp[1],trial[1])
    T1_DNA_int_data<-combine_sample_sets(path_to_ex[1],samples,nuc_data_types[5],"tsv",exp[1],trial[1])
    T1_geo_int_2D_data<-combine_sample_sets_2d(path_to_ex[1],samples,nuc_data_types_2d[1],exp[1],trial[1])
    
    T1_EDU_int_data<-combine_sample_sets(path_to_ex[1],samples,nuc_data_types[6],"tsv",exp[1],trial[1])
    
    T1_2d_EDU_int_data<-combine_sample_sets_2d(path_to_ex[1],samples,nuc_data_types_2d[28],exp[1],trial[1])
    
    
    T1_geo_int_2D_data$Label<-substring(T1_geo_int_2D_data$Label,5, nchar(T1_geo_int_2D_data$Label)-4)
    T1_geometrical_data$Label<-substring(T1_geometrical_data$Label,0, nchar(T1_geometrical_data$Label)-1)
    T1_DNA_int_data$Label<-substring(T1_DNA_int_data$Label,0, nchar(T1_DNA_int_data$Label)-1)
    T1_2d_EDU_int_data$Label<-substring(T1_2d_EDU_int_data$Label,5, nchar(T1_2d_EDU_int_data$Label)-4)
    T1_EDU_int_data$Label<-substring(T1_EDU_int_data$Label,0, nchar(T1_EDU_int_data$Label)-1)
    
    T1_combined<-merge(T1_geometrical_data,T1_geo_int_2D_data[,2:36], by="Label")
    T1_combined<-merge(T1_combined,T1_DNA_int_data[,c(3:5,12:18)], by="Label")
    
    rm(T1_geometrical_data,T1_geo_int_2D_data)
  }
  #T2
  {
    T2_geometrical_data<-combine_sample_sets(path_to_ex[2],samples,nuc_data_types[1],"tsv",exp[2],trial[2])
    T2_DNA_int_data<-combine_sample_sets(path_to_ex[2],samples,nuc_data_types[5],"tsv",exp[2],trial[2])
    T2_geo_int_2D_data<-combine_sample_sets_2d(path_to_ex[2],samples,nuc_data_types_2d[1],exp[2],trial[2])
    
    T2_EDU_int_data<-combine_sample_sets(path_to_ex[2],samples,nuc_data_types[6],"tsv",exp[2],trial[2])
    
    T2_2d_EDU_int_data<-combine_sample_sets_2d(path_to_ex[2],samples,nuc_data_types_2d[28],exp[2],trial[2])
    
    T2_geo_int_2D_data$Label<-substring(T2_geo_int_2D_data$Label,5, nchar(T2_geo_int_2D_data$Label)-4)
    T2_geometrical_data$Label<-substring(T2_geometrical_data$Label,0, nchar(T2_geometrical_data$Label)-1)
    T2_2d_EDU_int_data$Label<-substring(T2_2d_EDU_int_data$Label,5, nchar(T2_2d_EDU_int_data$Label)-4)
    T2_DNA_int_data$Label<-substring(T2_DNA_int_data$Label,0, nchar(T2_DNA_int_data$Label)-1)
    T2_EDU_int_data$Label<-substring(T2_EDU_int_data$Label,0, nchar(T2_EDU_int_data$Label)-1)
    
    T2_combined<-merge(T2_geometrical_data,T2_geo_int_2D_data[,2:36], by="Label")
    T2_combined<-merge(T2_combined,T2_DNA_int_data[,c(3:5,12:18)], by="Label")
    
    rm(T2_geometrical_data,T2_geo_int_2D_data)
  }
  #T3
  {
    T3_geometrical_data<-combine_sample_sets(path_to_ex[3],samples,nuc_data_types[1],"tsv",exp[3],trial[3])
    T3_DNA_int_data<-combine_sample_sets(path_to_ex[3],samples,nuc_data_types[5],"tsv",exp[3],trial[3])
    T3_geo_int_2D_data<-combine_sample_sets_2d(path_to_ex[3],samples,nuc_data_types_2d[1],exp[3],trial[3])
    
    T3_EDU_int_data<-combine_sample_sets(path_to_ex[3],samples,nuc_data_types[6],"tsv",exp[3],trial[3])
    
    T3_2d_EDU_int_data<-combine_sample_sets_2d(path_to_ex[3],samples,nuc_data_types_2d[28],exp[3],trial[3])
    
    T3_geo_int_2D_data$Label<-substring(T3_geo_int_2D_data$Label,5, nchar(T3_geo_int_2D_data$Label)-4)
    T3_geometrical_data$Label<-substring(T3_geometrical_data$Label,0, nchar(T3_geometrical_data$Label)-1)
    T3_2d_EDU_int_data$Label<-substring(T3_2d_EDU_int_data$Label,5, nchar(T3_2d_EDU_int_data$Label)-4)
    T3_DNA_int_data$Label<-substring(T3_DNA_int_data$Label,0, nchar(T3_DNA_int_data$Label)-1)
    T3_EDU_int_data$Label<-substring(T3_EDU_int_data$Label,0, nchar(T3_EDU_int_data$Label)-1)
    
    T3_combined<-merge(T3_geometrical_data,T3_geo_int_2D_data[,2:36], by="Label")
    T3_combined<-merge(T3_combined,T3_DNA_int_data[,c(3:5,12:18)], by="Label")
    
    rm(T3_geometrical_data,T3_geo_int_2D_data)
  }
  
  nuclei_2dint_data_EDU<-rbind(T1_2d_EDU_int_data,T2_2d_EDU_int_data,T3_2d_EDU_int_data)
  nuclei_2dint_data_EDU$Spheroid_id<-paste(nuclei_2dint_data_EDU$trial,nuclei_2dint_data_EDU$sample,
                                           substring(nuclei_2dint_data_EDU$Label,1,2),sep="_")
  spheroid_nuclei_2dint_data_EDU<-nuc_median_spheroid_levels(nuclei_2dint_data_EDU,parameters_2d_int_data)
  
  nuclei_int_data_EDU<-rbind(T1_EDU_int_data,T2_EDU_int_data,T3_EDU_int_data)
  nuclei_int_data_EDU$Spheroid_id<-paste(nuclei_int_data_EDU$trial,nuclei_int_data_EDU$sample,
                                         substring(nuclei_int_data_EDU$Label,1,2),sep="_")
  spheroid_nuclei_int_data_EDU<-nuc_median_spheroid_levels(nuclei_int_data_EDU,parameters_int_data)
  
  nuclei_int_data_dna<-rbind(T1_DNA_int_data,T2_DNA_int_data,T3_DNA_int_data)
  nuclei_int_data_dna$Spheroid_id<-paste(nuclei_int_data_dna$trial,nuclei_int_data_dna$sample,
                                         substring(nuclei_int_data_dna$Label,1,2),sep="_")
  spheroid_nuclei_int_data_dna<-nuc_median_spheroid_levels(nuclei_int_data_dna,parameters_int_data)
  
  
  nucleus_combined<-rbind(T1_combined,T2_combined,T3_combined)
  nucleus_combined$Spheroid_id<-paste(nucleus_combined$trial,nucleus_combined$sample,
                                      substring(nucleus_combined$Label,1,2),sep="_")
  spheroid_nucleus_combined<-nuc_median_spheroid_levels(nucleus_combined,parameters_geo_texture_data)
  
  
  nuclei_2dint_data<-list(EDU=nuclei_2dint_data_EDU)
  nuclei_3dint_data<-list(EDU=nuclei_int_data_EDU,dna=nuclei_int_data_dna)
  
  spheroid_nuclei_2dint_data<-list(EDU=spheroid_nuclei_2dint_data_EDU)
  spheroid_nuclei_3dint_data<-list(EDU=spheroid_nuclei_int_data_EDU,dna=spheroid_nuclei_int_data_dna)
  
  
  rm(T1_2d_EDU_int_data,T2_2d_EDU_int_data,T3_2d_EDU_int_data,
     T1_EDU_int_data,T2_EDU_int_data,T3_EDU_int_data,
     T1_DNA_int_data,T2_DNA_int_data,T3_DNA_int_data,
     T1_combined,T2_combined,T3_combined,
     nuclei_2dint_data_EDU, nuclei_int_data_EDU,
     nuclei_int_data_dna,
     spheroid_nuclei_2dint_data_EDU,spheroid_nuclei_int_data_EDU,spheroid_nuclei_int_data_dna)
  
}

#plotting
{
  dird<-"E:/HMF3A_reprogramming/R_analysis/EDU/black_bg_plots/"
  dir.create(dird)
  setwd(dird)
  cols<-colorRampPalette(c("red","yellow"))( 5)
  #Spheoid Proj_Int
  {
    
    dire<-paste(dird,"spheroid_proj_int/",sep="")
    dir.create(dire)
    setwd(dire)
    
    
    plot_black_bg_histogram_timeline(subset(spheroid_2dint_data$EDU,spheroid_2dint_data$EDU$trial=="T1"), cols,parameters_2d_int_data,paste("EDU ",labelling_2d_int_data,sep=""),"EDU_T1_")
    plot_black_bg_histogram_timeline(subset(spheroid_2dint_data$EDU,spheroid_2dint_data$EDU$trial=="T2"), cols,parameters_2d_int_data,paste("EDU ",labelling_2d_int_data,sep=""),"EDU_T2_")
    plot_black_bg_histogram_timeline(subset(spheroid_2dint_data$EDU,spheroid_2dint_data$EDU$trial=="T3"), cols,parameters_2d_int_data,paste("EDU ",labelling_2d_int_data,sep=""),"EDU_T3_")
    
    plot_black_bg_boxplot_timeline(subset(spheroid_2dint_data$EDU,spheroid_2dint_data$EDU$trial=="T1"), cols,parameters_2d_int_data,paste("EDU ",labelling_2d_int_data,sep=""),"EDU_T1_")
    plot_black_bg_boxplot_timeline(subset(spheroid_2dint_data$EDU,spheroid_2dint_data$EDU$trial=="T2"), cols,parameters_2d_int_data,paste("EDU ",labelling_2d_int_data,sep=""),"EDU_T2_")
    plot_black_bg_boxplot_timeline(subset(spheroid_2dint_data$EDU,spheroid_2dint_data$EDU$trial=="T3"), cols,parameters_2d_int_data,paste("EDU ",labelling_2d_int_data,sep=""),"EDU_T3_")
    
    plot_black_bg_barplot_timeline(spheroid_2dint_data$EDU, cols,parameters_2d_int_data,paste("EDU ",labelling_2d_int_data,sep=""),"EDU_")
    
    
  }
  #Spheoid 3D_Int
  {
    dire<-paste(dird,"spheroid_3d_int/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_black_bg_histogram_timeline(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T1_")
    plot_black_bg_histogram_timeline(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T2_")
    plot_black_bg_histogram_timeline(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T3_")
    
    plot_black_bg_boxplot_timeline(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T1_")
    plot_black_bg_boxplot_timeline(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T2_")
    plot_black_bg_boxplot_timeline(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T3_")
    
    plot_black_bg_barplot_timeline(spheroid_int_data$dna, cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_")
    
    plot_black_bg_histogram_timeline(subset(spheroid_int_data$EDU,spheroid_int_data$EDU$trial=="T1"), cols,parameters_int_data,paste("EDU ",labelling_int_data,sep=""),"EDU_T1_")
    plot_black_bg_histogram_timeline(subset(spheroid_int_data$EDU,spheroid_int_data$EDU$trial=="T2"), cols,parameters_int_data,paste("EDU ",labelling_int_data,sep=""),"EDU_T2_")
    plot_black_bg_histogram_timeline(subset(spheroid_int_data$EDU,spheroid_int_data$EDU$trial=="T3"), cols,parameters_int_data,paste("EDU ",labelling_int_data,sep=""),"EDU_T3_")
    
    plot_black_bg_boxplot_timeline(subset(spheroid_int_data$EDU,spheroid_int_data$EDU$trial=="T1"), cols,parameters_int_data,paste("EDU ",labelling_int_data,sep=""),"EDU_T1_")
    plot_black_bg_boxplot_timeline(subset(spheroid_int_data$EDU,spheroid_int_data$EDU$trial=="T2"), cols,parameters_int_data,paste("EDU ",labelling_int_data,sep=""),"EDU_T2_")
    plot_black_bg_boxplot_timeline(subset(spheroid_int_data$EDU,spheroid_int_data$EDU$trial=="T3"), cols,parameters_int_data,paste("EDU ",labelling_int_data,sep=""),"EDU_T3_")
    
    plot_black_bg_barplot_timeline(spheroid_int_data$EDU, cols,parameters_int_data,paste("EDU ",labelling_int_data,sep=""),"EDU_")
    
    
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
    
    
    plot_black_bg_histogram_timeline(subset(spheroid_angle$EDU,spheroid_angle$EDU$trial=="T1"), cols,parameters_angle_data,paste("EDU ",labelling_angle_data,sep=""),"EDU_T1_")
    plot_black_bg_histogram_timeline(subset(spheroid_angle$EDU,spheroid_angle$EDU$trial=="T2"), cols,parameters_angle_data,paste("EDU ",labelling_angle_data,sep=""),"EDU_T2_")
    plot_black_bg_histogram_timeline(subset(spheroid_angle$EDU,spheroid_angle$EDU$trial=="T3"), cols,parameters_angle_data,paste("EDU ",labelling_angle_data,sep=""),"EDU_T3_")
    
    plot_black_bg_boxplot_timeline(subset(spheroid_angle$EDU,spheroid_angle$EDU$trial=="T1"), cols,parameters_angle_data,paste("EDU ",labelling_angle_data,sep=""),"EDU_T1_")
    plot_black_bg_boxplot_timeline(subset(spheroid_angle$EDU,spheroid_angle$EDU$trial=="T2"), cols,parameters_angle_data,paste("EDU ",labelling_angle_data,sep=""),"EDU_T2_")
    plot_black_bg_boxplot_timeline(subset(spheroid_angle$EDU,spheroid_angle$EDU$trial=="T3"), cols,parameters_angle_data,paste("EDU ",labelling_angle_data,sep=""),"EDU_T3_")
    
    plot_black_bg_barplot_timeline(spheroid_angle$EDU, cols,parameters_angle_data,paste("EDU ",labelling_angle_data,sep=""),"EDU_")
    
    
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
    
    plot_black_bg_histogram_timeline(subset(spheroid_axial$EDU,spheroid_axial$EDU$trial=="T1"), cols,parameters_axial_data,paste("EDU ",labelling_axial_data,sep=""),"EDU_T1_")
    plot_black_bg_histogram_timeline(subset(spheroid_axial$EDU,spheroid_axial$EDU$trial=="T2"), cols,parameters_axial_data,paste("EDU ",labelling_axial_data,sep=""),"EDU_T2_")
    plot_black_bg_histogram_timeline(subset(spheroid_axial$EDU,spheroid_axial$EDU$trial=="T3"), cols,parameters_axial_data,paste("EDU ",labelling_axial_data,sep=""),"EDU_T3_")
    
    plot_black_bg_boxplot_timeline(subset(spheroid_axial$EDU,spheroid_axial$EDU$trial=="T1"), cols,parameters_axial_data,paste("EDU ",labelling_axial_data,sep=""),"EDU_T1_")
    plot_black_bg_boxplot_timeline(subset(spheroid_axial$EDU,spheroid_axial$EDU$trial=="T2"), cols,parameters_axial_data,paste("EDU ",labelling_axial_data,sep=""),"EDU_T2_")
    plot_black_bg_boxplot_timeline(subset(spheroid_axial$EDU,spheroid_axial$EDU$trial=="T3"), cols,parameters_axial_data,paste("EDU ",labelling_axial_data,sep=""),"EDU_T3_")
    
    plot_black_bg_barplot_timeline(spheroid_axial$EDU, cols,parameters_axial_data,paste("EDU ",labelling_axial_data,sep=""),"EDU_")
    
    
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
    
    plot_black_bg_histogram_timeline(subset(spheroid_radial$EDU,spheroid_radial$EDU$trial=="T1"), cols,parameters_radial_data,paste("EDU ",labelling_radial_data,sep=""),"EDU_T1_")
    plot_black_bg_histogram_timeline(subset(spheroid_radial$EDU,spheroid_radial$EDU$trial=="T2"), cols,parameters_radial_data,paste("EDU ",labelling_radial_data,sep=""),"EDU_T2_")
    plot_black_bg_histogram_timeline(subset(spheroid_radial$EDU,spheroid_radial$EDU$trial=="T3"), cols,parameters_radial_data,paste("EDU ",labelling_radial_data,sep=""),"EDU_T3_")
    
    plot_black_bg_boxplot_timeline(subset(spheroid_radial$EDU,spheroid_radial$EDU$trial=="T1"), cols,parameters_radial_data,paste("EDU ",labelling_radial_data,sep=""),"EDU_T1_")
    plot_black_bg_boxplot_timeline(subset(spheroid_radial$EDU,spheroid_radial$EDU$trial=="T2"), cols,parameters_radial_data,paste("EDU ",labelling_radial_data,sep=""),"EDU_T2_")
    plot_black_bg_boxplot_timeline(subset(spheroid_radial$EDU,spheroid_radial$EDU$trial=="T3"), cols,parameters_radial_data,paste("EDU ",labelling_radial_data,sep=""),"EDU_T3_")
    
    plot_black_bg_barplot_timeline(spheroid_radial$EDU, cols,parameters_radial_data,paste("EDU ",labelling_radial_data,sep=""),"EDU_")
    
    
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
    
    plot_black_bg_histogram_timeline(subset(nuclei_2dint_data$EDU,nuclei_2dint_data$EDU$trial=="T1"), cols,parameters_2d_int_data,paste("EDU ",nuc_labelling_2d_int_data,sep=""),"EDU_T1_")
    plot_black_bg_histogram_timeline(subset(nuclei_2dint_data$EDU,nuclei_2dint_data$EDU$trial=="T2"), cols,parameters_2d_int_data,paste("EDU ",nuc_labelling_2d_int_data,sep=""),"EDU_T2_")
    plot_black_bg_histogram_timeline(subset(nuclei_2dint_data$EDU,nuclei_2dint_data$EDU$trial=="T3"), cols,parameters_2d_int_data,paste("EDU ",nuc_labelling_2d_int_data,sep=""),"EDU_T3_")
    
    plot_black_bg_boxplot_timeline(subset(nuclei_2dint_data$EDU,nuclei_2dint_data$EDU$trial=="T1"), cols,parameters_2d_int_data,paste("F-EDU ",nuc_labelling_2d_int_data,sep=""),"EDU_T1_")
    plot_black_bg_boxplot_timeline(subset(nuclei_2dint_data$EDU,nuclei_2dint_data$EDU$trial=="T2"), cols,parameters_2d_int_data,paste("F-EDU ",nuc_labelling_2d_int_data,sep=""),"EDU_T2_")
    plot_black_bg_boxplot_timeline(subset(nuclei_2dint_data$EDU,nuclei_2dint_data$EDU$trial=="T3"), cols,parameters_2d_int_data,paste("F-EDU ",nuc_labelling_2d_int_data,sep=""),"EDU_T3_")
    
    plot_black_bg_barplot_timeline(nuclei_2dint_data$EDU, cols,parameters_2d_int_data,paste("EDU ",nuc_labelling_2d_int_data,sep=""),"EDU_")
    
    
    
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
    
    plot_black_bg_histogram_timeline(subset(nuclei_3dint_data$EDU,nuclei_3dint_data$EDU$trial=="T1"), cols,parameters_int_data,paste("EDU ",nuclabelling_int_data,sep=""),"EDU_T1_")
    plot_black_bg_histogram_timeline(subset(nuclei_3dint_data$EDU,nuclei_3dint_data$EDU$trial=="T2"), cols,parameters_int_data,paste("EDU ",nuclabelling_int_data,sep=""),"EDU_T2_")
    plot_black_bg_histogram_timeline(subset(nuclei_3dint_data$EDU,nuclei_3dint_data$EDU$trial=="T3"), cols,parameters_int_data,paste("EDU ",nuclabelling_int_data,sep=""),"EDU_T3_")
    
    plot_black_bg_boxplot_timeline(subset(nuclei_3dint_data$EDU,nuclei_3dint_data$EDU$trial=="T1"), cols,parameters_int_data,paste("EDU ",nuclabelling_int_data,sep=""),"EDU_T1_")
    plot_black_bg_boxplot_timeline(subset(nuclei_3dint_data$EDU,nuclei_3dint_data$EDU$trial=="T2"), cols,parameters_int_data,paste("EDU ",nuclabelling_int_data,sep=""),"EDU_T2_")
    plot_black_bg_boxplot_timeline(subset(nuclei_3dint_data$EDU,nuclei_3dint_data$EDU$trial=="T3"), cols,parameters_int_data,paste("EDU ",nuclabelling_int_data,sep=""),"EDU_T3_")
    
    plot_black_bg_barplot_timeline(nuclei_3dint_data$EDU, cols,parameters_int_data,paste("EDU ",nuclabelling_int_data,sep=""),"EDU_")
    
    
    
  }
  
  #Spheroid_median_nuclear_properties
  {
    dire<-paste(dird,"median_nuclear_properties/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_black_bg_histogram_timeline(subset(spheroid_nucleus_combined,spheroid_nucleus_combined$trial=="T1"), cols,parameters_geo_texture_data,paste("Nuclear ",labelling_geo_texture_data,sep=""),"geo_T1_")
    plot_black_bg_histogram_timeline(subset(spheroid_nucleus_combined,spheroid_nucleus_combined$trial=="T2"), cols,parameters_geo_texture_data,paste("Nuclear ",labelling_geo_texture_data,sep=""),"geo_T2_")
    plot_black_bg_histogram_timeline(subset(spheroid_nucleus_combined,spheroid_nucleus_combined$trial=="T3"), cols,parameters_geo_texture_data,paste("Nuclear ",labelling_geo_texture_data,sep=""),"geo_T3_")
    
    plot_black_bg_boxplot_timeline(subset(spheroid_nucleus_combined,spheroid_nucleus_combined$trial=="T1"), cols,parameters_geo_texture_data,labelling_geo_texture_data,"geo_T1_")
    plot_black_bg_boxplot_timeline(subset(spheroid_nucleus_combined,spheroid_nucleus_combined$trial=="T2"), cols,parameters_geo_texture_data,labelling_geo_texture_data,"geo_T2_")
    plot_black_bg_boxplot_timeline(subset(spheroid_nucleus_combined,spheroid_nucleus_combined$trial=="T3"), cols,parameters_geo_texture_data,labelling_geo_texture_data,"geo_T3_")
    
    plot_black_bg_barplot_timeline(spheroid_nucleus_combined, cols,parameters_geo_texture_data,paste("Nuclear ",labelling_geo_texture_data,sep=""),"geo_")
    
  }
  #Spheroid_Nuclear Proj_Int
  {
    
    dire<-paste(dird,"median_nuclear_proj_int/",sep="")
    dir.create(dire)
    setwd(dire)
    
    plot_black_bg_histogram_timeline(subset(spheroid_nuclei_2dint_data$EDU,spheroid_nuclei_2dint_data$EDU$trial=="T1"), cols,parameters_2d_int_data,paste("EDU ",nuc_labelling_2d_int_data,sep=""),"EDU_T1_")
    plot_black_bg_histogram_timeline(subset(spheroid_nuclei_2dint_data$EDU,spheroid_nuclei_2dint_data$EDU$trial=="T2"), cols,parameters_2d_int_data,paste("EDU ",nuc_labelling_2d_int_data,sep=""),"EDU_T2_")
    plot_black_bg_histogram_timeline(subset(spheroid_nuclei_2dint_data$EDU,spheroid_nuclei_2dint_data$EDU$trial=="T3"), cols,parameters_2d_int_data,paste("EDU ",nuc_labelling_2d_int_data,sep=""),"EDU_T3_")
    
    plot_black_bg_boxplot_timeline(subset(spheroid_nuclei_2dint_data$EDU,spheroid_nuclei_2dint_data$EDU$trial=="T1"), cols,parameters_2d_int_data,paste("EDU ",nuc_labelling_2d_int_data,sep=""),"EDU_T1_")
    plot_black_bg_boxplot_timeline(subset(spheroid_nuclei_2dint_data$EDU,spheroid_nuclei_2dint_data$EDU$trial=="T2"), cols,parameters_2d_int_data,paste("EDU ",nuc_labelling_2d_int_data,sep=""),"EDU_T2_")
    plot_black_bg_boxplot_timeline(subset(spheroid_nuclei_2dint_data$EDU,spheroid_nuclei_2dint_data$EDU$trial=="T3"), cols,parameters_2d_int_data,paste("EDU ",nuc_labelling_2d_int_data,sep=""),"EDU_T3_")
    
    plot_black_bg_barplot_timeline(spheroid_nuclei_2dint_data$EDU, cols,parameters_2d_int_data,paste("EDU ",nuc_labelling_2d_int_data,sep=""),"EDU_")
    
    
    
  }
  #spheroid_median_Nuclear 3D_Int
  {
    dire<-paste(dird,"median_nuclear_3d_int/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_black_bg_histogram_timeline(subset(spheroid_nuclei_3dint_data$dna,spheroid_nuclei_3dint_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T1_")
    plot_black_bg_histogram_timeline(subset(spheroid_nuclei_3dint_data$dna,spheroid_nuclei_3dint_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T2_")
    plot_black_bg_histogram_timeline(subset(spheroid_nuclei_3dint_data$dna,spheroid_nuclei_3dint_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T3_")
    
    plot_black_bg_boxplot_timeline(subset(spheroid_nuclei_3dint_data$dna,spheroid_nuclei_3dint_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T1_")
    plot_black_bg_boxplot_timeline(subset(spheroid_nuclei_3dint_data$dna,spheroid_nuclei_3dint_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T2_")
    plot_black_bg_boxplot_timeline(subset(spheroid_nuclei_3dint_data$dna,spheroid_nuclei_3dint_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T3_")
    
    plot_black_bg_barplot_timeline(spheroid_nuclei_3dint_data$dna, cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_")
    
    plot_black_bg_histogram_timeline(subset(spheroid_nuclei_3dint_data$EDU,spheroid_nuclei_3dint_data$EDU$trial=="T1"), cols,parameters_int_data,paste("EDU ",nuclabelling_int_data,sep=""),"EDU_T1_")
    plot_black_bg_histogram_timeline(subset(spheroid_nuclei_3dint_data$EDU,spheroid_nuclei_3dint_data$EDU$trial=="T2"), cols,parameters_int_data,paste("EDU ",nuclabelling_int_data,sep=""),"EDU_T2_")
    plot_black_bg_histogram_timeline(subset(spheroid_nuclei_3dint_data$EDU,spheroid_nuclei_3dint_data$EDU$trial=="T3"), cols,parameters_int_data,paste("EDU ",nuclabelling_int_data,sep=""),"EDU_T3_")
    
    plot_black_bg_boxplot_timeline(subset(spheroid_nuclei_3dint_data$EDU,spheroid_nuclei_3dint_data$EDU$trial=="T1"), cols,parameters_int_data,paste("EDU ",nuclabelling_int_data,sep=""),"EDU_T1_")
    plot_black_bg_boxplot_timeline(subset(spheroid_nuclei_3dint_data$EDU,spheroid_nuclei_3dint_data$EDU$trial=="T2"), cols,parameters_int_data,paste("EDU ",nuclabelling_int_data,sep=""),"EDU_T2_")
    plot_black_bg_boxplot_timeline(subset(spheroid_nuclei_3dint_data$EDU,spheroid_nuclei_3dint_data$EDU$trial=="T3"), cols,parameters_int_data,paste("EDU ",nuclabelling_int_data,sep=""),"EDU_T3_")
    
    plot_black_bg_barplot_timeline(spheroid_nuclei_3dint_data$EDU, cols,parameters_int_data,paste("EDU ",nuclabelling_int_data,sep=""),"EDU_")
    
    
    
  }
  
  
  dird<-"E:/HMF3A_reprogramming/R_analysis/EDU/white_bg_plots/"
  dir.create(dird)
  setwd(dird)
  cols<-colorRampPalette(c("red","yellow"))( 5) 
  #Spheoid Proj_Int
  {
    
    dire<-paste(dird,"spheroid_proj_int/",sep="")
    dir.create(dire)
    setwd(dire)
    
    
    plot_white_bg_histogram_timeline(subset(spheroid_2dint_data$EDU,spheroid_2dint_data$EDU$trial=="T1"), cols,parameters_2d_int_data,paste("EDU ",labelling_2d_int_data,sep=""),"EDU_T1_")
    plot_white_bg_histogram_timeline(subset(spheroid_2dint_data$EDU,spheroid_2dint_data$EDU$trial=="T2"), cols,parameters_2d_int_data,paste("EDU ",labelling_2d_int_data,sep=""),"EDU_T2_")
    plot_white_bg_histogram_timeline(subset(spheroid_2dint_data$EDU,spheroid_2dint_data$EDU$trial=="T3"), cols,parameters_2d_int_data,paste("EDU ",labelling_2d_int_data,sep=""),"EDU_T3_")
    
    plot_white_bg_boxplot_timeline(subset(spheroid_2dint_data$EDU,spheroid_2dint_data$EDU$trial=="T1"), cols,parameters_2d_int_data,paste("EDU ",labelling_2d_int_data,sep=""),"EDU_T1_")
    plot_white_bg_boxplot_timeline(subset(spheroid_2dint_data$EDU,spheroid_2dint_data$EDU$trial=="T2"), cols,parameters_2d_int_data,paste("EDU ",labelling_2d_int_data,sep=""),"EDU_T2_")
    plot_white_bg_boxplot_timeline(subset(spheroid_2dint_data$EDU,spheroid_2dint_data$EDU$trial=="T3"), cols,parameters_2d_int_data,paste("EDU ",labelling_2d_int_data,sep=""),"EDU_T3_")
    
    plot_white_bg_barplot_timeline(spheroid_2dint_data$EDU, cols,parameters_2d_int_data,paste("EDU ",labelling_2d_int_data,sep=""),"EDU_")
    
    
  }
  #Spheoid 3D_Int
  {
    dire<-paste(dird,"spheroid_3d_int/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_white_bg_histogram_timeline(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T1_")
    plot_white_bg_histogram_timeline(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T2_")
    plot_white_bg_histogram_timeline(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T3_")
    
    plot_white_bg_boxplot_timeline(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T1_")
    plot_white_bg_boxplot_timeline(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T2_")
    plot_white_bg_boxplot_timeline(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T3_")
    
    plot_white_bg_barplot_timeline(spheroid_int_data$dna, cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_")
    
    plot_white_bg_histogram_timeline(subset(spheroid_int_data$EDU,spheroid_int_data$EDU$trial=="T1"), cols,parameters_int_data,paste("EDU ",labelling_int_data,sep=""),"EDU_T1_")
    plot_white_bg_histogram_timeline(subset(spheroid_int_data$EDU,spheroid_int_data$EDU$trial=="T2"), cols,parameters_int_data,paste("EDU ",labelling_int_data,sep=""),"EDU_T2_")
    plot_white_bg_histogram_timeline(subset(spheroid_int_data$EDU,spheroid_int_data$EDU$trial=="T3"), cols,parameters_int_data,paste("EDU ",labelling_int_data,sep=""),"EDU_T3_")
    
    plot_white_bg_boxplot_timeline(subset(spheroid_int_data$EDU,spheroid_int_data$EDU$trial=="T1"), cols,parameters_int_data,paste("EDU ",labelling_int_data,sep=""),"EDU_T1_")
    plot_white_bg_boxplot_timeline(subset(spheroid_int_data$EDU,spheroid_int_data$EDU$trial=="T2"), cols,parameters_int_data,paste("EDU ",labelling_int_data,sep=""),"EDU_T2_")
    plot_white_bg_boxplot_timeline(subset(spheroid_int_data$EDU,spheroid_int_data$EDU$trial=="T3"), cols,parameters_int_data,paste("EDU ",labelling_int_data,sep=""),"EDU_T3_")
    
    plot_white_bg_barplot_timeline(spheroid_int_data$EDU, cols,parameters_int_data,paste("EDU ",labelling_int_data,sep=""),"EDU_")
    
    
  }
  #Spheoid angles
  {
    dire<-paste(dird,"spheroid_angles/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_white_bg_histogram_timeline(subset(spheroid_angle$dna,spheroid_angle$dna$trial=="T1"), cols,parameters_angle_data,paste("DNA ",labelling_angle_data,sep=""),"DNA_T1_")
    plot_white_bg_histogram_timeline(subset(spheroid_angle$dna,spheroid_angle$dna$trial=="T2"), cols,parameters_angle_data,paste("DNA ",labelling_angle_data,sep=""),"DNA_T2_")
    plot_white_bg_histogram_timeline(subset(spheroid_angle$dna,spheroid_angle$dna$trial=="T3"), cols,parameters_angle_data,paste("DNA ",labelling_angle_data,sep=""),"DNA_T3_")
    
    plot_white_bg_boxplot_timeline(subset(spheroid_angle$dna,spheroid_angle$dna$trial=="T1"), cols,parameters_angle_data,paste("DNA ",labelling_angle_data,sep=""),"DNA_T1_")
    plot_white_bg_boxplot_timeline(subset(spheroid_angle$dna,spheroid_angle$dna$trial=="T2"), cols,parameters_angle_data,paste("DNA ",labelling_angle_data,sep=""),"DNA_T2_")
    plot_white_bg_boxplot_timeline(subset(spheroid_angle$dna,spheroid_angle$dna$trial=="T3"), cols,parameters_angle_data,paste("DNA ",labelling_angle_data,sep=""),"DNA_T3_")
    
    
    plot_white_bg_histogram_timeline(subset(spheroid_angle$EDU,spheroid_angle$EDU$trial=="T1"), cols,parameters_angle_data,paste("EDU ",labelling_angle_data,sep=""),"EDU_T1_")
    plot_white_bg_histogram_timeline(subset(spheroid_angle$EDU,spheroid_angle$EDU$trial=="T2"), cols,parameters_angle_data,paste("EDU ",labelling_angle_data,sep=""),"EDU_T2_")
    plot_white_bg_histogram_timeline(subset(spheroid_angle$EDU,spheroid_angle$EDU$trial=="T3"), cols,parameters_angle_data,paste("EDU ",labelling_angle_data,sep=""),"EDU_T3_")
    
    plot_white_bg_boxplot_timeline(subset(spheroid_angle$EDU,spheroid_angle$EDU$trial=="T1"), cols,parameters_angle_data,paste("EDU ",labelling_angle_data,sep=""),"EDU_T1_")
    plot_white_bg_boxplot_timeline(subset(spheroid_angle$EDU,spheroid_angle$EDU$trial=="T2"), cols,parameters_angle_data,paste("EDU ",labelling_angle_data,sep=""),"EDU_T2_")
    plot_white_bg_boxplot_timeline(subset(spheroid_angle$EDU,spheroid_angle$EDU$trial=="T3"), cols,parameters_angle_data,paste("EDU ",labelling_angle_data,sep=""),"EDU_T3_")
    
    plot_white_bg_barplot_timeline(spheroid_angle$EDU, cols,parameters_angle_data,paste("EDU ",labelling_angle_data,sep=""),"EDU_")
    
    
  }
  #Spheoid axial distribution
  {
    dire<-paste(dird,"spheroid_axial_distribution/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_white_bg_histogram_timeline(subset(spheroid_axial$dna,spheroid_axial$dna$trial=="T1"), cols,parameters_axial_data,paste("DNA ",labelling_axial_data,sep=""),"DNA_T1_")
    plot_white_bg_histogram_timeline(subset(spheroid_axial$dna,spheroid_axial$dna$trial=="T2"), cols,parameters_axial_data,paste("DNA ",labelling_axial_data,sep=""),"DNA_T2_")
    plot_white_bg_histogram_timeline(subset(spheroid_axial$dna,spheroid_axial$dna$trial=="T3"), cols,parameters_axial_data,paste("DNA ",labelling_axial_data,sep=""),"DNA_T3_")
    
    plot_white_bg_boxplot_timeline(subset(spheroid_axial$dna,spheroid_axial$dna$trial=="T1"), cols,parameters_axial_data,paste("DNA ",labelling_axial_data,sep=""),"DNA_T1_")
    plot_white_bg_boxplot_timeline(subset(spheroid_axial$dna,spheroid_axial$dna$trial=="T2"), cols,parameters_axial_data,paste("DNA ",labelling_axial_data,sep=""),"DNA_T2_")
    plot_white_bg_boxplot_timeline(subset(spheroid_axial$dna,spheroid_axial$dna$trial=="T3"), cols,parameters_axial_data,paste("DNA ",labelling_axial_data,sep=""),"DNA_T3_")
    
    plot_white_bg_barplot_timeline(spheroid_axial$dna, cols,parameters_axial_data,paste("DNA ",labelling_axial_data,sep=""),"DNA_")
    
    plot_white_bg_histogram_timeline(subset(spheroid_axial$EDU,spheroid_axial$EDU$trial=="T1"), cols,parameters_axial_data,paste("EDU ",labelling_axial_data,sep=""),"EDU_T1_")
    plot_white_bg_histogram_timeline(subset(spheroid_axial$EDU,spheroid_axial$EDU$trial=="T2"), cols,parameters_axial_data,paste("EDU ",labelling_axial_data,sep=""),"EDU_T2_")
    plot_white_bg_histogram_timeline(subset(spheroid_axial$EDU,spheroid_axial$EDU$trial=="T3"), cols,parameters_axial_data,paste("EDU ",labelling_axial_data,sep=""),"EDU_T3_")
    
    plot_white_bg_boxplot_timeline(subset(spheroid_axial$EDU,spheroid_axial$EDU$trial=="T1"), cols,parameters_axial_data,paste("EDU ",labelling_axial_data,sep=""),"EDU_T1_")
    plot_white_bg_boxplot_timeline(subset(spheroid_axial$EDU,spheroid_axial$EDU$trial=="T2"), cols,parameters_axial_data,paste("EDU ",labelling_axial_data,sep=""),"EDU_T2_")
    plot_white_bg_boxplot_timeline(subset(spheroid_axial$EDU,spheroid_axial$EDU$trial=="T3"), cols,parameters_axial_data,paste("EDU ",labelling_axial_data,sep=""),"EDU_T3_")
    
    plot_white_bg_barplot_timeline(spheroid_axial$EDU, cols,parameters_axial_data,paste("EDU ",labelling_axial_data,sep=""),"EDU_")
    
    
  }
  #Spheoid radial distribution
  {
    dire<-paste(dird,"spheroid_radial_distribution/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_white_bg_histogram_timeline(subset(spheroid_radial$dna,spheroid_radial$dna$trial=="T1"), cols,parameters_radial_data,paste("DNA ",labelling_radial_data,sep=""),"DNA_T1_")
    plot_white_bg_histogram_timeline(subset(spheroid_radial$dna,spheroid_radial$dna$trial=="T2"), cols,parameters_radial_data,paste("DNA ",labelling_radial_data,sep=""),"DNA_T2_")
    plot_white_bg_histogram_timeline(subset(spheroid_radial$dna,spheroid_radial$dna$trial=="T3"), cols,parameters_radial_data,paste("DNA ",labelling_radial_data,sep=""),"DNA_T3_")
    
    plot_white_bg_boxplot_timeline(subset(spheroid_radial$dna,spheroid_radial$dna$trial=="T1"), cols,parameters_radial_data,paste("DNA ",labelling_radial_data,sep=""),"DNA_T1_")
    plot_white_bg_boxplot_timeline(subset(spheroid_radial$dna,spheroid_radial$dna$trial=="T2"), cols,parameters_radial_data,paste("DNA ",labelling_radial_data,sep=""),"DNA_T2_")
    plot_white_bg_boxplot_timeline(subset(spheroid_radial$dna,spheroid_radial$dna$trial=="T3"), cols,parameters_radial_data,paste("DNA ",labelling_radial_data,sep=""),"DNA_T3_")
    
    plot_white_bg_barplot_timeline(spheroid_radial$dna, cols,parameters_radial_data,paste("DNA ",labelling_radial_data,sep=""),"DNA_")
    
    plot_white_bg_histogram_timeline(subset(spheroid_radial$EDU,spheroid_radial$EDU$trial=="T1"), cols,parameters_radial_data,paste("EDU ",labelling_radial_data,sep=""),"EDU_T1_")
    plot_white_bg_histogram_timeline(subset(spheroid_radial$EDU,spheroid_radial$EDU$trial=="T2"), cols,parameters_radial_data,paste("EDU ",labelling_radial_data,sep=""),"EDU_T2_")
    plot_white_bg_histogram_timeline(subset(spheroid_radial$EDU,spheroid_radial$EDU$trial=="T3"), cols,parameters_radial_data,paste("EDU ",labelling_radial_data,sep=""),"EDU_T3_")
    
    plot_white_bg_boxplot_timeline(subset(spheroid_radial$EDU,spheroid_radial$EDU$trial=="T1"), cols,parameters_radial_data,paste("EDU ",labelling_radial_data,sep=""),"EDU_T1_")
    plot_white_bg_boxplot_timeline(subset(spheroid_radial$EDU,spheroid_radial$EDU$trial=="T2"), cols,parameters_radial_data,paste("EDU ",labelling_radial_data,sep=""),"EDU_T2_")
    plot_white_bg_boxplot_timeline(subset(spheroid_radial$EDU,spheroid_radial$EDU$trial=="T3"), cols,parameters_radial_data,paste("EDU ",labelling_radial_data,sep=""),"EDU_T3_")
    
    plot_white_bg_barplot_timeline(spheroid_radial$EDU, cols,parameters_radial_data,paste("EDU ",labelling_radial_data,sep=""),"EDU_")
    
    
  }
  #Spheoid geometeric properties
  {
    dire<-paste(dird,"spheroid_properties/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_white_bg_histogram_timeline(subset(spheroid_combined,spheroid_combined$trial=="T1"), cols,parameters_geometrical_data,paste("Spheroid ",labelling_geometrical_data,sep=""),"geo_T1_")
    plot_white_bg_histogram_timeline(subset(spheroid_combined,spheroid_combined$trial=="T2"), cols,parameters_geometrical_data,paste("Spheroid ",labelling_geometrical_data,sep=""),"geo_T2_")
    plot_white_bg_histogram_timeline(subset(spheroid_combined,spheroid_combined$trial=="T3"), cols,parameters_geometrical_data,paste("Spheroid ",labelling_geometrical_data,sep=""),"geo_T3_")
    
    plot_white_bg_boxplot_timeline(subset(spheroid_combined,spheroid_combined$trial=="T1"), cols,parameters_geometrical_data,labelling_geometrical_data,"geo_T1_")
    plot_white_bg_boxplot_timeline(subset(spheroid_combined,spheroid_combined$trial=="T2"), cols,parameters_geometrical_data,labelling_geometrical_data,"geo_T2_")
    plot_white_bg_boxplot_timeline(subset(spheroid_combined,spheroid_combined$trial=="T3"), cols,parameters_geometrical_data,labelling_geometrical_data,"geo_T3_")
    
    plot_white_bg_barplot_timeline(spheroid_combined, cols,parameters_geometrical_data,paste("Spheroid ",labelling_geometrical_data,sep=""),"geo_")
    
  }
  #Nuclear geometeric properties
  {
    dire<-paste(dird,"nuclear_properties/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_white_bg_histogram_timeline(subset(nucleus_combined,nucleus_combined$trial=="T1"), cols,parameters_geo_texture_data,paste("Nuclear ",labelling_geo_texture_data,sep=""),"geo_T1_")
    plot_white_bg_histogram_timeline(subset(nucleus_combined,nucleus_combined$trial=="T2"), cols,parameters_geo_texture_data,paste("Nuclear ",labelling_geo_texture_data,sep=""),"geo_T2_")
    plot_white_bg_histogram_timeline(subset(nucleus_combined,nucleus_combined$trial=="T3"), cols,parameters_geo_texture_data,paste("Nuclear ",labelling_geo_texture_data,sep=""),"geo_T3_")
    
    plot_white_bg_boxplot_timeline(subset(nucleus_combined,nucleus_combined$trial=="T1"), cols,parameters_geo_texture_data,labelling_geo_texture_data,"geo_T1_")
    plot_white_bg_boxplot_timeline(subset(nucleus_combined,nucleus_combined$trial=="T2"), cols,parameters_geo_texture_data,labelling_geo_texture_data,"geo_T2_")
    plot_white_bg_boxplot_timeline(subset(nucleus_combined,nucleus_combined$trial=="T3"), cols,parameters_geo_texture_data,labelling_geo_texture_data,"geo_T3_")
    
    plot_white_bg_barplot_timeline(nucleus_combined, cols,parameters_geo_texture_data,paste("Nuclear ",labelling_geo_texture_data,sep=""),"geo_")
    
  }
  #Nuclear Proj_Int
  {
    
    dire<-paste(dird,"nuclear_proj_int/",sep="")
    dir.create(dire)
    setwd(dire)
    
    plot_white_bg_histogram_timeline(subset(nuclei_2dint_data$EDU,nuclei_2dint_data$EDU$trial=="T1"), cols,parameters_2d_int_data,paste("EDU ",nuc_labelling_2d_int_data,sep=""),"EDU_T1_")
    plot_white_bg_histogram_timeline(subset(nuclei_2dint_data$EDU,nuclei_2dint_data$EDU$trial=="T2"), cols,parameters_2d_int_data,paste("EDU ",nuc_labelling_2d_int_data,sep=""),"EDU_T2_")
    plot_white_bg_histogram_timeline(subset(nuclei_2dint_data$EDU,nuclei_2dint_data$EDU$trial=="T3"), cols,parameters_2d_int_data,paste("EDU ",nuc_labelling_2d_int_data,sep=""),"EDU_T3_")
    
    plot_white_bg_boxplot_timeline(subset(nuclei_2dint_data$EDU,nuclei_2dint_data$EDU$trial=="T1"), cols,parameters_2d_int_data,paste("F-EDU ",nuc_labelling_2d_int_data,sep=""),"EDU_T1_")
    plot_white_bg_boxplot_timeline(subset(nuclei_2dint_data$EDU,nuclei_2dint_data$EDU$trial=="T2"), cols,parameters_2d_int_data,paste("F-EDU ",nuc_labelling_2d_int_data,sep=""),"EDU_T2_")
    plot_white_bg_boxplot_timeline(subset(nuclei_2dint_data$EDU,nuclei_2dint_data$EDU$trial=="T3"), cols,parameters_2d_int_data,paste("F-EDU ",nuc_labelling_2d_int_data,sep=""),"EDU_T3_")
    
    plot_white_bg_barplot_timeline(nuclei_2dint_data$EDU, cols,parameters_2d_int_data,paste("EDU ",nuc_labelling_2d_int_data,sep=""),"EDU_")
    
    
    
  }
  #Nuclear 3D_Int
  {
    dire<-paste(dird,"nuclear_3d_int/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_white_bg_histogram_timeline(subset(nuclei_3dint_data$dna,nuclei_3dint_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T1_")
    plot_white_bg_histogram_timeline(subset(nuclei_3dint_data$dna,nuclei_3dint_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T2_")
    plot_white_bg_histogram_timeline(subset(nuclei_3dint_data$dna,nuclei_3dint_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T3_")
    
    plot_white_bg_boxplot_timeline(subset(nuclei_3dint_data$dna,nuclei_3dint_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T1_")
    plot_white_bg_boxplot_timeline(subset(nuclei_3dint_data$dna,nuclei_3dint_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T2_")
    plot_white_bg_boxplot_timeline(subset(nuclei_3dint_data$dna,nuclei_3dint_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T3_")
    
    plot_white_bg_barplot_timeline(nuclei_3dint_data$dna, cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_")
    
    plot_white_bg_histogram_timeline(subset(nuclei_3dint_data$EDU,nuclei_3dint_data$EDU$trial=="T1"), cols,parameters_int_data,paste("EDU ",nuclabelling_int_data,sep=""),"EDU_T1_")
    plot_white_bg_histogram_timeline(subset(nuclei_3dint_data$EDU,nuclei_3dint_data$EDU$trial=="T2"), cols,parameters_int_data,paste("EDU ",nuclabelling_int_data,sep=""),"EDU_T2_")
    plot_white_bg_histogram_timeline(subset(nuclei_3dint_data$EDU,nuclei_3dint_data$EDU$trial=="T3"), cols,parameters_int_data,paste("EDU ",nuclabelling_int_data,sep=""),"EDU_T3_")
    
    plot_white_bg_boxplot_timeline(subset(nuclei_3dint_data$EDU,nuclei_3dint_data$EDU$trial=="T1"), cols,parameters_int_data,paste("EDU ",nuclabelling_int_data,sep=""),"EDU_T1_")
    plot_white_bg_boxplot_timeline(subset(nuclei_3dint_data$EDU,nuclei_3dint_data$EDU$trial=="T2"), cols,parameters_int_data,paste("EDU ",nuclabelling_int_data,sep=""),"EDU_T2_")
    plot_white_bg_boxplot_timeline(subset(nuclei_3dint_data$EDU,nuclei_3dint_data$EDU$trial=="T3"), cols,parameters_int_data,paste("EDU ",nuclabelling_int_data,sep=""),"EDU_T3_")
    
    plot_white_bg_barplot_timeline(nuclei_3dint_data$EDU, cols,parameters_int_data,paste("EDU ",nuclabelling_int_data,sep=""),"EDU_")
    
    
    
  }
  
  #Spheroid_median_nuclear_properties
  {
    dire<-paste(dird,"median_nuclear_properties/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_white_bg_histogram_timeline(subset(spheroid_nucleus_combined,spheroid_nucleus_combined$trial=="T1"), cols,parameters_geo_texture_data,paste("Median Spheroid ",labelling_geo_texture_data,sep=""),"geo_T1_")
    plot_white_bg_histogram_timeline(subset(spheroid_nucleus_combined,spheroid_nucleus_combined$trial=="T2"), cols,parameters_geo_texture_data,paste("Median Spheroid ",labelling_geo_texture_data,sep=""),"geo_T2_")
    plot_white_bg_histogram_timeline(subset(spheroid_nucleus_combined,spheroid_nucleus_combined$trial=="T3"), cols,parameters_geo_texture_data,paste("Median Spheroid ",labelling_geo_texture_data,sep=""),"geo_T3_")
    
    plot_white_bg_boxplot_timeline(subset(spheroid_nucleus_combined,spheroid_nucleus_combined$trial=="T1"), cols,parameters_geo_texture_data,labelling_geo_texture_data,"geo_T1_")
    plot_white_bg_boxplot_timeline(subset(spheroid_nucleus_combined,spheroid_nucleus_combined$trial=="T2"), cols,parameters_geo_texture_data,labelling_geo_texture_data,"geo_T2_")
    plot_white_bg_boxplot_timeline(subset(spheroid_nucleus_combined,spheroid_nucleus_combined$trial=="T3"), cols,parameters_geo_texture_data,labelling_geo_texture_data,"geo_T3_")
    
    plot_white_bg_barplot_timeline(spheroid_nucleus_combined, cols,parameters_geo_texture_data,paste("Median Spheroid ",labelling_geo_texture_data,sep=""),"geo_")
    
  }
  #Spheroid_Nuclear Proj_Int
  {
    
    dire<-paste(dird,"median_nuclear_proj_int/",sep="")
    dir.create(dire)
    setwd(dire)
    
    plot_white_bg_histogram_timeline(subset(spheroid_nuclei_2dint_data$EDU,spheroid_nuclei_2dint_data$EDU$trial=="T1"), cols,parameters_2d_int_data,paste("EDU ",nuc_labelling_2d_int_data,sep=""),"EDU_T1_")
    plot_white_bg_histogram_timeline(subset(spheroid_nuclei_2dint_data$EDU,spheroid_nuclei_2dint_data$EDU$trial=="T2"), cols,parameters_2d_int_data,paste("EDU ",nuc_labelling_2d_int_data,sep=""),"EDU_T2_")
    plot_white_bg_histogram_timeline(subset(spheroid_nuclei_2dint_data$EDU,spheroid_nuclei_2dint_data$EDU$trial=="T3"), cols,parameters_2d_int_data,paste("EDU ",nuc_labelling_2d_int_data,sep=""),"EDU_T3_")
    
    plot_white_bg_boxplot_timeline(subset(spheroid_nuclei_2dint_data$EDU,spheroid_nuclei_2dint_data$EDU$trial=="T1"), cols,parameters_2d_int_data,paste("EDU ",nuc_labelling_2d_int_data,sep=""),"EDU_T1_")
    plot_white_bg_boxplot_timeline(subset(spheroid_nuclei_2dint_data$EDU,spheroid_nuclei_2dint_data$EDU$trial=="T2"), cols,parameters_2d_int_data,paste("EDU ",nuc_labelling_2d_int_data,sep=""),"EDU_T2_")
    plot_white_bg_boxplot_timeline(subset(spheroid_nuclei_2dint_data$EDU,spheroid_nuclei_2dint_data$EDU$trial=="T3"), cols,parameters_2d_int_data,paste("EDU ",nuc_labelling_2d_int_data,sep=""),"EDU_T3_")
    
    plot_white_bg_barplot_timeline(spheroid_nuclei_2dint_data$EDU, cols,parameters_2d_int_data,paste("EDU ",nuc_labelling_2d_int_data,sep=""),"EDU_")
    
    
    
  }
  #spheroid_median_Nuclear 3D_Int
  {
    dire<-paste(dird,"median_nuclear_3d_int/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_white_bg_histogram_timeline(subset(spheroid_nuclei_3dint_data$dna,spheroid_nuclei_3dint_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T1_")
    plot_white_bg_histogram_timeline(subset(spheroid_nuclei_3dint_data$dna,spheroid_nuclei_3dint_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T2_")
    plot_white_bg_histogram_timeline(subset(spheroid_nuclei_3dint_data$dna,spheroid_nuclei_3dint_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T3_")
    
    plot_white_bg_boxplot_timeline(subset(spheroid_nuclei_3dint_data$dna,spheroid_nuclei_3dint_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T1_")
    plot_white_bg_boxplot_timeline(subset(spheroid_nuclei_3dint_data$dna,spheroid_nuclei_3dint_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T2_")
    plot_white_bg_boxplot_timeline(subset(spheroid_nuclei_3dint_data$dna,spheroid_nuclei_3dint_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T3_")
    
    plot_white_bg_barplot_timeline(spheroid_nuclei_3dint_data$dna, cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_")
    
    plot_white_bg_histogram_timeline(subset(spheroid_nuclei_3dint_data$EDU,spheroid_nuclei_3dint_data$EDU$trial=="T1"), cols,parameters_int_data,paste("EDU ",nuclabelling_int_data,sep=""),"EDU_T1_")
    plot_white_bg_histogram_timeline(subset(spheroid_nuclei_3dint_data$EDU,spheroid_nuclei_3dint_data$EDU$trial=="T2"), cols,parameters_int_data,paste("EDU ",nuclabelling_int_data,sep=""),"EDU_T2_")
    plot_white_bg_histogram_timeline(subset(spheroid_nuclei_3dint_data$EDU,spheroid_nuclei_3dint_data$EDU$trial=="T3"), cols,parameters_int_data,paste("EDU ",nuclabelling_int_data,sep=""),"EDU_T3_")
    
    plot_white_bg_boxplot_timeline(subset(spheroid_nuclei_3dint_data$EDU,spheroid_nuclei_3dint_data$EDU$trial=="T1"), cols,parameters_int_data,paste("EDU ",nuclabelling_int_data,sep=""),"EDU_T1_")
    plot_white_bg_boxplot_timeline(subset(spheroid_nuclei_3dint_data$EDU,spheroid_nuclei_3dint_data$EDU$trial=="T2"), cols,parameters_int_data,paste("EDU ",nuclabelling_int_data,sep=""),"EDU_T2_")
    plot_white_bg_boxplot_timeline(subset(spheroid_nuclei_3dint_data$EDU,spheroid_nuclei_3dint_data$EDU$trial=="T3"), cols,parameters_int_data,paste("EDU ",nuclabelling_int_data,sep=""),"EDU_T3_")
    
    plot_white_bg_barplot_timeline(spheroid_nuclei_3dint_data$EDU, cols,parameters_int_data,paste("EDU ",nuclabelling_int_data,sep=""),"EDU_")
    
    
    
  }
  
  
}

#write the collected data
{
  dir.create("E:/HMF3A_reprogramming/R_analysis/EDU/combined_data/")
  setwd("E:/HMF3A_reprogramming/R_analysis/EDU/combined_data/")
  
  write.csv(spheroid_combined, file="spheroid_combined.csv")
  write.csv(nucleus_combined, file="nucleus_combined.csv")
  write.csv(spheroid_nucleus_combined, file="spheroid_median_nucleus_combined.csv")
  
  write.csv(spheroid_radial$EDU, file="spheroid_radial_dist_EDU.csv")
  
  write.csv(spheroid_2dint_data$EDU, file="spheroid_2d_int_EDU.csv")
  write.csv(nuclei_2dint_data$EDU, file="nuclear_2d_int_EDU.csv")
  write.csv(spheroid_nuclei_2dint_data$EDU, file="spheroid_median_nuclear_2d_int_EDU.csv")
  
  
  write.csv(nuclei_3dint_data$EDU, file="nuclear_int_EDU.csv")
  write.csv(spheroid_nuclei_3dint_data$EDU, file="spheroid_median_nuclear_int_EDU.csv")
  write.csv(nuclei_3dint_data$dna, file="nuclear_int_dna.csv")
  write.csv(spheroid_nuclei_3dint_data$dna, file="spheroid_median_nuclear_int_dna.csv")
  
  write.csv(spheroid_axial$EDU, file="spheroid_axial_dist_EDU.csv")
  write.csv(spheroid_axial$dna, file="spheroid_axial_dist_dna.csv")
  
  write.csv(spheroid_angle$EDU, file="spheroid_angle_dist_EDU.csv")
  write.csv(spheroid_angle$dna, file="spheroid_angle_dist_dna.csv")
  
  write.csv(spheroid_int_data$EDU, file="spheroid_int_EDU.csv")
  write.csv(spheroid_int_data$dna, file="spheroid_int_dna.csv")
  
  
}

#EDU positive cells calculations
{
  #fraction of EDU positive cells with time
  {
    T1<-subset(nuclei_3dint_data$EDU,nuclei_3dint_data$EDU$trial=="T1")
    T1$Mean<-T1$Mean/max(T1$Mean)
    T1$status<-ifelse(T1$Mean>0.1,1,0)
    T2<-subset(nuclei_3dint_data$EDU,nuclei_3dint_data$EDU$trial=="T2")
    T2$Mean<-T2$Mean/max(T2$Mean)
    T2$status<-ifelse(T2$Mean>0.1,1,0)
    T3<-subset(nuclei_3dint_data$EDU,nuclei_3dint_data$EDU$trial=="T3")
    T3$Mean<-T3$Mean/max(T3$Mean)
    T3$status<-ifelse(T3$Mean>0.1,1,0)
    {
      T2_D0<-nrow(subset(T2,T2$status==1 & T2$sample=="D0"))/nrow(subset(T2, T2$sample=="D0"))
      T2_D2<-nrow(subset(T2,T2$status==1 & T2$sample=="D2"))/nrow(subset(T2, T2$sample=="D2"))
      T2_D4<-nrow(subset(T2,T2$status==1 & T2$sample=="D4"))/nrow(subset(T2, T2$sample=="D4"))
      T2_D6<-nrow(subset(T2,T2$status==1 & T2$sample=="D6"))/nrow(subset(T2, T2$sample=="D6"))
      T2_D8<-nrow(subset(T2,T2$status==1 & T2$sample=="D8"))/nrow(subset(T2, T2$sample=="D8"))
      
      T1_D0<-nrow(subset(T1,T1$status==1 & T1$sample=="D0"))/nrow(subset(T1, T1$sample=="D0"))
      T1_D2<-nrow(subset(T1,T1$status==1 & T1$sample=="D2"))/nrow(subset(T1, T1$sample=="D2"))
      T1_D4<-nrow(subset(T1,T1$status==1 & T1$sample=="D4"))/nrow(subset(T1, T1$sample=="D4"))
      T1_D6<-nrow(subset(T1,T1$status==1 & T1$sample=="D6"))/nrow(subset(T1, T1$sample=="D6"))
      T1_D8<-nrow(subset(T1,T1$status==1 & T1$sample=="D8"))/nrow(subset(T1, T1$sample=="D8"))
      
      T3_D0<-nrow(subset(T3,T3$status==1 & T3$sample=="D0"))/nrow(subset(T3, T3$sample=="D0"))
      T3_D2<-nrow(subset(T3,T3$status==1 & T3$sample=="D2"))/nrow(subset(T3, T3$sample=="D2"))
      T3_D4<-nrow(subset(T3,T3$status==1 & T3$sample=="D4"))/nrow(subset(T3, T3$sample=="D4"))
      T3_D6<-nrow(subset(T3,T3$status==1 & T3$sample=="D6"))/nrow(subset(T3, T3$sample=="D6"))
      T3_D8<-nrow(subset(T3,T3$status==1 & T3$sample=="D8"))/nrow(subset(T3, T3$sample=="D8"))
      
    }
    pop_means<-c(mean(c(T2_D0,T1_D0,T3_D0),na.rm=T),mean(c(T2_D2,T1_D2,T3_D2)),mean(c(T2_D4,T1_D4,T3_D4)),
                 mean(c(T1_D6,T2_D6,T3_D6)),mean(c(T1_D8,T2_D8,T3_D8)))
    pop_sd<-c(sd(c(T2_D0,T1_D0,T3_D0),na.rm=T),sd(c(T2_D2,T1_D2,T3_D2)),sd(c(T2_D4,T1_D4,T3_D4)),
              sd(c(T1_D6,T2_D6,T3_D6)),sd(c(T1_D8,T2_D8,T3_D8)))
    ymax<-max((pop_means+pop_sd)[-1])*1.1
    
    dird<-"E:/HMF3A_reprogramming/R_analysis/EDU/black_bg_plots/"
    setwd(dird)
    png(filename="fraction_EDU_positive_cells.png", units="in",width=2, height=2, pointsize=7,  res=1200)
    par(font.axis = 2,font.lab=2,mar=c(4,5,1,1), font=2,bg='black',fg='white',col.lab='white',col.axis='white')
    foo<-barplot(pop_means,las=1,ylim=c(0,ymax),names=c("D0","D2","D4","D6","D8"),ylab="Fraction of EDU\npositive cells %")
    box()
    segments(foo,pop_means-pop_sd,foo,pop_means+pop_sd)
    dev.off()
    
    dird<-"E:/HMF3A_reprogramming/R_analysis/EDU/white_bg_plots/"
    setwd(dird)
    png(filename="fraction_EDU_positive_cells.png", units="in",width=2, height=2, pointsize=7,  res=1200)
    par(font.axis = 2,font.lab=2,mar=c(4,5,1,1), font=2,bg='white',fg='black',col.lab='black',col.axis='black')
    foo<-barplot(pop_means,las=1,ylim=c(0,ymax),names=c("D0","D2","D4","D6","D8"),ylab="Fraction of EDU\npositive cells %")
    box()
    segments(foo,pop_means-pop_sd,foo,pop_means+pop_sd)
    dev.off()
    
    EDU_3d<-rbind(T1,T2,T3)
  }
  #fraction of EDU postive cells vs cell number
  {
    fields<-levels(as.factor(EDU_3d$Spheroid_id))
    field_EDU<-as.data.frame(matrix(nrow=length(fields),ncol=8))
    colnames(field_EDU)<-c("Spheroid_id","Positive_fraction","sample","trial","num_cells","status","Volume","atleast_60")
    for(index in 1:length(fields)){
      field_EDU[index,1]<-fields[index]
      temp<-subset(EDU_3d,EDU_3d$Spheroid_id==fields[index])
      field_EDU[index,2]<-nrow(subset(temp,temp$status==1))/nrow(temp)
      field_EDU[index,3]<-temp$sample[1]
      field_EDU[index,4]<-temp$trial[1]
      field_EDU[index,5]<-nrow(temp)
      field_EDU[index,6]<-ifelse(nrow(subset(temp,temp$status==1))>=1,1,0)
      field_EDU[index,8]<-ifelse(field_EDU[index,2]>=0.6,1,0)
      
      temp<-subset(nucleus_combined,nucleus_combined$Spheroid_id==fields[index])
      field_EDU[index,7]<-sum(temp$Vol..unit.)
      
    }
    
    D8<-subset(field_EDU,field_EDU$sample=="D8")
    D6<-subset(field_EDU,field_EDU$sample=="D6")
    D4<-subset(field_EDU,field_EDU$sample=="D4")
    D2<-subset(field_EDU,field_EDU$sample=="D2")
    ymax=max(as.numeric(as.character(field_EDU$num_cells)))
    
    dird<-"E:/HMF3A_reprogramming/R_analysis/EDU/black_bg_plots/"
    setwd(dird)
    png(filename="Edu_fraction_vs_numb_cells_scatter.png.png", units="in",width=2, height=2, pointsize=7,  res=1200)
    par(font.axis = 2,font.lab=2,mar=c(4,4.5,1,1), font=2,bg='black',fg='white',col.lab='white',col.axis='white',cex.lab=0.8)
    plot(field_EDU$Positive_fraction~field_EDU$num_cells,ylim=c(0,1),xlim=c(0,ymax),
         xlab="Number of Cells per island",ylab="% Edu positive cells ",las=1)
    points(D6$Positive_fraction~D6$num_cells,pch=19,col=cols[4])
    points(D4$Positive_fraction~D4$num_cells,pch=19,col=cols[3])
    points(D2$Positive_fraction~D2$num_cells,pch=19,col=cols[2])
    points(D8$Positive_fraction~D8$num_cells,pch=19,col=cols[5])
    dev.off()
    
    png(filename="Edu_fraction_vs_numb_cells_scatter_D8.png.png", units="in",width=2, height=2, pointsize=7,  res=1200)
    par(font.axis = 2,font.lab=2,mar=c(4,4.5,1,1), font=2,bg='black',fg='white',col.lab='white',col.axis='white',cex.lab=0.8)
    plot(field_EDU$Positive_fraction~field_EDU$num_cells,ylim=c(0,1),xlim=c(0,ymax),
         xlab="Number of Cells per island",ylab="% Edu positive cells ",las=1,col=NA)
    points(D8$Positive_fraction~D8$num_cells,pch=19,col=alpha.col("gold",alpha = 0.5))
    dev.off()
    
    
    dird<-"E:/HMF3A_reprogramming/R_analysis/EDU/white_bg_plots/"
    setwd(dird)
    png(filename="Edu_fraction_vs_numb_cells_scatter.png", units="in",width=2, height=2, pointsize=7,  res=1200)
    par(font.axis = 2,font.lab=2,mar=c(4,4.5,1,1), font=2,bg='white',fg='black',col.lab='black',col.axis='black',cex.lab=0.8)
    plot(field_EDU$Positive_fraction~field_EDU$num_cells,ylim=c(0,1),xlim=c(0,ymax),
         xlab="Number of Cells per island",ylab="% Edu positive cells ",las=1)
    points(D6$Positive_fraction~D6$num_cells,pch=19,col=cols[4])
    points(D4$Positive_fraction~D4$num_cells,pch=19,col=cols[3])
    points(D2$Positive_fraction~D2$num_cells,pch=19,col=cols[2])
    points(D8$Positive_fraction~D8$num_cells,pch=19,col=cols[5])
    dev.off()
    
    png(filename="Edu_fraction_vs_numb_cells_scatter_D8.png.png", units="in",width=2, height=2, pointsize=7,  res=1200)
    par(font.axis = 2,font.lab=2,mar=c(4,4.5,1,1), font=2,bg='white',fg='black',col.lab='black',col.axis='black',cex.lab=0.8)
    plot(field_EDU$Positive_fraction~field_EDU$num_cells,ylim=c(0,1),xlim=c(0,ymax),
         xlab="Number of Cells per island",ylab="% Edu positive cells ",las=1,col=NA)
    points(D8$Positive_fraction~D8$num_cells,pch=19,col=alpha.col("grey",alpha = 0.5))
    dev.off()
  }
  #fraction of EDU postive islands
  {
    positive<-subset(field_EDU,field_EDU$status==1)
    negative<-subset(field_EDU,field_EDU$status==0)
    
    pos_D8<-subset(positive,positive$sample=="D8")
    neg_D8<-subset(negative,negative$sample=="D8")
    pos_D6<-subset(positive,positive$sample=="D6")
    neg_D6<-subset(negative,negative$sample=="D6")
    pos_D4<-subset(positive,positive$sample=="D4")
    neg_D4<-subset(negative,negative$sample=="D4")
    pos_D2<-subset(positive,positive$sample=="D2")
    neg_D2<-subset(negative,negative$sample=="D2")
    
    dird<-"E:/HMF3A_reprogramming/R_analysis/EDU/black_bg_plots/"
    setwd(dird)
    png(filename="Edu_fraction_vs_numb_cells.png", units="in",width=2, height=2, pointsize=7,  res=1200)
    par(font.axis = 2,font.lab=2,mar=c(4,4.5,1,1), font=2,bg='black',fg='white',col.lab='white',col.axis='white',cex.lab=0.8)
    boxplot(pos_D2$num_cells,neg_D2$num_cells,
            pos_D4$num_cells,neg_D4$num_cells,
            pos_D6$num_cells,neg_D6$num_cells,
            pos_D8$num_cells,neg_D8$num_cells,
            las=1,lty=1,pch=19,col=c("cyan","magenta"), ylab="Number of cells",cex=0.8,lwd=0.6)
    axis(1, at = c(1.5,3.5,5.5,7.5), labels=c("D2","D4",'D6',"D8"),tick =F)
    legend("topleft",legend=c("EDU positive","EDU negative"),fill=c("cyan","magenta"), horiz=F,cex=0.8)
    dev.off()
    
    
    
    T1_pos<-c(nrow(subset(pos_D2,pos_D2$trial=="T1"))/nrow(subset(D2,D2$trial=="T1")),
              nrow(subset(pos_D4,pos_D4$trial=="T1"))/nrow(subset(D4,D4$trial=="T1")),
              nrow(subset(pos_D6,pos_D6$trial=="T1"))/nrow(subset(D6,D6$trial=="T1")),
              nrow(subset(pos_D8,pos_D8$trial=="T1"))/nrow(subset(D8,D8$trial=="T1")))
    T2_pos<-c(nrow(subset(pos_D2,pos_D2$trial=="T2"))/nrow(subset(D2,D2$trial=="T2")),
              nrow(subset(pos_D4,pos_D4$trial=="T2"))/nrow(subset(D4,D4$trial=="T2")),
              nrow(subset(pos_D6,pos_D6$trial=="T2"))/nrow(subset(D6,D6$trial=="T2")),
              nrow(subset(pos_D8,pos_D8$trial=="T2"))/nrow(subset(D8,D8$trial=="T2")))
    T3_pos<-c(nrow(subset(pos_D2,pos_D2$trial=="T3"))/nrow(subset(D2,D2$trial=="T3")),
              nrow(subset(pos_D4,pos_D4$trial=="T3"))/nrow(subset(D4,D4$trial=="T3")),
              nrow(subset(pos_D6,pos_D6$trial=="T3"))/nrow(subset(D6,D6$trial=="T3")),
              nrow(subset(pos_D8,pos_D8$trial=="T3"))/nrow(subset(D8,D8$trial=="T3")))
    pop_means<-c(0,apply(rbind(T1_pos,T2_pos,T3_pos),2,mean))
    pop_sd<-c(0,apply(rbind(T1_pos,T2_pos,T3_pos),2,sd))
    
    png(filename="fraction_EDU_positive_island.png", units="in",width=2, height=2, pointsize=7,  res=1200)
    par(font.axis = 2,font.lab=2,mar=c(4,5,1,1), font=2,bg='black',fg='white',col.lab='white',col.axis='white')
    foo<-barplot(pop_means,las=1,ylim=c(0,1),names=c("D0","D2","D4","D6","D8"),ylab="Fraction of islands\nwith EDU %")
    box()
    segments(foo,pop_means-pop_sd,foo,pop_means+pop_sd)
    dev.off()
    
    
    dird<-"E:/HMF3A_reprogramming/R_analysis/EDU/white_bg_plots/"
    setwd(dird)
    png(filename="Edu_fraction_vs_numb_cells.png", units="in",width=2, height=2, pointsize=7,  res=1200)
    par(font.axis = 2,font.lab=2,mar=c(4,4.5,1,1), font=2,bg='white',fg='black',col.lab='black',col.axis='black',cex.lab=0.8)
    boxplot(pos_D2$num_cells,neg_D2$num_cells,
            pos_D4$num_cells,neg_D4$num_cells,
            pos_D6$num_cells,neg_D6$num_cells,
            pos_D8$num_cells,neg_D8$num_cells,
            las=1,lty=1,pch=19,col=c("cyan","magenta"), ylab="Number of cells",cex=0.8,lwd=0.6)
    axis(1, at = c(1.5,3.5,5.5,7.5), labels=c("D2","D4",'D6',"D8"),tick =F)
    legend("topleft",legend=c("EDU positive","EDU negative"),fill=c("cyan","magenta"), horiz=F,cex=0.8)
    dev.off()
    
    png(filename="fraction_EDU_positive_island.png", units="in",width=2, height=2, pointsize=7,  res=1200)
    par(font.axis = 2,font.lab=2,mar=c(4,5,1,1), font=2,bg='white',fg='black',col.lab='black',col.axis='black')
    foo<-barplot(pop_means,las=1,ylim=c(0,1),names=c("D0","D2","D4","D6","D8"),ylab="Fraction of islands\nwith EDU %")
    box()
    segments(foo,pop_means-pop_sd,foo,pop_means+pop_sd)
    dev.off()
  }
  
}
