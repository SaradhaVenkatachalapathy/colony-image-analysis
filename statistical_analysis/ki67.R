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
                                "DC_Min.","DC_Max","DC_Mean.","DC_SD","FeretAngle",
                                "CX","CY","CZ",
                                "Proj. Area", "Proj. Perimeter","Proj. Major Axis", "Proj. Minor Axis", "Proj. Circularity",
                                "Proj. Feret","Proj. Min Feret","Proj. A.R","Proj. Roundness","Proj. Solidity",
                                
                                "HCcontent","ECcontent","HCvolume","ECvolume","HC_EC_content","HC_EC_volume",
                                "EFD1","EFD2","EFD3","EFD4","EFD5","EFD6","EFD7","EFD8","EFD9","EFD10",                                  
                                
                                "Deg_0_Angular_2nd\n_Moment _step5","Deg_0\nContrast _step5","Deg_0_Correlation\n _step5",
                                "Deg_0_Inverse\nDifference_Moment _step5","Deg_0_Entropy\n_step5",
                                "Deg_90_Angular_2nd\n_Moment _step5","Deg_90\nContrast _step5","Deg_90_Correlation\n _step5",
                                "Deg_90_Inverse\nDifference_Moment _step5","Deg_90_Entropy\n_step5",
                                "Deg_180_Angular_2nd\n_Moment _step5","Deg_180\nContrast _step5", "Deg_180_Correlation\n _step5",
                                "Deg_180_Inverse\nDifference_Moment _step5","Deg_180_Entropy\n_step5",                   
                                "Deg_270_Angular_2nd\n_Moment _step5","Deg_270\nContrast _step5","Deg_270_Correlation\n _step5",
                                "Deg_270_Inverse\nDifference_Moment _step5","Deg_270_Entropy\n_step5",
                                
                                "Deg_0_Angular_2nd\n_Moment _step15","Deg_0\nContrast _step15","Deg_0_Correlation\n _step15",
                                "Deg_0_Inverse\nDifference_Moment _step15" ,"Deg_0_Entropy\n_step15",
                                "Deg_90_Angular_2nd\n_Moment _step15","Deg_90\nContrast _step15","Deg_90_Correlation\n _step15",
                                "Deg_90_Inverse\nDifference_Moment _step15","Deg_90_Entropy\n_step15",                   
                                "Deg_180_Angular_2nd\n_Moment _step15","Deg_180\nContrast _step15","Deg_180_Correlation\n _step15",              
                                "Deg_180_Inverse\nDifference_Moment _step15","Deg_180_Entropy\n_step15",
                                "Deg_270_Angular_2nd\n_Moment _step15","Deg_270\nContrast _step15","Deg_270_Correlation\n _step15",
                                "Deg_270_Inverse\nDifference_Moment _step15","Deg_270_Entropy\n_step15",
                                
                                "Deg_0_Angular_2nd\n_Moment _step25","Deg_0\nContrast _step25","Deg_0_Correlation\n _step25",
                                "Deg_0_Inverse\nDifference_Moment _step25","Deg_0_Entropy\n_step25",
                                "Deg_90_Angular_2nd\n_Moment _step25","Deg_90\nContrast _step25","Deg_90_Correlation\n _step25",               
                                "Deg_90_Inverse\nDifference_Moment _step25","Deg_90_Entropy\n_step25",
                                "Deg_180_Angular_2nd\n_Moment _step25","Deg_180\nContrast _step25","Deg_180_Correlation\n _step25",
                                "Deg_180_Inverse\nDifference_Moment _step25","Deg_180_Entropy\n_step25",
                                "Deg_270_Angular_2nd\n_Moment _step25","Deg_270\nContrast _step25",               
                                "Deg_270_Correlation\n _step25","Deg_270_Inverse\nDifference_Moment _step25","Deg_270_Entropy\n_step25",
                                
                                "Deg_0_Angular_2nd\n_Moment _step35","Deg_0\nContrast _step35","Deg_0_Correlation\n _step35",                
                                "Deg_0_Inverse\nDifference_Moment _step35","Deg_0_Entropy\n_step35",
                                "Deg_90_Angular_2nd\n_Moment _step35","Deg_90\nContrast _step35","Deg_90_Correlation\n _step35",  
                                "Deg_90_Inverse\nDifference_Moment _step35","Deg_90_Entropy\n_step35",
                                "Deg_180_Angular_2nd\n_Moment _step35","Deg_180\nContrast _step35","Deg_180_Correlation\n _step35",              
                                "Deg_180_Inverse\nDifference_Moment _step35","Deg_180_Entropy\n_step35",
                                "Deg_270_Angular_2nd\n_Moment _step35","Deg_270\nContrast _step35","Deg_270_Correlation\n _step35",
                                "Deg_270_Inverse\nDifference_Moment _step35","Deg_270_Entropy\n_step35",
                                
                                "Deg_0_Angular_2nd\n_Moment _step45","Deg_0\nContrast _step45","Deg_0_Correlation\n _step45",
                                "Deg_0_Inverse\nDifference_Moment _step45","Deg_0_Entropy\n_step45",
                                "Deg_90_Angular_2nd\n_Moment _step45","Deg_90\nContrast _step45","Deg_90_Correlation\n _step45",
                                "Deg_90_Inverse\nDifference_Moment _step45","Deg_90_Entropy\n_step45",                   
                                "Deg_180_Angular_2nd\n_Moment _step45","Deg_180\nContrast _step45","Deg_180_Correlation\n _step45",               
                                "Deg_180_Inverse\nDifference_Moment _step45","Deg_180_Entropy\n_step45",
                                "Deg_270_Angular_2nd\n_Moment _step45","Deg_270\nContrast _step45","Deg_270_Correlation\n _step45",
                                "Deg_270_Inverse\nDifference_Moment _step45","Deg_270_Entropy\n_step45",
                                
                                "Deg_0_Angular_2nd\n_Moment _step100","Deg_0\nContrast _step100","Deg_0_Correlation\n _step100",
                                "Deg_0_Inverse\nDifference_Moment _step100","Deg_0_Entropy\n_step100",                                     
                                "Deg_90_Angular_2nd\n_Moment _step100","Deg_90\nContrast _step100","Deg_90_Correlation\n _step100",
                                "Deg_90_Inverse\nDifference_Moment _step100" ,"Deg_90_Entropy\n_step100",
                                "Deg_180_Angular_2nd\n_Moment _step100","Deg_180\nContrast _step100","Deg_180_Correlation\n _step100",
                                "Deg_180_Inverse\nDifference_Moment _step100","Deg_180_Entropy\n_step100",               
                                "Deg_270_Angular_2nd\n_Moment _step100","Deg_270\nContrast _step100","Deg_270_Correlation\n _step100",            
                                "Deg_270_Inverse\nDifference_Moment _step100","Deg_270_Entropy\n_step100")
  
}

source("E:/HMF3A_reprogramming/R_analysis/functions.R")
# initialise the experiment details and the subdirectories required
{
  path_to_ex<-c("E:/HMF3A_reprogramming/Timeline_Ki67/20191214_s14_hmf3a_reprogramming_ki67_488/",
                "E:/HMF3A_reprogramming/Timeline_Ki67/20191223_hmf3a_reprogramming_s16_ki67_488_DAPI/",
                "E:/HMF3A_reprogramming/Timeline_Ki67/20191227_hmf_s17_oct4_488_ki67_568_h3k9ac/")
  
  exp<-c("20191214_s14_hmf3a_reprogramming_ki67_488","20191223_hmf3a_reprogramming_s16_ki67_488_DAPI","20191227_hmf_s17_oct4_488_ki67_568_h3k9ac")
  
  data_types<-c("/3D geometrical data spheroid/",
                "/3D int_data spheroid/DNA/","/3D int_data spheroid/KI67/")
  data_types_2d<-c("/2D_measures_spheroid/2D_spheroid_DNA.csv","/2D_measures_spheroid/2D_spheroid_KI67.csv")
  nuc_data_types<-c("/3D geometrical data/","/3D ellipsoid/","/3D geometerical_simple/","/3D shape measure/",
                    "/3D int_data/DNA/","/3D int_data/KI67/",
                    "/cell_2microns_measure/KI67/","/cell_2microns_measure/DNA/")
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
                       "/2D_measures_nuclei/2D_nucleusKI67.csv")
 
  samples<-c("D0","D2","D4","D6","D8")
  trial<-c("T1","T2","T3")
  
  dird_apo<-"E:/HMF3A_reprogramming/R_analysis/ki67/"
  dir.create(dird_apo)
  dird_ki67<-paste(dird_apo,"spheroid_shells_ki67", sep="")
  dird_dna<-paste(dird_apo,"spheroid_shells_dna", sep="")
  dir.create(dird_ki67)
  dir.create(dird_dna)
}

# read in all the required features for spheroids
{
  plot(0:1,0:1)
  #T1
  {
    geometrical_data<-combine_sample_sets(path_to_ex[1],samples,data_types[1],"tsv",exp[1],trial[1])
    T1DNA_int_data<-combine_sample_sets(path_to_ex[1],samples,data_types[2],"tsv",exp[1],trial[1])
    T1KI67_int_data<-combine_sample_sets(path_to_ex[1],samples,data_types[3],"tsv",exp[1],trial[1])
    geo_int_2D_data<-combine_sample_sets_2d(path_to_ex[1],samples,data_types_2d[1],exp[1],trial[1])
    
    T1_2d_data_KI67<-combine_sample_sets_2d(path_to_ex[1],samples,data_types_2d[2],exp[1],trial[1])
    
    shell_Ki67_T1<-shell_analsysis(path_to_ex[1],samples,"/Shells/KI67",exp[1],trial[1],dird_ki67)
    shell_dna_T1<-shell_analsysis(path_to_ex[1],samples,"/Shells/DNA",exp[1],trial[1],dird_dna)
    
    T1_dna_axial<-axial_analysis(path_to_ex[1],samples,exp[1],trial[1],"/DNA/")
    T1_Ki67_axial<-axial_analysis(path_to_ex[1],samples,exp[1],trial[1],"/KI67/")
    
    T1_Ki67_angle<-angle_analsysis(path_to_ex[1],samples,"/KI67",exp[1],trial[1])
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
    T2KI67_int_data<-combine_sample_sets(path_to_ex[2],samples,data_types[3],"tsv",exp[2],trial[2])
    geo_int_2D_data<-combine_sample_sets_2d(path_to_ex[2],samples,data_types_2d[1],exp[2],trial[2])
    
    T2_2d_data_KI67<-combine_sample_sets_2d(path_to_ex[2],samples,data_types_2d[2],exp[2],trial[2])
    
    shell_Ki67_T2<-shell_analsysis(path_to_ex[2],samples,"/Shells/KI67",exp[2],trial[2],dird_ki67)
    shell_dna_T2<-shell_analsysis(path_to_ex[2],samples,"/Shells/DNA",exp[2],trial[2],dird_dna)
    
    T2_dna_axial<-axial_analysis(path_to_ex[2],samples,exp[2],trial[2],"/DNA/")
    T2_Ki67_axial<-axial_analysis(path_to_ex[2],samples,exp[2],trial[2],"/KI67/")
    
    T2_Ki67_angle<-angle_analsysis(path_to_ex[2],samples,"/KI67",exp[2],trial[2])
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
    T3KI67_int_data<-combine_sample_sets(path_to_ex[3],samples,data_types[3],"tsv",exp[3],trial[3])
    geo_int_2D_data<-combine_sample_sets_2d(path_to_ex[3],samples,data_types_2d[1],exp[3],trial[3])
    
    T3_2d_data_KI67<-combine_sample_sets_2d(path_to_ex[3],samples,data_types_2d[2],exp[3],trial[3])
    
    shell_Ki67_T3<-shell_analsysis(path_to_ex[3],samples,"/Shells/KI67",exp[3],trial[3],dird_ki67)
    shell_dna_T3<-shell_analsysis(path_to_ex[3],samples,"/Shells/DNA",exp[3],trial[3],dird_dna)
    
    T3_dna_axial<-axial_analysis(path_to_ex[3],samples,exp[3],trial[3],"/DNA/")
    T3_Ki67_axial<-axial_analysis(path_to_ex[3],samples,exp[3],trial[3],"/KI67/")
    
    T3_Ki67_angle<-angle_analsysis(path_to_ex[3],samples,"/KI67",exp[3],trial[3])
    T3_dna_angle<-angle_analsysis(path_to_ex[3],samples,"/DNA",exp[3],trial[3])
    
    geo_int_2D_data$Label<-substring(geo_int_2D_data$Label,5, nchar(geo_int_2D_data$Label)-4)
    geometrical_data$Label<-substring(geometrical_data$Label,0, nchar(geometrical_data$Label)-1)
    T3combined<-merge(geometrical_data,geo_int_2D_data, by="Label")
    
    rm(geometrical_data,geo_int_2D_data)
    
  }
  
  spheroid_2dint_data_Ki67<-rbind(T1_2d_data_KI67,T2_2d_data_KI67,T3_2d_data_KI67)
  spheroid_2dint_data_Ki67$Label<-substring(spheroid_2dint_data_Ki67$Label,5,(nchar(spheroid_2dint_data_Ki67$Label)-4))
  spheroid_2dint_data_Ki67$Spheroid_id<-paste(spheroid_2dint_data_Ki67$trial,spheroid_2dint_data_Ki67$sample,
                                              substring(spheroid_2dint_data_Ki67$Label,1,2),sep="_")
  
  spheroid_int_data_Ki67<-rbind(T1KI67_int_data,T2KI67_int_data,T3KI67_int_data)
  spheroid_int_data_Ki67$Label<-substring(spheroid_int_data_Ki67$Label,1,(nchar(spheroid_int_data_Ki67$Label)-1))
  spheroid_int_data_Ki67$Spheroid_id<-paste(spheroid_int_data_Ki67$trial,spheroid_int_data_Ki67$sample,
                                              substring(spheroid_int_data_Ki67$Label,1,2),sep="_")
  
  spheroid_int_data_dna<-rbind(T1DNA_int_data,T2DNA_int_data,T3DNA_int_data)
  spheroid_int_data_dna$Label<-substring(spheroid_int_data_dna$Label,1,(nchar(spheroid_int_data_dna$Label)-1))
  spheroid_int_data_dna$Spheroid_id<-paste(spheroid_int_data_dna$trial,spheroid_int_data_dna$sample,
                                            substring(spheroid_int_data_dna$Label,1,2),sep="_")
  
  spheroid_radial_Ki67<-rbind(shell_Ki67_T1,shell_Ki67_T2,shell_Ki67_T3)
  spheroid_radial_Ki67$Spheroid_id<-paste(spheroid_radial_Ki67$trial,spheroid_radial_Ki67$sample,
                                           substring(spheroid_radial_Ki67$Label,1,2),sep="_")
  
  spheroid_radial_dna<-rbind(shell_dna_T1,shell_dna_T2,shell_dna_T3)
  spheroid_radial_dna$Spheroid_id<-paste(spheroid_radial_dna$trial,spheroid_radial_dna$sample,
                                          substring(spheroid_radial_dna$Label,1,2),sep="_")
  
  spheroid_axial_dna<-rbind(T1_dna_axial,T2_dna_axial,T3_dna_axial)
  spheroid_axial_dna$Spheroid_id<-paste(spheroid_axial_dna$trial,spheroid_axial_dna$sample,
                                         substring(spheroid_axial_dna$Label,1,2),sep="_")
  
  spheroid_axial_Ki67<-rbind(T1_Ki67_axial,T2_Ki67_axial,T3_Ki67_axial)
  spheroid_axial_Ki67$Spheroid_id<-paste(spheroid_axial_Ki67$trial,spheroid_axial_Ki67$sample,
                                        substring(spheroid_axial_Ki67$Label,1,2),sep="_")
  
  spheroid_angle_Ki67<-rbind(T1_Ki67_angle,T2_Ki67_angle,T3_Ki67_angle)
  spheroid_angle_Ki67$Spheroid_id<-paste(spheroid_angle_Ki67$trial,spheroid_angle_Ki67$sample,
                                         substring(spheroid_angle_Ki67$Label,1,2),sep="_")
  
  spheroid_angle_dna<-rbind(T1_dna_angle,T2_dna_angle,T3_dna_angle)
  spheroid_angle_dna$Spheroid_id<-paste(spheroid_angle_dna$trial,spheroid_angle_dna$sample,
                                         substring(spheroid_angle_dna$Label,1,2),sep="_")
  
  spheroid_combined<-rbind(T1combined,T2combined,T3combined)
  colnames(spheroid_combined)[86:87]<-c("sample","trial")
  spheroid_combined$Spheroid_id<-paste(spheroid_combined$trial,spheroid_combined$sample,
                                      substring(spheroid_combined$Label,1,2),sep="_")
  
  spheroid_int_data<-list( Ki67=spheroid_int_data_Ki67,
                          dna=spheroid_int_data_dna)
  spheroid_radial<-list(Ki67=spheroid_radial_Ki67,
                        dna=spheroid_radial_dna)
  spheroid_axial<-list( Ki67=spheroid_axial_Ki67,
                       dna=spheroid_axial_dna)
  spheroid_angle<-list(Ki67=spheroid_angle_Ki67,
                       dna=spheroid_angle_dna)
  spheroid_2dint_data<-list( Ki67=spheroid_2dint_data_Ki67)

  
  rm(T1combined,T2combined,T3combined,
     T1KI67_int_data,T2KI67_int_data,T3KI67_int_data,
     shell_Ki67_T1,shell_Ki67_T2,shell_Ki67_T3,
     shell_dna_T1,shell_dna_T2,shell_dna_T3,
     T1_dna_axial,T2_dna_axial,T3_dna_axial,
     T1_Ki67_axial,T2_Ki67_axial,T3_Ki67_axial,
      T1_Ki67_angle,T2_Ki67_angle,T3_Ki67_angle,
     T1_dna_angle,T2_dna_angle,T3_dna_angle,
     T1DNA_int_data,T2DNA_int_data,T3DNA_int_data,
    spheroid_angle_Ki67,spheroid_angle_dna,
      spheroid_axial_Ki67,spheroid_axial_dna,
      spheroid_radial_Ki67,spheroid_radial_dna,
      spheroid_int_data_Ki67,spheroid_int_data_dna,
     T1_2d_data_KI67,T2_2d_data_KI67,T3_2d_data_KI67,
      spheroid_2dint_data_Ki67)
  
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
    
    T1_KI67_int_data<-combine_sample_sets(path_to_ex[1],samples,nuc_data_types[6],"tsv",exp[1],trial[1])
    
    T1_2d_KI67_int_data<-combine_sample_sets_2d(path_to_ex[1],samples,nuc_data_types_2d[28],exp[1],trial[1])
    
    T1_cell_KI67_int_data<-combine_sample_sets(path_to_ex[1],samples,nuc_data_types[7],"tsv",exp[1],trial[1])
    T1_cell_DNA_int_data<-combine_sample_sets(path_to_ex[1],samples,nuc_data_types[8],"tsv",exp[1],trial[1])
    
    colnames(T1_compact_2D_data)[2]<-"Label"
    T1_compact_2D_data$Label<-substring(T1_compact_2D_data$Label,0, nchar(T1_compact_2D_data$Label)-4)
    T1_geo_int_2D_data$Label<-substring(T1_geo_int_2D_data$Label,5, nchar(T1_geo_int_2D_data$Label)-4)
    T1_geometrical_data$Label<-substring(T1_geometrical_data$Label,0, nchar(T1_geometrical_data$Label)-1)
    T1_DNA_int_data$Label<-substring(T1_DNA_int_data$Label,0, nchar(T1_DNA_int_data$Label)-1)
    T1_2d_KI67_int_data$Label<-substring(T1_2d_KI67_int_data$Label,5, nchar(T1_2d_KI67_int_data$Label)-4)
    T1_KI67_int_data$Label<-substring(T1_KI67_int_data$Label,0, nchar(T1_KI67_int_data$Label)-1)
    T1_cell_KI67_int_data$Label<-substring(T1_cell_KI67_int_data$Label,0, nchar(T1_cell_KI67_int_data$Label)-1)
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
    
    T2_KI67_int_data<-combine_sample_sets(path_to_ex[2],samples,nuc_data_types[6],"tsv",exp[2],trial[2])
    
    T2_2d_KI67_int_data<-combine_sample_sets_2d(path_to_ex[2],samples,nuc_data_types_2d[28],exp[2],trial[2])
    
    T2_cell_KI67_int_data<-combine_sample_sets(path_to_ex[2],samples,nuc_data_types[7],"tsv",exp[2],trial[2])
    T2_cell_DNA_int_data<-combine_sample_sets(path_to_ex[2],samples,nuc_data_types[8],"tsv",exp[2],trial[2])
    
    colnames(T2_compact_2D_data)[2]<-"Label"
    T2_compact_2D_data$Label<-substring(T2_compact_2D_data$Label,0, nchar(T2_compact_2D_data$Label)-4)
    T2_geo_int_2D_data$Label<-substring(T2_geo_int_2D_data$Label,5, nchar(T2_geo_int_2D_data$Label)-4)
    T2_geometrical_data$Label<-substring(T2_geometrical_data$Label,0, nchar(T2_geometrical_data$Label)-1)
    T2_2d_KI67_int_data$Label<-substring(T2_2d_KI67_int_data$Label,5, nchar(T2_2d_KI67_int_data$Label)-4)
    T2_DNA_int_data$Label<-substring(T2_DNA_int_data$Label,0, nchar(T2_DNA_int_data$Label)-1)
    T2_KI67_int_data$Label<-substring(T2_KI67_int_data$Label,0, nchar(T2_KI67_int_data$Label)-1)
    T2_cell_KI67_int_data$Label<-substring(T2_cell_KI67_int_data$Label,0, nchar(T2_cell_KI67_int_data$Label)-1)
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
    
    T3_KI67_int_data<-combine_sample_sets(path_to_ex[3],samples,nuc_data_types[6],"tsv",exp[3],trial[3])
    
    T3_2d_KI67_int_data<-combine_sample_sets_2d(path_to_ex[3],samples,nuc_data_types_2d[28],exp[3],trial[3])
    
    T3_cell_KI67_int_data<-combine_sample_sets(path_to_ex[3],samples,nuc_data_types[7],"tsv",exp[3],trial[3])
    T3_cell_DNA_int_data<-combine_sample_sets(path_to_ex[3],samples,nuc_data_types[8],"tsv",exp[3],trial[3])
    
    colnames(T3_compact_2D_data)[2]<-"Label"
    T3_compact_2D_data$Label<-substring(T3_compact_2D_data$Label,0, nchar(T3_compact_2D_data$Label)-4)
    T3_geo_int_2D_data$Label<-substring(T3_geo_int_2D_data$Label,5, nchar(T3_geo_int_2D_data$Label)-4)
    T3_geometrical_data$Label<-substring(T3_geometrical_data$Label,0, nchar(T3_geometrical_data$Label)-1)
    T3_2d_KI67_int_data$Label<-substring(T3_2d_KI67_int_data$Label,5, nchar(T3_2d_KI67_int_data$Label)-4)
    T3_DNA_int_data$Label<-substring(T3_DNA_int_data$Label,0, nchar(T3_DNA_int_data$Label)-1)
    T3_KI67_int_data$Label<-substring(T3_KI67_int_data$Label,0, nchar(T3_KI67_int_data$Label)-1)
    T3_cell_KI67_int_data$Label<-substring(T3_cell_KI67_int_data$Label,0, nchar(T3_cell_KI67_int_data$Label)-1)
    T3_cell_DNA_int_data$Label<-substring(T3_cell_DNA_int_data$Label,0, nchar(T3_cell_DNA_int_data$Label)-1)
    
    T3_combined<-merge(T3_geometrical_data,T3_geo_int_2D_data[,2:36], by="Label")
    T3_combined<-merge(T3_combined,T3_compact_2D_data[,2:8], by="Label")
    T3_combined<-merge(T3_combined,T3_DNA_int_data[,c(3:5,12:18)], by="Label")
    T3_combined<-merge(T3_combined,T3_edf_2D_data[,1:11], by="Label")
    T3_combined<-cbind(T3_combined,T3_gclm[,-c(1,128:131)])
    
    rm(T3_geometrical_data,T3_geo_int_2D_data,T3_compact_2D_data,T3_edf_2D_data,T3_gclm)
  }
  
  nuclei_2dint_data_KI67<-rbind(T1_2d_KI67_int_data,T2_2d_KI67_int_data,T3_2d_KI67_int_data)
  nuclei_2dint_data_KI67$Spheroid_id<-paste(nuclei_2dint_data_KI67$trial,nuclei_2dint_data_KI67$sample,
                                              substring(nuclei_2dint_data_KI67$Label,1,2),sep="_")
  spheroid_nuclei_2dint_data_KI67<-nuc_median_spheroid_levels(nuclei_2dint_data_KI67,parameters_2d_int_data)
  
  nuclei_int_data_KI67<-rbind(T1_KI67_int_data,T2_KI67_int_data,T3_KI67_int_data)
  nuclei_int_data_KI67$Spheroid_id<-paste(nuclei_int_data_KI67$trial,nuclei_int_data_KI67$sample,
                                            substring(nuclei_int_data_KI67$Label,1,2),sep="_")
  spheroid_nuclei_int_data_KI67<-nuc_median_spheroid_levels(nuclei_int_data_KI67,parameters_int_data)
  
  nuclei_int_data_dna<-rbind(T1_DNA_int_data,T2_DNA_int_data,T3_DNA_int_data)
  nuclei_int_data_dna$Spheroid_id<-paste(nuclei_int_data_dna$trial,nuclei_int_data_dna$sample,
                                          substring(nuclei_int_data_dna$Label,1,2),sep="_")
  spheroid_nuclei_int_data_dna<-nuc_median_spheroid_levels(nuclei_int_data_dna,parameters_int_data)
  
  cell_int_data_KI67<-rbind(T1_cell_KI67_int_data,T2_cell_KI67_int_data,T3_cell_KI67_int_data)
  cell_int_data_KI67$Spheroid_id<-paste(cell_int_data_KI67$trial,cell_int_data_KI67$sample,
                                         substring(cell_int_data_KI67$Label,1,2),sep="_")
  spheroid_cell_int_data_KI67<-nuc_median_spheroid_levels(cell_int_data_KI67,parameters_int_data)
  
  cell_int_data_dna<-rbind(T1_cell_DNA_int_data,T2_cell_DNA_int_data,T3_cell_DNA_int_data)
  cell_int_data_dna$Spheroid_id<-paste(cell_int_data_dna$trial,cell_int_data_dna$sample,
                                         substring(cell_int_data_dna$Label,1,2),sep="_")
  spheroid_cell_int_data_dna<-nuc_median_spheroid_levels(cell_int_data_dna,parameters_int_data)
  
  
  nucleus_combined<-rbind(T1_combined,T2_combined,T3_combined)
  nucleus_combined$Spheroid_id<-paste(nucleus_combined$trial,nucleus_combined$sample,
                                       substring(nucleus_combined$Label,1,2),sep="_")
  spheroid_nucleus_combined<-nuc_median_spheroid_levels(nucleus_combined,parameters_geo_texture_data)
  
  
  nuclei_2dint_data<-list(KI67=nuclei_2dint_data_KI67)
  nuclei_3dint_data<-list(KI67=nuclei_int_data_KI67,dna=nuclei_int_data_dna)
  cell_3dint_data<-list(KI67=cell_int_data_KI67, dna=cell_int_data_dna)
  
  spheroid_nuclei_2dint_data<-list(KI67=spheroid_nuclei_2dint_data_KI67)
  spheroid_nuclei_3dint_data<-list(KI67=spheroid_nuclei_int_data_KI67,dna=spheroid_nuclei_int_data_dna)
  spheroid_cell_3dint_data<-list(KI67=spheroid_cell_int_data_KI67, dna=spheroid_cell_int_data_dna)
  
  
  rm(T1_2d_KI67_int_data,T2_2d_KI67_int_data,T3_2d_KI67_int_data,
     T1_KI67_int_data,T2_KI67_int_data,T3_KI67_int_data,
     T1_DNA_int_data,T2_DNA_int_data,T3_DNA_int_data,
     T1_cell_KI67_int_data,T2_cell_KI67_int_data,T3_cell_KI67_int_data,
     T1_cell_DNA_int_data,T2_cell_DNA_int_data,T3_cell_DNA_int_data,
     T1_combined,T2_combined,T3_combined,
     nuclei_2dint_data_KI67, nuclei_int_data_KI67,
     nuclei_int_data_dna,cell_int_data_KI67,cell_int_data_dna,
     spheroid_nuclei_2dint_data_KI67,spheroid_nuclei_int_data_KI67,spheroid_nuclei_int_data_dna,
     spheroid_cell_int_data_dna,spheroid_cell_int_data_KI67)
  
}


#plotting
{
  dird<-"E:/HMF3A_reprogramming/R_analysis/Ki67/black_bg_plots/"
  dir.create(dird)
  setwd(dird)
  cols<-colorRampPalette(c("red","yellow"))( 5) 
  #Spheoid Proj_Int
  {
    
    dire<-paste(dird,"spheroid_proj_int/",sep="")
    dir.create(dire)
    setwd(dire)
    
      
    plot_black_bg_histogram_timeline(subset(spheroid_2dint_data$Ki67,spheroid_2dint_data$Ki67$trial=="T1"), cols,parameters_2d_int_data,paste("Ki67 ",labelling_2d_int_data,sep=""),"Ki67_T1_")
    plot_black_bg_histogram_timeline(subset(spheroid_2dint_data$Ki67,spheroid_2dint_data$Ki67$trial=="T2"), cols,parameters_2d_int_data,paste("Ki67 ",labelling_2d_int_data,sep=""),"Ki67_T2_")
    plot_black_bg_histogram_timeline(subset(spheroid_2dint_data$Ki67,spheroid_2dint_data$Ki67$trial=="T3"), cols,parameters_2d_int_data,paste("Ki67 ",labelling_2d_int_data,sep=""),"Ki67_T3_")
    
    plot_black_bg_boxplot_timeline(subset(spheroid_2dint_data$Ki67,spheroid_2dint_data$Ki67$trial=="T1"), cols,parameters_2d_int_data,paste("Ki67 ",labelling_2d_int_data,sep=""),"Ki67_T1_")
    plot_black_bg_boxplot_timeline(subset(spheroid_2dint_data$Ki67,spheroid_2dint_data$Ki67$trial=="T2"), cols,parameters_2d_int_data,paste("Ki67 ",labelling_2d_int_data,sep=""),"Ki67_T2_")
    plot_black_bg_boxplot_timeline(subset(spheroid_2dint_data$Ki67,spheroid_2dint_data$Ki67$trial=="T3"), cols,parameters_2d_int_data,paste("Ki67 ",labelling_2d_int_data,sep=""),"Ki67_T3_")
    
    plot_black_bg_barplot_timeline(spheroid_2dint_data$Ki67, cols,parameters_2d_int_data,paste("Ki67 ",labelling_2d_int_data,sep=""),"Ki67_")
    
    
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
    
     plot_black_bg_histogram_timeline(subset(spheroid_int_data$Ki67,spheroid_int_data$Ki67$trial=="T1"), cols,parameters_int_data,paste("Ki67 ",labelling_int_data,sep=""),"Ki67_T1_")
    plot_black_bg_histogram_timeline(subset(spheroid_int_data$Ki67,spheroid_int_data$Ki67$trial=="T2"), cols,parameters_int_data,paste("Ki67 ",labelling_int_data,sep=""),"Ki67_T2_")
    plot_black_bg_histogram_timeline(subset(spheroid_int_data$Ki67,spheroid_int_data$Ki67$trial=="T3"), cols,parameters_int_data,paste("Ki67 ",labelling_int_data,sep=""),"Ki67_T3_")
    
    plot_black_bg_boxplot_timeline(subset(spheroid_int_data$Ki67,spheroid_int_data$Ki67$trial=="T1"), cols,parameters_int_data,paste("Ki67 ",labelling_int_data,sep=""),"Ki67_T1_")
    plot_black_bg_boxplot_timeline(subset(spheroid_int_data$Ki67,spheroid_int_data$Ki67$trial=="T2"), cols,parameters_int_data,paste("Ki67 ",labelling_int_data,sep=""),"Ki67_T2_")
    plot_black_bg_boxplot_timeline(subset(spheroid_int_data$Ki67,spheroid_int_data$Ki67$trial=="T3"), cols,parameters_int_data,paste("Ki67 ",labelling_int_data,sep=""),"Ki67_T3_")
    
    plot_black_bg_barplot_timeline(spheroid_int_data$Ki67, cols,parameters_int_data,paste("Ki67 ",labelling_int_data,sep=""),"Ki67_")
    
    
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
    
    
    plot_black_bg_histogram_timeline(subset(spheroid_angle$Ki67,spheroid_angle$Ki67$trial=="T1"), cols,parameters_angle_data,paste("Ki67 ",labelling_angle_data,sep=""),"Ki67_T1_")
    plot_black_bg_histogram_timeline(subset(spheroid_angle$Ki67,spheroid_angle$Ki67$trial=="T2"), cols,parameters_angle_data,paste("Ki67 ",labelling_angle_data,sep=""),"Ki67_T2_")
    plot_black_bg_histogram_timeline(subset(spheroid_angle$Ki67,spheroid_angle$Ki67$trial=="T3"), cols,parameters_angle_data,paste("Ki67 ",labelling_angle_data,sep=""),"Ki67_T3_")
    
    plot_black_bg_boxplot_timeline(subset(spheroid_angle$Ki67,spheroid_angle$Ki67$trial=="T1"), cols,parameters_angle_data,paste("Ki67 ",labelling_angle_data,sep=""),"Ki67_T1_")
    plot_black_bg_boxplot_timeline(subset(spheroid_angle$Ki67,spheroid_angle$Ki67$trial=="T2"), cols,parameters_angle_data,paste("Ki67 ",labelling_angle_data,sep=""),"Ki67_T2_")
    plot_black_bg_boxplot_timeline(subset(spheroid_angle$Ki67,spheroid_angle$Ki67$trial=="T3"), cols,parameters_angle_data,paste("Ki67 ",labelling_angle_data,sep=""),"Ki67_T3_")
    
    plot_black_bg_barplot_timeline(spheroid_angle$Ki67, cols,parameters_angle_data,paste("Ki67 ",labelling_angle_data,sep=""),"Ki67_")
    
    
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
    
    plot_black_bg_histogram_timeline(subset(spheroid_axial$Ki67,spheroid_axial$Ki67$trial=="T1"), cols,parameters_axial_data,paste("Ki67 ",labelling_axial_data,sep=""),"Ki67_T1_")
    plot_black_bg_histogram_timeline(subset(spheroid_axial$Ki67,spheroid_axial$Ki67$trial=="T2"), cols,parameters_axial_data,paste("Ki67 ",labelling_axial_data,sep=""),"Ki67_T2_")
    plot_black_bg_histogram_timeline(subset(spheroid_axial$Ki67,spheroid_axial$Ki67$trial=="T3"), cols,parameters_axial_data,paste("Ki67 ",labelling_axial_data,sep=""),"Ki67_T3_")
    
    plot_black_bg_boxplot_timeline(subset(spheroid_axial$Ki67,spheroid_axial$Ki67$trial=="T1"), cols,parameters_axial_data,paste("Ki67 ",labelling_axial_data,sep=""),"Ki67_T1_")
    plot_black_bg_boxplot_timeline(subset(spheroid_axial$Ki67,spheroid_axial$Ki67$trial=="T2"), cols,parameters_axial_data,paste("Ki67 ",labelling_axial_data,sep=""),"Ki67_T2_")
    plot_black_bg_boxplot_timeline(subset(spheroid_axial$Ki67,spheroid_axial$Ki67$trial=="T3"), cols,parameters_axial_data,paste("Ki67 ",labelling_axial_data,sep=""),"Ki67_T3_")
    
    plot_black_bg_barplot_timeline(spheroid_axial$Ki67, cols,parameters_axial_data,paste("Ki67 ",labelling_axial_data,sep=""),"Ki67_")
    
    
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
    
    plot_black_bg_histogram_timeline(subset(spheroid_radial$Ki67,spheroid_radial$Ki67$trial=="T1"), cols,parameters_radial_data,paste("Ki67 ",labelling_radial_data,sep=""),"Ki67_T1_")
    plot_black_bg_histogram_timeline(subset(spheroid_radial$Ki67,spheroid_radial$Ki67$trial=="T2"), cols,parameters_radial_data,paste("Ki67 ",labelling_radial_data,sep=""),"Ki67_T2_")
    plot_black_bg_histogram_timeline(subset(spheroid_radial$Ki67,spheroid_radial$Ki67$trial=="T3"), cols,parameters_radial_data,paste("Ki67 ",labelling_radial_data,sep=""),"Ki67_T3_")
    
    plot_black_bg_boxplot_timeline(subset(spheroid_radial$Ki67,spheroid_radial$Ki67$trial=="T1"), cols,parameters_radial_data,paste("Ki67 ",labelling_radial_data,sep=""),"Ki67_T1_")
    plot_black_bg_boxplot_timeline(subset(spheroid_radial$Ki67,spheroid_radial$Ki67$trial=="T2"), cols,parameters_radial_data,paste("Ki67 ",labelling_radial_data,sep=""),"Ki67_T2_")
    plot_black_bg_boxplot_timeline(subset(spheroid_radial$Ki67,spheroid_radial$Ki67$trial=="T3"), cols,parameters_radial_data,paste("Ki67 ",labelling_radial_data,sep=""),"Ki67_T3_")
    
    plot_black_bg_barplot_timeline(spheroid_radial$Ki67, cols,parameters_radial_data,paste("Ki67 ",labelling_radial_data,sep=""),"Ki67_")
    
    
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
    
    plot_black_bg_histogram_timeline(subset(nuclei_2dint_data$KI67,nuclei_2dint_data$KI67$trial=="T1"), cols,parameters_2d_int_data,paste("KI67 ",nuc_labelling_2d_int_data,sep=""),"KI67_T1_")
    plot_black_bg_histogram_timeline(subset(nuclei_2dint_data$KI67,nuclei_2dint_data$KI67$trial=="T2"), cols,parameters_2d_int_data,paste("KI67 ",nuc_labelling_2d_int_data,sep=""),"KI67_T2_")
    plot_black_bg_histogram_timeline(subset(nuclei_2dint_data$KI67,nuclei_2dint_data$KI67$trial=="T3"), cols,parameters_2d_int_data,paste("KI67 ",nuc_labelling_2d_int_data,sep=""),"KI67_T3_")
    
    plot_black_bg_boxplot_timeline(subset(nuclei_2dint_data$KI67,nuclei_2dint_data$KI67$trial=="T1"), cols,parameters_2d_int_data,paste("F-KI67 ",nuc_labelling_2d_int_data,sep=""),"KI67_T1_")
    plot_black_bg_boxplot_timeline(subset(nuclei_2dint_data$KI67,nuclei_2dint_data$KI67$trial=="T2"), cols,parameters_2d_int_data,paste("F-KI67 ",nuc_labelling_2d_int_data,sep=""),"KI67_T2_")
    plot_black_bg_boxplot_timeline(subset(nuclei_2dint_data$KI67,nuclei_2dint_data$KI67$trial=="T3"), cols,parameters_2d_int_data,paste("F-KI67 ",nuc_labelling_2d_int_data,sep=""),"KI67_T3_")
    
    plot_black_bg_barplot_timeline(nuclei_2dint_data$KI67, cols,parameters_2d_int_data,paste("KI67 ",nuc_labelling_2d_int_data,sep=""),"KI67_")
    
    
    
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
    
    plot_black_bg_histogram_timeline(subset(nuclei_3dint_data$KI67,nuclei_3dint_data$KI67$trial=="T1"), cols,parameters_int_data,paste("KI67 ",nuclabelling_int_data,sep=""),"KI67_T1_")
    plot_black_bg_histogram_timeline(subset(nuclei_3dint_data$KI67,nuclei_3dint_data$KI67$trial=="T2"), cols,parameters_int_data,paste("KI67 ",nuclabelling_int_data,sep=""),"KI67_T2_")
    plot_black_bg_histogram_timeline(subset(nuclei_3dint_data$KI67,nuclei_3dint_data$KI67$trial=="T3"), cols,parameters_int_data,paste("KI67 ",nuclabelling_int_data,sep=""),"KI67_T3_")
    
    plot_black_bg_boxplot_timeline(subset(nuclei_3dint_data$KI67,nuclei_3dint_data$KI67$trial=="T1"), cols,parameters_int_data,paste("KI67 ",nuclabelling_int_data,sep=""),"KI67_T1_")
    plot_black_bg_boxplot_timeline(subset(nuclei_3dint_data$KI67,nuclei_3dint_data$KI67$trial=="T2"), cols,parameters_int_data,paste("KI67 ",nuclabelling_int_data,sep=""),"KI67_T2_")
    plot_black_bg_boxplot_timeline(subset(nuclei_3dint_data$KI67,nuclei_3dint_data$KI67$trial=="T3"), cols,parameters_int_data,paste("KI67 ",nuclabelling_int_data,sep=""),"KI67_T3_")
    
    plot_black_bg_barplot_timeline(nuclei_3dint_data$KI67, cols,parameters_int_data,paste("KI67 ",nuclabelling_int_data,sep=""),"KI67_")
    
    
    
  }
  #Cellular 3D_Int
  {
    dire<-paste(dird,"cellular_3d_int/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_black_bg_histogram_timeline(subset(cell_3dint_data$dna,cell_3dint_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T1_")
    plot_black_bg_histogram_timeline(subset(cell_3dint_data$dna,cell_3dint_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T2_")
    plot_black_bg_histogram_timeline(subset(cell_3dint_data$dna,cell_3dint_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T3_")
    
    plot_black_bg_boxplot_timeline(subset(cell_3dint_data$dna,cell_3dint_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T1_")
    plot_black_bg_boxplot_timeline(subset(cell_3dint_data$dna,cell_3dint_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T2_")
    plot_black_bg_boxplot_timeline(subset(cell_3dint_data$dna,cell_3dint_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T3_")
    
    plot_black_bg_barplot_timeline(cell_3dint_data$dna, cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_")
    
    plot_black_bg_histogram_timeline(subset(cell_3dint_data$KI67,cell_3dint_data$KI67$trial=="T1"), cols,parameters_int_data,paste("KI67 ",celllabelling_int_data,sep=""),"KI67_T1_")
    plot_black_bg_histogram_timeline(subset(cell_3dint_data$KI67,cell_3dint_data$KI67$trial=="T2"), cols,parameters_int_data,paste("KI67 ",celllabelling_int_data,sep=""),"KI67_T2_")
    plot_black_bg_histogram_timeline(subset(cell_3dint_data$KI67,cell_3dint_data$KI67$trial=="T3"), cols,parameters_int_data,paste("KI67 ",celllabelling_int_data,sep=""),"KI67_T3_")
    
    plot_black_bg_boxplot_timeline(subset(cell_3dint_data$KI67,cell_3dint_data$KI67$trial=="T1"), cols,parameters_int_data,paste("KI67 ",celllabelling_int_data,sep=""),"KI67_T1_")
    plot_black_bg_boxplot_timeline(subset(cell_3dint_data$KI67,cell_3dint_data$KI67$trial=="T2"), cols,parameters_int_data,paste("KI67 ",celllabelling_int_data,sep=""),"KI67_T2_")
    plot_black_bg_boxplot_timeline(subset(cell_3dint_data$KI67,cell_3dint_data$KI67$trial=="T3"), cols,parameters_int_data,paste("KI67 ",celllabelling_int_data,sep=""),"KI67_T3_")
    
    plot_black_bg_barplot_timeline(cell_3dint_data$KI67, cols,parameters_int_data,paste("KI67 ",celllabelling_int_data,sep=""),"KI67_")
    
    
    
    
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
    
    plot_black_bg_histogram_timeline(subset(spheroid_nuclei_2dint_data$KI67,spheroid_nuclei_2dint_data$KI67$trial=="T1"), cols,parameters_2d_int_data,paste("KI67 ",nuc_labelling_2d_int_data,sep=""),"KI67_T1_")
    plot_black_bg_histogram_timeline(subset(spheroid_nuclei_2dint_data$KI67,spheroid_nuclei_2dint_data$KI67$trial=="T2"), cols,parameters_2d_int_data,paste("KI67 ",nuc_labelling_2d_int_data,sep=""),"KI67_T2_")
    plot_black_bg_histogram_timeline(subset(spheroid_nuclei_2dint_data$KI67,spheroid_nuclei_2dint_data$KI67$trial=="T3"), cols,parameters_2d_int_data,paste("KI67 ",nuc_labelling_2d_int_data,sep=""),"KI67_T3_")
    
    plot_black_bg_boxplot_timeline(subset(spheroid_nuclei_2dint_data$KI67,spheroid_nuclei_2dint_data$KI67$trial=="T1"), cols,parameters_2d_int_data,paste("KI67 ",nuc_labelling_2d_int_data,sep=""),"KI67_T1_")
    plot_black_bg_boxplot_timeline(subset(spheroid_nuclei_2dint_data$KI67,spheroid_nuclei_2dint_data$KI67$trial=="T2"), cols,parameters_2d_int_data,paste("KI67 ",nuc_labelling_2d_int_data,sep=""),"KI67_T2_")
    plot_black_bg_boxplot_timeline(subset(spheroid_nuclei_2dint_data$KI67,spheroid_nuclei_2dint_data$KI67$trial=="T3"), cols,parameters_2d_int_data,paste("KI67 ",nuc_labelling_2d_int_data,sep=""),"KI67_T3_")
    
    plot_black_bg_barplot_timeline(spheroid_nuclei_2dint_data$KI67, cols,parameters_2d_int_data,paste("KI67 ",nuc_labelling_2d_int_data,sep=""),"KI67_")
    
    
    
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
    
    plot_black_bg_histogram_timeline(subset(spheroid_nuclei_3dint_data$KI67,spheroid_nuclei_3dint_data$KI67$trial=="T1"), cols,parameters_int_data,paste("KI67 ",nuclabelling_int_data,sep=""),"KI67_T1_")
    plot_black_bg_histogram_timeline(subset(spheroid_nuclei_3dint_data$KI67,spheroid_nuclei_3dint_data$KI67$trial=="T2"), cols,parameters_int_data,paste("KI67 ",nuclabelling_int_data,sep=""),"KI67_T2_")
    plot_black_bg_histogram_timeline(subset(spheroid_nuclei_3dint_data$KI67,spheroid_nuclei_3dint_data$KI67$trial=="T3"), cols,parameters_int_data,paste("KI67 ",nuclabelling_int_data,sep=""),"KI67_T3_")
    
    plot_black_bg_boxplot_timeline(subset(spheroid_nuclei_3dint_data$KI67,spheroid_nuclei_3dint_data$KI67$trial=="T1"), cols,parameters_int_data,paste("KI67 ",nuclabelling_int_data,sep=""),"KI67_T1_")
    plot_black_bg_boxplot_timeline(subset(spheroid_nuclei_3dint_data$KI67,spheroid_nuclei_3dint_data$KI67$trial=="T2"), cols,parameters_int_data,paste("KI67 ",nuclabelling_int_data,sep=""),"KI67_T2_")
    plot_black_bg_boxplot_timeline(subset(spheroid_nuclei_3dint_data$KI67,spheroid_nuclei_3dint_data$KI67$trial=="T3"), cols,parameters_int_data,paste("KI67 ",nuclabelling_int_data,sep=""),"KI67_T3_")
    
    plot_black_bg_barplot_timeline(spheroid_nuclei_3dint_data$KI67, cols,parameters_int_data,paste("KI67 ",nuclabelling_int_data,sep=""),"KI67_")
    
    
    
  }
  #Spheroid_Cellular 3D_Int
  {
    dire<-paste(dird,"median_cellular_3d_int/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_black_bg_histogram_timeline(subset(spheroid_cell_3dint_data$dna,spheroid_cell_3dint_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T1_")
    plot_black_bg_histogram_timeline(subset(spheroid_cell_3dint_data$dna,spheroid_cell_3dint_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T2_")
    plot_black_bg_histogram_timeline(subset(spheroid_cell_3dint_data$dna,spheroid_cell_3dint_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T3_")
    
    plot_black_bg_boxplot_timeline(subset(spheroid_cell_3dint_data$dna,spheroid_cell_3dint_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T1_")
    plot_black_bg_boxplot_timeline(subset(spheroid_cell_3dint_data$dna,spheroid_cell_3dint_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T2_")
    plot_black_bg_boxplot_timeline(subset(spheroid_cell_3dint_data$dna,spheroid_cell_3dint_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T3_")
    
    plot_black_bg_barplot_timeline(spheroid_cell_3dint_data$dna, cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_")
    
    plot_black_bg_histogram_timeline(subset(spheroid_cell_3dint_data$KI67,spheroid_cell_3dint_data$KI67$trial=="T1"), cols,parameters_int_data,paste("KI67 ",celllabelling_int_data,sep=""),"KI67_T1_")
    plot_black_bg_histogram_timeline(subset(spheroid_cell_3dint_data$KI67,spheroid_cell_3dint_data$KI67$trial=="T2"), cols,parameters_int_data,paste("KI67 ",celllabelling_int_data,sep=""),"KI67_T2_")
    plot_black_bg_histogram_timeline(subset(spheroid_cell_3dint_data$KI67,spheroid_cell_3dint_data$KI67$trial=="T3"), cols,parameters_int_data,paste("KI67 ",celllabelling_int_data,sep=""),"KI67_T3_")
    
    plot_black_bg_boxplot_timeline(subset(spheroid_cell_3dint_data$KI67,spheroid_cell_3dint_data$KI67$trial=="T1"), cols,parameters_int_data,paste("KI67 ",celllabelling_int_data,sep=""),"KI67_T1_")
    plot_black_bg_boxplot_timeline(subset(spheroid_cell_3dint_data$KI67,spheroid_cell_3dint_data$KI67$trial=="T2"), cols,parameters_int_data,paste("KI67 ",celllabelling_int_data,sep=""),"KI67_T2_")
    plot_black_bg_boxplot_timeline(subset(spheroid_cell_3dint_data$KI67,spheroid_cell_3dint_data$KI67$trial=="T3"), cols,parameters_int_data,paste("KI67 ",celllabelling_int_data,sep=""),"KI67_T3_")
    
    plot_black_bg_barplot_timeline(spheroid_cell_3dint_data$KI67, cols,parameters_int_data,paste("KI67 ",celllabelling_int_data,sep=""),"KI67_")
    
    
    
    
  }
  
  dird<-"E:/HMF3A_reprogramming/R_analysis/Ki67/white_bg_plots/"
  dir.create(dird)
  setwd(dird)
  cols<-colorRampPalette(c("red","yellow"))( 5) 
  #Spheoid Proj_Int
  {
    
    dire<-paste(dird,"spheroid_proj_int/",sep="")
    dir.create(dire)
    setwd(dire)
    
    
    plot_white_bg_histogram_timeline(subset(spheroid_2dint_data$Ki67,spheroid_2dint_data$Ki67$trial=="T1"), cols,parameters_2d_int_data,paste("Ki67 ",labelling_2d_int_data,sep=""),"Ki67_T1_")
    plot_white_bg_histogram_timeline(subset(spheroid_2dint_data$Ki67,spheroid_2dint_data$Ki67$trial=="T2"), cols,parameters_2d_int_data,paste("Ki67 ",labelling_2d_int_data,sep=""),"Ki67_T2_")
    plot_white_bg_histogram_timeline(subset(spheroid_2dint_data$Ki67,spheroid_2dint_data$Ki67$trial=="T3"), cols,parameters_2d_int_data,paste("Ki67 ",labelling_2d_int_data,sep=""),"Ki67_T3_")
    
    plot_white_bg_boxplot_timeline(subset(spheroid_2dint_data$Ki67,spheroid_2dint_data$Ki67$trial=="T1"), cols,parameters_2d_int_data,paste("Ki67 ",labelling_2d_int_data,sep=""),"Ki67_T1_")
    plot_white_bg_boxplot_timeline(subset(spheroid_2dint_data$Ki67,spheroid_2dint_data$Ki67$trial=="T2"), cols,parameters_2d_int_data,paste("Ki67 ",labelling_2d_int_data,sep=""),"Ki67_T2_")
    plot_white_bg_boxplot_timeline(subset(spheroid_2dint_data$Ki67,spheroid_2dint_data$Ki67$trial=="T3"), cols,parameters_2d_int_data,paste("Ki67 ",labelling_2d_int_data,sep=""),"Ki67_T3_")
    
    plot_white_bg_barplot_timeline(spheroid_2dint_data$Ki67, cols,parameters_2d_int_data,paste("Ki67 ",labelling_2d_int_data,sep=""),"Ki67_")
    
    
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
    
    plot_white_bg_histogram_timeline(subset(spheroid_int_data$Ki67,spheroid_int_data$Ki67$trial=="T1"), cols,parameters_int_data,paste("Ki67 ",labelling_int_data,sep=""),"Ki67_T1_")
    plot_white_bg_histogram_timeline(subset(spheroid_int_data$Ki67,spheroid_int_data$Ki67$trial=="T2"), cols,parameters_int_data,paste("Ki67 ",labelling_int_data,sep=""),"Ki67_T2_")
    plot_white_bg_histogram_timeline(subset(spheroid_int_data$Ki67,spheroid_int_data$Ki67$trial=="T3"), cols,parameters_int_data,paste("Ki67 ",labelling_int_data,sep=""),"Ki67_T3_")
    
    plot_white_bg_boxplot_timeline(subset(spheroid_int_data$Ki67,spheroid_int_data$Ki67$trial=="T1"), cols,parameters_int_data,paste("Ki67 ",labelling_int_data,sep=""),"Ki67_T1_")
    plot_white_bg_boxplot_timeline(subset(spheroid_int_data$Ki67,spheroid_int_data$Ki67$trial=="T2"), cols,parameters_int_data,paste("Ki67 ",labelling_int_data,sep=""),"Ki67_T2_")
    plot_white_bg_boxplot_timeline(subset(spheroid_int_data$Ki67,spheroid_int_data$Ki67$trial=="T3"), cols,parameters_int_data,paste("Ki67 ",labelling_int_data,sep=""),"Ki67_T3_")
    
    plot_white_bg_barplot_timeline(spheroid_int_data$Ki67, cols,parameters_int_data,paste("Ki67 ",labelling_int_data,sep=""),"Ki67_")
    
    
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
    
    
    plot_white_bg_histogram_timeline(subset(spheroid_angle$Ki67,spheroid_angle$Ki67$trial=="T1"), cols,parameters_angle_data,paste("Ki67 ",labelling_angle_data,sep=""),"Ki67_T1_")
    plot_white_bg_histogram_timeline(subset(spheroid_angle$Ki67,spheroid_angle$Ki67$trial=="T2"), cols,parameters_angle_data,paste("Ki67 ",labelling_angle_data,sep=""),"Ki67_T2_")
    plot_white_bg_histogram_timeline(subset(spheroid_angle$Ki67,spheroid_angle$Ki67$trial=="T3"), cols,parameters_angle_data,paste("Ki67 ",labelling_angle_data,sep=""),"Ki67_T3_")
    
    plot_white_bg_boxplot_timeline(subset(spheroid_angle$Ki67,spheroid_angle$Ki67$trial=="T1"), cols,parameters_angle_data,paste("Ki67 ",labelling_angle_data,sep=""),"Ki67_T1_")
    plot_white_bg_boxplot_timeline(subset(spheroid_angle$Ki67,spheroid_angle$Ki67$trial=="T2"), cols,parameters_angle_data,paste("Ki67 ",labelling_angle_data,sep=""),"Ki67_T2_")
    plot_white_bg_boxplot_timeline(subset(spheroid_angle$Ki67,spheroid_angle$Ki67$trial=="T3"), cols,parameters_angle_data,paste("Ki67 ",labelling_angle_data,sep=""),"Ki67_T3_")
    
    plot_white_bg_barplot_timeline(spheroid_angle$Ki67, cols,parameters_angle_data,paste("Ki67 ",labelling_angle_data,sep=""),"Ki67_")
    
    
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
    
    plot_white_bg_histogram_timeline(subset(spheroid_axial$Ki67,spheroid_axial$Ki67$trial=="T1"), cols,parameters_axial_data,paste("Ki67 ",labelling_axial_data,sep=""),"Ki67_T1_")
    plot_white_bg_histogram_timeline(subset(spheroid_axial$Ki67,spheroid_axial$Ki67$trial=="T2"), cols,parameters_axial_data,paste("Ki67 ",labelling_axial_data,sep=""),"Ki67_T2_")
    plot_white_bg_histogram_timeline(subset(spheroid_axial$Ki67,spheroid_axial$Ki67$trial=="T3"), cols,parameters_axial_data,paste("Ki67 ",labelling_axial_data,sep=""),"Ki67_T3_")
    
    plot_white_bg_boxplot_timeline(subset(spheroid_axial$Ki67,spheroid_axial$Ki67$trial=="T1"), cols,parameters_axial_data,paste("Ki67 ",labelling_axial_data,sep=""),"Ki67_T1_")
    plot_white_bg_boxplot_timeline(subset(spheroid_axial$Ki67,spheroid_axial$Ki67$trial=="T2"), cols,parameters_axial_data,paste("Ki67 ",labelling_axial_data,sep=""),"Ki67_T2_")
    plot_white_bg_boxplot_timeline(subset(spheroid_axial$Ki67,spheroid_axial$Ki67$trial=="T3"), cols,parameters_axial_data,paste("Ki67 ",labelling_axial_data,sep=""),"Ki67_T3_")
    
    plot_white_bg_barplot_timeline(spheroid_axial$Ki67, cols,parameters_axial_data,paste("Ki67 ",labelling_axial_data,sep=""),"Ki67_")
    
    
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
    
    plot_white_bg_histogram_timeline(subset(spheroid_radial$Ki67,spheroid_radial$Ki67$trial=="T1"), cols,parameters_radial_data,paste("Ki67 ",labelling_radial_data,sep=""),"Ki67_T1_")
    plot_white_bg_histogram_timeline(subset(spheroid_radial$Ki67,spheroid_radial$Ki67$trial=="T2"), cols,parameters_radial_data,paste("Ki67 ",labelling_radial_data,sep=""),"Ki67_T2_")
    plot_white_bg_histogram_timeline(subset(spheroid_radial$Ki67,spheroid_radial$Ki67$trial=="T3"), cols,parameters_radial_data,paste("Ki67 ",labelling_radial_data,sep=""),"Ki67_T3_")
    
    plot_white_bg_boxplot_timeline(subset(spheroid_radial$Ki67,spheroid_radial$Ki67$trial=="T1"), cols,parameters_radial_data,paste("Ki67 ",labelling_radial_data,sep=""),"Ki67_T1_")
    plot_white_bg_boxplot_timeline(subset(spheroid_radial$Ki67,spheroid_radial$Ki67$trial=="T2"), cols,parameters_radial_data,paste("Ki67 ",labelling_radial_data,sep=""),"Ki67_T2_")
    plot_white_bg_boxplot_timeline(subset(spheroid_radial$Ki67,spheroid_radial$Ki67$trial=="T3"), cols,parameters_radial_data,paste("Ki67 ",labelling_radial_data,sep=""),"Ki67_T3_")
    
    plot_white_bg_barplot_timeline(spheroid_radial$Ki67, cols,parameters_radial_data,paste("Ki67 ",labelling_radial_data,sep=""),"Ki67_")
    
    
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
    
    plot_white_bg_histogram_timeline(subset(nuclei_2dint_data$KI67,nuclei_2dint_data$KI67$trial=="T1"), cols,parameters_2d_int_data,paste("KI67 ",nuc_labelling_2d_int_data,sep=""),"KI67_T1_")
    plot_white_bg_histogram_timeline(subset(nuclei_2dint_data$KI67,nuclei_2dint_data$KI67$trial=="T2"), cols,parameters_2d_int_data,paste("KI67 ",nuc_labelling_2d_int_data,sep=""),"KI67_T2_")
    plot_white_bg_histogram_timeline(subset(nuclei_2dint_data$KI67,nuclei_2dint_data$KI67$trial=="T3"), cols,parameters_2d_int_data,paste("KI67 ",nuc_labelling_2d_int_data,sep=""),"KI67_T3_")
    
    plot_white_bg_boxplot_timeline(subset(nuclei_2dint_data$KI67,nuclei_2dint_data$KI67$trial=="T1"), cols,parameters_2d_int_data,paste("F-KI67 ",nuc_labelling_2d_int_data,sep=""),"KI67_T1_")
    plot_white_bg_boxplot_timeline(subset(nuclei_2dint_data$KI67,nuclei_2dint_data$KI67$trial=="T2"), cols,parameters_2d_int_data,paste("F-KI67 ",nuc_labelling_2d_int_data,sep=""),"KI67_T2_")
    plot_white_bg_boxplot_timeline(subset(nuclei_2dint_data$KI67,nuclei_2dint_data$KI67$trial=="T3"), cols,parameters_2d_int_data,paste("F-KI67 ",nuc_labelling_2d_int_data,sep=""),"KI67_T3_")
    
    plot_white_bg_barplot_timeline(nuclei_2dint_data$KI67, cols,parameters_2d_int_data,paste("KI67 ",nuc_labelling_2d_int_data,sep=""),"KI67_")
    
    
    
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
    
    plot_white_bg_histogram_timeline(subset(nuclei_3dint_data$KI67,nuclei_3dint_data$KI67$trial=="T1"), cols,parameters_int_data,paste("KI67 ",nuclabelling_int_data,sep=""),"KI67_T1_")
    plot_white_bg_histogram_timeline(subset(nuclei_3dint_data$KI67,nuclei_3dint_data$KI67$trial=="T2"), cols,parameters_int_data,paste("KI67 ",nuclabelling_int_data,sep=""),"KI67_T2_")
    plot_white_bg_histogram_timeline(subset(nuclei_3dint_data$KI67,nuclei_3dint_data$KI67$trial=="T3"), cols,parameters_int_data,paste("KI67 ",nuclabelling_int_data,sep=""),"KI67_T3_")
    
    plot_white_bg_boxplot_timeline(subset(nuclei_3dint_data$KI67,nuclei_3dint_data$KI67$trial=="T1"), cols,parameters_int_data,paste("KI67 ",nuclabelling_int_data,sep=""),"KI67_T1_")
    plot_white_bg_boxplot_timeline(subset(nuclei_3dint_data$KI67,nuclei_3dint_data$KI67$trial=="T2"), cols,parameters_int_data,paste("KI67 ",nuclabelling_int_data,sep=""),"KI67_T2_")
    plot_white_bg_boxplot_timeline(subset(nuclei_3dint_data$KI67,nuclei_3dint_data$KI67$trial=="T3"), cols,parameters_int_data,paste("KI67 ",nuclabelling_int_data,sep=""),"KI67_T3_")
    
    plot_white_bg_barplot_timeline(nuclei_3dint_data$KI67, cols,parameters_int_data,paste("KI67 ",nuclabelling_int_data,sep=""),"KI67_")
    
    
    
  }
  #Cellular 3D_Int
  {
    dire<-paste(dird,"cellular_3d_int/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_white_bg_histogram_timeline(subset(cell_3dint_data$dna,cell_3dint_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T1_")
    plot_white_bg_histogram_timeline(subset(cell_3dint_data$dna,cell_3dint_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T2_")
    plot_white_bg_histogram_timeline(subset(cell_3dint_data$dna,cell_3dint_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T3_")
    
    plot_white_bg_boxplot_timeline(subset(cell_3dint_data$dna,cell_3dint_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T1_")
    plot_white_bg_boxplot_timeline(subset(cell_3dint_data$dna,cell_3dint_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T2_")
    plot_white_bg_boxplot_timeline(subset(cell_3dint_data$dna,cell_3dint_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T3_")
    
    plot_white_bg_barplot_timeline(cell_3dint_data$dna, cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_")
    
    plot_white_bg_histogram_timeline(subset(cell_3dint_data$KI67,cell_3dint_data$KI67$trial=="T1"), cols,parameters_int_data,paste("KI67 ",celllabelling_int_data,sep=""),"KI67_T1_")
    plot_white_bg_histogram_timeline(subset(cell_3dint_data$KI67,cell_3dint_data$KI67$trial=="T2"), cols,parameters_int_data,paste("KI67 ",celllabelling_int_data,sep=""),"KI67_T2_")
    plot_white_bg_histogram_timeline(subset(cell_3dint_data$KI67,cell_3dint_data$KI67$trial=="T3"), cols,parameters_int_data,paste("KI67 ",celllabelling_int_data,sep=""),"KI67_T3_")
    
    plot_white_bg_boxplot_timeline(subset(cell_3dint_data$KI67,cell_3dint_data$KI67$trial=="T1"), cols,parameters_int_data,paste("KI67 ",celllabelling_int_data,sep=""),"KI67_T1_")
    plot_white_bg_boxplot_timeline(subset(cell_3dint_data$KI67,cell_3dint_data$KI67$trial=="T2"), cols,parameters_int_data,paste("KI67 ",celllabelling_int_data,sep=""),"KI67_T2_")
    plot_white_bg_boxplot_timeline(subset(cell_3dint_data$KI67,cell_3dint_data$KI67$trial=="T3"), cols,parameters_int_data,paste("KI67 ",celllabelling_int_data,sep=""),"KI67_T3_")
    
    plot_white_bg_barplot_timeline(cell_3dint_data$KI67, cols,parameters_int_data,paste("KI67 ",celllabelling_int_data,sep=""),"KI67_")
    
    
    
    
  }
  #Spheroid_median_nuclear_properties
  {
    dire<-paste(dird,"median_nuclear_properties/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_white_bg_histogram_timeline(subset(spheroid_nucleus_combined,spheroid_nucleus_combined$trial=="T1"), cols,parameters_geo_texture_data,paste("Nuclear ",labelling_geo_texture_data,sep=""),"geo_T1_")
    plot_white_bg_histogram_timeline(subset(spheroid_nucleus_combined,spheroid_nucleus_combined$trial=="T2"), cols,parameters_geo_texture_data,paste("Nuclear ",labelling_geo_texture_data,sep=""),"geo_T2_")
    plot_white_bg_histogram_timeline(subset(spheroid_nucleus_combined,spheroid_nucleus_combined$trial=="T3"), cols,parameters_geo_texture_data,paste("Nuclear ",labelling_geo_texture_data,sep=""),"geo_T3_")
    
    plot_white_bg_boxplot_timeline(subset(spheroid_nucleus_combined,spheroid_nucleus_combined$trial=="T1"), cols,parameters_geo_texture_data,labelling_geo_texture_data,"geo_T1_")
    plot_white_bg_boxplot_timeline(subset(spheroid_nucleus_combined,spheroid_nucleus_combined$trial=="T2"), cols,parameters_geo_texture_data,labelling_geo_texture_data,"geo_T2_")
    plot_white_bg_boxplot_timeline(subset(spheroid_nucleus_combined,spheroid_nucleus_combined$trial=="T3"), cols,parameters_geo_texture_data,labelling_geo_texture_data,"geo_T3_")
    
    plot_white_bg_barplot_timeline(spheroid_nucleus_combined, cols,parameters_geo_texture_data,paste("Nuclear ",labelling_geo_texture_data,sep=""),"geo_")
    
  }
  #Spheroid_Nuclear Proj_Int
  {
    
    dire<-paste(dird,"median_nuclear_proj_int/",sep="")
    dir.create(dire)
    setwd(dire)
    
    plot_white_bg_histogram_timeline(subset(spheroid_nuclei_2dint_data$KI67,spheroid_nuclei_2dint_data$KI67$trial=="T1"), cols,parameters_2d_int_data,paste("KI67 ",nuc_labelling_2d_int_data,sep=""),"KI67_T1_")
    plot_white_bg_histogram_timeline(subset(spheroid_nuclei_2dint_data$KI67,spheroid_nuclei_2dint_data$KI67$trial=="T2"), cols,parameters_2d_int_data,paste("KI67 ",nuc_labelling_2d_int_data,sep=""),"KI67_T2_")
    plot_white_bg_histogram_timeline(subset(spheroid_nuclei_2dint_data$KI67,spheroid_nuclei_2dint_data$KI67$trial=="T3"), cols,parameters_2d_int_data,paste("KI67 ",nuc_labelling_2d_int_data,sep=""),"KI67_T3_")
    
    plot_white_bg_boxplot_timeline(subset(spheroid_nuclei_2dint_data$KI67,spheroid_nuclei_2dint_data$KI67$trial=="T1"), cols,parameters_2d_int_data,paste("KI67 ",nuc_labelling_2d_int_data,sep=""),"KI67_T1_")
    plot_white_bg_boxplot_timeline(subset(spheroid_nuclei_2dint_data$KI67,spheroid_nuclei_2dint_data$KI67$trial=="T2"), cols,parameters_2d_int_data,paste("KI67 ",nuc_labelling_2d_int_data,sep=""),"KI67_T2_")
    plot_white_bg_boxplot_timeline(subset(spheroid_nuclei_2dint_data$KI67,spheroid_nuclei_2dint_data$KI67$trial=="T3"), cols,parameters_2d_int_data,paste("KI67 ",nuc_labelling_2d_int_data,sep=""),"KI67_T3_")
    
    plot_white_bg_barplot_timeline(spheroid_nuclei_2dint_data$KI67, cols,parameters_2d_int_data,paste("KI67 ",nuc_labelling_2d_int_data,sep=""),"KI67_")
    
    
    
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
    
    plot_white_bg_histogram_timeline(subset(spheroid_nuclei_3dint_data$KI67,spheroid_nuclei_3dint_data$KI67$trial=="T1"), cols,parameters_int_data,paste("KI67 ",nuclabelling_int_data,sep=""),"KI67_T1_")
    plot_white_bg_histogram_timeline(subset(spheroid_nuclei_3dint_data$KI67,spheroid_nuclei_3dint_data$KI67$trial=="T2"), cols,parameters_int_data,paste("KI67 ",nuclabelling_int_data,sep=""),"KI67_T2_")
    plot_white_bg_histogram_timeline(subset(spheroid_nuclei_3dint_data$KI67,spheroid_nuclei_3dint_data$KI67$trial=="T3"), cols,parameters_int_data,paste("KI67 ",nuclabelling_int_data,sep=""),"KI67_T3_")
    
    plot_white_bg_boxplot_timeline(subset(spheroid_nuclei_3dint_data$KI67,spheroid_nuclei_3dint_data$KI67$trial=="T1"), cols,parameters_int_data,paste("KI67 ",nuclabelling_int_data,sep=""),"KI67_T1_")
    plot_white_bg_boxplot_timeline(subset(spheroid_nuclei_3dint_data$KI67,spheroid_nuclei_3dint_data$KI67$trial=="T2"), cols,parameters_int_data,paste("KI67 ",nuclabelling_int_data,sep=""),"KI67_T2_")
    plot_white_bg_boxplot_timeline(subset(spheroid_nuclei_3dint_data$KI67,spheroid_nuclei_3dint_data$KI67$trial=="T3"), cols,parameters_int_data,paste("KI67 ",nuclabelling_int_data,sep=""),"KI67_T3_")
    
    plot_white_bg_barplot_timeline(spheroid_nuclei_3dint_data$KI67, cols,parameters_int_data,paste("KI67 ",nuclabelling_int_data,sep=""),"KI67_")
    
    
    
  }
  #Spheroid_Cellular 3D_Int
  {
    dire<-paste(dird,"median_cellular_3d_int/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_white_bg_histogram_timeline(subset(spheroid_cell_3dint_data$dna,spheroid_cell_3dint_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T1_")
    plot_white_bg_histogram_timeline(subset(spheroid_cell_3dint_data$dna,spheroid_cell_3dint_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T2_")
    plot_white_bg_histogram_timeline(subset(spheroid_cell_3dint_data$dna,spheroid_cell_3dint_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T3_")
    
    plot_white_bg_boxplot_timeline(subset(spheroid_cell_3dint_data$dna,spheroid_cell_3dint_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T1_")
    plot_white_bg_boxplot_timeline(subset(spheroid_cell_3dint_data$dna,spheroid_cell_3dint_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T2_")
    plot_white_bg_boxplot_timeline(subset(spheroid_cell_3dint_data$dna,spheroid_cell_3dint_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T3_")
    
    plot_white_bg_barplot_timeline(spheroid_cell_3dint_data$dna, cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_")
    
    plot_white_bg_histogram_timeline(subset(spheroid_cell_3dint_data$KI67,spheroid_cell_3dint_data$KI67$trial=="T1"), cols,parameters_int_data,paste("KI67 ",celllabelling_int_data,sep=""),"KI67_T1_")
    plot_white_bg_histogram_timeline(subset(spheroid_cell_3dint_data$KI67,spheroid_cell_3dint_data$KI67$trial=="T2"), cols,parameters_int_data,paste("KI67 ",celllabelling_int_data,sep=""),"KI67_T2_")
    plot_white_bg_histogram_timeline(subset(spheroid_cell_3dint_data$KI67,spheroid_cell_3dint_data$KI67$trial=="T3"), cols,parameters_int_data,paste("KI67 ",celllabelling_int_data,sep=""),"KI67_T3_")
    
    plot_white_bg_boxplot_timeline(subset(spheroid_cell_3dint_data$KI67,spheroid_cell_3dint_data$KI67$trial=="T1"), cols,parameters_int_data,paste("KI67 ",celllabelling_int_data,sep=""),"KI67_T1_")
    plot_white_bg_boxplot_timeline(subset(spheroid_cell_3dint_data$KI67,spheroid_cell_3dint_data$KI67$trial=="T2"), cols,parameters_int_data,paste("KI67 ",celllabelling_int_data,sep=""),"KI67_T2_")
    plot_white_bg_boxplot_timeline(subset(spheroid_cell_3dint_data$KI67,spheroid_cell_3dint_data$KI67$trial=="T3"), cols,parameters_int_data,paste("KI67 ",celllabelling_int_data,sep=""),"KI67_T3_")
    
    plot_white_bg_barplot_timeline(spheroid_cell_3dint_data$KI67, cols,parameters_int_data,paste("KI67 ",celllabelling_int_data,sep=""),"KI67_")
    
    
    
    
  }
  
}

#write the collected data
{
  dir.create("E:/HMF3A_reprogramming/R_analysis/Ki67/combined_data/")
  setwd("E:/HMF3A_reprogramming/R_analysis/Ki67/combined_data/")
  
  write.csv(spheroid_combined, file="spheroid_combined.csv")
  write.csv(nucleus_combined, file="nucleus_combined.csv")
  write.csv(spheroid_nucleus_combined, file="spheroid_median_nucleus_combined.csv")
  
  write.csv(spheroid_radial$Ki67, file="spheroid_radial_dist_KI67.csv")
  
  write.csv(spheroid_2dint_data$Ki67, file="spheroid_2d_int_KI67.csv")
  write.csv(nuclei_2dint_data$KI67, file="nuclear_2d_int_KI67.csv")
  write.csv(spheroid_nuclei_2dint_data$KI67, file="spheroid_median_nuclear_2d_int_KI67.csv")
  
  
  write.csv(nuclei_3dint_data$KI67, file="nuclear_int_KI67.csv")
  write.csv(spheroid_nuclei_3dint_data$KI67, file="spheroid_median_nuclear_int_KI67.csv")
  write.csv(nuclei_3dint_data$dna, file="nuclear_int_dna.csv")
  write.csv(spheroid_nuclei_3dint_data$dna, file="spheroid_median_nuclear_int_dna.csv")
  
  write.csv(cell_3dint_data$KI67, file="cellular_int_KI67.csv")
  write.csv(spheroid_cell_3dint_data$KI67, file="spheroid_median_cellular_int_KI67.csv")
  write.csv(cell_3dint_data$dna, file="cellular_int_dna.csv")
  write.csv(spheroid_cell_3dint_data$dna, file="spheroid_median_cellular_int_dna.csv")
  
  write.csv(spheroid_axial$Ki67, file="spheroid_axial_dist_KI67.csv")
  write.csv(spheroid_axial$dna, file="spheroid_axial_dist_dna.csv")
  
  write.csv(spheroid_angle$Ki67, file="spheroid_angle_dist_KI67.csv")
  write.csv(spheroid_angle$dna, file="spheroid_angle_dist_dna.csv")
  
  write.csv(spheroid_int_data$Ki67, file="spheroid_int_KI67.csv")
  write.csv(spheroid_int_data$dna, file="spheroid_int_dna.csv")
  
  
}


#fraction of Ki67 positive cells with time
{
  T1<-subset(nuclei_3dint_data$KI6,nuclei_3dint_data$KI67$trial=="T1")
  T1$Max<-T1$Max/max(T1$Max)
  T1$status<-ifelse(T1$Max>0.5,1,0)
  T2<-subset(nuclei_3dint_data$KI6,nuclei_3dint_data$KI67$trial=="T2")
  T2$Max<-T2$Max/max(T2$Max)
  T2$status<-ifelse(T2$Max>0.5,1,0)
  T3<-subset(nuclei_3dint_data$KI6,nuclei_3dint_data$KI67$trial=="T3")
  T3$Max<-T3$Max/max(T3$Max)
  T3$status<-ifelse(T3$Max>0.4,1,0)
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
  
  pop_means<-c(mean(c(T2_D0,T1_D0,T3_D0)),mean(c(T2_D2,T1_D2,T3_D2)),mean(c(T2_D4,T1_D4,T3_D4)),
               mean(c(T1_D6,T2_D6,T3_D6)),mean(c(T1_D8,T2_D8,T3_D8)))
  pop_sd<-c(sd(c(T2_D0,T1_D0,T3_D0)),sd(c(T2_D2,T1_D2,T3_D2)),sd(c(T2_D4,T1_D4,T3_D4)),
            sd(c(T1_D6,T2_D6,T3_D6)),sd(c(T1_D8,T2_D8,T3_D8)))
  ymax<-max(pop_means+pop_sd)*1.1
  
  dird<-"E:/HMF3A_reprogramming/R_analysis/Ki67/black_bg_plots/"
  setwd(dird)
  png(filename="fraction_ki67_positive_cells.png", units="in",width=2, height=2, pointsize=7,  res=1200)
  par(font.axis = 2,font.lab=2,mar=c(4,5,1,1), font=2,bg='black',fg='white',col.lab='white',col.axis='white')
  foo<-barplot(pop_means,las=1,ylim=c(0,ymax),names=c("D0","D2","D4","D6","D8"),ylab="Fraction of KI67\npositive cells %")
  box()
  segments(foo,pop_means-pop_sd,foo,pop_means+pop_sd)
  dev.off()
  
  dird<-"E:/HMF3A_reprogramming/R_analysis/Ki67/white_bg_plots/"
  setwd(dird)
  png(filename="fraction_ki67_positive_cells.png", units="in",width=2, height=2, pointsize=7,  res=1200)
  par(font.axis = 2,font.lab=2,mar=c(4,5,1,1), font=2,bg='white',fg='black',col.lab='black',col.axis='black')
  foo<-barplot(pop_means,las=1,ylim=c(0,ymax),names=c("D0","D2","D4","D6","D8"),ylab="Fraction of KI67\npositive cells %")
  box()
  segments(foo,pop_means-pop_sd,foo,pop_means+pop_sd)
  dev.off()
  
  Ki67_3d<-rbind(T1,T2,T3)
  fields<-levels(as.factor(Ki67_3d$Spheroid_id))
  field_Ki67<-as.data.frame(matrix(nrow=length(fields),ncol=6))
  colnames(field_Ki67)<-c("Spheroid_id","Positive_fraction","sample","trial","num_cells","status")
  for(index in 1:length(fields)){
    field_Ki67[index,1]<-fields[index]
    temp<-subset(Ki67_3d,Ki67_3d$Spheroid_id==fields[index])
    field_Ki67[index,2]<-nrow(subset(temp,temp$status==1))/nrow(temp)
    field_Ki67[index,3]<-temp$sample[1]
    field_Ki67[index,4]<-temp$trial[1]
    field_Ki67[index,5]<-nrow(temp)
    field_Ki67[index,6]<-ifelse(nrow(subset(temp,temp$status==1))>1,1,0)
    
  }
  
  merged_data<-merge(cbind(Spheroid_id=spheroid_combined$Spheroid_id,Volume=spheroid_combined$Vol..unit.),
                     field_Ki67,by="Spheroid_id",stringsAsFactors=FALSE)
  
  D8<-subset(merged_data,merged_data$sample=="D8")
  plot(merged_data$Positive_fraction~merged_data$num_cells)
  plot(merged_data$Positive_fraction~as.numeric(as.character(merged_data$Volume)))
  
  plot(D8$status~D8$num_cells)
}


