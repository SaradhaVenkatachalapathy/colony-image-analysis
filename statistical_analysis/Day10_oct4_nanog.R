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
  path_to_ex<-c("E:/HMF3A_reprogramming/D10_actin_nanog_oct4/20200129_hmf3a_D10_Nanog_568_oct4_647_actin_488_DAPI/",
                "E:/HMF3A_reprogramming/D10_actin_nanog_oct4/20200203_hmf3a_D10_Nanog_568_oct4_647_actin_488_DAPI/",
                "E:/HMF3A_reprogramming/D10_actin_nanog_oct4/20200204_hmf3a_D10_Nanog_568_oct4_647_actin_488_DAPI/")
  
  exp<-c("20200129_hmf3a_D10_Nanog_568_oct4_647_actin_488_DAPI","20200203_hmf3a_D10_Nanog_568_oct4_647_actin_488_DAPI","20200204_hmf3a_D10_Nanog_568_oct4_647_actin_488_DAPI")
  
  data_types<-c("/3D geometrical data spheroid/",
                "/3D int_data spheroid/DNA/","/3D int_data spheroid/NANOG/","/3D int_data spheroid/ACTIN/","/3D int_data spheroid/OCT4/")
  
  data_types_2d<-c("/2D_measures_spheroid/2D_spheroid_DNA.csv","/2D_measures_spheroid/2D_spheroid_NANOG.csv",
                   "/2D_measures_spheroid/2D_spheroid_ACTIN.csv","/2D_measures_spheroid/2D_spheroid_OCT4.csv")
  
  nuc_data_types<-c("/3D geometrical data/","/3D ellipsoid/","/3D geometerical_simple/","/3D shape measure/",
                    "/3D int_data/DNA/","/3D int_data/NANOG/","/3D int_data/ACTIN/","/3D int_data/OCT4/")
  
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
                       "/2D_measures_nuclei/2D_nucleus_NANOG.csv", "/2D_measures_nuclei/2D_nucleus_OCT4.csv",
                       "/2D_measures_nuclei/2D_nucleusACTIN.csv")
  
  samples<-c("D10","D10")
  trial<-c("T1","T2","T3")
  
  dird_apo<-"E:/HMF3A_reprogramming/R_analysis/Day10_oct4_nanog/"
  dir.create(dird_apo)
  dird_OCT4<-paste(dird_apo,"spheroid_shells_OCT4", sep="")
  dird_ACTIN<-paste(dird_apo,"spheroid_shells_ACTIN", sep="")
  dird_NANOG<-paste(dird_apo,"spheroid_shells_NANOG", sep="")
  dird_dna<-paste(dird_apo,"spheroid_shells_dna", sep="")
  dir.create(dird_NANOG)
  dir.create(dird_dna)
  dir.create(dird_ACTIN)
  dir.create(dird_OCT4)
  
}

# read in all the required features for spheroids
{
  plot(0:1,0:1)
  #T1
  {
    geometrical_data<-combine_sample_sets(path_to_ex[1],samples,data_types[1],"tsv",exp[1],trial[1])
    T1DNA_int_data<-combine_sample_sets(path_to_ex[1],samples,data_types[2],"tsv",exp[1],trial[1])
    T1NANOG_int_data<-combine_sample_sets(path_to_ex[1],samples,data_types[3],"tsv",exp[1],trial[1])
    T1ACTIN_int_data<-combine_sample_sets(path_to_ex[1],samples,data_types[4],"tsv",exp[1],trial[1])
    T1OCT4_int_data<-combine_sample_sets(path_to_ex[1],samples,data_types[5],"tsv",exp[1],trial[1])
    geo_int_2D_data<-combine_sample_sets_2d(path_to_ex[1],samples,data_types_2d[1],exp[1],trial[1])
    
    T1_2d_data_NANOG<-combine_sample_sets_2d(path_to_ex[1],samples,data_types_2d[2],exp[1],trial[1])
    T1_2d_data_ACTIN<-combine_sample_sets_2d(path_to_ex[1],samples,data_types_2d[3],exp[1],trial[1])
    T1_2d_data_OCT4 <-combine_sample_sets_2d(path_to_ex[1],samples,data_types_2d[4],exp[1],trial[1])
    
    shell_NANOG_T1<-shell_analsysis(path_to_ex[1],samples,"/Shells/NANOG",exp[1],trial[1],dird_NANOG)
    shell_ACTIN_T1<-shell_analsysis(path_to_ex[1],samples,"/Shells/ACTIN",exp[1],trial[1],dird_ACTIN)
    shell_OCT4_T1<-shell_analsysis(path_to_ex[1],samples,"/Shells/OCT4",exp[1],trial[1],dird_OCT4)
    shell_dna_T1<-shell_analsysis(path_to_ex[1],samples,"/Shells/DNA",exp[1],trial[1],dird_dna)
    
    T1_dna_axial<-axial_analysis(path_to_ex[1],samples,exp[1],trial[1],"/DNA/")
    T1_NANOG_axial<-axial_analysis(path_to_ex[1],samples,exp[1],trial[1],"/NANOG/")
    T1_ACTIN_axial<-axial_analysis(path_to_ex[1],samples,exp[1],trial[1],"/ACTIN/")
    T1_OCT4_axial<-axial_analysis(path_to_ex[1],samples,exp[1],trial[1],"/OCT4/")
    
    T1_NANOG_angle<-angle_analsysis(path_to_ex[1],samples,"/NANOG",exp[1],trial[1])
    T1_ACTIN_angle<-angle_analsysis(path_to_ex[1],samples,"/ACTIN",exp[1],trial[1])
    T1_OCT4_angle<-angle_analsysis(path_to_ex[1],samples,"/OCT4",exp[1],trial[1])
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
    T2NANOG_int_data<-combine_sample_sets(path_to_ex[2],samples,data_types[3],"tsv",exp[2],trial[2])
    T2ACTIN_int_data<-combine_sample_sets(path_to_ex[2],samples,data_types[4],"tsv",exp[2],trial[2])
    T2OCT4_int_data<-combine_sample_sets(path_to_ex[2],samples,data_types[5],"tsv",exp[2],trial[2])
    geo_int_2D_data<-combine_sample_sets_2d(path_to_ex[2],samples,data_types_2d[1],exp[2],trial[2])
    
    T2_2d_data_NANOG<-combine_sample_sets_2d(path_to_ex[2],samples,data_types_2d[2],exp[2],trial[2])
    T2_2d_data_ACTIN<-combine_sample_sets_2d(path_to_ex[2],samples,data_types_2d[3],exp[2],trial[2])
    T2_2d_data_OCT4 <-combine_sample_sets_2d(path_to_ex[2],samples,data_types_2d[4],exp[2],trial[2])
    
    shell_NANOG_T2<-shell_analsysis(path_to_ex[2],samples,"/Shells/NANOG",exp[2],trial[2],dird_NANOG)
    shell_ACTIN_T2<-shell_analsysis(path_to_ex[2],samples,"/Shells/ACTIN",exp[2],trial[2],dird_ACTIN)
    shell_OCT4_T2<-shell_analsysis(path_to_ex[2],samples,"/Shells/OCT4",exp[2],trial[2],dird_OCT4)
    shell_dna_T2<-shell_analsysis(path_to_ex[2],samples,"/Shells/DNA",exp[2],trial[2],dird_dna)
    
    T2_dna_axial<-axial_analysis(path_to_ex[2],samples,exp[2],trial[2],"/DNA/")
    T2_NANOG_axial<-axial_analysis(path_to_ex[2],samples,exp[2],trial[2],"/NANOG/")
    T2_ACTIN_axial<-axial_analysis(path_to_ex[2],samples,exp[2],trial[2],"/ACTIN/")
    T2_OCT4_axial<-axial_analysis(path_to_ex[2],samples,exp[2],trial[2],"/OCT4/")
    
    T2_NANOG_angle<-angle_analsysis(path_to_ex[2],samples,"/NANOG",exp[2],trial[2])
    T2_ACTIN_angle<-angle_analsysis(path_to_ex[2],samples,"/ACTIN",exp[2],trial[2])
    T2_OCT4_angle<-angle_analsysis(path_to_ex[2],samples,"/OCT4",exp[2],trial[2])
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
    T3NANOG_int_data<-combine_sample_sets(path_to_ex[3],samples,data_types[3],"tsv",exp[3],trial[3])
    T3ACTIN_int_data<-combine_sample_sets(path_to_ex[3],samples,data_types[4],"tsv",exp[3],trial[3])
    T3OCT4_int_data<-combine_sample_sets(path_to_ex[3],samples,data_types[5],"tsv",exp[3],trial[3])
    geo_int_2D_data<-combine_sample_sets_2d(path_to_ex[3],samples,data_types_2d[1],exp[3],trial[3])
    
    T3_2d_data_NANOG<-combine_sample_sets_2d(path_to_ex[3],samples,data_types_2d[2],exp[3],trial[3])
    T3_2d_data_ACTIN<-combine_sample_sets_2d(path_to_ex[3],samples,data_types_2d[3],exp[3],trial[3])
    T3_2d_data_OCT4 <-combine_sample_sets_2d(path_to_ex[3],samples,data_types_2d[4],exp[3],trial[3])
    
    shell_NANOG_T3<-shell_analsysis(path_to_ex[3],samples,"/Shells/NANOG",exp[3],trial[3],dird_NANOG)
    shell_ACTIN_T3<-shell_analsysis(path_to_ex[3],samples,"/Shells/ACTIN",exp[3],trial[3],dird_ACTIN)
    shell_OCT4_T3<-shell_analsysis(path_to_ex[3],samples,"/Shells/OCT4",exp[3],trial[3],dird_OCT4)
    shell_dna_T3<-shell_analsysis(path_to_ex[3],samples,"/Shells/DNA",exp[3],trial[3],dird_dna)
    
    T3_dna_axial<-axial_analysis(path_to_ex[3],samples,exp[3],trial[3],"/DNA/")
    T3_NANOG_axial<-axial_analysis(path_to_ex[3],samples,exp[3],trial[3],"/NANOG/")
    T3_ACTIN_axial<-axial_analysis(path_to_ex[3],samples,exp[3],trial[3],"/ACTIN/")
    T3_OCT4_axial<-axial_analysis(path_to_ex[3],samples,exp[3],trial[3],"/OCT4/")
    
    T3_NANOG_angle<-angle_analsysis(path_to_ex[3],samples,"/NANOG",exp[3],trial[3])
    T3_ACTIN_angle<-angle_analsysis(path_to_ex[3],samples,"/ACTIN",exp[3],trial[3])
    T3_OCT4_angle<-angle_analsysis(path_to_ex[3],samples,"/OCT4",exp[3],trial[3])
    T3_dna_angle<-angle_analsysis(path_to_ex[3],samples,"/DNA",exp[3],trial[3])
    
    geo_int_2D_data$Label<-substring(geo_int_2D_data$Label,5, nchar(geo_int_2D_data$Label)-4)
    geometrical_data$Label<-substring(geometrical_data$Label,0, nchar(geometrical_data$Label)-1)
    
    T3combined<-merge(geometrical_data,geo_int_2D_data, by="Label")
    
    rm(geometrical_data,geo_int_2D_data)
    
  }
  
  spheroid_2dint_data_NANOG<-rbind(T1_2d_data_NANOG,T2_2d_data_NANOG,T3_2d_data_NANOG)
  spheroid_2dint_data_NANOG$Label<-substring(spheroid_2dint_data_NANOG$Label,5,(nchar(spheroid_2dint_data_NANOG$Label)-4))
  spheroid_2dint_data_NANOG$Spheroid_id<-paste(spheroid_2dint_data_NANOG$trial,spheroid_2dint_data_NANOG$sample,
                                               substring(spheroid_2dint_data_NANOG$Label,1,2),sep="_")
  
  spheroid_2dint_data_ACTIN<-rbind(T1_2d_data_ACTIN,T2_2d_data_ACTIN,T3_2d_data_ACTIN)
  spheroid_2dint_data_ACTIN$Label<-substring(spheroid_2dint_data_ACTIN$Label,5,(nchar(spheroid_2dint_data_ACTIN$Label)-4))
  spheroid_2dint_data_ACTIN$Spheroid_id<-paste(spheroid_2dint_data_ACTIN$trial,spheroid_2dint_data_ACTIN$sample,
                                               substring(spheroid_2dint_data_ACTIN$Label,1,2),sep="_")
  
  spheroid_2dint_data_OCT4<-rbind(T1_2d_data_OCT4,T2_2d_data_OCT4,T3_2d_data_OCT4)
  spheroid_2dint_data_OCT4$Label<-substring(spheroid_2dint_data_OCT4$Label,5,(nchar(spheroid_2dint_data_OCT4$Label)-4))
  spheroid_2dint_data_OCT4$Spheroid_id<-paste(spheroid_2dint_data_OCT4$trial,spheroid_2dint_data_OCT4$sample,
                                              substring(spheroid_2dint_data_OCT4$Label,1,2),sep="_")
  ###
  spheroid_int_data_NANOG<-rbind(T1NANOG_int_data,T2NANOG_int_data,T3NANOG_int_data)
  spheroid_int_data_NANOG$Label<-substring(spheroid_int_data_NANOG$Label,1,(nchar(spheroid_int_data_NANOG$Label)-1))
  spheroid_int_data_NANOG$Spheroid_id<-paste(spheroid_int_data_NANOG$trial,spheroid_int_data_NANOG$sample,
                                             substring(spheroid_int_data_NANOG$Label,1,2),sep="_")
  
  spheroid_int_data_ACTIN<-rbind(T1ACTIN_int_data,T2ACTIN_int_data,T3ACTIN_int_data)
  spheroid_int_data_ACTIN$Label<-substring(spheroid_int_data_ACTIN$Label,1,(nchar(spheroid_int_data_ACTIN$Label)-1))
  spheroid_int_data_ACTIN$Spheroid_id<-paste(spheroid_int_data_ACTIN$trial,spheroid_int_data_ACTIN$sample,
                                             substring(spheroid_int_data_ACTIN$Label,1,2),sep="_")
  
  spheroid_int_data_OCT4<-rbind(T1OCT4_int_data,T2OCT4_int_data,T3OCT4_int_data)
  spheroid_int_data_OCT4$Label<-substring(spheroid_int_data_OCT4$Label,1,(nchar(spheroid_int_data_OCT4$Label)-1))
  spheroid_int_data_OCT4$Spheroid_id<-paste(spheroid_int_data_OCT4$trial,spheroid_int_data_OCT4$sample,
                                            substring(spheroid_int_data_OCT4$Label,1,2),sep="_")
  
  spheroid_int_data_dna<-rbind(T1DNA_int_data,T2DNA_int_data,T3DNA_int_data)
  spheroid_int_data_dna$Label<-substring(spheroid_int_data_dna$Label,1,(nchar(spheroid_int_data_dna$Label)-1))
  spheroid_int_data_dna$Spheroid_id<-paste(spheroid_int_data_dna$trial,spheroid_int_data_dna$sample,
                                           substring(spheroid_int_data_dna$Label,1,2),sep="_")
  ###
  spheroid_radial_NANOG<-rbind(shell_NANOG_T1,shell_NANOG_T2,shell_NANOG_T3)
  spheroid_radial_NANOG$Spheroid_id<-paste(spheroid_radial_NANOG$trial,spheroid_radial_NANOG$sample,
                                           substring(spheroid_radial_NANOG$Label,1,2),sep="_")
  
  spheroid_radial_ACTIN<-rbind(shell_ACTIN_T1,shell_ACTIN_T2,shell_ACTIN_T3)
  spheroid_radial_ACTIN$Spheroid_id<-paste(spheroid_radial_ACTIN$trial,spheroid_radial_ACTIN$sample,
                                           substring(spheroid_radial_ACTIN$Label,1,2),sep="_")
  
  spheroid_radial_OCT4<-rbind(shell_OCT4_T1,shell_OCT4_T2,shell_OCT4_T3)
  spheroid_radial_OCT4$Spheroid_id<-paste(spheroid_radial_OCT4$trial,spheroid_radial_OCT4$sample,
                                          substring(spheroid_radial_OCT4$Label,1,2),sep="_")
  
  spheroid_radial_dna<-rbind(shell_dna_T1,shell_dna_T2,shell_dna_T3)
  spheroid_radial_dna$Spheroid_id<-paste(spheroid_radial_dna$trial,spheroid_radial_dna$sample,
                                         substring(spheroid_radial_dna$Label,1,2),sep="_")
  ####
  spheroid_axial_dna<-rbind(T1_dna_axial,T2_dna_axial,T3_dna_axial)
  spheroid_axial_dna$Spheroid_id<-paste(spheroid_axial_dna$trial,spheroid_axial_dna$sample,
                                        substring(spheroid_axial_dna$Label,1,2),sep="_")
  
  spheroid_axial_NANOG<-rbind(T1_NANOG_axial,T2_NANOG_axial,T3_NANOG_axial)
  spheroid_axial_NANOG$Spheroid_id<-paste(spheroid_axial_NANOG$trial,spheroid_axial_NANOG$sample,
                                          substring(spheroid_axial_NANOG$Label,1,2),sep="_")
  
  spheroid_axial_ACTIN<-rbind(T1_ACTIN_axial,T2_ACTIN_axial,T3_ACTIN_axial)
  spheroid_axial_ACTIN$Spheroid_id<-paste(spheroid_axial_ACTIN$trial,spheroid_axial_ACTIN$sample,
                                          substring(spheroid_axial_ACTIN$Label,1,2),sep="_")
  
  spheroid_axial_OCT4<-rbind(T1_OCT4_axial,T2_OCT4_axial,T3_OCT4_axial)
  spheroid_axial_OCT4$Spheroid_id<-paste(spheroid_axial_OCT4$trial,spheroid_axial_OCT4$sample,
                                         substring(spheroid_axial_OCT4$Label,1,2),sep="_")
  
  
  #### 
  spheroid_angle_NANOG<-rbind(T1_NANOG_angle,T2_NANOG_angle,T3_NANOG_angle)
  spheroid_angle_NANOG$Spheroid_id<-paste(spheroid_angle_NANOG$trial,spheroid_angle_NANOG$sample,
                                          substring(spheroid_angle_NANOG$Label,1,2),sep="_")
  
  spheroid_angle_ACTIN<-rbind(T1_ACTIN_angle,T2_ACTIN_angle,T3_ACTIN_angle)
  spheroid_angle_ACTIN$Spheroid_id<-paste(spheroid_angle_ACTIN$trial,spheroid_angle_ACTIN$sample,
                                          substring(spheroid_angle_ACTIN$Label,1,2),sep="_")
  
  spheroid_angle_OCT4<-rbind(T1_OCT4_angle,T2_OCT4_angle,T3_OCT4_angle)
  spheroid_angle_OCT4$Spheroid_id<-paste(spheroid_angle_OCT4$trial,spheroid_angle_OCT4$sample,
                                         substring(spheroid_angle_OCT4$Label,1,2),sep="_")
  
  spheroid_angle_dna<-rbind(T1_dna_angle,T2_dna_angle,T3_dna_angle)
  spheroid_angle_dna$Spheroid_id<-paste(spheroid_angle_dna$trial,spheroid_angle_dna$sample,
                                        substring(spheroid_angle_dna$Label,1,2),sep="_")
  
  ####
  spheroid_combined<-rbind(T1combined,T2combined,T3combined)
  colnames(spheroid_combined)[86:87]<-c("sample","trial")
  spheroid_combined$Spheroid_id<-paste(spheroid_combined$trial,spheroid_combined$sample,
                                       substring(spheroid_combined$Label,1,2),sep="_")
  
  spheroid_int_data<-list( NANOG=spheroid_int_data_NANOG,
                           dna=spheroid_int_data_dna,
                           OCT4=spheroid_int_data_OCT4,
                           ACTIN=spheroid_int_data_ACTIN)
  
  spheroid_radial<-list(NANOG=spheroid_radial_NANOG,
                        dna=spheroid_radial_dna,
                        OCT4=spheroid_radial_OCT4,
                        ACTIN=spheroid_radial_ACTIN)
  
  spheroid_axial<-list( NANOG=spheroid_axial_NANOG,
                        dna=spheroid_axial_dna,
                        OCT4=spheroid_axial_OCT4,
                        ACTIN=spheroid_axial_ACTIN)
  
  spheroid_angle<-list(NANOG=spheroid_angle_NANOG,
                       dna=spheroid_angle_dna,
                       OCT4=spheroid_angle_OCT4,
                       ACTIN=spheroid_angle_ACTIN)
  
  spheroid_2dint_data<-list( NANOG=spheroid_2dint_data_NANOG,
                             OCT4=spheroid_2dint_data_OCT4,
                             ACTIN=spheroid_2dint_data_ACTIN)
  
  
  rm(T1_2d_data_NANOG,T2_2d_data_NANOG,T3_2d_data_NANOG,T1_2d_data_ACTIN,T2_2d_data_ACTIN,T3_2d_data_ACTIN,
     T1_2d_data_OCT4,T2_2d_data_OCT4,T3_2d_data_OCT4,T1NANOG_int_data,T2NANOG_int_data,T3NANOG_int_data,
     T1ACTIN_int_data,T2ACTIN_int_data,T3ACTIN_int_data,T1OCT4_int_data,T2OCT4_int_data,T3OCT4_int_data,
     T1DNA_int_data,T2DNA_int_data,T3DNA_int_data,
     shell_NANOG_T1,shell_NANOG_T2,shell_NANOG_T3,shell_ACTIN_T1,shell_ACTIN_T2,shell_ACTIN_T3,
     shell_OCT4_T1,shell_OCT4_T2,shell_OCT4_T3,shell_dna_T1,shell_dna_T2,shell_dna_T3,
     T1_dna_axial,T2_dna_axial,T3_dna_axial,T1_NANOG_axial,T2_NANOG_axial,T3_NANOG_axial,
     T1_ACTIN_axial,T2_ACTIN_axial,T3_ACTIN_axial,T1_OCT4_axial,T2_OCT4_axial,T3_OCT4_axial,
     T1_NANOG_angle,T2_NANOG_angle,T3_NANOG_angle,T1_ACTIN_angle,T2_ACTIN_angle,T3_ACTIN_angle,
     T1_OCT4_angle,T2_OCT4_angle,T3_OCT4_angle,T1_dna_angle,T2_dna_angle,T3_dna_angle,
     T1combined,T2combined,T3combined,
     spheroid_int_data_NANOG,spheroid_int_data_dna,spheroid_int_data_OCT4,spheroid_int_data_ACTIN,
     spheroid_radial_NANOG,spheroid_radial_dna,spheroid_radial_OCT4,spheroid_radial_ACTIN,
     spheroid_axial_NANOG,spheroid_axial_dna,spheroid_axial_OCT4,spheroid_axial_ACTIN,
     spheroid_angle_NANOG,spheroid_angle_dna,spheroid_angle_OCT4,spheroid_angle_ACTIN,
     spheroid_2dint_data_NANOG,spheroid_2dint_data_OCT4,spheroid_2dint_data_ACTIN)
  
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
    
    T1_NANOG_int_data<-combine_sample_sets(path_to_ex[1],samples,nuc_data_types[6],"tsv",exp[1],trial[1])
    T1_ACTIN_int_data<-combine_sample_sets(path_to_ex[1],samples,nuc_data_types[7],"tsv",exp[1],trial[1])
    T1_OCT4_int_data<-combine_sample_sets(path_to_ex[1],samples,nuc_data_types[8],"tsv",exp[1],trial[1])
    
    T1_2d_NANOG_int_data<-combine_sample_sets_2d(path_to_ex[1],samples,nuc_data_types_2d[28],exp[1],trial[1])
    T1_2d_OCT4_int_data<-combine_sample_sets_2d(path_to_ex[1],samples,nuc_data_types_2d[29],exp[1],trial[1])
    T1_2d_ACTIN_int_data<-combine_sample_sets_2d(path_to_ex[1],samples,nuc_data_types_2d[30],exp[1],trial[1])
    
    
    colnames(T1_compact_2D_data)[2]<-"Label"
    T1_compact_2D_data$Label<-substring(T1_compact_2D_data$Label,0, nchar(T1_compact_2D_data$Label)-4)
    T1_geo_int_2D_data$Label<-substring(T1_geo_int_2D_data$Label,5, nchar(T1_geo_int_2D_data$Label)-4)
    T1_geometrical_data$Label<-substring(T1_geometrical_data$Label,0, nchar(T1_geometrical_data$Label)-1)
    T1_DNA_int_data$Label<-substring(T1_DNA_int_data$Label,0, nchar(T1_DNA_int_data$Label)-1)
    T1_2d_NANOG_int_data$Label<-substring(T1_2d_NANOG_int_data$Label,5, nchar(T1_2d_NANOG_int_data$Label)-4)
    T1_2d_OCT4_int_data$Label<-substring(T1_2d_OCT4_int_data$Label,5, nchar(T1_2d_OCT4_int_data$Label)-4)
    T1_2d_ACTIN_int_data$Label<-substring(T1_2d_ACTIN_int_data$Label,5, nchar(T1_2d_ACTIN_int_data$Label)-4)
    T1_NANOG_int_data$Label<-substring(T1_NANOG_int_data$Label,0, nchar(T1_NANOG_int_data$Label)-1)
    T1_ACTIN_int_data$Label<-substring(T1_ACTIN_int_data$Label,0, nchar(T1_ACTIN_int_data$Label)-1)
    T1_OCT4_int_data$Label<-substring(T1_OCT4_int_data$Label,0, nchar(T1_OCT4_int_data$Label)-1)
    
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
    
    T2_NANOG_int_data<-combine_sample_sets(path_to_ex[2],samples,nuc_data_types[6],"tsv",exp[2],trial[2])
    T2_ACTIN_int_data<-combine_sample_sets(path_to_ex[2],samples,nuc_data_types[7],"tsv",exp[2],trial[2])
    T2_OCT4_int_data<-combine_sample_sets(path_to_ex[2],samples,nuc_data_types[8],"tsv",exp[2],trial[2])
    
    T2_2d_NANOG_int_data<-combine_sample_sets_2d(path_to_ex[2],samples,nuc_data_types_2d[28],exp[2],trial[2])
    T2_2d_OCT4_int_data<-combine_sample_sets_2d(path_to_ex[2],samples,nuc_data_types_2d[29],exp[2],trial[2])
    T2_2d_ACTIN_int_data<-combine_sample_sets_2d(path_to_ex[2],samples,nuc_data_types_2d[30],exp[2],trial[2])
    
    
    colnames(T2_compact_2D_data)[2]<-"Label"
    T2_compact_2D_data$Label<-substring(T2_compact_2D_data$Label,0, nchar(T2_compact_2D_data$Label)-4)
    T2_geo_int_2D_data$Label<-substring(T2_geo_int_2D_data$Label,5, nchar(T2_geo_int_2D_data$Label)-4)
    T2_geometrical_data$Label<-substring(T2_geometrical_data$Label,0, nchar(T2_geometrical_data$Label)-1)
    T2_DNA_int_data$Label<-substring(T2_DNA_int_data$Label,0, nchar(T2_DNA_int_data$Label)-1)
    T2_2d_NANOG_int_data$Label<-substring(T2_2d_NANOG_int_data$Label,5, nchar(T2_2d_NANOG_int_data$Label)-4)
    T2_2d_OCT4_int_data$Label<-substring(T2_2d_OCT4_int_data$Label,5, nchar(T2_2d_OCT4_int_data$Label)-4)
    T2_2d_ACTIN_int_data$Label<-substring(T2_2d_ACTIN_int_data$Label,5, nchar(T2_2d_ACTIN_int_data$Label)-4)
    T2_NANOG_int_data$Label<-substring(T2_NANOG_int_data$Label,0, nchar(T2_NANOG_int_data$Label)-1)
    T2_ACTIN_int_data$Label<-substring(T2_ACTIN_int_data$Label,0, nchar(T2_ACTIN_int_data$Label)-1)
    T2_OCT4_int_data$Label<-substring(T2_OCT4_int_data$Label,0, nchar(T2_OCT4_int_data$Label)-1)
    
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
    
    T3_NANOG_int_data<-combine_sample_sets(path_to_ex[3],samples,nuc_data_types[6],"tsv",exp[3],trial[3])
    T3_ACTIN_int_data<-combine_sample_sets(path_to_ex[3],samples,nuc_data_types[7],"tsv",exp[3],trial[3])
    T3_OCT4_int_data<-combine_sample_sets(path_to_ex[3],samples,nuc_data_types[8],"tsv",exp[3],trial[3])
    
    T3_2d_NANOG_int_data<-combine_sample_sets_2d(path_to_ex[3],samples,nuc_data_types_2d[28],exp[3],trial[3])
    T3_2d_OCT4_int_data<-combine_sample_sets_2d(path_to_ex[3],samples,nuc_data_types_2d[29],exp[3],trial[3])
    T3_2d_ACTIN_int_data<-combine_sample_sets_2d(path_to_ex[3],samples,nuc_data_types_2d[30],exp[3],trial[3])
    
    
    colnames(T3_compact_2D_data)[2]<-"Label"
    T3_compact_2D_data$Label<-substring(T3_compact_2D_data$Label,0, nchar(T3_compact_2D_data$Label)-4)
    T3_geo_int_2D_data$Label<-substring(T3_geo_int_2D_data$Label,5, nchar(T3_geo_int_2D_data$Label)-4)
    T3_geometrical_data$Label<-substring(T3_geometrical_data$Label,0, nchar(T3_geometrical_data$Label)-1)
    T3_DNA_int_data$Label<-substring(T3_DNA_int_data$Label,0, nchar(T3_DNA_int_data$Label)-1)
    T3_2d_NANOG_int_data$Label<-substring(T3_2d_NANOG_int_data$Label,5, nchar(T3_2d_NANOG_int_data$Label)-4)
    T3_2d_OCT4_int_data$Label<-substring(T3_2d_OCT4_int_data$Label,5, nchar(T3_2d_OCT4_int_data$Label)-4)
    T3_2d_ACTIN_int_data$Label<-substring(T3_2d_ACTIN_int_data$Label,5, nchar(T3_2d_ACTIN_int_data$Label)-4)
    T3_NANOG_int_data$Label<-substring(T3_NANOG_int_data$Label,0, nchar(T3_NANOG_int_data$Label)-1)
    T3_ACTIN_int_data$Label<-substring(T3_ACTIN_int_data$Label,0, nchar(T3_ACTIN_int_data$Label)-1)
    T3_OCT4_int_data$Label<-substring(T3_OCT4_int_data$Label,0, nchar(T3_OCT4_int_data$Label)-1)
    
    T3_combined<-merge(T3_geometrical_data,T3_geo_int_2D_data[,2:36], by="Label")
    T3_combined<-merge(T3_combined,T3_compact_2D_data[,2:8], by="Label")
    T3_combined<-merge(T3_combined,T3_DNA_int_data[,c(3:5,12:18)], by="Label")
    T3_combined<-merge(T3_combined,T3_edf_2D_data[,1:11], by="Label")
    T3_combined<-cbind(T3_combined,T3_gclm[,-c(1,128:131)])
    
    rm(T3_geometrical_data,T3_geo_int_2D_data,T3_compact_2D_data,T3_edf_2D_data,T3_gclm)
  }
  
  nuclei_2dint_data_NANOG<-rbind(T1_2d_NANOG_int_data,T2_2d_NANOG_int_data,T3_2d_NANOG_int_data)
  nuclei_2dint_data_NANOG$Spheroid_id<-paste(nuclei_2dint_data_NANOG$trial,nuclei_2dint_data_NANOG$sample,
                                             substring(nuclei_2dint_data_NANOG$Label,1,2),sep="_")
  spheroid_nuclei_2dint_data_NANOG<-nuc_median_spheroid_levels(nuclei_2dint_data_NANOG,parameters_2d_int_data)
  
  nuclei_2dint_data_ACTIN<-rbind(T1_2d_ACTIN_int_data,T2_2d_ACTIN_int_data,T3_2d_ACTIN_int_data)
  nuclei_2dint_data_ACTIN$Spheroid_id<-paste(nuclei_2dint_data_ACTIN$trial,nuclei_2dint_data_ACTIN$sample,
                                             substring(nuclei_2dint_data_ACTIN$Label,1,2),sep="_")
  spheroid_nuclei_2dint_data_ACTIN<-nuc_median_spheroid_levels(nuclei_2dint_data_ACTIN,parameters_2d_int_data)
  
  nuclei_2dint_data_OCT4<-rbind(T1_2d_OCT4_int_data,T2_2d_OCT4_int_data,T3_2d_OCT4_int_data)
  nuclei_2dint_data_OCT4$Spheroid_id<-paste(nuclei_2dint_data_OCT4$trial,nuclei_2dint_data_OCT4$sample,
                                            substring(nuclei_2dint_data_OCT4$Label,1,2),sep="_")
  spheroid_nuclei_2dint_data_OCT4<-nuc_median_spheroid_levels(nuclei_2dint_data_OCT4,parameters_2d_int_data)
  ###
  nuclei_int_data_NANOG<-rbind(T1_NANOG_int_data,T2_NANOG_int_data,T3_NANOG_int_data)
  nuclei_int_data_NANOG$Spheroid_id<-paste(nuclei_int_data_NANOG$trial,nuclei_int_data_NANOG$sample,
                                           substring(nuclei_int_data_NANOG$Label,1,2),sep="_")
  spheroid_nuclei_int_data_NANOG<-nuc_median_spheroid_levels(nuclei_int_data_NANOG,parameters_int_data)
  
  nuclei_int_data_OCT4<-rbind(T1_OCT4_int_data,T2_OCT4_int_data,T3_OCT4_int_data)
  nuclei_int_data_OCT4$Spheroid_id<-paste(nuclei_int_data_OCT4$trial,nuclei_int_data_OCT4$sample,
                                          substring(nuclei_int_data_OCT4$Label,1,2),sep="_")
  spheroid_nuclei_int_data_OCT4<-nuc_median_spheroid_levels(nuclei_int_data_OCT4,parameters_int_data)
  
  nuclei_int_data_ACTIN<-rbind(T1_ACTIN_int_data,T2_ACTIN_int_data,T3_ACTIN_int_data)
  nuclei_int_data_ACTIN$Spheroid_id<-paste(nuclei_int_data_ACTIN$trial,nuclei_int_data_ACTIN$sample,
                                           substring(nuclei_int_data_ACTIN$Label,1,2),sep="_")
  spheroid_nuclei_int_data_ACTIN<-nuc_median_spheroid_levels(nuclei_int_data_ACTIN,parameters_int_data)
  
  nuclei_int_data_dna<-rbind(T1_DNA_int_data,T2_DNA_int_data,T3_DNA_int_data)
  nuclei_int_data_dna$Spheroid_id<-paste(nuclei_int_data_dna$trial,nuclei_int_data_dna$sample,
                                         substring(nuclei_int_data_dna$Label,1,2),sep="_")
  spheroid_nuclei_int_data_dna<-nuc_median_spheroid_levels(nuclei_int_data_dna,parameters_int_data)
  ####
  
  nucleus_combined<-rbind(T1_combined,T2_combined,T3_combined)
  nucleus_combined$Spheroid_id<-paste(nucleus_combined$trial,nucleus_combined$sample,
                                      substring(nucleus_combined$Label,1,2),sep="_")
  spheroid_nucleus_combined<-nuc_median_spheroid_levels(nucleus_combined,parameters_geo_texture_data)
  ####
  
  nuclei_2dint_data<-list(NANOG=nuclei_2dint_data_NANOG,ACTIN=nuclei_2dint_data_ACTIN,OCT4=nuclei_2dint_data_OCT4)
  nuclei_3dint_data<-list(NANOG=nuclei_int_data_NANOG,dna=nuclei_int_data_dna,ACTIN=nuclei_int_data_ACTIN,OCT4=nuclei_int_data_OCT4)
  
  spheroid_nuclei_2dint_data<-list(NANOG=spheroid_nuclei_2dint_data_NANOG, ACTIN=spheroid_nuclei_2dint_data_ACTIN,
                                   OCT4=spheroid_nuclei_2dint_data_OCT4)
  spheroid_nuclei_3dint_data<-list(NANOG=spheroid_nuclei_int_data_NANOG,dna=spheroid_nuclei_int_data_dna,
                                   ACTIN=spheroid_nuclei_int_data_ACTIN,OCT4=spheroid_nuclei_int_data_OCT4)
  ####
  
  
  rm(T1_2d_NANOG_int_data,T2_2d_NANOG_int_data,T3_2d_NANOG_int_data,T1_2d_ACTIN_int_data,T2_2d_ACTIN_int_data,T3_2d_ACTIN_int_data,
     T1_2d_OCT4_int_data,T2_2d_OCT4_int_data,T3_2d_OCT4_int_data,
     T1_NANOG_int_data,T2_NANOG_int_data,T3_NANOG_int_data,T1_OCT4_int_data,T2_OCT4_int_data,T3_OCT4_int_data,
     T1_ACTIN_int_data,T2_ACTIN_int_data,T3_ACTIN_int_data,T1_DNA_int_data,T2_DNA_int_data,T3_DNA_int_data,
     T1_combined,T2_combined,T3_combined,nuclei_2dint_data_NANOG,nuclei_2dint_data_ACTIN,nuclei_2dint_data_OCT4,
     nuclei_int_data_NANOG,nuclei_int_data_dna,nuclei_int_data_ACTIN,nuclei_int_data_OCT4,
     spheroid_nuclei_2dint_data_NANOG,spheroid_nuclei_2dint_data_ACTIN, spheroid_nuclei_2dint_data_OCT4,
     spheroid_nuclei_int_data_NANOG,spheroid_nuclei_int_data_dna,spheroid_nuclei_int_data_ACTIN,spheroid_nuclei_int_data_OCT4)
  
}

#plotting
{
  dird<-"E:/HMF3A_reprogramming/R_analysis/NANOG/black_bg_plots/"
  dir.create(dird)
  setwd(dird)
  cols<-colorRampPalette(c("red","yellow"))( 5) 
  #Spheoid Proj_Int
  {
    dire<-paste(dird,"spheroid_proj_int/",sep="")
    dir.create(dire)
    setwd(dire)
    
    plot_black_bg_histogram(subset(spheroid_2dint_data$NANOG,spheroid_2dint_data$NANOG$trial=="T1"), cols,parameters_2d_int_data,paste("NANOG ",labelling_2d_int_data,sep=""),"NANOG_T1_")
    plot_black_bg_histogram(subset(spheroid_2dint_data$NANOG,spheroid_2dint_data$NANOG$trial=="T2"), cols,parameters_2d_int_data,paste("NANOG ",labelling_2d_int_data,sep=""),"NANOG_T2_")
    plot_black_bg_histogram(subset(spheroid_2dint_data$NANOG,spheroid_2dint_data$NANOG$trial=="T3"), cols,parameters_2d_int_data,paste("NANOG ",labelling_2d_int_data,sep=""),"NANOG_T3_")
    
    plot_black_bg_histogram(subset(spheroid_2dint_data$OCT4,spheroid_2dint_data$OCT4$trial=="T1"), cols,parameters_2d_int_data,paste("OCT4 ",labelling_2d_int_data,sep=""),"OCT4_T1_")
    plot_black_bg_histogram(subset(spheroid_2dint_data$OCT4,spheroid_2dint_data$OCT4$trial=="T2"), cols,parameters_2d_int_data,paste("OCT4 ",labelling_2d_int_data,sep=""),"OCT4_T2_")
    plot_black_bg_histogram(subset(spheroid_2dint_data$OCT4,spheroid_2dint_data$OCT4$trial=="T3"), cols,parameters_2d_int_data,paste("OCT4 ",labelling_2d_int_data,sep=""),"OCT4_T3_")
    
    plot_black_bg_histogram(subset(spheroid_2dint_data$ACTIN,spheroid_2dint_data$ACTIN$trial=="T1"), cols,parameters_2d_int_data,paste("ACTIN ",labelling_2d_int_data,sep=""),"ACTIN_T1_")
    plot_black_bg_histogram(subset(spheroid_2dint_data$ACTIN,spheroid_2dint_data$ACTIN$trial=="T2"), cols,parameters_2d_int_data,paste("ACTIN ",labelling_2d_int_data,sep=""),"ACTIN_T2_")
    plot_black_bg_histogram(subset(spheroid_2dint_data$ACTIN,spheroid_2dint_data$ACTIN$trial=="T3"), cols,parameters_2d_int_data,paste("ACTIN ",labelling_2d_int_data,sep=""),"ACTIN_T3_")
  }
  #Spheoid 3D_Int
  {
    dire<-paste(dird,"spheroid_3d_int/",sep="")
    dir.create(dire)
    setwd(dire)
    
    plot_black_bg_histogram(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T1_")
    plot_black_bg_histogram(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T2_")
    plot_black_bg_histogram(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T3_")
    
    
    plot_black_bg_histogram(subset(spheroid_int_data$NANOG,spheroid_int_data$NANOG$trial=="T1"), cols,parameters_int_data,paste("NANOG ",labelling_int_data,sep=""),"NANOG_T1_")
    plot_black_bg_histogram(subset(spheroid_int_data$NANOG,spheroid_int_data$NANOG$trial=="T2"), cols,parameters_int_data,paste("NANOG ",labelling_int_data,sep=""),"NANOG_T2_")
    plot_black_bg_histogram(subset(spheroid_int_data$NANOG,spheroid_int_data$NANOG$trial=="T3"), cols,parameters_int_data,paste("NANOG ",labelling_int_data,sep=""),"NANOG_T3_")
    
    plot_black_bg_histogram(subset(spheroid_int_data$OCT4,spheroid_int_data$OCT4$trial=="T1"), cols,parameters_int_data,paste("OCT4 ",labelling_int_data,sep=""),"OCT4_T1_")
    plot_black_bg_histogram(subset(spheroid_int_data$OCT4,spheroid_int_data$OCT4$trial=="T2"), cols,parameters_int_data,paste("OCT4 ",labelling_int_data,sep=""),"OCT4_T2_")
    plot_black_bg_histogram(subset(spheroid_int_data$OCT4,spheroid_int_data$OCT4$trial=="T3"), cols,parameters_int_data,paste("OCT4 ",labelling_int_data,sep=""),"OCT4_T3_")
    
    
    plot_black_bg_histogram(subset(spheroid_int_data$ACTIN,spheroid_int_data$ACTIN$trial=="T1"), cols,parameters_int_data,paste("ACTIN ",labelling_int_data,sep=""),"ACTIN_T1_")
    plot_black_bg_histogram(subset(spheroid_int_data$ACTIN,spheroid_int_data$ACTIN$trial=="T2"), cols,parameters_int_data,paste("ACTIN ",labelling_int_data,sep=""),"ACTIN_T2_")
    plot_black_bg_histogram(subset(spheroid_int_data$ACTIN,spheroid_int_data$ACTIN$trial=="T3"), cols,parameters_int_data,paste("ACTIN ",labelling_int_data,sep=""),"ACTIN_T3_")
    
  }
  #Spheoid angles
  {
    dire<-paste(dird,"spheroid_angles/",sep="")
    dir.create(dire)
    setwd(dire)
    
    plot_black_bg_histogram(subset(spheroid_angle$dna,spheroid_angle$dna$trial=="T1"), cols,parameters_angle_data,paste("DNA ",labelling_angle_data,sep=""),"DNA_T1_")
    plot_black_bg_histogram(subset(spheroid_angle$dna,spheroid_angle$dna$trial=="T2"), cols,parameters_angle_data,paste("DNA ",labelling_angle_data,sep=""),"DNA_T2_")
    plot_black_bg_histogram(subset(spheroid_angle$dna,spheroid_angle$dna$trial=="T3"), cols,parameters_angle_data,paste("DNA ",labelling_angle_data,sep=""),"DNA_T3_")
    
    plot_black_bg_histogram(subset(spheroid_angle$NANOG,spheroid_angle$NANOG$trial=="T1"), cols,parameters_angle_data,paste("NANOG ",labelling_angle_data,sep=""),"NANOG_T1_")
    plot_black_bg_histogram(subset(spheroid_angle$NANOG,spheroid_angle$NANOG$trial=="T2"), cols,parameters_angle_data,paste("NANOG ",labelling_angle_data,sep=""),"NANOG_T2_")
    plot_black_bg_histogram(subset(spheroid_angle$NANOG,spheroid_angle$NANOG$trial=="T3"), cols,parameters_angle_data,paste("NANOG ",labelling_angle_data,sep=""),"NANOG_T3_")
    
    plot_black_bg_histogram(subset(spheroid_angle$OCT4,spheroid_angle$OCT4$trial=="T1"), cols,parameters_angle_data,paste("OCT4 ",labelling_angle_data,sep=""),"OCT4_T1_")
    plot_black_bg_histogram(subset(spheroid_angle$OCT4,spheroid_angle$OCT4$trial=="T2"), cols,parameters_angle_data,paste("OCT4 ",labelling_angle_data,sep=""),"OCT4_T2_")
    plot_black_bg_histogram(subset(spheroid_angle$OCT4,spheroid_angle$OCT4$trial=="T3"), cols,parameters_angle_data,paste("OCT4 ",labelling_angle_data,sep=""),"OCT4_T3_")
    
    plot_black_bg_histogram(subset(spheroid_angle$ACTIN,spheroid_angle$ACTIN$trial=="T1"), cols,parameters_angle_data,paste("ACTIN ",labelling_angle_data,sep=""),"ACTIN_T1_")
    plot_black_bg_histogram(subset(spheroid_angle$ACTIN,spheroid_angle$ACTIN$trial=="T2"), cols,parameters_angle_data,paste("ACTIN ",labelling_angle_data,sep=""),"ACTIN_T2_")
    plot_black_bg_histogram(subset(spheroid_angle$ACTIN,spheroid_angle$ACTIN$trial=="T3"), cols,parameters_angle_data,paste("ACTIN ",labelling_angle_data,sep=""),"ACTIN_T3_")
    
  }
  #Spheoid axial distribution
  {
    dire<-paste(dird,"spheroid_axial_distribution/",sep="")
    dir.create(dire)
    setwd(dire)
    
    plot_black_bg_histogram(subset(spheroid_axial$dna,spheroid_axial$dna$trial=="T1"), cols,parameters_axial_data,paste("DNA ",labelling_axial_data,sep=""),"DNA_T1_")
    plot_black_bg_histogram(subset(spheroid_axial$dna,spheroid_axial$dna$trial=="T2"), cols,parameters_axial_data,paste("DNA ",labelling_axial_data,sep=""),"DNA_T2_")
    plot_black_bg_histogram(subset(spheroid_axial$dna,spheroid_axial$dna$trial=="T3"), cols,parameters_axial_data,paste("DNA ",labelling_axial_data,sep=""),"DNA_T3_")
    
    plot_black_bg_histogram(subset(spheroid_axial$NANOG,spheroid_axial$NANOG$trial=="T1"), cols,parameters_axial_data,paste("NANOG ",labelling_axial_data,sep=""),"NANOG_T1_")
    plot_black_bg_histogram(subset(spheroid_axial$NANOG,spheroid_axial$NANOG$trial=="T2"), cols,parameters_axial_data,paste("NANOG ",labelling_axial_data,sep=""),"NANOG_T2_")
    plot_black_bg_histogram(subset(spheroid_axial$NANOG,spheroid_axial$NANOG$trial=="T3"), cols,parameters_axial_data,paste("NANOG ",labelling_axial_data,sep=""),"NANOG_T3_")
    
    plot_black_bg_histogram(subset(spheroid_axial$ACTIN,spheroid_axial$ACTIN$trial=="T1"), cols,parameters_axial_data,paste("ACTIN ",labelling_axial_data,sep=""),"ACTIN_T1_")
    plot_black_bg_histogram(subset(spheroid_axial$ACTIN,spheroid_axial$ACTIN$trial=="T2"), cols,parameters_axial_data,paste("ACTIN ",labelling_axial_data,sep=""),"ACTIN_T2_")
    plot_black_bg_histogram(subset(spheroid_axial$ACTIN,spheroid_axial$ACTIN$trial=="T3"), cols,parameters_axial_data,paste("ACTIN ",labelling_axial_data,sep=""),"ACTIN_T3_")
    
    plot_black_bg_histogram(subset(spheroid_axial$OCT4,spheroid_axial$OCT4$trial=="T1"), cols,parameters_axial_data,paste("OCT4 ",labelling_axial_data,sep=""),"OCT4_T1_")
    plot_black_bg_histogram(subset(spheroid_axial$OCT4,spheroid_axial$OCT4$trial=="T2"), cols,parameters_axial_data,paste("OCT4 ",labelling_axial_data,sep=""),"OCT4_T2_")
    plot_black_bg_histogram(subset(spheroid_axial$OCT4,spheroid_axial$OCT4$trial=="T3"), cols,parameters_axial_data,paste("OCT4 ",labelling_axial_data,sep=""),"OCT4_T3_")
    
    
    
  }
  #Spheoid radial distribution
  {
    dire<-paste(dird,"spheroid_radial_distribution/",sep="")
    dir.create(dire)
    setwd(dire)
    
    plot_black_bg_histogram(subset(spheroid_radial$dna,spheroid_radial$dna$trial=="T1"), cols,parameters_radial_data,paste("DNA ",labelling_radial_data,sep=""),"DNA_T1_")
    plot_black_bg_histogram(subset(spheroid_radial$dna,spheroid_radial$dna$trial=="T2"), cols,parameters_radial_data,paste("DNA ",labelling_radial_data,sep=""),"DNA_T2_")
    plot_black_bg_histogram(subset(spheroid_radial$dna,spheroid_radial$dna$trial=="T3"), cols,parameters_radial_data,paste("DNA ",labelling_radial_data,sep=""),"DNA_T3_")
    
    plot_black_bg_histogram(subset(spheroid_radial$NANOG,spheroid_radial$NANOG$trial=="T1"), cols,parameters_radial_data,paste("NANOG ",labelling_radial_data,sep=""),"NANOG_T1_")
    plot_black_bg_histogram(subset(spheroid_radial$NANOG,spheroid_radial$NANOG$trial=="T2"), cols,parameters_radial_data,paste("NANOG ",labelling_radial_data,sep=""),"NANOG_T2_")
    plot_black_bg_histogram(subset(spheroid_radial$NANOG,spheroid_radial$NANOG$trial=="T3"), cols,parameters_radial_data,paste("NANOG ",labelling_radial_data,sep=""),"NANOG_T3_")
    
    plot_black_bg_histogram(subset(spheroid_radial$OCT4,spheroid_radial$OCT4$trial=="T1"), cols,parameters_radial_data,paste("OCT4 ",labelling_radial_data,sep=""),"OCT4_T1_")
    plot_black_bg_histogram(subset(spheroid_radial$OCT4,spheroid_radial$OCT4$trial=="T2"), cols,parameters_radial_data,paste("OCT4 ",labelling_radial_data,sep=""),"OCT4_T2_")
    plot_black_bg_histogram(subset(spheroid_radial$OCT4,spheroid_radial$OCT4$trial=="T3"), cols,parameters_radial_data,paste("OCT4 ",labelling_radial_data,sep=""),"OCT4_T3_")
    
    plot_black_bg_histogram(subset(spheroid_radial$ACTIN,spheroid_radial$ACTIN$trial=="T1"), cols,parameters_radial_data,paste("ACTIN ",labelling_radial_data,sep=""),"ACTIN_T1_")
    plot_black_bg_histogram(subset(spheroid_radial$ACTIN,spheroid_radial$ACTIN$trial=="T2"), cols,parameters_radial_data,paste("ACTIN ",labelling_radial_data,sep=""),"ACTIN_T2_")
    plot_black_bg_histogram(subset(spheroid_radial$ACTIN,spheroid_radial$ACTIN$trial=="T3"), cols,parameters_radial_data,paste("ACTIN ",labelling_radial_data,sep=""),"ACTIN_T3_")
    
    
  }
  #Spheoid geometeric properties
  {
    dire<-paste(dird,"spheroid_properties/",sep="")
    dir.create(dire)
    setwd(dire)
    
    plot_black_bg_histogram(subset(spheroid_combined,spheroid_combined$trial=="T1"), cols,parameters_geometrical_data,paste("Spheroid ",labelling_geometrical_data,sep=""),"geo_T1_")
    plot_black_bg_histogram(subset(spheroid_combined,spheroid_combined$trial=="T2"), cols,parameters_geometrical_data,paste("Spheroid ",labelling_geometrical_data,sep=""),"geo_T2_")
    plot_black_bg_histogram(subset(spheroid_combined,spheroid_combined$trial=="T3"), cols,parameters_geometrical_data,paste("Spheroid ",labelling_geometrical_data,sep=""),"geo_T3_")
    
  }
  #Nuclear geometeric properties
  {
    dire<-paste(dird,"nuclear_properties/",sep="")
    dir.create(dire)
    setwd(dire)
    
    plot_black_bg_histogram(subset(nucleus_combined,nucleus_combined$trial=="T1"), cols,parameters_geo_texture_data,paste("Nuclear ",labelling_geo_texture_data,sep=""),"geo_T1_")
    plot_black_bg_histogram(subset(nucleus_combined,nucleus_combined$trial=="T2"), cols,parameters_geo_texture_data,paste("Nuclear ",labelling_geo_texture_data,sep=""),"geo_T2_")
    plot_black_bg_histogram(subset(nucleus_combined,nucleus_combined$trial=="T3"), cols,parameters_geo_texture_data,paste("Nuclear ",labelling_geo_texture_data,sep=""),"geo_T3_")
    
  }
  #Nuclear Proj_Int
  {
    
    dire<-paste(dird,"nuclear_proj_int/",sep="")
    dir.create(dire)
    setwd(dire)
    
    plot_black_bg_histogram(subset(nuclei_2dint_data$NANOG,nuclei_2dint_data$NANOG$trial=="T1"), cols,parameters_2d_int_data,paste("NANOG ",nuc_labelling_2d_int_data,sep=""),"NANOG_T1_")
    plot_black_bg_histogram(subset(nuclei_2dint_data$NANOG,nuclei_2dint_data$NANOG$trial=="T2"), cols,parameters_2d_int_data,paste("NANOG ",nuc_labelling_2d_int_data,sep=""),"NANOG_T2_")
    plot_black_bg_histogram(subset(nuclei_2dint_data$NANOG,nuclei_2dint_data$NANOG$trial=="T3"), cols,parameters_2d_int_data,paste("NANOG ",nuc_labelling_2d_int_data,sep=""),"NANOG_T3_")
    
    plot_black_bg_histogram(subset(nuclei_2dint_data$OCT4,nuclei_2dint_data$OCT4$trial=="T1"), cols,parameters_2d_int_data,paste("OCT4 ",nuc_labelling_2d_int_data,sep=""),"OCT4_T1_")
    plot_black_bg_histogram(subset(nuclei_2dint_data$OCT4,nuclei_2dint_data$OCT4$trial=="T2"), cols,parameters_2d_int_data,paste("OCT4 ",nuc_labelling_2d_int_data,sep=""),"OCT4_T2_")
    plot_black_bg_histogram(subset(nuclei_2dint_data$OCT4,nuclei_2dint_data$OCT4$trial=="T3"), cols,parameters_2d_int_data,paste("OCT4 ",nuc_labelling_2d_int_data,sep=""),"OCT4_T3_")
    
    plot_black_bg_histogram(subset(nuclei_2dint_data$ACTIN,nuclei_2dint_data$ACTIN$trial=="T1"), cols,parameters_2d_int_data,paste("ACTIN ",nuc_labelling_2d_int_data,sep=""),"ACTIN_T1_")
    plot_black_bg_histogram(subset(nuclei_2dint_data$ACTIN,nuclei_2dint_data$ACTIN$trial=="T2"), cols,parameters_2d_int_data,paste("ACTIN ",nuc_labelling_2d_int_data,sep=""),"ACTIN_T2_")
    plot_black_bg_histogram(subset(nuclei_2dint_data$ACTIN,nuclei_2dint_data$ACTIN$trial=="T3"), cols,parameters_2d_int_data,paste("ACTIN ",nuc_labelling_2d_int_data,sep=""),"ACTIN_T3_")
    
  }
  #Nuclear 3D_Int
  {
    dire<-paste(dird,"nuclear_3d_int/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_black_bg_histogram(subset(nuclei_3dint_data$dna,nuclei_3dint_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T1_")
    plot_black_bg_histogram(subset(nuclei_3dint_data$dna,nuclei_3dint_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T2_")
    plot_black_bg_histogram(subset(nuclei_3dint_data$dna,nuclei_3dint_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T3_")
    
    
    plot_black_bg_histogram(subset(nuclei_3dint_data$NANOG,nuclei_3dint_data$NANOG$trial=="T1"), cols,parameters_int_data,paste("NANOG ",nuclabelling_int_data,sep=""),"NANOG_T1_")
    plot_black_bg_histogram(subset(nuclei_3dint_data$NANOG,nuclei_3dint_data$NANOG$trial=="T2"), cols,parameters_int_data,paste("NANOG ",nuclabelling_int_data,sep=""),"NANOG_T2_")
    plot_black_bg_histogram(subset(nuclei_3dint_data$NANOG,nuclei_3dint_data$NANOG$trial=="T3"), cols,parameters_int_data,paste("NANOG ",nuclabelling_int_data,sep=""),"NANOG_T3_")
    
    plot_black_bg_histogram(subset(nuclei_3dint_data$OCT4,nuclei_3dint_data$OCT4$trial=="T1"), cols,parameters_int_data,paste("OCT4 ",nuclabelling_int_data,sep=""),"OCT4_T1_")
    plot_black_bg_histogram(subset(nuclei_3dint_data$OCT4,nuclei_3dint_data$OCT4$trial=="T2"), cols,parameters_int_data,paste("OCT4 ",nuclabelling_int_data,sep=""),"OCT4_T2_")
    plot_black_bg_histogram(subset(nuclei_3dint_data$OCT4,nuclei_3dint_data$OCT4$trial=="T3"), cols,parameters_int_data,paste("OCT4 ",nuclabelling_int_data,sep=""),"OCT4_T3_")
    
    
    plot_black_bg_histogram(subset(nuclei_3dint_data$ACTIN,nuclei_3dint_data$ACTIN$trial=="T1"), cols,parameters_int_data,paste("ACTIN ",nuclabelling_int_data,sep=""),"ACTIN_T1_")
    plot_black_bg_histogram(subset(nuclei_3dint_data$ACTIN,nuclei_3dint_data$ACTIN$trial=="T2"), cols,parameters_int_data,paste("ACTIN ",nuclabelling_int_data,sep=""),"ACTIN_T2_")
    plot_black_bg_histogram(subset(nuclei_3dint_data$ACTIN,nuclei_3dint_data$ACTIN$trial=="T3"), cols,parameters_int_data,paste("ACTIN ",nuclabelling_int_data,sep=""),"ACTIN_T3_")
    
  }
  
  #Spheroid_median_nuclear_properties
  {
    dire<-paste(dird,"median_nuclear_properties/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_black_bg_histogram(subset(spheroid_nucleus_combined,spheroid_nucleus_combined$trial=="T1"), cols,parameters_geo_texture_data,paste("Nuclear ",labelling_geo_texture_data,sep=""),"geo_T1_")
    plot_black_bg_histogram(subset(spheroid_nucleus_combined,spheroid_nucleus_combined$trial=="T2"), cols,parameters_geo_texture_data,paste("Nuclear ",labelling_geo_texture_data,sep=""),"geo_T2_")
    plot_black_bg_histogram(subset(spheroid_nucleus_combined,spheroid_nucleus_combined$trial=="T3"), cols,parameters_geo_texture_data,paste("Nuclear ",labelling_geo_texture_data,sep=""),"geo_T3_")
    
    
  }
  #Spheroid_Nuclear Proj_Int
  {
    
    dire<-paste(dird,"median_nuclear_proj_int/",sep="")
    dir.create(dire)
    setwd(dire)
    
    plot_black_bg_histogram(subset(spheroid_nuclei_2dint_data$ACTIN,spheroid_nuclei_2dint_data$ACTIN$trial=="T1"), cols,parameters_2d_int_data,paste("ACTIN ",nuc_labelling_2d_int_data,sep=""),"ACTIN_T1_")
    plot_black_bg_histogram(subset(spheroid_nuclei_2dint_data$ACTIN,spheroid_nuclei_2dint_data$ACTIN$trial=="T2"), cols,parameters_2d_int_data,paste("ACTIN ",nuc_labelling_2d_int_data,sep=""),"ACTIN_T2_")
    plot_black_bg_histogram(subset(spheroid_nuclei_2dint_data$ACTIN,spheroid_nuclei_2dint_data$ACTIN$trial=="T3"), cols,parameters_2d_int_data,paste("ACTIN ",nuc_labelling_2d_int_data,sep=""),"ACTIN_T3_")
    
    plot_black_bg_histogram(subset(spheroid_nuclei_2dint_data$NANOG,spheroid_nuclei_2dint_data$NANOG$trial=="T1"), cols,parameters_2d_int_data,paste("NANOG ",nuc_labelling_2d_int_data,sep=""),"NANOG_T1_")
    plot_black_bg_histogram(subset(spheroid_nuclei_2dint_data$NANOG,spheroid_nuclei_2dint_data$NANOG$trial=="T2"), cols,parameters_2d_int_data,paste("NANOG ",nuc_labelling_2d_int_data,sep=""),"NANOG_T2_")
    plot_black_bg_histogram(subset(spheroid_nuclei_2dint_data$NANOG,spheroid_nuclei_2dint_data$NANOG$trial=="T3"), cols,parameters_2d_int_data,paste("NANOG ",nuc_labelling_2d_int_data,sep=""),"NANOG_T3_")
    
    plot_black_bg_histogram(subset(spheroid_nuclei_2dint_data$OCT4,spheroid_nuclei_2dint_data$OCT4$trial=="T1"), cols,parameters_2d_int_data,paste("OCT4 ",nuc_labelling_2d_int_data,sep=""),"OCT4_T1_")
    plot_black_bg_histogram(subset(spheroid_nuclei_2dint_data$OCT4,spheroid_nuclei_2dint_data$OCT4$trial=="T2"), cols,parameters_2d_int_data,paste("OCT4 ",nuc_labelling_2d_int_data,sep=""),"OCT4_T2_")
    plot_black_bg_histogram(subset(spheroid_nuclei_2dint_data$OCT4,spheroid_nuclei_2dint_data$OCT4$trial=="T3"), cols,parameters_2d_int_data,paste("OCT4 ",nuc_labelling_2d_int_data,sep=""),"OCT4_T3_")
    
  }
  #spheroid_median_Nuclear 3D_Int
  {
    dire<-paste(dird,"median_nuclear_3d_int/",sep="")
    dir.create(dire)
    setwd(dire)
    
    plot_black_bg_histogram(subset(spheroid_nuclei_3dint_data$dna,spheroid_nuclei_3dint_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T1_")
    plot_black_bg_histogram(subset(spheroid_nuclei_3dint_data$dna,spheroid_nuclei_3dint_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T2_")
    plot_black_bg_histogram(subset(spheroid_nuclei_3dint_data$dna,spheroid_nuclei_3dint_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T3_")
    
    plot_black_bg_histogram(subset(spheroid_nuclei_3dint_data$OCT4,spheroid_nuclei_3dint_data$OCT4$trial=="T1"), cols,parameters_int_data,paste("OCT4 ",nuclabelling_int_data,sep=""),"OCT4_T1_")
    plot_black_bg_histogram(subset(spheroid_nuclei_3dint_data$OCT4,spheroid_nuclei_3dint_data$OCT4$trial=="T2"), cols,parameters_int_data,paste("OCT4 ",nuclabelling_int_data,sep=""),"OCT4_T2_")
    plot_black_bg_histogram(subset(spheroid_nuclei_3dint_data$OCT4,spheroid_nuclei_3dint_data$OCT4$trial=="T3"), cols,parameters_int_data,paste("OCT4 ",nuclabelling_int_data,sep=""),"OCT4_T3_")
    
    plot_black_bg_histogram(subset(spheroid_nuclei_3dint_data$NANOG,spheroid_nuclei_3dint_data$NANOG$trial=="T1"), cols,parameters_int_data,paste("NANOG ",nuclabelling_int_data,sep=""),"NANOG_T1_")
    plot_black_bg_histogram(subset(spheroid_nuclei_3dint_data$NANOG,spheroid_nuclei_3dint_data$NANOG$trial=="T2"), cols,parameters_int_data,paste("NANOG ",nuclabelling_int_data,sep=""),"NANOG_T2_")
    plot_black_bg_histogram(subset(spheroid_nuclei_3dint_data$NANOG,spheroid_nuclei_3dint_data$NANOG$trial=="T3"), cols,parameters_int_data,paste("NANOG ",nuclabelling_int_data,sep=""),"NANOG_T3_")
    
    plot_black_bg_histogram(subset(spheroid_nuclei_3dint_data$ACTIN,spheroid_nuclei_3dint_data$ACTIN$trial=="T1"), cols,parameters_int_data,paste("ACTIN ",nuclabelling_int_data,sep=""),"ACTIN_T1_")
    plot_black_bg_histogram(subset(spheroid_nuclei_3dint_data$ACTIN,spheroid_nuclei_3dint_data$ACTIN$trial=="T2"), cols,parameters_int_data,paste("ACTIN ",nuclabelling_int_data,sep=""),"ACTIN_T2_")
    plot_black_bg_histogram(subset(spheroid_nuclei_3dint_data$ACTIN,spheroid_nuclei_3dint_data$ACTIN$trial=="T3"), cols,parameters_int_data,paste("ACTIN ",nuclabelling_int_data,sep=""),"ACTIN_T3_")
    
    
    
  }
  
  
  dird<-"E:/HMF3A_reprogramming/R_analysis/NANOG/white_bg_plots/"
  dir.create(dird)
  setwd(dird)
  cols<-colorRampPalette(c("red","yellow"))( 5) 
  
  #Spheoid Proj_Int
  {
    dire<-paste(dird,"spheroid_proj_int/",sep="")
    dir.create(dire)
    setwd(dire)
    
    plot_white_bg_histogram(subset(spheroid_2dint_data$NANOG,spheroid_2dint_data$NANOG$trial=="T1"), cols,parameters_2d_int_data,paste("NANOG ",labelling_2d_int_data,sep=""),"NANOG_T1_")
    plot_white_bg_histogram(subset(spheroid_2dint_data$NANOG,spheroid_2dint_data$NANOG$trial=="T2"), cols,parameters_2d_int_data,paste("NANOG ",labelling_2d_int_data,sep=""),"NANOG_T2_")
    plot_white_bg_histogram(subset(spheroid_2dint_data$NANOG,spheroid_2dint_data$NANOG$trial=="T3"), cols,parameters_2d_int_data,paste("NANOG ",labelling_2d_int_data,sep=""),"NANOG_T3_")
    
    plot_white_bg_histogram(subset(spheroid_2dint_data$OCT4,spheroid_2dint_data$OCT4$trial=="T1"), cols,parameters_2d_int_data,paste("OCT4 ",labelling_2d_int_data,sep=""),"OCT4_T1_")
    plot_white_bg_histogram(subset(spheroid_2dint_data$OCT4,spheroid_2dint_data$OCT4$trial=="T2"), cols,parameters_2d_int_data,paste("OCT4 ",labelling_2d_int_data,sep=""),"OCT4_T2_")
    plot_white_bg_histogram(subset(spheroid_2dint_data$OCT4,spheroid_2dint_data$OCT4$trial=="T3"), cols,parameters_2d_int_data,paste("OCT4 ",labelling_2d_int_data,sep=""),"OCT4_T3_")
    
    plot_white_bg_histogram(subset(spheroid_2dint_data$ACTIN,spheroid_2dint_data$ACTIN$trial=="T1"), cols,parameters_2d_int_data,paste("ACTIN ",labelling_2d_int_data,sep=""),"ACTIN_T1_")
    plot_white_bg_histogram(subset(spheroid_2dint_data$ACTIN,spheroid_2dint_data$ACTIN$trial=="T2"), cols,parameters_2d_int_data,paste("ACTIN ",labelling_2d_int_data,sep=""),"ACTIN_T2_")
    plot_white_bg_histogram(subset(spheroid_2dint_data$ACTIN,spheroid_2dint_data$ACTIN$trial=="T3"), cols,parameters_2d_int_data,paste("ACTIN ",labelling_2d_int_data,sep=""),"ACTIN_T3_")
  }
  #Spheoid 3D_Int
  {
    dire<-paste(dird,"spheroid_3d_int/",sep="")
    dir.create(dire)
    setwd(dire)
    
    plot_white_bg_histogram(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T1_")
    plot_white_bg_histogram(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T2_")
    plot_white_bg_histogram(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T3_")
    
    
    plot_white_bg_histogram(subset(spheroid_int_data$NANOG,spheroid_int_data$NANOG$trial=="T1"), cols,parameters_int_data,paste("NANOG ",labelling_int_data,sep=""),"NANOG_T1_")
    plot_white_bg_histogram(subset(spheroid_int_data$NANOG,spheroid_int_data$NANOG$trial=="T2"), cols,parameters_int_data,paste("NANOG ",labelling_int_data,sep=""),"NANOG_T2_")
    plot_white_bg_histogram(subset(spheroid_int_data$NANOG,spheroid_int_data$NANOG$trial=="T3"), cols,parameters_int_data,paste("NANOG ",labelling_int_data,sep=""),"NANOG_T3_")
    
    plot_white_bg_histogram(subset(spheroid_int_data$OCT4,spheroid_int_data$OCT4$trial=="T1"), cols,parameters_int_data,paste("OCT4 ",labelling_int_data,sep=""),"OCT4_T1_")
    plot_white_bg_histogram(subset(spheroid_int_data$OCT4,spheroid_int_data$OCT4$trial=="T2"), cols,parameters_int_data,paste("OCT4 ",labelling_int_data,sep=""),"OCT4_T2_")
    plot_white_bg_histogram(subset(spheroid_int_data$OCT4,spheroid_int_data$OCT4$trial=="T3"), cols,parameters_int_data,paste("OCT4 ",labelling_int_data,sep=""),"OCT4_T3_")
    
    
    plot_white_bg_histogram(subset(spheroid_int_data$ACTIN,spheroid_int_data$ACTIN$trial=="T1"), cols,parameters_int_data,paste("ACTIN ",labelling_int_data,sep=""),"ACTIN_T1_")
    plot_white_bg_histogram(subset(spheroid_int_data$ACTIN,spheroid_int_data$ACTIN$trial=="T2"), cols,parameters_int_data,paste("ACTIN ",labelling_int_data,sep=""),"ACTIN_T2_")
    plot_white_bg_histogram(subset(spheroid_int_data$ACTIN,spheroid_int_data$ACTIN$trial=="T3"), cols,parameters_int_data,paste("ACTIN ",labelling_int_data,sep=""),"ACTIN_T3_")
    
  }
  #Spheoid angles
  {
    dire<-paste(dird,"spheroid_angles/",sep="")
    dir.create(dire)
    setwd(dire)
    
    plot_white_bg_histogram(subset(spheroid_angle$dna,spheroid_angle$dna$trial=="T1"), cols,parameters_angle_data,paste("DNA ",labelling_angle_data,sep=""),"DNA_T1_")
    plot_white_bg_histogram(subset(spheroid_angle$dna,spheroid_angle$dna$trial=="T2"), cols,parameters_angle_data,paste("DNA ",labelling_angle_data,sep=""),"DNA_T2_")
    plot_white_bg_histogram(subset(spheroid_angle$dna,spheroid_angle$dna$trial=="T3"), cols,parameters_angle_data,paste("DNA ",labelling_angle_data,sep=""),"DNA_T3_")
    
    plot_white_bg_histogram(subset(spheroid_angle$NANOG,spheroid_angle$NANOG$trial=="T1"), cols,parameters_angle_data,paste("NANOG ",labelling_angle_data,sep=""),"NANOG_T1_")
    plot_white_bg_histogram(subset(spheroid_angle$NANOG,spheroid_angle$NANOG$trial=="T2"), cols,parameters_angle_data,paste("NANOG ",labelling_angle_data,sep=""),"NANOG_T2_")
    plot_white_bg_histogram(subset(spheroid_angle$NANOG,spheroid_angle$NANOG$trial=="T3"), cols,parameters_angle_data,paste("NANOG ",labelling_angle_data,sep=""),"NANOG_T3_")
    
    plot_white_bg_histogram(subset(spheroid_angle$OCT4,spheroid_angle$OCT4$trial=="T1"), cols,parameters_angle_data,paste("OCT4 ",labelling_angle_data,sep=""),"OCT4_T1_")
    plot_white_bg_histogram(subset(spheroid_angle$OCT4,spheroid_angle$OCT4$trial=="T2"), cols,parameters_angle_data,paste("OCT4 ",labelling_angle_data,sep=""),"OCT4_T2_")
    plot_white_bg_histogram(subset(spheroid_angle$OCT4,spheroid_angle$OCT4$trial=="T3"), cols,parameters_angle_data,paste("OCT4 ",labelling_angle_data,sep=""),"OCT4_T3_")
    
    plot_white_bg_histogram(subset(spheroid_angle$ACTIN,spheroid_angle$ACTIN$trial=="T1"), cols,parameters_angle_data,paste("ACTIN ",labelling_angle_data,sep=""),"ACTIN_T1_")
    plot_white_bg_histogram(subset(spheroid_angle$ACTIN,spheroid_angle$ACTIN$trial=="T2"), cols,parameters_angle_data,paste("ACTIN ",labelling_angle_data,sep=""),"ACTIN_T2_")
    plot_white_bg_histogram(subset(spheroid_angle$ACTIN,spheroid_angle$ACTIN$trial=="T3"), cols,parameters_angle_data,paste("ACTIN ",labelling_angle_data,sep=""),"ACTIN_T3_")
    
  }
  #Spheoid axial distribution
  {
    dire<-paste(dird,"spheroid_axial_distribution/",sep="")
    dir.create(dire)
    setwd(dire)
    
    plot_white_bg_histogram(subset(spheroid_axial$dna,spheroid_axial$dna$trial=="T1"), cols,parameters_axial_data,paste("DNA ",labelling_axial_data,sep=""),"DNA_T1_")
    plot_white_bg_histogram(subset(spheroid_axial$dna,spheroid_axial$dna$trial=="T2"), cols,parameters_axial_data,paste("DNA ",labelling_axial_data,sep=""),"DNA_T2_")
    plot_white_bg_histogram(subset(spheroid_axial$dna,spheroid_axial$dna$trial=="T3"), cols,parameters_axial_data,paste("DNA ",labelling_axial_data,sep=""),"DNA_T3_")
    
    plot_white_bg_histogram(subset(spheroid_axial$NANOG,spheroid_axial$NANOG$trial=="T1"), cols,parameters_axial_data,paste("NANOG ",labelling_axial_data,sep=""),"NANOG_T1_")
    plot_white_bg_histogram(subset(spheroid_axial$NANOG,spheroid_axial$NANOG$trial=="T2"), cols,parameters_axial_data,paste("NANOG ",labelling_axial_data,sep=""),"NANOG_T2_")
    plot_white_bg_histogram(subset(spheroid_axial$NANOG,spheroid_axial$NANOG$trial=="T3"), cols,parameters_axial_data,paste("NANOG ",labelling_axial_data,sep=""),"NANOG_T3_")
    
    plot_white_bg_histogram(subset(spheroid_axial$ACTIN,spheroid_axial$ACTIN$trial=="T1"), cols,parameters_axial_data,paste("ACTIN ",labelling_axial_data,sep=""),"ACTIN_T1_")
    plot_white_bg_histogram(subset(spheroid_axial$ACTIN,spheroid_axial$ACTIN$trial=="T2"), cols,parameters_axial_data,paste("ACTIN ",labelling_axial_data,sep=""),"ACTIN_T2_")
    plot_white_bg_histogram(subset(spheroid_axial$ACTIN,spheroid_axial$ACTIN$trial=="T3"), cols,parameters_axial_data,paste("ACTIN ",labelling_axial_data,sep=""),"ACTIN_T3_")
    
    plot_white_bg_histogram(subset(spheroid_axial$OCT4,spheroid_axial$OCT4$trial=="T1"), cols,parameters_axial_data,paste("OCT4 ",labelling_axial_data,sep=""),"OCT4_T1_")
    plot_white_bg_histogram(subset(spheroid_axial$OCT4,spheroid_axial$OCT4$trial=="T2"), cols,parameters_axial_data,paste("OCT4 ",labelling_axial_data,sep=""),"OCT4_T2_")
    plot_white_bg_histogram(subset(spheroid_axial$OCT4,spheroid_axial$OCT4$trial=="T3"), cols,parameters_axial_data,paste("OCT4 ",labelling_axial_data,sep=""),"OCT4_T3_")
    
    
    
  }
  #Spheoid radial distribution
  {
    dire<-paste(dird,"spheroid_radial_distribution/",sep="")
    dir.create(dire)
    setwd(dire)
    
    plot_white_bg_histogram(subset(spheroid_radial$dna,spheroid_radial$dna$trial=="T1"), cols,parameters_radial_data,paste("DNA ",labelling_radial_data,sep=""),"DNA_T1_")
    plot_white_bg_histogram(subset(spheroid_radial$dna,spheroid_radial$dna$trial=="T2"), cols,parameters_radial_data,paste("DNA ",labelling_radial_data,sep=""),"DNA_T2_")
    plot_white_bg_histogram(subset(spheroid_radial$dna,spheroid_radial$dna$trial=="T3"), cols,parameters_radial_data,paste("DNA ",labelling_radial_data,sep=""),"DNA_T3_")
    
    plot_white_bg_histogram(subset(spheroid_radial$NANOG,spheroid_radial$NANOG$trial=="T1"), cols,parameters_radial_data,paste("NANOG ",labelling_radial_data,sep=""),"NANOG_T1_")
    plot_white_bg_histogram(subset(spheroid_radial$NANOG,spheroid_radial$NANOG$trial=="T2"), cols,parameters_radial_data,paste("NANOG ",labelling_radial_data,sep=""),"NANOG_T2_")
    plot_white_bg_histogram(subset(spheroid_radial$NANOG,spheroid_radial$NANOG$trial=="T3"), cols,parameters_radial_data,paste("NANOG ",labelling_radial_data,sep=""),"NANOG_T3_")
    
    plot_white_bg_histogram(subset(spheroid_radial$OCT4,spheroid_radial$OCT4$trial=="T1"), cols,parameters_radial_data,paste("OCT4 ",labelling_radial_data,sep=""),"OCT4_T1_")
    plot_white_bg_histogram(subset(spheroid_radial$OCT4,spheroid_radial$OCT4$trial=="T2"), cols,parameters_radial_data,paste("OCT4 ",labelling_radial_data,sep=""),"OCT4_T2_")
    plot_white_bg_histogram(subset(spheroid_radial$OCT4,spheroid_radial$OCT4$trial=="T3"), cols,parameters_radial_data,paste("OCT4 ",labelling_radial_data,sep=""),"OCT4_T3_")
    
    plot_white_bg_histogram(subset(spheroid_radial$ACTIN,spheroid_radial$ACTIN$trial=="T1"), cols,parameters_radial_data,paste("ACTIN ",labelling_radial_data,sep=""),"ACTIN_T1_")
    plot_white_bg_histogram(subset(spheroid_radial$ACTIN,spheroid_radial$ACTIN$trial=="T2"), cols,parameters_radial_data,paste("ACTIN ",labelling_radial_data,sep=""),"ACTIN_T2_")
    plot_white_bg_histogram(subset(spheroid_radial$ACTIN,spheroid_radial$ACTIN$trial=="T3"), cols,parameters_radial_data,paste("ACTIN ",labelling_radial_data,sep=""),"ACTIN_T3_")
    
    
  }
  #Spheoid geometeric properties
  {
    dire<-paste(dird,"spheroid_properties/",sep="")
    dir.create(dire)
    setwd(dire)
    
    plot_white_bg_histogram(subset(spheroid_combined,spheroid_combined$trial=="T1"), cols,parameters_geometrical_data,paste("Spheroid ",labelling_geometrical_data,sep=""),"geo_T1_")
    plot_white_bg_histogram(subset(spheroid_combined,spheroid_combined$trial=="T2"), cols,parameters_geometrical_data,paste("Spheroid ",labelling_geometrical_data,sep=""),"geo_T2_")
    plot_white_bg_histogram(subset(spheroid_combined,spheroid_combined$trial=="T3"), cols,parameters_geometrical_data,paste("Spheroid ",labelling_geometrical_data,sep=""),"geo_T3_")
    
  }
  #Nuclear geometeric properties
  {
    dire<-paste(dird,"nuclear_properties/",sep="")
    dir.create(dire)
    setwd(dire)
    
    plot_white_bg_histogram(subset(nucleus_combined,nucleus_combined$trial=="T1"), cols,parameters_geo_texture_data,paste("Nuclear ",labelling_geo_texture_data,sep=""),"geo_T1_")
    plot_white_bg_histogram(subset(nucleus_combined,nucleus_combined$trial=="T2"), cols,parameters_geo_texture_data,paste("Nuclear ",labelling_geo_texture_data,sep=""),"geo_T2_")
    plot_white_bg_histogram(subset(nucleus_combined,nucleus_combined$trial=="T3"), cols,parameters_geo_texture_data,paste("Nuclear ",labelling_geo_texture_data,sep=""),"geo_T3_")
    
  }
  #Nuclear Proj_Int
  {
    
    dire<-paste(dird,"nuclear_proj_int/",sep="")
    dir.create(dire)
    setwd(dire)
    
    plot_white_bg_histogram(subset(nuclei_2dint_data$NANOG,nuclei_2dint_data$NANOG$trial=="T1"), cols,parameters_2d_int_data,paste("NANOG ",nuc_labelling_2d_int_data,sep=""),"NANOG_T1_")
    plot_white_bg_histogram(subset(nuclei_2dint_data$NANOG,nuclei_2dint_data$NANOG$trial=="T2"), cols,parameters_2d_int_data,paste("NANOG ",nuc_labelling_2d_int_data,sep=""),"NANOG_T2_")
    plot_white_bg_histogram(subset(nuclei_2dint_data$NANOG,nuclei_2dint_data$NANOG$trial=="T3"), cols,parameters_2d_int_data,paste("NANOG ",nuc_labelling_2d_int_data,sep=""),"NANOG_T3_")
    
    plot_white_bg_histogram(subset(nuclei_2dint_data$OCT4,nuclei_2dint_data$OCT4$trial=="T1"), cols,parameters_2d_int_data,paste("OCT4 ",nuc_labelling_2d_int_data,sep=""),"OCT4_T1_")
    plot_white_bg_histogram(subset(nuclei_2dint_data$OCT4,nuclei_2dint_data$OCT4$trial=="T2"), cols,parameters_2d_int_data,paste("OCT4 ",nuc_labelling_2d_int_data,sep=""),"OCT4_T2_")
    plot_white_bg_histogram(subset(nuclei_2dint_data$OCT4,nuclei_2dint_data$OCT4$trial=="T3"), cols,parameters_2d_int_data,paste("OCT4 ",nuc_labelling_2d_int_data,sep=""),"OCT4_T3_")
    
    plot_white_bg_histogram(subset(nuclei_2dint_data$ACTIN,nuclei_2dint_data$ACTIN$trial=="T1"), cols,parameters_2d_int_data,paste("ACTIN ",nuc_labelling_2d_int_data,sep=""),"ACTIN_T1_")
    plot_white_bg_histogram(subset(nuclei_2dint_data$ACTIN,nuclei_2dint_data$ACTIN$trial=="T2"), cols,parameters_2d_int_data,paste("ACTIN ",nuc_labelling_2d_int_data,sep=""),"ACTIN_T2_")
    plot_white_bg_histogram(subset(nuclei_2dint_data$ACTIN,nuclei_2dint_data$ACTIN$trial=="T3"), cols,parameters_2d_int_data,paste("ACTIN ",nuc_labelling_2d_int_data,sep=""),"ACTIN_T3_")
    
  }
  #Nuclear 3D_Int
  {
    dire<-paste(dird,"nuclear_3d_int/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_white_bg_histogram(subset(nuclei_3dint_data$dna,nuclei_3dint_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T1_")
    plot_white_bg_histogram(subset(nuclei_3dint_data$dna,nuclei_3dint_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T2_")
    plot_white_bg_histogram(subset(nuclei_3dint_data$dna,nuclei_3dint_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T3_")
    
    
    plot_white_bg_histogram(subset(nuclei_3dint_data$NANOG,nuclei_3dint_data$NANOG$trial=="T1"), cols,parameters_int_data,paste("NANOG ",nuclabelling_int_data,sep=""),"NANOG_T1_")
    plot_white_bg_histogram(subset(nuclei_3dint_data$NANOG,nuclei_3dint_data$NANOG$trial=="T2"), cols,parameters_int_data,paste("NANOG ",nuclabelling_int_data,sep=""),"NANOG_T2_")
    plot_white_bg_histogram(subset(nuclei_3dint_data$NANOG,nuclei_3dint_data$NANOG$trial=="T3"), cols,parameters_int_data,paste("NANOG ",nuclabelling_int_data,sep=""),"NANOG_T3_")
    
    plot_white_bg_histogram(subset(nuclei_3dint_data$OCT4,nuclei_3dint_data$OCT4$trial=="T1"), cols,parameters_int_data,paste("OCT4 ",nuclabelling_int_data,sep=""),"OCT4_T1_")
    plot_white_bg_histogram(subset(nuclei_3dint_data$OCT4,nuclei_3dint_data$OCT4$trial=="T2"), cols,parameters_int_data,paste("OCT4 ",nuclabelling_int_data,sep=""),"OCT4_T2_")
    plot_white_bg_histogram(subset(nuclei_3dint_data$OCT4,nuclei_3dint_data$OCT4$trial=="T3"), cols,parameters_int_data,paste("OCT4 ",nuclabelling_int_data,sep=""),"OCT4_T3_")
    
    
    plot_white_bg_histogram(subset(nuclei_3dint_data$ACTIN,nuclei_3dint_data$ACTIN$trial=="T1"), cols,parameters_int_data,paste("ACTIN ",nuclabelling_int_data,sep=""),"ACTIN_T1_")
    plot_white_bg_histogram(subset(nuclei_3dint_data$ACTIN,nuclei_3dint_data$ACTIN$trial=="T2"), cols,parameters_int_data,paste("ACTIN ",nuclabelling_int_data,sep=""),"ACTIN_T2_")
    plot_white_bg_histogram(subset(nuclei_3dint_data$ACTIN,nuclei_3dint_data$ACTIN$trial=="T3"), cols,parameters_int_data,paste("ACTIN ",nuclabelling_int_data,sep=""),"ACTIN_T3_")
    
  }
  
  #Spheroid_median_nuclear_properties
  {
    dire<-paste(dird,"median_nuclear_properties/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_white_bg_histogram(subset(spheroid_nucleus_combined,spheroid_nucleus_combined$trial=="T1"), cols,parameters_geo_texture_data,paste("Nuclear ",labelling_geo_texture_data,sep=""),"geo_T1_")
    plot_white_bg_histogram(subset(spheroid_nucleus_combined,spheroid_nucleus_combined$trial=="T2"), cols,parameters_geo_texture_data,paste("Nuclear ",labelling_geo_texture_data,sep=""),"geo_T2_")
    plot_white_bg_histogram(subset(spheroid_nucleus_combined,spheroid_nucleus_combined$trial=="T3"), cols,parameters_geo_texture_data,paste("Nuclear ",labelling_geo_texture_data,sep=""),"geo_T3_")
    
    
  }
  #Spheroid_Nuclear Proj_Int
  {
    
    dire<-paste(dird,"median_nuclear_proj_int/",sep="")
    dir.create(dire)
    setwd(dire)
    
    plot_white_bg_histogram(subset(spheroid_nuclei_2dint_data$ACTIN,spheroid_nuclei_2dint_data$ACTIN$trial=="T1"), cols,parameters_2d_int_data,paste("ACTIN ",nuc_labelling_2d_int_data,sep=""),"ACTIN_T1_")
    plot_white_bg_histogram(subset(spheroid_nuclei_2dint_data$ACTIN,spheroid_nuclei_2dint_data$ACTIN$trial=="T2"), cols,parameters_2d_int_data,paste("ACTIN ",nuc_labelling_2d_int_data,sep=""),"ACTIN_T2_")
    plot_white_bg_histogram(subset(spheroid_nuclei_2dint_data$ACTIN,spheroid_nuclei_2dint_data$ACTIN$trial=="T3"), cols,parameters_2d_int_data,paste("ACTIN ",nuc_labelling_2d_int_data,sep=""),"ACTIN_T3_")
    
    plot_white_bg_histogram(subset(spheroid_nuclei_2dint_data$NANOG,spheroid_nuclei_2dint_data$NANOG$trial=="T1"), cols,parameters_2d_int_data,paste("NANOG ",nuc_labelling_2d_int_data,sep=""),"NANOG_T1_")
    plot_white_bg_histogram(subset(spheroid_nuclei_2dint_data$NANOG,spheroid_nuclei_2dint_data$NANOG$trial=="T2"), cols,parameters_2d_int_data,paste("NANOG ",nuc_labelling_2d_int_data,sep=""),"NANOG_T2_")
    plot_white_bg_histogram(subset(spheroid_nuclei_2dint_data$NANOG,spheroid_nuclei_2dint_data$NANOG$trial=="T3"), cols,parameters_2d_int_data,paste("NANOG ",nuc_labelling_2d_int_data,sep=""),"NANOG_T3_")
    
    plot_white_bg_histogram(subset(spheroid_nuclei_2dint_data$OCT4,spheroid_nuclei_2dint_data$OCT4$trial=="T1"), cols,parameters_2d_int_data,paste("OCT4 ",nuc_labelling_2d_int_data,sep=""),"OCT4_T1_")
    plot_white_bg_histogram(subset(spheroid_nuclei_2dint_data$OCT4,spheroid_nuclei_2dint_data$OCT4$trial=="T2"), cols,parameters_2d_int_data,paste("OCT4 ",nuc_labelling_2d_int_data,sep=""),"OCT4_T2_")
    plot_white_bg_histogram(subset(spheroid_nuclei_2dint_data$OCT4,spheroid_nuclei_2dint_data$OCT4$trial=="T3"), cols,parameters_2d_int_data,paste("OCT4 ",nuc_labelling_2d_int_data,sep=""),"OCT4_T3_")
    
  }
  #spheroid_median_Nuclear 3D_Int
  {
    dire<-paste(dird,"median_nuclear_3d_int/",sep="")
    dir.create(dire)
    setwd(dire)
    
    plot_white_bg_histogram(subset(spheroid_nuclei_3dint_data$dna,spheroid_nuclei_3dint_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T1_")
    plot_white_bg_histogram(subset(spheroid_nuclei_3dint_data$dna,spheroid_nuclei_3dint_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T2_")
    plot_white_bg_histogram(subset(spheroid_nuclei_3dint_data$dna,spheroid_nuclei_3dint_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T3_")
    
    plot_white_bg_histogram(subset(spheroid_nuclei_3dint_data$OCT4,spheroid_nuclei_3dint_data$OCT4$trial=="T1"), cols,parameters_int_data,paste("OCT4 ",nuclabelling_int_data,sep=""),"OCT4_T1_")
    plot_white_bg_histogram(subset(spheroid_nuclei_3dint_data$OCT4,spheroid_nuclei_3dint_data$OCT4$trial=="T2"), cols,parameters_int_data,paste("OCT4 ",nuclabelling_int_data,sep=""),"OCT4_T2_")
    plot_white_bg_histogram(subset(spheroid_nuclei_3dint_data$OCT4,spheroid_nuclei_3dint_data$OCT4$trial=="T3"), cols,parameters_int_data,paste("OCT4 ",nuclabelling_int_data,sep=""),"OCT4_T3_")
    
    plot_white_bg_histogram(subset(spheroid_nuclei_3dint_data$NANOG,spheroid_nuclei_3dint_data$NANOG$trial=="T1"), cols,parameters_int_data,paste("NANOG ",nuclabelling_int_data,sep=""),"NANOG_T1_")
    plot_white_bg_histogram(subset(spheroid_nuclei_3dint_data$NANOG,spheroid_nuclei_3dint_data$NANOG$trial=="T2"), cols,parameters_int_data,paste("NANOG ",nuclabelling_int_data,sep=""),"NANOG_T2_")
    plot_white_bg_histogram(subset(spheroid_nuclei_3dint_data$NANOG,spheroid_nuclei_3dint_data$NANOG$trial=="T3"), cols,parameters_int_data,paste("NANOG ",nuclabelling_int_data,sep=""),"NANOG_T3_")
    
    plot_white_bg_histogram(subset(spheroid_nuclei_3dint_data$ACTIN,spheroid_nuclei_3dint_data$ACTIN$trial=="T1"), cols,parameters_int_data,paste("ACTIN ",nuclabelling_int_data,sep=""),"ACTIN_T1_")
    plot_white_bg_histogram(subset(spheroid_nuclei_3dint_data$ACTIN,spheroid_nuclei_3dint_data$ACTIN$trial=="T2"), cols,parameters_int_data,paste("ACTIN ",nuclabelling_int_data,sep=""),"ACTIN_T2_")
    plot_white_bg_histogram(subset(spheroid_nuclei_3dint_data$ACTIN,spheroid_nuclei_3dint_data$ACTIN$trial=="T3"), cols,parameters_int_data,paste("ACTIN ",nuclabelling_int_data,sep=""),"ACTIN_T3_")
    
    
    
  }
}


#write the collected data
{
  dir.create("E:/HMF3A_reprogramming/R_analysis/NANOG/combined_data/")
  setwd("E:/HMF3A_reprogramming/R_analysis/NANOG/combined_data/")
  
  write.csv(spheroid_combined, file="spheroid_combined.csv")
  write.csv(nucleus_combined, file="nucleus_combined.csv")
  write.csv(spheroid_nucleus_combined, file="spheroid_median_nucleus_combined.csv")
  
  write.csv(spheroid_radial$NANOG, file="spheroid_radial_dist_NANOG.csv")
  write.csv(spheroid_radial$ACTIN, file="spheroid_radial_dist_ACTIN.csv")
  write.csv(spheroid_radial$OCT4, file="spheroid_radial_dist_OCT4.csv")
  write.csv(spheroid_radial$dna, file="spheroid_radial_dist_DNA.csv")
  
  write.csv(spheroid_2dint_data$NANOG, file="spheroid_2d_int_NANOG.csv")
  write.csv(nuclei_2dint_data$NANOG, file="nuclear_2d_int_NANOG.csv")
  write.csv(spheroid_nuclei_2dint_data$NANOG, file="spheroid_median_nuclear_2d_int_NANOG.csv")
  
  write.csv(spheroid_2dint_data$ACTIN, file="spheroid_2d_int_ACTIN.csv")
  write.csv(nuclei_2dint_data$ACTIN, file="nuclear_2d_int_ACTIN.csv")
  write.csv(spheroid_nuclei_2dint_data$ACTIN, file="spheroid_median_nuclear_2d_int_ACTIN.csv")
  
  write.csv(spheroid_2dint_data$OCT4, file="spheroid_2d_int_OCT4.csv")
  write.csv(nuclei_2dint_data$OCT4, file="nuclear_2d_int_OCT4.csv")
  write.csv(spheroid_nuclei_2dint_data$OCT4, file="spheroid_median_nuclear_2d_int_OCT4.csv")
  
  write.csv(nuclei_3dint_data$NANOG, file="nuclear_int_NANOG.csv")
  write.csv(spheroid_nuclei_3dint_data$NANOG, file="spheroid_median_nuclear_int_NANOG.csv")
  write.csv(nuclei_3dint_data$OCT4, file="nuclear_int_OCT4.csv")
  write.csv(spheroid_nuclei_3dint_data$OCT4, file="spheroid_median_nuclear_int_OCT4.csv")
  write.csv(nuclei_3dint_data$ACTIN, file="nuclear_int_ACTIN.csv")
  write.csv(spheroid_nuclei_3dint_data$ACTIN, file="spheroid_median_nuclear_int_ACTIN.csv")
  write.csv(nuclei_3dint_data$dna, file="nuclear_int_dna.csv")
  write.csv(spheroid_nuclei_3dint_data$dna, file="spheroid_median_nuclear_int_dna.csv")
  
  write.csv(spheroid_axial$ACTIN, file="spheroid_axial_dist_ACTIN.csv")
  write.csv(spheroid_axial$OCT4, file="spheroid_axial_dist_OCT4.csv")
  write.csv(spheroid_axial$NANOG, file="spheroid_axial_dist_NANOG.csv")
  write.csv(spheroid_axial$dna, file="spheroid_axial_dist_dna.csv")
  
  write.csv(spheroid_angle$ACTIN, file="spheroid_angle_dist_ACTIN.csv")
  write.csv(spheroid_angle$OCT4, file="spheroid_angle_dist_OCT4.csv")
  write.csv(spheroid_angle$NANOG, file="spheroid_angle_dist_NANOG.csv")
  write.csv(spheroid_angle$dna, file="spheroid_angle_dist_dna.csv")
  
  write.csv(spheroid_int_data$NANOG, file="spheroid_int_NANOG.csv")
  write.csv(spheroid_int_data$dna, file="spheroid_int_dna.csv")
  write.csv(spheroid_int_data$OCT4, file="spheroid_int_OCT4.csv")
  write.csv(spheroid_int_data$ACTIN, file="spheroid_int_ACTIN.csv")
  
  
}

