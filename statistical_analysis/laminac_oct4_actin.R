#functions
source("E:/HMF3A_reprogramming/R_analysis/functions.R")

# initialise the experiment details and the subdirectories required
{
  path_to_ex<-c("E:/HMF3A_reprogramming/Timeline_laminac_oct4/20191113_hmf_s7_repg_oct4_488_laminac_568_actin_647/",
                "E:/HMF3A_reprogramming/Timeline_laminac_oct4/20191127_hmf_s8_repg_oct4_488_laminac_568_actin_647/",
                "E:/HMF3A_reprogramming/Timeline_laminac_oct4/20191210_s13_hmf3a_reprogramming_oct4_488_laminAC_568_actin_647/",
                "E:/HMF3A_reprogramming/Timeline_laminac_oct4/20191213_s14_hmf3a_reprogramming_oct4_488_laminAC_568_actin_647/")
  
  exp<-c("20191113_hmf_s7_repg_oct4_488_laminac_568_actin_647","20191127_hmf_s8_repg_oct4_488_laminac_568_actin_647",
         "20191210_s13_hmf3a_reprogramming_oct4_488_laminAC_568_actin_647","20191213_s14_hmf3a_reprogramming_oct4_488_laminAC_568_actin_647")
  
  data_types<-c("/3D geometrical data spheroid/",
                "/3D int_data spheroid/DNA/","/3D int_data spheroid/ACTIN/","/3D int_data spheroid/LAMINAC/",
                "/3D int_data spheroid/OCT4")
  data_types_2d<-c("/2D_measures_spheroid/2D_spheroid_DNA.csv","/2D_measures_spheroid/2D_spheroid_ACTIN.csv",
                   "/2D_measures_spheroid/2D_spheroid_LAMINAC.csv","/2D_measures_spheroid/2D_spheroid_OCT4.csv")
  nuc_data_types<-c("/3D geometrical data/","/3D ellipsoid/","/3D geometerical_simple/","/3D shape measure/",
                    "/3D int_data/DNA/","/3D int_data/ACTIN/","/3D int_data/LAMINAC/","/3D int_data/OCT4/",
                    "/cell_2microns_measure/ACTIN/","/cell_2microns_measure/LAMINAC/","/cell_2microns_measure/OCT4/","/cell_2microns_measure/DNA/")
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
                       "/2D_measures_nuclei/2D_nucleusOCT4.csv","/2D_measures_nuclei/2D_nucleus_ACTIN.csv","/2D_measures_nuclei/2D_nucleus_LAMINAC.csv")
  
  
  samples<-c("D0","D2","D4","D6","D8")
  trial<-c("T1","T2","T3","T4")
  
  dird_apo<-"E:/HMF3A_reprogramming/R_analysis/actin_laminac_oct4/"
  dir.create(dird_apo)
  dird_oct4<-paste(dird_apo,"spheroid_shells_oct4", sep="")
  dird_actin<-paste(dird_apo,"spheroid_shells_actin", sep="")
  dird_LAMINAC<-paste(dird_apo,"spheroid_shells_LAMINAC", sep="")
  dird_dna<-paste(dird_apo,"spheroid_shells_dna", sep="")
  
  dir.create(dird_oct4)
  dir.create(dird_actin)
  dir.create(dird_LAMINAC)
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
    T1LAMINAC_int_data<-combine_sample_sets(path_to_ex[1],samples,data_types[4],"tsv",exp[1],trial[1])
    T1OCT4_int_data<-combine_sample_sets(path_to_ex[1],samples,data_types[5],"tsv",exp[1],trial[1])
    geo_int_2D_data<-combine_sample_sets_2d(path_to_ex[1],samples,data_types_2d[1],exp[1],trial[1])
    
    T1_2d_data_ACTIN<-combine_sample_sets_2d(path_to_ex[1],samples,data_types_2d[2],exp[1],trial[1])
    T1_2d_data_LAMINAC<-combine_sample_sets_2d(path_to_ex[1],samples,data_types_2d[3],exp[1],trial[1])
    T1_2d_data_OCT4<-combine_sample_sets_2d(path_to_ex[1],samples,data_types_2d[4],exp[1],trial[1])
  
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
    T2LAMINAC_int_data<-combine_sample_sets(path_to_ex[2],samples,data_types[4],"tsv",exp[2],trial[2])
    T2OCT4_int_data<-combine_sample_sets(path_to_ex[2],samples,data_types[5],"tsv",exp[2],trial[2])
    geo_int_2D_data<-combine_sample_sets_2d(path_to_ex[2],samples,data_types_2d[1],exp[2],trial[2])
    
    T2_2d_data_ACTIN<-combine_sample_sets_2d(path_to_ex[2],samples,data_types_2d[2],exp[2],trial[2])
    T2_2d_data_LAMINAC<-combine_sample_sets_2d(path_to_ex[2],samples,data_types_2d[3],exp[2],trial[2])
    T2_2d_data_OCT4<-combine_sample_sets_2d(path_to_ex[2],samples,data_types_2d[4],exp[2],trial[2])
    
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
    T3LAMINAC_int_data<-combine_sample_sets(path_to_ex[3],samples,data_types[4],"tsv",exp[3],trial[3])
    T3OCT4_int_data<-combine_sample_sets(path_to_ex[3],samples,data_types[5],"tsv",exp[3],trial[3])
    geo_int_2D_data<-combine_sample_sets_2d(path_to_ex[3],samples,data_types_2d[1],exp[3],trial[3])
    
    T3_2d_data_ACTIN<-combine_sample_sets_2d(path_to_ex[3],samples,data_types_2d[2],exp[3],trial[3])
    T3_2d_data_LAMINAC<-combine_sample_sets_2d(path_to_ex[3],samples,data_types_2d[3],exp[3],trial[3])
    T3_2d_data_OCT4<-combine_sample_sets_2d(path_to_ex[3],samples,data_types_2d[4],exp[3],trial[3])
    
    geo_int_2D_data$Label<-substring(geo_int_2D_data$Label,5, nchar(geo_int_2D_data$Label)-4)
    geometrical_data$Label<-substring(geometrical_data$Label,0, nchar(geometrical_data$Label)-1)
    T3combined<-merge(geometrical_data,geo_int_2D_data, by="Label")
    
    rm(geometrical_data,geo_int_2D_data)
    
  }
  
  #T4
  {
    geometrical_data<-combine_sample_sets(path_to_ex[4],samples,data_types[1],"tsv",exp[4],trial[4])
    T4DNA_int_data<-combine_sample_sets(path_to_ex[4],samples,data_types[2],"tsv",exp[4],trial[4])
    T4ACTIN_int_data<-combine_sample_sets(path_to_ex[4],samples,data_types[3],"tsv",exp[4],trial[4])
    T4LAMINAC_int_data<-combine_sample_sets(path_to_ex[4],samples,data_types[4],"tsv",exp[4],trial[4])
    T4OCT4_int_data<-combine_sample_sets(path_to_ex[4],samples,data_types[5],"tsv",exp[4],trial[4])
    geo_int_2D_data<-combine_sample_sets_2d(path_to_ex[4],samples,data_types_2d[1],exp[4],trial[4])
    
    T4_2d_data_ACTIN<-combine_sample_sets_2d(path_to_ex[4],samples,data_types_2d[2],exp[4],trial[4])
    T4_2d_data_LAMINAC<-combine_sample_sets_2d(path_to_ex[4],samples,data_types_2d[3],exp[4],trial[4])
    T4_2d_data_OCT4<-combine_sample_sets_2d(path_to_ex[4],samples,data_types_2d[4],exp[4],trial[4])
  
    geo_int_2D_data$Label<-substring(geo_int_2D_data$Label,5, nchar(geo_int_2D_data$Label)-4)
    geometrical_data$Label<-substring(geometrical_data$Label,0, nchar(geometrical_data$Label)-1)
    T4combined<-merge(geometrical_data,geo_int_2D_data, by="Label")
    
    rm(geometrical_data,geo_int_2D_data)
    
  }
  spheroid_2dint_data_actin<-rbind(T1_2d_data_ACTIN,T2_2d_data_ACTIN,T3_2d_data_ACTIN,T4_2d_data_ACTIN)
  spheroid_2dint_data_oct4<-rbind(T1_2d_data_OCT4,T2_2d_data_OCT4,T3_2d_data_OCT4,T4_2d_data_OCT4)
  spheroid_2dint_data_LAMINAC<-rbind(T1_2d_data_LAMINAC,T2_2d_data_LAMINAC,T3_2d_data_LAMINAC,T4_2d_data_LAMINAC)
  
  spheroid_int_data_actin<-rbind(T1ACTIN_int_data,T2ACTIN_int_data,T3ACTIN_int_data,T4ACTIN_int_data)
  spheroid_int_data_oct4<-rbind(T1OCT4_int_data,T2OCT4_int_data,T3OCT4_int_data,T4OCT4_int_data)
  spheroid_int_data_LAMINAC<-rbind(T1LAMINAC_int_data,T2LAMINAC_int_data,T3LAMINAC_int_data,T4LAMINAC_int_data)
  spheroid_int_data_dna<-rbind(T1DNA_int_data,T2DNA_int_data,T3DNA_int_data,T4DNA_int_data)
  
  
  spheroid_int_data<-list(actin=spheroid_int_data_actin, oct4=spheroid_int_data_oct4,
                          LAMINAC=spheroid_int_data_LAMINAC,dna=spheroid_int_data_dna)
  spheroid_2dint_data<-list(actin=spheroid_2dint_data_actin, oct4=spheroid_2dint_data_oct4,
                            LAMINAC=spheroid_2dint_data_LAMINAC)
  
  
  spheroid_combined<-rbind(T1combined,T2combined,T3combined,T4combined)
  colnames(spheroid_combined)[86:87]<-c("sample","trial")
  rm(T1combined,T2combined,T3combined,T4combined,
     T1ACTIN_int_data,T2ACTIN_int_data,T3ACTIN_int_data,
     T1OCT4_int_data,T2OCT4_int_data,T3OCT4_int_data,T1LAMINAC_int_data,T2LAMINAC_int_data,T3LAMINAC_int_data,
     T1DNA_int_data,T2DNA_int_data,T3DNA_int_data,
     spheroid_int_data_actin, spheroid_int_data_oct4,spheroid_int_data_LAMINAC,spheroid_int_data_dna,
     T1_2d_data_ACTIN,T2_2d_data_ACTIN,T3_2d_data_ACTIN,T1_2d_data_OCT4,T2_2d_data_OCT4,T3_2d_data_OCT4,
     T1_2d_data_LAMINAC,T2_2d_data_LAMINAC,T3_2d_data_LAMINAC,spheroid_2dint_data_actin, spheroid_2dint_data_oct4,
     spheroid_2dint_data_LAMINAC,T4_2d_data_ACTIN,T4_2d_data_OCT4,T4_2d_data_LAMINAC,T4ACTIN_int_data,T4OCT4_int_data,
     T4LAMINAC_int_data,T4DNA_int_data)
  
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
    T1_LAMINAC_int_data<-combine_sample_sets(path_to_ex[1],samples,nuc_data_types[7],"tsv",exp[1],trial[1])
    T1_OCT4_int_data<-combine_sample_sets(path_to_ex[1],samples,nuc_data_types[8],"tsv",exp[1],trial[1])
    
    T1_2d_ACTIN_int_data<-combine_sample_sets_2d(path_to_ex[1],samples,nuc_data_types_2d[29],exp[1],trial[1])
    T1_2d_LAMINAC_int_data<-combine_sample_sets_2d(path_to_ex[1],samples,nuc_data_types_2d[30],exp[1],trial[1])
    T1_2d_OCT4_int_data<-combine_sample_sets_2d(path_to_ex[1],samples,nuc_data_types_2d[28],exp[1],trial[1])
    
    T1_cell_ACTIN_int_data<-combine_sample_sets(path_to_ex[1],samples,nuc_data_types[9],"tsv",exp[1],trial[1])
    T1_cell_LAMINAC_int_data<-combine_sample_sets(path_to_ex[1],samples,nuc_data_types[10],"tsv",exp[1],trial[1])
    T1_cell_OCT4_int_data<-combine_sample_sets(path_to_ex[1],samples,nuc_data_types[11],"tsv",exp[1],trial[1])
    T1_cell_DNA_int_data<-combine_sample_sets(path_to_ex[1],samples,nuc_data_types[12],"tsv",exp[1],trial[1])
    
    colnames(T1_compact_2D_data)[2]<-"Label"
    T1_compact_2D_data$Label<-substring(T1_compact_2D_data$Label,0, nchar(T1_compact_2D_data$Label)-4)
    T1_geo_int_2D_data$Label<-substring(T1_geo_int_2D_data$Label,5, nchar(T1_geo_int_2D_data$Label)-4)
    T1_geometrical_data$Label<-substring(T1_geometrical_data$Label,0, nchar(T1_geometrical_data$Label)-1)
    T1_DNA_int_data$Label<-substring(T1_DNA_int_data$Label,0, nchar(T1_DNA_int_data$Label)-1)
    T1_2d_LAMINAC_int_data$Label<-substring(T1_2d_LAMINAC_int_data$Label,5, nchar(T1_2d_LAMINAC_int_data$Label)-4)
    T1_2d_ACTIN_int_data$Label<-substring(T1_2d_ACTIN_int_data$Label,5, nchar(T1_2d_ACTIN_int_data$Label)-4)
    T1_2d_OCT4_int_data$Label<-substring(T1_2d_OCT4_int_data$Label,5, nchar(T1_2d_OCT4_int_data$Label)-4)
    T1_ACTIN_int_data$Label<-substring(T1_ACTIN_int_data$Label,0, nchar(T1_ACTIN_int_data$Label)-1)
    T1_LAMINAC_int_data$Label<-substring(T1_LAMINAC_int_data$Label,0, nchar(T1_LAMINAC_int_data$Label)-1)
    T1_OCT4_int_data$Label<-substring(T1_OCT4_int_data$Label,0, nchar(T1_OCT4_int_data$Label)-1)
    T1_cell_ACTIN_int_data$Label<-substring(T1_cell_ACTIN_int_data$Label,0, nchar(T1_cell_ACTIN_int_data$Label)-1)
    T1_cell_LAMINAC_int_data$Label<-substring(T1_cell_LAMINAC_int_data$Label,0, nchar(T1_cell_LAMINAC_int_data$Label)-1)
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
    T2_LAMINAC_int_data<-combine_sample_sets(path_to_ex[2],samples,nuc_data_types[7],"tsv",exp[2],trial[2])
    T2_OCT4_int_data<-combine_sample_sets(path_to_ex[2],samples,nuc_data_types[8],"tsv",exp[2],trial[2])
    
    T2_2d_ACTIN_int_data<-combine_sample_sets_2d(path_to_ex[2],samples,nuc_data_types_2d[29],exp[2],trial[2])
    T2_2d_LAMINAC_int_data<-combine_sample_sets_2d(path_to_ex[2],samples,nuc_data_types_2d[30],exp[2],trial[2])
    T2_2d_OCT4_int_data<-combine_sample_sets_2d(path_to_ex[2],samples,nuc_data_types_2d[28],exp[2],trial[2])
    
    T2_cell_ACTIN_int_data<-combine_sample_sets(path_to_ex[2],samples,nuc_data_types[9],"tsv",exp[2],trial[2])
    T2_cell_LAMINAC_int_data<-combine_sample_sets(path_to_ex[2],samples,nuc_data_types[10],"tsv",exp[2],trial[2])
    T2_cell_OCT4_int_data<-combine_sample_sets(path_to_ex[2],samples,nuc_data_types[11],"tsv",exp[2],trial[2])
    T2_cell_DNA_int_data<-combine_sample_sets(path_to_ex[2],samples,nuc_data_types[12],"tsv",exp[2],trial[2])
    
    colnames(T2_compact_2D_data)[2]<-"Label"
    T2_compact_2D_data$Label<-substring(T2_compact_2D_data$Label,0, nchar(T2_compact_2D_data$Label)-4)
    T2_geo_int_2D_data$Label<-substring(T2_geo_int_2D_data$Label,5, nchar(T2_geo_int_2D_data$Label)-4)
    T2_geometrical_data$Label<-substring(T2_geometrical_data$Label,0, nchar(T2_geometrical_data$Label)-1)
    T2_2d_LAMINAC_int_data$Label<-substring(T2_2d_LAMINAC_int_data$Label,5, nchar(T2_2d_LAMINAC_int_data$Label)-4)
    T2_2d_ACTIN_int_data$Label<-substring(T2_2d_ACTIN_int_data$Label,5, nchar(T2_2d_ACTIN_int_data$Label)-4)
    T2_2d_OCT4_int_data$Label<-substring(T2_2d_OCT4_int_data$Label,5, nchar(T2_2d_OCT4_int_data$Label)-4)
    T2_DNA_int_data$Label<-substring(T2_DNA_int_data$Label,0, nchar(T2_DNA_int_data$Label)-1)
    T2_ACTIN_int_data$Label<-substring(T2_ACTIN_int_data$Label,0, nchar(T2_ACTIN_int_data$Label)-1)
    T2_LAMINAC_int_data$Label<-substring(T2_LAMINAC_int_data$Label,0, nchar(T2_LAMINAC_int_data$Label)-1)
    T2_OCT4_int_data$Label<-substring(T2_OCT4_int_data$Label,0, nchar(T2_OCT4_int_data$Label)-1)
    T2_cell_ACTIN_int_data$Label<-substring(T2_cell_ACTIN_int_data$Label,0, nchar(T2_cell_ACTIN_int_data$Label)-1)
    T2_cell_LAMINAC_int_data$Label<-substring(T2_cell_LAMINAC_int_data$Label,0, nchar(T2_cell_LAMINAC_int_data$Label)-1)
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
    T3_LAMINAC_int_data<-combine_sample_sets(path_to_ex[3],samples,nuc_data_types[7],"tsv",exp[3],trial[3])
    T3_OCT4_int_data<-combine_sample_sets(path_to_ex[3],samples,nuc_data_types[8],"tsv",exp[3],trial[3])
    
    T3_2d_ACTIN_int_data<-combine_sample_sets_2d(path_to_ex[3],samples,nuc_data_types_2d[29],exp[3],trial[3])
    T3_2d_LAMINAC_int_data<-combine_sample_sets_2d(path_to_ex[3],samples,nuc_data_types_2d[30],exp[3],trial[3])
    T3_2d_OCT4_int_data<-combine_sample_sets_2d(path_to_ex[3],samples,nuc_data_types_2d[28],exp[3],trial[3])
    
    T3_cell_ACTIN_int_data<-combine_sample_sets(path_to_ex[3],samples,nuc_data_types[9],"tsv",exp[3],trial[3])
    T3_cell_LAMINAC_int_data<-combine_sample_sets(path_to_ex[3],samples,nuc_data_types[10],"tsv",exp[3],trial[3])
    T3_cell_OCT4_int_data<-combine_sample_sets(path_to_ex[3],samples,nuc_data_types[11],"tsv",exp[3],trial[3])
    T3_cell_DNA_int_data<-combine_sample_sets(path_to_ex[3],samples,nuc_data_types[12],"tsv",exp[3],trial[3])
    
    colnames(T3_compact_2D_data)[2]<-"Label"
    T3_compact_2D_data$Label<-substring(T3_compact_2D_data$Label,0, nchar(T3_compact_2D_data$Label)-4)
    T3_geo_int_2D_data$Label<-substring(T3_geo_int_2D_data$Label,5, nchar(T3_geo_int_2D_data$Label)-4)
    T3_geometrical_data$Label<-substring(T3_geometrical_data$Label,0, nchar(T3_geometrical_data$Label)-1)
    T3_2d_LAMINAC_int_data$Label<-substring(T3_2d_LAMINAC_int_data$Label,5, nchar(T3_2d_LAMINAC_int_data$Label)-4)
    T3_2d_ACTIN_int_data$Label<-substring(T3_2d_ACTIN_int_data$Label,5, nchar(T3_2d_ACTIN_int_data$Label)-4)
    T3_2d_OCT4_int_data$Label<-substring(T3_2d_OCT4_int_data$Label,5, nchar(T3_2d_OCT4_int_data$Label)-4)
    T3_DNA_int_data$Label<-substring(T3_DNA_int_data$Label,0, nchar(T3_DNA_int_data$Label)-1)
    T3_ACTIN_int_data$Label<-substring(T3_ACTIN_int_data$Label,0, nchar(T3_ACTIN_int_data$Label)-1)
    T3_LAMINAC_int_data$Label<-substring(T3_LAMINAC_int_data$Label,0, nchar(T3_LAMINAC_int_data$Label)-1)
    T3_OCT4_int_data$Label<-substring(T3_OCT4_int_data$Label,0, nchar(T3_OCT4_int_data$Label)-1)
    T3_cell_ACTIN_int_data$Label<-substring(T3_cell_ACTIN_int_data$Label,0, nchar(T3_cell_ACTIN_int_data$Label)-1)
    T3_cell_LAMINAC_int_data$Label<-substring(T3_cell_LAMINAC_int_data$Label,0, nchar(T3_cell_LAMINAC_int_data$Label)-1)
    T3_cell_OCT4_int_data$Label<-substring(T3_cell_OCT4_int_data$Label,0, nchar(T3_cell_OCT4_int_data$Label)-1)
    T3_cell_DNA_int_data$Label<-substring(T3_cell_DNA_int_data$Label,0, nchar(T3_cell_DNA_int_data$Label)-1)
    
    T3_combined<-merge(T3_geometrical_data,T3_geo_int_2D_data[,2:36], by="Label")
    T3_combined<-merge(T3_combined,T3_compact_2D_data[,2:8], by="Label")
    T3_combined<-merge(T3_combined,T3_DNA_int_data[,c(3:5,12:18)], by="Label")
    T3_combined<-merge(T3_combined,T3_edf_2D_data[,1:11], by="Label")
    T3_combined<-cbind(T3_combined,T3_gclm[,-c(1,128:131)])
    
    rm(T3_geometrical_data,T3_geo_int_2D_data,T3_compact_2D_data,T3_edf_2D_data,T3_gclm)
  }
  #T4
  {
    T4_geometrical_data<-combine_sample_sets(path_to_ex[4],samples,nuc_data_types[1],"tsv",exp[4],trial[4])
    T4_DNA_int_data<-combine_sample_sets(path_to_ex[4],samples,nuc_data_types[5],"tsv",exp[4],trial[4])
    T4_geo_int_2D_data<-combine_sample_sets_2d(path_to_ex[4],samples,nuc_data_types_2d[1],exp[4],trial[4])
    T4_compact_2D_data<-combine_sample_sets_2d(path_to_ex[4],samples,nuc_data_types_2d[2],exp[4],trial[4])
    T4_edf_2D_data<-combine_sample_sets_2d(path_to_ex[4],samples,nuc_data_types_2d[3],exp[4],trial[4])
    T4_gclm<-read_glcm(path_to_ex[4],samples,nuc_data_types_2d,exp[4],trial[4])
    
    T4_ACTIN_int_data<-combine_sample_sets(path_to_ex[4],samples,nuc_data_types[6],"tsv",exp[4],trial[4])
    T4_LAMINAC_int_data<-combine_sample_sets(path_to_ex[4],samples,nuc_data_types[7],"tsv",exp[4],trial[4])
    T4_OCT4_int_data<-combine_sample_sets(path_to_ex[4],samples,nuc_data_types[8],"tsv",exp[4],trial[4])
    
    T4_2d_ACTIN_int_data<-combine_sample_sets_2d(path_to_ex[4],samples,nuc_data_types_2d[29],exp[4],trial[4])
    T4_2d_LAMINAC_int_data<-combine_sample_sets_2d(path_to_ex[4],samples,nuc_data_types_2d[30],exp[4],trial[4])
    T4_2d_OCT4_int_data<-combine_sample_sets_2d(path_to_ex[4],samples,nuc_data_types_2d[28],exp[4],trial[4])
    
    T4_cell_ACTIN_int_data<-combine_sample_sets(path_to_ex[4],samples,nuc_data_types[9],"tsv",exp[4],trial[4])
    T4_cell_LAMINAC_int_data<-combine_sample_sets(path_to_ex[4],samples,nuc_data_types[10],"tsv",exp[4],trial[4])
    T4_cell_OCT4_int_data<-combine_sample_sets(path_to_ex[4],samples,nuc_data_types[11],"tsv",exp[4],trial[4])
    T4_cell_DNA_int_data<-combine_sample_sets(path_to_ex[4],samples,nuc_data_types[12],"tsv",exp[4],trial[4])
    
    colnames(T4_compact_2D_data)[2]<-"Label"
    T4_compact_2D_data$Label<-substring(T4_compact_2D_data$Label,0, nchar(T4_compact_2D_data$Label)-4)
    T4_geo_int_2D_data$Label<-substring(T4_geo_int_2D_data$Label,5, nchar(T4_geo_int_2D_data$Label)-4)
    T4_geometrical_data$Label<-substring(T4_geometrical_data$Label,0, nchar(T4_geometrical_data$Label)-1)
    T4_2d_LAMINAC_int_data$Label<-substring(T4_2d_LAMINAC_int_data$Label,5, nchar(T4_2d_LAMINAC_int_data$Label)-4)
    T4_2d_ACTIN_int_data$Label<-substring(T4_2d_ACTIN_int_data$Label,5, nchar(T4_2d_ACTIN_int_data$Label)-4)
    T4_2d_OCT4_int_data$Label<-substring(T4_2d_OCT4_int_data$Label,5, nchar(T4_2d_OCT4_int_data$Label)-4)
    T4_DNA_int_data$Label<-substring(T4_DNA_int_data$Label,0, nchar(T4_DNA_int_data$Label)-1)
    T4_ACTIN_int_data$Label<-substring(T4_ACTIN_int_data$Label,0, nchar(T4_ACTIN_int_data$Label)-1)
    T4_LAMINAC_int_data$Label<-substring(T4_LAMINAC_int_data$Label,0, nchar(T4_LAMINAC_int_data$Label)-1)
    T4_OCT4_int_data$Label<-substring(T4_OCT4_int_data$Label,0, nchar(T4_OCT4_int_data$Label)-1)
    T4_cell_ACTIN_int_data$Label<-substring(T4_cell_ACTIN_int_data$Label,0, nchar(T4_cell_ACTIN_int_data$Label)-1)
    T4_cell_LAMINAC_int_data$Label<-substring(T4_cell_LAMINAC_int_data$Label,0, nchar(T4_cell_LAMINAC_int_data$Label)-1)
    T4_cell_OCT4_int_data$Label<-substring(T4_cell_OCT4_int_data$Label,0, nchar(T4_cell_OCT4_int_data$Label)-1)
    T4_cell_DNA_int_data$Label<-substring(T4_cell_DNA_int_data$Label,0, nchar(T4_cell_DNA_int_data$Label)-1)
    
    T4_combined<-merge(T4_geometrical_data,T4_geo_int_2D_data[,2:36], by="Label")
    T4_combined<-merge(T4_combined,T4_compact_2D_data[,2:8], by="Label")
    T4_combined<-merge(T4_combined,T4_DNA_int_data[,c(3:5,12:18)], by="Label")
    T4_combined<-merge(T4_combined,T4_edf_2D_data[,1:11], by="Label")
    T4_combined<-cbind(T4_combined,T4_gclm[,-c(1,128:131)])
    
    rm(T4_geometrical_data,T4_geo_int_2D_data,T4_compact_2D_data,T4_edf_2D_data,T4_gclm)
  }
  nuclei_2dint_data_actin<-rbind(T1_2d_ACTIN_int_data,T2_2d_ACTIN_int_data,T3_2d_ACTIN_int_data,T4_2d_ACTIN_int_data)
  nuclei_2dint_data_oct4<-rbind(T1_2d_OCT4_int_data,T2_2d_OCT4_int_data,T3_2d_OCT4_int_data,T4_2d_OCT4_int_data)
  nuclei_2dint_data_LAMINAC<-rbind(T1_2d_LAMINAC_int_data,T2_2d_LAMINAC_int_data,T3_2d_LAMINAC_int_data,T4_2d_LAMINAC_int_data)
  
  nuclei_int_data_actin<-rbind(T1_ACTIN_int_data,T2_ACTIN_int_data,T3_ACTIN_int_data,T4_ACTIN_int_data)
  nuclei_int_data_oct4<-rbind(T1_OCT4_int_data,T2_OCT4_int_data,T3_OCT4_int_data,T4_OCT4_int_data)
  nuclei_int_data_LAMINAC<-rbind(T1_LAMINAC_int_data,T2_LAMINAC_int_data,T3_LAMINAC_int_data,T4_LAMINAC_int_data)
  nuclei_int_data_dna<-rbind(T1_DNA_int_data,T2_DNA_int_data,T3_DNA_int_data,T4_DNA_int_data)
  
  cell_int_data_actin<-rbind(T1_cell_ACTIN_int_data,T2_cell_ACTIN_int_data,T3_cell_ACTIN_int_data,T4_cell_ACTIN_int_data)
  cell_int_data_oct4<-rbind(T1_cell_OCT4_int_data,T2_cell_OCT4_int_data,T3_cell_OCT4_int_data,T4_cell_OCT4_int_data)
  cell_int_data_LAMINAC<-rbind(T1_cell_LAMINAC_int_data,T2_cell_LAMINAC_int_data,T3_cell_LAMINAC_int_data,T4_cell_LAMINAC_int_data)
  cell_int_data_dna<-rbind(T1_cell_DNA_int_data,T2_cell_DNA_int_data,T3_cell_DNA_int_data,T4_cell_DNA_int_data)
  
  nuclei_2dint_data<-list(actin=nuclei_2dint_data_actin, oct4=nuclei_2dint_data_oct4,
                          LAMINAC=nuclei_2dint_data_LAMINAC)
  nuclei_3dint_data<-list(actin=nuclei_int_data_actin, oct4=nuclei_int_data_oct4,
                          LAMINAC=nuclei_int_data_LAMINAC,dna=nuclei_int_data_dna)
  cell_3dint_data<-list(actin=cell_int_data_actin, oct4=cell_int_data_oct4,
                        LAMINAC=cell_int_data_LAMINAC,dna=cell_int_data_dna)
  
  nucleus_combined<-rbind(T1_combined,T2_combined,T3_combined,T4_combined)
  
  rm(T1_2d_ACTIN_int_data,T2_2d_ACTIN_int_data,T3_2d_ACTIN_int_data,T4_2d_ACTIN_int_data,
     T1_2d_OCT4_int_data,T2_2d_OCT4_int_data,T3_2d_OCT4_int_data,T4_2d_OCT4_int_data,
     T1_2d_LAMINAC_int_data,T2_2d_LAMINAC_int_data,T3_2d_LAMINAC_int_data,T4_2d_LAMINAC_int_data,
     T1_ACTIN_int_data,T2_ACTIN_int_data,T3_ACTIN_int_data,T4_ACTIN_int_data,
     T1_OCT4_int_data,T2_OCT4_int_data,T3_OCT4_int_data,T4_OCT4_int_data,
     T1_LAMINAC_int_data,T2_LAMINAC_int_data,T3_LAMINAC_int_data,T4_LAMINAC_int_data,
     T1_DNA_int_data,T2_DNA_int_data,T3_DNA_int_data,T4_DNA_int_data,
     T1_cell_ACTIN_int_data,T2_cell_ACTIN_int_data,T3_cell_ACTIN_int_data,T4_cell_ACTIN_int_data,
     T1_cell_OCT4_int_data,T2_cell_OCT4_int_data,T3_cell_OCT4_int_data,T4_cell_OCT4_int_data,
     T1_cell_LAMINAC_int_data,T2_cell_LAMINAC_int_data,T3_cell_LAMINAC_int_data,T4_cell_LAMINAC_int_data,
     T1_cell_DNA_int_data,T2_cell_DNA_int_data,T3_cell_DNA_int_data,T4_cell_DNA_int_data,
     T1_combined,T2_combined,T3_combined,T4_combined,
     nuclei_2dint_data_actin, nuclei_2dint_data_oct4,nuclei_2dint_data_LAMINAC,nuclei_int_data_actin, nuclei_int_data_oct4,
     nuclei_int_data_LAMINAC,nuclei_int_data_dna,cell_int_data_actin, cell_int_data_oct4,cell_int_data_LAMINAC,cell_int_data_dna)
  
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
  dird<-"E:/HMF3A_reprogramming/R_analysis/actin_LAMINAC_oct4/black_bg_plots/"
  dir.create(dird)
  setwd(dird)
  cols<-colorRampPalette(c("red","yellow"))( 5) 
  
  #Spheoid Proj_Int
  {
    
    dire<-paste(dird,"spheroid_proj_int/",sep="")
    dir.create(dire)
    setwd(dire)
    
    plot_black_bg_histogram_timeline(subset(spheroid_2dint_data$actin,spheroid_2dint_data$actin$trial=="T1"), cols,parameters_2d_int_data,paste("F-Actin ",labelling_2d_int_data,sep=""),"ACTIN_T1_")
    plot_black_bg_histogram_timeline(subset(spheroid_2dint_data$actin,spheroid_2dint_data$actin$trial=="T2"), cols,parameters_2d_int_data,paste("F-Actin ",labelling_2d_int_data,sep=""),"ACTIN_T2_")
    plot_black_bg_histogram_timeline(subset(spheroid_2dint_data$actin,spheroid_2dint_data$actin$trial=="T3"), cols,parameters_2d_int_data,paste("F-Actin ",labelling_2d_int_data,sep=""),"ACTIN_T3_")
    plot_black_bg_histogram_timeline(subset(spheroid_2dint_data$actin,spheroid_2dint_data$actin$trial=="T4"), cols,parameters_2d_int_data,paste("F-Actin ",labelling_2d_int_data,sep=""),"ACTIN_T4_")
    
    plot_black_bg_boxplot_timeline(subset(spheroid_2dint_data$actin,spheroid_2dint_data$actin$trial=="T1"), cols,parameters_2d_int_data,paste("F-Actin ",labelling_2d_int_data,sep=""),"ACTIN_T1_")
    plot_black_bg_boxplot_timeline(subset(spheroid_2dint_data$actin,spheroid_2dint_data$actin$trial=="T2"), cols,parameters_2d_int_data,paste("F-Actin ",labelling_2d_int_data,sep=""),"ACTIN_T2_")
    plot_black_bg_boxplot_timeline(subset(spheroid_2dint_data$actin,spheroid_2dint_data$actin$trial=="T3"), cols,parameters_2d_int_data,paste("F-Actin ",labelling_2d_int_data,sep=""),"ACTIN_T3_")
    plot_black_bg_boxplot_timeline(subset(spheroid_2dint_data$actin,spheroid_2dint_data$actin$trial=="T4"), cols,parameters_2d_int_data,paste("F-Actin ",labelling_2d_int_data,sep=""),"ACTIN_T4_")
    
    plot_black_bg_barplot_timeline_4_trials(spheroid_2dint_data$actin, cols,parameters_2d_int_data,paste("F-Actin ",labelling_2d_int_data,sep=""),"ACTIN_")
    
    
    plot_black_bg_histogram_timeline(subset(spheroid_2dint_data$LAMINAC,spheroid_2dint_data$LAMINAC$trial=="T1"), cols,parameters_2d_int_data,paste("LAMINAC ",labelling_2d_int_data,sep=""),"LAMINAC_T1_")
    plot_black_bg_histogram_timeline(subset(spheroid_2dint_data$LAMINAC,spheroid_2dint_data$LAMINAC$trial=="T2"), cols,parameters_2d_int_data,paste("LAMINAC ",labelling_2d_int_data,sep=""),"LAMINAC_T2_")
    plot_black_bg_histogram_timeline(subset(spheroid_2dint_data$LAMINAC,spheroid_2dint_data$LAMINAC$trial=="T3"), cols,parameters_2d_int_data,paste("LAMINAC ",labelling_2d_int_data,sep=""),"LAMINAC_T3_")
    plot_black_bg_histogram_timeline(subset(spheroid_2dint_data$LAMINAC,spheroid_2dint_data$LAMINAC$trial=="T4"), cols,parameters_2d_int_data,paste("LAMINAC ",labelling_2d_int_data,sep=""),"LAMINAC_T3_")
    
    plot_black_bg_boxplot_timeline(subset(spheroid_2dint_data$LAMINAC,spheroid_2dint_data$LAMINAC$trial=="T1"), cols,parameters_2d_int_data,paste("LAMINAC ",labelling_2d_int_data,sep=""),"LAMINAC_T1_")
    plot_black_bg_boxplot_timeline(subset(spheroid_2dint_data$LAMINAC,spheroid_2dint_data$LAMINAC$trial=="T2"), cols,parameters_2d_int_data,paste("LAMINAC ",labelling_2d_int_data,sep=""),"LAMINAC_T2_")
    plot_black_bg_boxplot_timeline(subset(spheroid_2dint_data$LAMINAC,spheroid_2dint_data$LAMINAC$trial=="T3"), cols,parameters_2d_int_data,paste("LAMINAC ",labelling_2d_int_data,sep=""),"LAMINAC_T3_")
    plot_black_bg_boxplot_timeline(subset(spheroid_2dint_data$LAMINAC,spheroid_2dint_data$LAMINAC$trial=="T4"), cols,parameters_2d_int_data,paste("LAMINAC ",labelling_2d_int_data,sep=""),"LAMINAC_T4_")
    
    plot_black_bg_barplot_timeline_4_trials(spheroid_2dint_data$LAMINAC, cols,parameters_2d_int_data,paste("LAMINAC ",labelling_2d_int_data,sep=""),"LAMINAC_")
    
    plot_black_bg_histogram_timeline(subset(spheroid_2dint_data$oct4,spheroid_2dint_data$oct4$trial=="T1"), cols,parameters_2d_int_data,paste("Oct4 ",labelling_2d_int_data,sep=""),"Oct4_T1_")
    plot_black_bg_histogram_timeline(subset(spheroid_2dint_data$oct4,spheroid_2dint_data$oct4$trial=="T2"), cols,parameters_2d_int_data,paste("Oct4 ",labelling_2d_int_data,sep=""),"Oct4_T2_")
    plot_black_bg_histogram_timeline(subset(spheroid_2dint_data$oct4,spheroid_2dint_data$oct4$trial=="T3"), cols,parameters_2d_int_data,paste("Oct4 ",labelling_2d_int_data,sep=""),"Oct4_T3_")
    plot_black_bg_histogram_timeline(subset(spheroid_2dint_data$oct4,spheroid_2dint_data$oct4$trial=="T4"), cols,parameters_2d_int_data,paste("Oct4 ",labelling_2d_int_data,sep=""),"Oct4_T4_")
    
    plot_black_bg_boxplot_timeline(subset(spheroid_2dint_data$oct4,spheroid_2dint_data$oct4$trial=="T1"), cols,parameters_2d_int_data,paste("Oct4 ",labelling_2d_int_data,sep=""),"Oct4_T1_")
    plot_black_bg_boxplot_timeline(subset(spheroid_2dint_data$oct4,spheroid_2dint_data$oct4$trial=="T2"), cols,parameters_2d_int_data,paste("Oct4 ",labelling_2d_int_data,sep=""),"Oct4_T2_")
    plot_black_bg_boxplot_timeline(subset(spheroid_2dint_data$oct4,spheroid_2dint_data$oct4$trial=="T3"), cols,parameters_2d_int_data,paste("Oct4 ",labelling_2d_int_data,sep=""),"Oct4_T3_")
    plot_black_bg_boxplot_timeline(subset(spheroid_2dint_data$oct4,spheroid_2dint_data$oct4$trial=="T4"), cols,parameters_2d_int_data,paste("Oct4 ",labelling_2d_int_data,sep=""),"Oct4_T4_")
    
    plot_black_bg_barplot_timeline_4_trials(spheroid_2dint_data$oct4, cols,parameters_2d_int_data,paste("Oct4 ",labelling_2d_int_data,sep=""),"Oct4_")
    
    
  }
  #Spheoid 3D_Int
  {
    dire<-paste(dird,"spheroid_3d_int/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_black_bg_histogram_timeline(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T1_")
    plot_black_bg_histogram_timeline(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T2_")
    plot_black_bg_histogram_timeline(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T3_")
    plot_black_bg_histogram_timeline(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T4"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T4_")
    
    plot_black_bg_boxplot_timeline(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T1_")
    plot_black_bg_boxplot_timeline(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T2_")
    plot_black_bg_boxplot_timeline(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T3_")
    plot_black_bg_boxplot_timeline(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T4"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T4_")
    
    plot_black_bg_barplot_timeline_4_trials(spheroid_int_data$dna, cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_")
    
    plot_black_bg_histogram_timeline(subset(spheroid_int_data$actin,spheroid_int_data$actin$trial=="T1"), cols,parameters_int_data,paste("F-Actin ",labelling_int_data,sep=""),"ACTIN_T1_")
    plot_black_bg_histogram_timeline(subset(spheroid_int_data$actin,spheroid_int_data$actin$trial=="T2"), cols,parameters_int_data,paste("F-Actin ",labelling_int_data,sep=""),"ACTIN_T2_")
    plot_black_bg_histogram_timeline(subset(spheroid_int_data$actin,spheroid_int_data$actin$trial=="T3"), cols,parameters_int_data,paste("F-Actin ",labelling_int_data,sep=""),"ACTIN_T3_")
    plot_black_bg_histogram_timeline(subset(spheroid_int_data$actin,spheroid_int_data$actin$trial=="T4"), cols,parameters_int_data,paste("F-Actin ",labelling_int_data,sep=""),"ACTIN_T4_")
    
    plot_black_bg_boxplot_timeline(subset(spheroid_int_data$actin,spheroid_int_data$actin$trial=="T1"), cols,parameters_int_data,paste("F-Actin ",labelling_int_data,sep=""),"ACTIN_T1_")
    plot_black_bg_boxplot_timeline(subset(spheroid_int_data$actin,spheroid_int_data$actin$trial=="T2"), cols,parameters_int_data,paste("F-Actin ",labelling_int_data,sep=""),"ACTIN_T2_")
    plot_black_bg_boxplot_timeline(subset(spheroid_int_data$actin,spheroid_int_data$actin$trial=="T3"), cols,parameters_int_data,paste("F-Actin ",labelling_int_data,sep=""),"ACTIN_T3_")
    plot_black_bg_boxplot_timeline(subset(spheroid_int_data$actin,spheroid_int_data$actin$trial=="T4"), cols,parameters_int_data,paste("F-Actin ",labelling_int_data,sep=""),"ACTIN_T4_")
    
    plot_black_bg_barplot_timeline_4_trials(spheroid_int_data$actin, cols,parameters_int_data,paste("F-Actin ",labelling_int_data,sep=""),"ACTIN_")
    
    
    plot_black_bg_histogram_timeline(subset(spheroid_int_data$LAMINAC,spheroid_int_data$LAMINAC$trial=="T1"), cols,parameters_int_data,paste("LAMINAC ",labelling_int_data,sep=""),"LAMINAC_T1_")
    plot_black_bg_histogram_timeline(subset(spheroid_int_data$LAMINAC,spheroid_int_data$LAMINAC$trial=="T2"), cols,parameters_int_data,paste("LAMINAC ",labelling_int_data,sep=""),"LAMINAC_T2_")
    plot_black_bg_histogram_timeline(subset(spheroid_int_data$LAMINAC,spheroid_int_data$LAMINAC$trial=="T3"), cols,parameters_int_data,paste("LAMINAC ",labelling_int_data,sep=""),"LAMINAC_T3_")
    plot_black_bg_histogram_timeline(subset(spheroid_int_data$LAMINAC,spheroid_int_data$LAMINAC$trial=="T4"), cols,parameters_int_data,paste("LAMINAC ",labelling_int_data,sep=""),"LAMINAC_T4_")
    
    plot_black_bg_boxplot_timeline(subset(spheroid_int_data$LAMINAC,spheroid_int_data$LAMINAC$trial=="T1"), cols,parameters_int_data,paste("LAMINAC ",labelling_int_data,sep=""),"LAMINAC_T1_")
    plot_black_bg_boxplot_timeline(subset(spheroid_int_data$LAMINAC,spheroid_int_data$LAMINAC$trial=="T2"), cols,parameters_int_data,paste("LAMINAC ",labelling_int_data,sep=""),"LAMINAC_T2_")
    plot_black_bg_boxplot_timeline(subset(spheroid_int_data$LAMINAC,spheroid_int_data$LAMINAC$trial=="T3"), cols,parameters_int_data,paste("LAMINAC ",labelling_int_data,sep=""),"LAMINAC_T3_")
    plot_black_bg_boxplot_timeline(subset(spheroid_int_data$LAMINAC,spheroid_int_data$LAMINAC$trial=="T4"), cols,parameters_int_data,paste("LAMINAC ",labelling_int_data,sep=""),"LAMINAC_T4_")
    
    plot_black_bg_barplot_timeline_4_trials(spheroid_int_data$LAMINAC, cols,parameters_int_data,paste("LAMINAC ",labelling_int_data,sep=""),"LAMINAC_")
    
    plot_black_bg_histogram_timeline(subset(spheroid_int_data$oct4,spheroid_int_data$oct4$trial=="T1"), cols,parameters_int_data,paste("Oct4 ",labelling_int_data,sep=""),"Oct4_T1_")
    plot_black_bg_histogram_timeline(subset(spheroid_int_data$oct4,spheroid_int_data$oct4$trial=="T2"), cols,parameters_int_data,paste("Oct4 ",labelling_int_data,sep=""),"Oct4_T2_")
    plot_black_bg_histogram_timeline(subset(spheroid_int_data$oct4,spheroid_int_data$oct4$trial=="T3"), cols,parameters_int_data,paste("Oct4 ",labelling_int_data,sep=""),"Oct4_T3_")
    plot_black_bg_histogram_timeline(subset(spheroid_int_data$oct4,spheroid_int_data$oct4$trial=="T4"), cols,parameters_int_data,paste("Oct4 ",labelling_int_data,sep=""),"Oct4_T4_")
    
    plot_black_bg_boxplot_timeline(subset(spheroid_int_data$oct4,spheroid_int_data$oct4$trial=="T1"), cols,parameters_int_data,paste("Oct4 ",labelling_int_data,sep=""),"Oct4_T1_")
    plot_black_bg_boxplot_timeline(subset(spheroid_int_data$oct4,spheroid_int_data$oct4$trial=="T2"), cols,parameters_int_data,paste("Oct4 ",labelling_int_data,sep=""),"Oct4_T2_")
    plot_black_bg_boxplot_timeline(subset(spheroid_int_data$oct4,spheroid_int_data$oct4$trial=="T3"), cols,parameters_int_data,paste("Oct4 ",labelling_int_data,sep=""),"Oct4_T3_")
    plot_black_bg_boxplot_timeline(subset(spheroid_int_data$oct4,spheroid_int_data$oct4$trial=="T4"), cols,parameters_int_data,paste("Oct4 ",labelling_int_data,sep=""),"Oct4_T4_")
    
    plot_black_bg_barplot_timeline_4_trials(spheroid_int_data$oct4, cols,parameters_int_data,paste("Oct4 ",labelling_int_data,sep=""),"Oct4_")
    
    
  }

  #Spheoid geometeric properties
  {
    dire<-paste(dird,"spheroid_properties/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_black_bg_histogram_timeline(subset(spheroid_combined,spheroid_combined$trial=="T1"), cols,parameters_geometrical_data,paste("Spheroid ",labelling_geometrical_data,sep=""),"geo_T1_")
    plot_black_bg_histogram_timeline(subset(spheroid_combined,spheroid_combined$trial=="T2"), cols,parameters_geometrical_data,paste("Spheroid ",labelling_geometrical_data,sep=""),"geo_T2_")
    plot_black_bg_histogram_timeline(subset(spheroid_combined,spheroid_combined$trial=="T3"), cols,parameters_geometrical_data,paste("Spheroid ",labelling_geometrical_data,sep=""),"geo_T3_")
    plot_black_bg_histogram_timeline(subset(spheroid_combined,spheroid_combined$trial=="T4"), cols,parameters_geometrical_data,paste("Spheroid ",labelling_geometrical_data,sep=""),"geo_T4_")
    
    plot_black_bg_boxplot_timeline(subset(spheroid_combined,spheroid_combined$trial=="T1"), cols,parameters_geometrical_data,labelling_geometrical_data,"geo_T1_")
    plot_black_bg_boxplot_timeline(subset(spheroid_combined,spheroid_combined$trial=="T2"), cols,parameters_geometrical_data,labelling_geometrical_data,"geo_T2_")
    plot_black_bg_boxplot_timeline(subset(spheroid_combined,spheroid_combined$trial=="T3"), cols,parameters_geometrical_data,labelling_geometrical_data,"geo_T3_")
    plot_black_bg_boxplot_timeline(subset(spheroid_combined,spheroid_combined$trial=="T4"), cols,parameters_geometrical_data,labelling_geometrical_data,"geo_T4_")
    
    plot_black_bg_barplot_timeline_4_trials(spheroid_combined, cols,parameters_geometrical_data,paste("Spheroid ",labelling_geometrical_data,sep=""),"geo_")
    
  }
  #Nuclear geometeric properties
  {
    dire<-paste(dird,"nuclear_properties/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_black_bg_histogram_timeline(subset(nucleus_combined,nucleus_combined$trial=="T1"), cols,parameters_geo_texture_data,paste("Nuclear ",labelling_geo_texture_data,sep=""),"geo_T1_")
    plot_black_bg_histogram_timeline(subset(nucleus_combined,nucleus_combined$trial=="T2"), cols,parameters_geo_texture_data,paste("Nuclear ",labelling_geo_texture_data,sep=""),"geo_T2_")
    plot_black_bg_histogram_timeline(subset(nucleus_combined,nucleus_combined$trial=="T3"), cols,parameters_geo_texture_data,paste("Nuclear ",labelling_geo_texture_data,sep=""),"geo_T3_")
    plot_black_bg_histogram_timeline(subset(nucleus_combined,nucleus_combined$trial=="T4"), cols,parameters_geo_texture_data,paste("Nuclear ",labelling_geo_texture_data,sep=""),"geo_T4_")
    
    plot_black_bg_boxplot_timeline(subset(nucleus_combined,nucleus_combined$trial=="T1"), cols,parameters_geo_texture_data,labelling_geo_texture_data,"geo_T1_")
    plot_black_bg_boxplot_timeline(subset(nucleus_combined,nucleus_combined$trial=="T2"), cols,parameters_geo_texture_data,labelling_geo_texture_data,"geo_T2_")
    plot_black_bg_boxplot_timeline(subset(nucleus_combined,nucleus_combined$trial=="T3"), cols,parameters_geo_texture_data,labelling_geo_texture_data,"geo_T3_")
    plot_black_bg_boxplot_timeline(subset(nucleus_combined,nucleus_combined$trial=="T4"), cols,parameters_geo_texture_data,labelling_geo_texture_data,"geo_T4_")
    
    plot_black_bg_barplot_timeline_4_trials(nucleus_combined, cols,parameters_geo_texture_data,paste("Nuclear ",labelling_geo_texture_data,sep=""),"geo_")
    
  }
  #Nuclear Proj_Int
  {
    
    dire<-paste(dird,"nuclear_proj_int/",sep="")
    dir.create(dire)
    setwd(dire)
    
    plot_black_bg_histogram_timeline(subset(nuclei_2dint_data$actin,nuclei_2dint_data$actin$trial=="T1"), cols,parameters_2d_int_data,paste("F-Actin ",nuc_labelling_2d_int_data,sep=""),"ACTIN_T1_")
    plot_black_bg_histogram_timeline(subset(nuclei_2dint_data$actin,nuclei_2dint_data$actin$trial=="T2"), cols,parameters_2d_int_data,paste("F-Actin ",nuc_labelling_2d_int_data,sep=""),"ACTIN_T2_")
    plot_black_bg_histogram_timeline(subset(nuclei_2dint_data$actin,nuclei_2dint_data$actin$trial=="T3"), cols,parameters_2d_int_data,paste("F-Actin ",nuc_labelling_2d_int_data,sep=""),"ACTIN_T3_")
    plot_black_bg_histogram_timeline(subset(nuclei_2dint_data$actin,nuclei_2dint_data$actin$trial=="T4"), cols,parameters_2d_int_data,paste("F-Actin ",nuc_labelling_2d_int_data,sep=""),"ACTIN_T4_")
    
    plot_black_bg_boxplot_timeline(subset(nuclei_2dint_data$actin,nuclei_2dint_data$actin$trial=="T1"), cols,parameters_2d_int_data,paste("F-Actin ",nuc_labelling_2d_int_data,sep=""),"ACTIN_T1_")
    plot_black_bg_boxplot_timeline(subset(nuclei_2dint_data$actin,nuclei_2dint_data$actin$trial=="T2"), cols,parameters_2d_int_data,paste("F-Actin ",nuc_labelling_2d_int_data,sep=""),"ACTIN_T2_")
    plot_black_bg_boxplot_timeline(subset(nuclei_2dint_data$actin,nuclei_2dint_data$actin$trial=="T3"), cols,parameters_2d_int_data,paste("F-Actin ",nuc_labelling_2d_int_data,sep=""),"ACTIN_T3_")
    plot_black_bg_boxplot_timeline(subset(nuclei_2dint_data$actin,nuclei_2dint_data$actin$trial=="T4"), cols,parameters_2d_int_data,paste("F-Actin ",nuc_labelling_2d_int_data,sep=""),"ACTIN_T4_")
    
    plot_black_bg_barplot_timeline_4_trials(nuclei_2dint_data$actin, cols,parameters_2d_int_data,paste("F-Actin ",nuc_labelling_2d_int_data,sep=""),"ACTIN_")
    
    
    plot_black_bg_histogram_timeline(subset(nuclei_2dint_data$LAMINAC,nuclei_2dint_data$LAMINAC$trial=="T1"), cols,parameters_2d_int_data,paste("LAMINAC ",nuc_labelling_2d_int_data,sep=""),"LAMINAC_T1_")
    plot_black_bg_histogram_timeline(subset(nuclei_2dint_data$LAMINAC,nuclei_2dint_data$LAMINAC$trial=="T2"), cols,parameters_2d_int_data,paste("LAMINAC ",nuc_labelling_2d_int_data,sep=""),"LAMINAC_T2_")
    plot_black_bg_histogram_timeline(subset(nuclei_2dint_data$LAMINAC,nuclei_2dint_data$LAMINAC$trial=="T3"), cols,parameters_2d_int_data,paste("LAMINAC ",nuc_labelling_2d_int_data,sep=""),"LAMINAC_T3_")
    plot_black_bg_histogram_timeline(subset(nuclei_2dint_data$LAMINAC,nuclei_2dint_data$LAMINAC$trial=="T4"), cols,parameters_2d_int_data,paste("LAMINAC ",nuc_labelling_2d_int_data,sep=""),"LAMINAC_T4_")
    
    plot_black_bg_boxplot_timeline(subset(nuclei_2dint_data$LAMINAC,nuclei_2dint_data$LAMINAC$trial=="T1"), cols,parameters_2d_int_data,paste("LAMINAC ",nuc_labelling_2d_int_data,sep=""),"LAMINAC_T1_")
    plot_black_bg_boxplot_timeline(subset(nuclei_2dint_data$LAMINAC,nuclei_2dint_data$LAMINAC$trial=="T2"), cols,parameters_2d_int_data,paste("LAMINAC ",nuc_labelling_2d_int_data,sep=""),"LAMINAC_T2_")
    plot_black_bg_boxplot_timeline(subset(nuclei_2dint_data$LAMINAC,nuclei_2dint_data$LAMINAC$trial=="T3"), cols,parameters_2d_int_data,paste("LAMINAC ",nuc_labelling_2d_int_data,sep=""),"LAMINAC_T3_")
    plot_black_bg_boxplot_timeline(subset(nuclei_2dint_data$LAMINAC,nuclei_2dint_data$LAMINAC$trial=="T4"), cols,parameters_2d_int_data,paste("LAMINAC ",nuc_labelling_2d_int_data,sep=""),"LAMINAC_T4_")
    
    plot_black_bg_barplot_timeline_4_trials(nuclei_2dint_data$LAMINAC, cols,parameters_2d_int_data,paste("LAMINAC ",nuc_labelling_2d_int_data,sep=""),"LAMINAC_")
    
    plot_black_bg_histogram_timeline(subset(nuclei_2dint_data$oct4,nuclei_2dint_data$oct4$trial=="T1"), cols,parameters_2d_int_data,paste("Oct4 ",nuc_labelling_2d_int_data,sep=""),"Oct4_T1_")
    plot_black_bg_histogram_timeline(subset(nuclei_2dint_data$oct4,nuclei_2dint_data$oct4$trial=="T2"), cols,parameters_2d_int_data,paste("Oct4 ",nuc_labelling_2d_int_data,sep=""),"Oct4_T2_")
    plot_black_bg_histogram_timeline(subset(nuclei_2dint_data$oct4,nuclei_2dint_data$oct4$trial=="T3"), cols,parameters_2d_int_data,paste("Oct4 ",nuc_labelling_2d_int_data,sep=""),"Oct4_T3_")
    plot_black_bg_histogram_timeline(subset(nuclei_2dint_data$oct4,nuclei_2dint_data$oct4$trial=="T4"), cols,parameters_2d_int_data,paste("Oct4 ",nuc_labelling_2d_int_data,sep=""),"Oct4_T4_")
    
    plot_black_bg_boxplot_timeline(subset(nuclei_2dint_data$oct4,nuclei_2dint_data$oct4$trial=="T1"), cols,parameters_2d_int_data,paste("Oct4 ",nuc_labelling_2d_int_data,sep=""),"Oct4_T1_")
    plot_black_bg_boxplot_timeline(subset(nuclei_2dint_data$oct4,nuclei_2dint_data$oct4$trial=="T2"), cols,parameters_2d_int_data,paste("Oct4 ",nuc_labelling_2d_int_data,sep=""),"Oct4_T2_")
    plot_black_bg_boxplot_timeline(subset(nuclei_2dint_data$oct4,nuclei_2dint_data$oct4$trial=="T3"), cols,parameters_2d_int_data,paste("Oct4 ",nuc_labelling_2d_int_data,sep=""),"Oct4_T3_")
    plot_black_bg_boxplot_timeline(subset(nuclei_2dint_data$oct4,nuclei_2dint_data$oct4$trial=="T4"), cols,parameters_2d_int_data,paste("Oct4 ",nuc_labelling_2d_int_data,sep=""),"Oct4_T4_")
    
    plot_black_bg_barplot_timeline_4_trials(nuclei_2dint_data$oct4, cols,parameters_2d_int_data,paste("Oct4 ",nuc_labelling_2d_int_data,sep=""),"Oct4_")
    
    
  }
  #Nuclear 3D_Int
  {
    dire<-paste(dird,"nuclear_3d_int/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_black_bg_histogram_timeline(subset(nuclei_3dint_data$dna,nuclei_3dint_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T1_")
    plot_black_bg_histogram_timeline(subset(nuclei_3dint_data$dna,nuclei_3dint_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T2_")
    plot_black_bg_histogram_timeline(subset(nuclei_3dint_data$dna,nuclei_3dint_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T3_")
    plot_black_bg_histogram_timeline(subset(nuclei_3dint_data$dna,nuclei_3dint_data$dna$trial=="T4"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T4_")
    
    plot_black_bg_boxplot_timeline(subset(nuclei_3dint_data$dna,nuclei_3dint_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T1_")
    plot_black_bg_boxplot_timeline(subset(nuclei_3dint_data$dna,nuclei_3dint_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T2_")
    plot_black_bg_boxplot_timeline(subset(nuclei_3dint_data$dna,nuclei_3dint_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T3_")
    plot_black_bg_boxplot_timeline(subset(nuclei_3dint_data$dna,nuclei_3dint_data$dna$trial=="T4"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T4_")
    
    plot_black_bg_barplot_timeline_4_trials(nuclei_3dint_data$dna, cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_")
    
    plot_black_bg_histogram_timeline(subset(nuclei_3dint_data$actin,nuclei_3dint_data$actin$trial=="T1"), cols,parameters_int_data,paste("F-Actin ",nuclabelling_int_data,sep=""),"ACTIN_T1_")
    plot_black_bg_histogram_timeline(subset(nuclei_3dint_data$actin,nuclei_3dint_data$actin$trial=="T2"), cols,parameters_int_data,paste("F-Actin ",nuclabelling_int_data,sep=""),"ACTIN_T2_")
    plot_black_bg_histogram_timeline(subset(nuclei_3dint_data$actin,nuclei_3dint_data$actin$trial=="T3"), cols,parameters_int_data,paste("F-Actin ",nuclabelling_int_data,sep=""),"ACTIN_T3_")
    plot_black_bg_histogram_timeline(subset(nuclei_3dint_data$actin,nuclei_3dint_data$actin$trial=="T4"), cols,parameters_int_data,paste("F-Actin ",nuclabelling_int_data,sep=""),"ACTIN_T4_")
    
    plot_black_bg_boxplot_timeline(subset(nuclei_3dint_data$actin,nuclei_3dint_data$actin$trial=="T1"), cols,parameters_int_data,paste("F-Actin ",nuclabelling_int_data,sep=""),"ACTIN_T1_")
    plot_black_bg_boxplot_timeline(subset(nuclei_3dint_data$actin,nuclei_3dint_data$actin$trial=="T2"), cols,parameters_int_data,paste("F-Actin ",nuclabelling_int_data,sep=""),"ACTIN_T2_")
    plot_black_bg_boxplot_timeline(subset(nuclei_3dint_data$actin,nuclei_3dint_data$actin$trial=="T3"), cols,parameters_int_data,paste("F-Actin ",nuclabelling_int_data,sep=""),"ACTIN_T3_")
    plot_black_bg_boxplot_timeline(subset(nuclei_3dint_data$actin,nuclei_3dint_data$actin$trial=="T4"), cols,parameters_int_data,paste("F-Actin ",nuclabelling_int_data,sep=""),"ACTIN_T4_")
    
    plot_black_bg_barplot_timeline_4_trials(nuclei_3dint_data$actin, cols,parameters_int_data,paste("F-Actin ",nuclabelling_int_data,sep=""),"ACTIN_")
    
    
    plot_black_bg_histogram_timeline(subset(nuclei_3dint_data$LAMINAC,nuclei_3dint_data$LAMINAC$trial=="T1"), cols,parameters_int_data,paste("LAMINAC ",nuclabelling_int_data,sep=""),"LAMINAC_T1_")
    plot_black_bg_histogram_timeline(subset(nuclei_3dint_data$LAMINAC,nuclei_3dint_data$LAMINAC$trial=="T2"), cols,parameters_int_data,paste("LAMINAC ",nuclabelling_int_data,sep=""),"LAMINAC_T2_")
    plot_black_bg_histogram_timeline(subset(nuclei_3dint_data$LAMINAC,nuclei_3dint_data$LAMINAC$trial=="T3"), cols,parameters_int_data,paste("LAMINAC ",nuclabelling_int_data,sep=""),"LAMINAC_T3_")
    plot_black_bg_histogram_timeline(subset(nuclei_3dint_data$LAMINAC,nuclei_3dint_data$LAMINAC$trial=="T4"), cols,parameters_int_data,paste("LAMINAC ",nuclabelling_int_data,sep=""),"LAMINAC_T4_")
    
    plot_black_bg_boxplot_timeline(subset(nuclei_3dint_data$LAMINAC,nuclei_3dint_data$LAMINAC$trial=="T1"), cols,parameters_int_data,paste("LAMINAC ",nuclabelling_int_data,sep=""),"LAMINAC_T1_")
    plot_black_bg_boxplot_timeline(subset(nuclei_3dint_data$LAMINAC,nuclei_3dint_data$LAMINAC$trial=="T2"), cols,parameters_int_data,paste("LAMINAC ",nuclabelling_int_data,sep=""),"LAMINAC_T2_")
    plot_black_bg_boxplot_timeline(subset(nuclei_3dint_data$LAMINAC,nuclei_3dint_data$LAMINAC$trial=="T3"), cols,parameters_int_data,paste("LAMINAC ",nuclabelling_int_data,sep=""),"LAMINAC_T3_")
    plot_black_bg_boxplot_timeline(subset(nuclei_3dint_data$LAMINAC,nuclei_3dint_data$LAMINAC$trial=="T4"), cols,parameters_int_data,paste("LAMINAC ",nuclabelling_int_data,sep=""),"LAMINAC_T4_")
    
    plot_black_bg_barplot_timeline_4_trials(nuclei_3dint_data$LAMINAC, cols,parameters_int_data,paste("LAMINAC ",nuclabelling_int_data,sep=""),"LAMINAC_")
    
    plot_black_bg_histogram_timeline(subset(nuclei_3dint_data$oct4,nuclei_3dint_data$oct4$trial=="T1"), cols,parameters_int_data,paste("Oct4 ",nuclabelling_int_data,sep=""),"Oct4_T1_")
    plot_black_bg_histogram_timeline(subset(nuclei_3dint_data$oct4,nuclei_3dint_data$oct4$trial=="T2"), cols,parameters_int_data,paste("Oct4 ",nuclabelling_int_data,sep=""),"Oct4_T2_")
    plot_black_bg_histogram_timeline(subset(nuclei_3dint_data$oct4,nuclei_3dint_data$oct4$trial=="T3"), cols,parameters_int_data,paste("Oct4 ",nuclabelling_int_data,sep=""),"Oct4_T3_")
    plot_black_bg_histogram_timeline(subset(nuclei_3dint_data$oct4,nuclei_3dint_data$oct4$trial=="T4"), cols,parameters_int_data,paste("Oct4 ",nuclabelling_int_data,sep=""),"Oct4_T4_")
    
    plot_black_bg_boxplot_timeline(subset(nuclei_3dint_data$oct4,nuclei_3dint_data$oct4$trial=="T1"), cols,parameters_int_data,paste("Oct4 ",nuclabelling_int_data,sep=""),"Oct4_T1_")
    plot_black_bg_boxplot_timeline(subset(nuclei_3dint_data$oct4,nuclei_3dint_data$oct4$trial=="T2"), cols,parameters_int_data,paste("Oct4 ",nuclabelling_int_data,sep=""),"Oct4_T2_")
    plot_black_bg_boxplot_timeline(subset(nuclei_3dint_data$oct4,nuclei_3dint_data$oct4$trial=="T3"), cols,parameters_int_data,paste("Oct4 ",nuclabelling_int_data,sep=""),"Oct4_T3_")
    plot_black_bg_boxplot_timeline(subset(nuclei_3dint_data$oct4,nuclei_3dint_data$oct4$trial=="T4"), cols,parameters_int_data,paste("Oct4 ",nuclabelling_int_data,sep=""),"Oct4_T4_")
    
    plot_black_bg_barplot_timeline_4_trials(nuclei_3dint_data$oct4, cols,parameters_int_data,paste("Oct4 ",nuclabelling_int_data,sep=""),"Oct4_")
    
    
  }
  #Cellular 3D_Int
  {
    dire<-paste(dird,"cellular_3d_int/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_black_bg_histogram_timeline(subset(cell_3dint_data$dna,cell_3dint_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T1_")
    plot_black_bg_histogram_timeline(subset(cell_3dint_data$dna,cell_3dint_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T2_")
    plot_black_bg_histogram_timeline(subset(cell_3dint_data$dna,cell_3dint_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T3_")
    plot_black_bg_histogram_timeline(subset(cell_3dint_data$dna,cell_3dint_data$dna$trial=="T4"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T4_")
    
    plot_black_bg_boxplot_timeline(subset(cell_3dint_data$dna,cell_3dint_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T1_")
    plot_black_bg_boxplot_timeline(subset(cell_3dint_data$dna,cell_3dint_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T2_")
    plot_black_bg_boxplot_timeline(subset(cell_3dint_data$dna,cell_3dint_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T3_")
    plot_black_bg_boxplot_timeline(subset(cell_3dint_data$dna,cell_3dint_data$dna$trial=="T4"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T4_")
    
    plot_black_bg_barplot_timeline_4_trials(cell_3dint_data$dna, cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_")
    
    plot_black_bg_histogram_timeline(subset(cell_3dint_data$actin,cell_3dint_data$actin$trial=="T1"), cols,parameters_int_data,paste("F-Actin ",celllabelling_int_data,sep=""),"ACTIN_T1_")
    plot_black_bg_histogram_timeline(subset(cell_3dint_data$actin,cell_3dint_data$actin$trial=="T2"), cols,parameters_int_data,paste("F-Actin ",celllabelling_int_data,sep=""),"ACTIN_T2_")
    plot_black_bg_histogram_timeline(subset(cell_3dint_data$actin,cell_3dint_data$actin$trial=="T3"), cols,parameters_int_data,paste("F-Actin ",celllabelling_int_data,sep=""),"ACTIN_T3_")
    plot_black_bg_histogram_timeline(subset(cell_3dint_data$actin,cell_3dint_data$actin$trial=="T4"), cols,parameters_int_data,paste("F-Actin ",celllabelling_int_data,sep=""),"ACTIN_T4_")
    
    plot_black_bg_boxplot_timeline(subset(cell_3dint_data$actin,cell_3dint_data$actin$trial=="T1"), cols,parameters_int_data,paste("F-Actin ",celllabelling_int_data,sep=""),"ACTIN_T1_")
    plot_black_bg_boxplot_timeline(subset(cell_3dint_data$actin,cell_3dint_data$actin$trial=="T2"), cols,parameters_int_data,paste("F-Actin ",celllabelling_int_data,sep=""),"ACTIN_T2_")
    plot_black_bg_boxplot_timeline(subset(cell_3dint_data$actin,cell_3dint_data$actin$trial=="T3"), cols,parameters_int_data,paste("F-Actin ",celllabelling_int_data,sep=""),"ACTIN_T3_")
    plot_black_bg_boxplot_timeline(subset(cell_3dint_data$actin,cell_3dint_data$actin$trial=="T4"), cols,parameters_int_data,paste("F-Actin ",celllabelling_int_data,sep=""),"ACTIN_T4_")
    
    plot_black_bg_barplot_timeline_4_trials(cell_3dint_data$actin, cols,parameters_int_data,paste("F-Actin ",celllabelling_int_data,sep=""),"ACTIN_")
    
    
    plot_black_bg_histogram_timeline(subset(cell_3dint_data$LAMINAC,cell_3dint_data$LAMINAC$trial=="T1"), cols,parameters_int_data,paste("LAMINAC ",celllabelling_int_data,sep=""),"LAMINAC_T1_")
    plot_black_bg_histogram_timeline(subset(cell_3dint_data$LAMINAC,cell_3dint_data$LAMINAC$trial=="T2"), cols,parameters_int_data,paste("LAMINAC ",celllabelling_int_data,sep=""),"LAMINAC_T2_")
    plot_black_bg_histogram_timeline(subset(cell_3dint_data$LAMINAC,cell_3dint_data$LAMINAC$trial=="T3"), cols,parameters_int_data,paste("LAMINAC ",celllabelling_int_data,sep=""),"LAMINAC_T3_")
    plot_black_bg_histogram_timeline(subset(cell_3dint_data$LAMINAC,cell_3dint_data$LAMINAC$trial=="T4"), cols,parameters_int_data,paste("LAMINAC ",celllabelling_int_data,sep=""),"LAMINAC_T4_")
    
    plot_black_bg_boxplot_timeline(subset(cell_3dint_data$LAMINAC,cell_3dint_data$LAMINAC$trial=="T1"), cols,parameters_int_data,paste("LAMINAC ",celllabelling_int_data,sep=""),"LAMINAC_T1_")
    plot_black_bg_boxplot_timeline(subset(cell_3dint_data$LAMINAC,cell_3dint_data$LAMINAC$trial=="T2"), cols,parameters_int_data,paste("LAMINAC ",celllabelling_int_data,sep=""),"LAMINAC_T2_")
    plot_black_bg_boxplot_timeline(subset(cell_3dint_data$LAMINAC,cell_3dint_data$LAMINAC$trial=="T3"), cols,parameters_int_data,paste("LAMINAC ",celllabelling_int_data,sep=""),"LAMINAC_T3_")
    plot_black_bg_boxplot_timeline(subset(cell_3dint_data$LAMINAC,cell_3dint_data$LAMINAC$trial=="T4"), cols,parameters_int_data,paste("LAMINAC ",celllabelling_int_data,sep=""),"LAMINAC_T4_")
    
    plot_black_bg_barplot_timeline_4_trials(cell_3dint_data$LAMINAC, cols,parameters_int_data,paste("LAMINAC ",celllabelling_int_data,sep=""),"LAMINAC_")
    
    plot_black_bg_histogram_timeline(subset(cell_3dint_data$oct4,cell_3dint_data$oct4$trial=="T1"), cols,parameters_int_data,paste("Oct4 ",celllabelling_int_data,sep=""),"Oct4_T1_")
    plot_black_bg_histogram_timeline(subset(cell_3dint_data$oct4,cell_3dint_data$oct4$trial=="T2"), cols,parameters_int_data,paste("Oct4 ",celllabelling_int_data,sep=""),"Oct4_T2_")
    plot_black_bg_histogram_timeline(subset(cell_3dint_data$oct4,cell_3dint_data$oct4$trial=="T3"), cols,parameters_int_data,paste("Oct4 ",celllabelling_int_data,sep=""),"Oct4_T3_")
    plot_black_bg_histogram_timeline(subset(cell_3dint_data$oct4,cell_3dint_data$oct4$trial=="T4"), cols,parameters_int_data,paste("Oct4 ",celllabelling_int_data,sep=""),"Oct4_T4_")
    
    plot_black_bg_boxplot_timeline(subset(cell_3dint_data$oct4,cell_3dint_data$oct4$trial=="T1"), cols,parameters_int_data,paste("Oct4 ",celllabelling_int_data,sep=""),"Oct4_T1_")
    plot_black_bg_boxplot_timeline(subset(cell_3dint_data$oct4,cell_3dint_data$oct4$trial=="T2"), cols,parameters_int_data,paste("Oct4 ",celllabelling_int_data,sep=""),"Oct4_T2_")
    plot_black_bg_boxplot_timeline(subset(cell_3dint_data$oct4,cell_3dint_data$oct4$trial=="T3"), cols,parameters_int_data,paste("Oct4 ",celllabelling_int_data,sep=""),"Oct4_T3_")
    plot_black_bg_boxplot_timeline(subset(cell_3dint_data$oct4,cell_3dint_data$oct4$trial=="T4"), cols,parameters_int_data,paste("Oct4 ",celllabelling_int_data,sep=""),"Oct4_T4_")
    
    plot_black_bg_barplot_timeline_4_trials(cell_3dint_data$oct4, cols,parameters_int_data,paste("Oct4 ",celllabelling_int_data,sep=""),"Oct4_")
    
    
  }
  
  dird<-"E:/HMF3A_reprogramming/R_analysis/actin_LAMINAC_oct4/white_bg_plots/"
  dir.create(dird)
  setwd(dird)
  cols<-colorRampPalette(c("red","yellow"))( 5) 
  
  #Spheoid Proj_Int
  {
    
    dire<-paste(dird,"spheroid_proj_int/",sep="")
    dir.create(dire)
    setwd(dire)
    
    plot_white_bg_histogram_timeline(subset(spheroid_2dint_data$actin,spheroid_2dint_data$actin$trial=="T1"), cols,parameters_2d_int_data,paste("F-Actin ",labelling_2d_int_data,sep=""),"ACTIN_T1_")
    plot_white_bg_histogram_timeline(subset(spheroid_2dint_data$actin,spheroid_2dint_data$actin$trial=="T2"), cols,parameters_2d_int_data,paste("F-Actin ",labelling_2d_int_data,sep=""),"ACTIN_T2_")
    plot_white_bg_histogram_timeline(subset(spheroid_2dint_data$actin,spheroid_2dint_data$actin$trial=="T3"), cols,parameters_2d_int_data,paste("F-Actin ",labelling_2d_int_data,sep=""),"ACTIN_T3_")
    plot_white_bg_histogram_timeline(subset(spheroid_2dint_data$actin,spheroid_2dint_data$actin$trial=="T4"), cols,parameters_2d_int_data,paste("F-Actin ",labelling_2d_int_data,sep=""),"ACTIN_T4_")
    
    plot_white_bg_boxplot_timeline(subset(spheroid_2dint_data$actin,spheroid_2dint_data$actin$trial=="T1"), cols,parameters_2d_int_data,paste("F-Actin ",labelling_2d_int_data,sep=""),"ACTIN_T1_")
    plot_white_bg_boxplot_timeline(subset(spheroid_2dint_data$actin,spheroid_2dint_data$actin$trial=="T2"), cols,parameters_2d_int_data,paste("F-Actin ",labelling_2d_int_data,sep=""),"ACTIN_T2_")
    plot_white_bg_boxplot_timeline(subset(spheroid_2dint_data$actin,spheroid_2dint_data$actin$trial=="T3"), cols,parameters_2d_int_data,paste("F-Actin ",labelling_2d_int_data,sep=""),"ACTIN_T3_")
    plot_white_bg_boxplot_timeline(subset(spheroid_2dint_data$actin,spheroid_2dint_data$actin$trial=="T4"), cols,parameters_2d_int_data,paste("F-Actin ",labelling_2d_int_data,sep=""),"ACTIN_T4_")
    
    plot_white_bg_barplot_timeline_4_trials(spheroid_2dint_data$actin, cols,parameters_2d_int_data,paste("F-Actin ",labelling_2d_int_data,sep=""),"ACTIN_")
    
    
    plot_white_bg_histogram_timeline(subset(spheroid_2dint_data$LAMINAC,spheroid_2dint_data$LAMINAC$trial=="T1"), cols,parameters_2d_int_data,paste("LAMINAC ",labelling_2d_int_data,sep=""),"LAMINAC_T1_")
    plot_white_bg_histogram_timeline(subset(spheroid_2dint_data$LAMINAC,spheroid_2dint_data$LAMINAC$trial=="T2"), cols,parameters_2d_int_data,paste("LAMINAC ",labelling_2d_int_data,sep=""),"LAMINAC_T2_")
    plot_white_bg_histogram_timeline(subset(spheroid_2dint_data$LAMINAC,spheroid_2dint_data$LAMINAC$trial=="T3"), cols,parameters_2d_int_data,paste("LAMINAC ",labelling_2d_int_data,sep=""),"LAMINAC_T3_")
    plot_white_bg_histogram_timeline(subset(spheroid_2dint_data$LAMINAC,spheroid_2dint_data$LAMINAC$trial=="T4"), cols,parameters_2d_int_data,paste("LAMINAC ",labelling_2d_int_data,sep=""),"LAMINAC_T3_")
    
    plot_white_bg_boxplot_timeline(subset(spheroid_2dint_data$LAMINAC,spheroid_2dint_data$LAMINAC$trial=="T1"), cols,parameters_2d_int_data,paste("LAMINAC ",labelling_2d_int_data,sep=""),"LAMINAC_T1_")
    plot_white_bg_boxplot_timeline(subset(spheroid_2dint_data$LAMINAC,spheroid_2dint_data$LAMINAC$trial=="T2"), cols,parameters_2d_int_data,paste("LAMINAC ",labelling_2d_int_data,sep=""),"LAMINAC_T2_")
    plot_white_bg_boxplot_timeline(subset(spheroid_2dint_data$LAMINAC,spheroid_2dint_data$LAMINAC$trial=="T3"), cols,parameters_2d_int_data,paste("LAMINAC ",labelling_2d_int_data,sep=""),"LAMINAC_T3_")
    plot_white_bg_boxplot_timeline(subset(spheroid_2dint_data$LAMINAC,spheroid_2dint_data$LAMINAC$trial=="T4"), cols,parameters_2d_int_data,paste("LAMINAC ",labelling_2d_int_data,sep=""),"LAMINAC_T4_")
    
    plot_white_bg_barplot_timeline_4_trials(spheroid_2dint_data$LAMINAC, cols,parameters_2d_int_data,paste("LAMINAC ",labelling_2d_int_data,sep=""),"LAMINAC_")
    
    plot_white_bg_histogram_timeline(subset(spheroid_2dint_data$oct4,spheroid_2dint_data$oct4$trial=="T1"), cols,parameters_2d_int_data,paste("Oct4 ",labelling_2d_int_data,sep=""),"Oct4_T1_")
    plot_white_bg_histogram_timeline(subset(spheroid_2dint_data$oct4,spheroid_2dint_data$oct4$trial=="T2"), cols,parameters_2d_int_data,paste("Oct4 ",labelling_2d_int_data,sep=""),"Oct4_T2_")
    plot_white_bg_histogram_timeline(subset(spheroid_2dint_data$oct4,spheroid_2dint_data$oct4$trial=="T3"), cols,parameters_2d_int_data,paste("Oct4 ",labelling_2d_int_data,sep=""),"Oct4_T3_")
    plot_white_bg_histogram_timeline(subset(spheroid_2dint_data$oct4,spheroid_2dint_data$oct4$trial=="T4"), cols,parameters_2d_int_data,paste("Oct4 ",labelling_2d_int_data,sep=""),"Oct4_T4_")
    
    plot_white_bg_boxplot_timeline(subset(spheroid_2dint_data$oct4,spheroid_2dint_data$oct4$trial=="T1"), cols,parameters_2d_int_data,paste("Oct4 ",labelling_2d_int_data,sep=""),"Oct4_T1_")
    plot_white_bg_boxplot_timeline(subset(spheroid_2dint_data$oct4,spheroid_2dint_data$oct4$trial=="T2"), cols,parameters_2d_int_data,paste("Oct4 ",labelling_2d_int_data,sep=""),"Oct4_T2_")
    plot_white_bg_boxplot_timeline(subset(spheroid_2dint_data$oct4,spheroid_2dint_data$oct4$trial=="T3"), cols,parameters_2d_int_data,paste("Oct4 ",labelling_2d_int_data,sep=""),"Oct4_T3_")
    plot_white_bg_boxplot_timeline(subset(spheroid_2dint_data$oct4,spheroid_2dint_data$oct4$trial=="T4"), cols,parameters_2d_int_data,paste("Oct4 ",labelling_2d_int_data,sep=""),"Oct4_T4_")
    
    plot_white_bg_barplot_timeline_4_trials(spheroid_2dint_data$oct4, cols,parameters_2d_int_data,paste("Oct4 ",labelling_2d_int_data,sep=""),"Oct4_")
    
    
  }
  #Spheoid 3D_Int
  {
    dire<-paste(dird,"spheroid_3d_int/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_white_bg_histogram_timeline(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T1_")
    plot_white_bg_histogram_timeline(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T2_")
    plot_white_bg_histogram_timeline(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T3_")
    plot_white_bg_histogram_timeline(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T4"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T4_")
    
    plot_white_bg_boxplot_timeline(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T1_")
    plot_white_bg_boxplot_timeline(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T2_")
    plot_white_bg_boxplot_timeline(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T3_")
    plot_white_bg_boxplot_timeline(subset(spheroid_int_data$dna,spheroid_int_data$dna$trial=="T4"), cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_T4_")
    
    plot_white_bg_barplot_timeline_4_trials(spheroid_int_data$dna, cols,parameters_int_data,paste("DNA ",labelling_int_data,sep=""),"DNA_")
    
    plot_white_bg_histogram_timeline(subset(spheroid_int_data$actin,spheroid_int_data$actin$trial=="T1"), cols,parameters_int_data,paste("F-Actin ",labelling_int_data,sep=""),"ACTIN_T1_")
    plot_white_bg_histogram_timeline(subset(spheroid_int_data$actin,spheroid_int_data$actin$trial=="T2"), cols,parameters_int_data,paste("F-Actin ",labelling_int_data,sep=""),"ACTIN_T2_")
    plot_white_bg_histogram_timeline(subset(spheroid_int_data$actin,spheroid_int_data$actin$trial=="T3"), cols,parameters_int_data,paste("F-Actin ",labelling_int_data,sep=""),"ACTIN_T3_")
    plot_white_bg_histogram_timeline(subset(spheroid_int_data$actin,spheroid_int_data$actin$trial=="T4"), cols,parameters_int_data,paste("F-Actin ",labelling_int_data,sep=""),"ACTIN_T4_")
    
    plot_white_bg_boxplot_timeline(subset(spheroid_int_data$actin,spheroid_int_data$actin$trial=="T1"), cols,parameters_int_data,paste("F-Actin ",labelling_int_data,sep=""),"ACTIN_T1_")
    plot_white_bg_boxplot_timeline(subset(spheroid_int_data$actin,spheroid_int_data$actin$trial=="T2"), cols,parameters_int_data,paste("F-Actin ",labelling_int_data,sep=""),"ACTIN_T2_")
    plot_white_bg_boxplot_timeline(subset(spheroid_int_data$actin,spheroid_int_data$actin$trial=="T3"), cols,parameters_int_data,paste("F-Actin ",labelling_int_data,sep=""),"ACTIN_T3_")
    plot_white_bg_boxplot_timeline(subset(spheroid_int_data$actin,spheroid_int_data$actin$trial=="T4"), cols,parameters_int_data,paste("F-Actin ",labelling_int_data,sep=""),"ACTIN_T4_")
    
    plot_white_bg_barplot_timeline_4_trials(spheroid_int_data$actin, cols,parameters_int_data,paste("F-Actin ",labelling_int_data,sep=""),"ACTIN_")
    
    
    plot_white_bg_histogram_timeline(subset(spheroid_int_data$LAMINAC,spheroid_int_data$LAMINAC$trial=="T1"), cols,parameters_int_data,paste("LAMINAC ",labelling_int_data,sep=""),"LAMINAC_T1_")
    plot_white_bg_histogram_timeline(subset(spheroid_int_data$LAMINAC,spheroid_int_data$LAMINAC$trial=="T2"), cols,parameters_int_data,paste("LAMINAC ",labelling_int_data,sep=""),"LAMINAC_T2_")
    plot_white_bg_histogram_timeline(subset(spheroid_int_data$LAMINAC,spheroid_int_data$LAMINAC$trial=="T3"), cols,parameters_int_data,paste("LAMINAC ",labelling_int_data,sep=""),"LAMINAC_T3_")
    plot_white_bg_histogram_timeline(subset(spheroid_int_data$LAMINAC,spheroid_int_data$LAMINAC$trial=="T4"), cols,parameters_int_data,paste("LAMINAC ",labelling_int_data,sep=""),"LAMINAC_T4_")
    
    plot_white_bg_boxplot_timeline(subset(spheroid_int_data$LAMINAC,spheroid_int_data$LAMINAC$trial=="T1"), cols,parameters_int_data,paste("LAMINAC ",labelling_int_data,sep=""),"LAMINAC_T1_")
    plot_white_bg_boxplot_timeline(subset(spheroid_int_data$LAMINAC,spheroid_int_data$LAMINAC$trial=="T2"), cols,parameters_int_data,paste("LAMINAC ",labelling_int_data,sep=""),"LAMINAC_T2_")
    plot_white_bg_boxplot_timeline(subset(spheroid_int_data$LAMINAC,spheroid_int_data$LAMINAC$trial=="T3"), cols,parameters_int_data,paste("LAMINAC ",labelling_int_data,sep=""),"LAMINAC_T3_")
    plot_white_bg_boxplot_timeline(subset(spheroid_int_data$LAMINAC,spheroid_int_data$LAMINAC$trial=="T4"), cols,parameters_int_data,paste("LAMINAC ",labelling_int_data,sep=""),"LAMINAC_T4_")
    
    plot_white_bg_barplot_timeline_4_trials(spheroid_int_data$LAMINAC, cols,parameters_int_data,paste("LAMINAC ",labelling_int_data,sep=""),"LAMINAC_")
    
    plot_white_bg_histogram_timeline(subset(spheroid_int_data$oct4,spheroid_int_data$oct4$trial=="T1"), cols,parameters_int_data,paste("Oct4 ",labelling_int_data,sep=""),"Oct4_T1_")
    plot_white_bg_histogram_timeline(subset(spheroid_int_data$oct4,spheroid_int_data$oct4$trial=="T2"), cols,parameters_int_data,paste("Oct4 ",labelling_int_data,sep=""),"Oct4_T2_")
    plot_white_bg_histogram_timeline(subset(spheroid_int_data$oct4,spheroid_int_data$oct4$trial=="T3"), cols,parameters_int_data,paste("Oct4 ",labelling_int_data,sep=""),"Oct4_T3_")
    plot_white_bg_histogram_timeline(subset(spheroid_int_data$oct4,spheroid_int_data$oct4$trial=="T4"), cols,parameters_int_data,paste("Oct4 ",labelling_int_data,sep=""),"Oct4_T4_")
    
    plot_white_bg_boxplot_timeline(subset(spheroid_int_data$oct4,spheroid_int_data$oct4$trial=="T1"), cols,parameters_int_data,paste("Oct4 ",labelling_int_data,sep=""),"Oct4_T1_")
    plot_white_bg_boxplot_timeline(subset(spheroid_int_data$oct4,spheroid_int_data$oct4$trial=="T2"), cols,parameters_int_data,paste("Oct4 ",labelling_int_data,sep=""),"Oct4_T2_")
    plot_white_bg_boxplot_timeline(subset(spheroid_int_data$oct4,spheroid_int_data$oct4$trial=="T3"), cols,parameters_int_data,paste("Oct4 ",labelling_int_data,sep=""),"Oct4_T3_")
    plot_white_bg_boxplot_timeline(subset(spheroid_int_data$oct4,spheroid_int_data$oct4$trial=="T4"), cols,parameters_int_data,paste("Oct4 ",labelling_int_data,sep=""),"Oct4_T4_")
    
    plot_white_bg_barplot_timeline_4_trials(spheroid_int_data$oct4, cols,parameters_int_data,paste("Oct4 ",labelling_int_data,sep=""),"Oct4_")
    
    
  }
  
  #Spheoid geometeric properties
  {
    dire<-paste(dird,"spheroid_properties/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_white_bg_histogram_timeline(subset(spheroid_combined,spheroid_combined$trial=="T1"), cols,parameters_geometrical_data,paste("Spheroid ",labelling_geometrical_data,sep=""),"geo_T1_")
    plot_white_bg_histogram_timeline(subset(spheroid_combined,spheroid_combined$trial=="T2"), cols,parameters_geometrical_data,paste("Spheroid ",labelling_geometrical_data,sep=""),"geo_T2_")
    plot_white_bg_histogram_timeline(subset(spheroid_combined,spheroid_combined$trial=="T3"), cols,parameters_geometrical_data,paste("Spheroid ",labelling_geometrical_data,sep=""),"geo_T3_")
    plot_white_bg_histogram_timeline(subset(spheroid_combined,spheroid_combined$trial=="T4"), cols,parameters_geometrical_data,paste("Spheroid ",labelling_geometrical_data,sep=""),"geo_T4_")
    
    plot_white_bg_boxplot_timeline(subset(spheroid_combined,spheroid_combined$trial=="T1"), cols,parameters_geometrical_data,labelling_geometrical_data,"geo_T1_")
    plot_white_bg_boxplot_timeline(subset(spheroid_combined,spheroid_combined$trial=="T2"), cols,parameters_geometrical_data,labelling_geometrical_data,"geo_T2_")
    plot_white_bg_boxplot_timeline(subset(spheroid_combined,spheroid_combined$trial=="T3"), cols,parameters_geometrical_data,labelling_geometrical_data,"geo_T3_")
    plot_white_bg_boxplot_timeline(subset(spheroid_combined,spheroid_combined$trial=="T4"), cols,parameters_geometrical_data,labelling_geometrical_data,"geo_T4_")
    
    plot_white_bg_barplot_timeline_4_trials(spheroid_combined, cols,parameters_geometrical_data,paste("Spheroid ",labelling_geometrical_data,sep=""),"geo_")
    
  }
  #Nuclear geometeric properties
  {
    dire<-paste(dird,"nuclear_properties/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_white_bg_histogram_timeline(subset(nucleus_combined,nucleus_combined$trial=="T1"), cols,parameters_geo_texture_data,paste("Nuclear ",labelling_geo_texture_data,sep=""),"geo_T1_")
    plot_white_bg_histogram_timeline(subset(nucleus_combined,nucleus_combined$trial=="T2"), cols,parameters_geo_texture_data,paste("Nuclear ",labelling_geo_texture_data,sep=""),"geo_T2_")
    plot_white_bg_histogram_timeline(subset(nucleus_combined,nucleus_combined$trial=="T3"), cols,parameters_geo_texture_data,paste("Nuclear ",labelling_geo_texture_data,sep=""),"geo_T3_")
    plot_white_bg_histogram_timeline(subset(nucleus_combined,nucleus_combined$trial=="T4"), cols,parameters_geo_texture_data,paste("Nuclear ",labelling_geo_texture_data,sep=""),"geo_T4_")
    
    plot_white_bg_boxplot_timeline(subset(nucleus_combined,nucleus_combined$trial=="T1"), cols,parameters_geo_texture_data,labelling_geo_texture_data,"geo_T1_")
    plot_white_bg_boxplot_timeline(subset(nucleus_combined,nucleus_combined$trial=="T2"), cols,parameters_geo_texture_data,labelling_geo_texture_data,"geo_T2_")
    plot_white_bg_boxplot_timeline(subset(nucleus_combined,nucleus_combined$trial=="T3"), cols,parameters_geo_texture_data,labelling_geo_texture_data,"geo_T3_")
    plot_white_bg_boxplot_timeline(subset(nucleus_combined,nucleus_combined$trial=="T4"), cols,parameters_geo_texture_data,labelling_geo_texture_data,"geo_T4_")
    
    plot_white_bg_barplot_timeline_4_trials(nucleus_combined, cols,parameters_geo_texture_data,paste("Nuclear ",labelling_geo_texture_data,sep=""),"geo_")
    
  }
  #Nuclear Proj_Int
  {
    
    dire<-paste(dird,"nuclear_proj_int/",sep="")
    dir.create(dire)
    setwd(dire)
    
    plot_white_bg_histogram_timeline(subset(nuclei_2dint_data$actin,nuclei_2dint_data$actin$trial=="T1"), cols,parameters_2d_int_data,paste("F-Actin ",nuc_labelling_2d_int_data,sep=""),"ACTIN_T1_")
    plot_white_bg_histogram_timeline(subset(nuclei_2dint_data$actin,nuclei_2dint_data$actin$trial=="T2"), cols,parameters_2d_int_data,paste("F-Actin ",nuc_labelling_2d_int_data,sep=""),"ACTIN_T2_")
    plot_white_bg_histogram_timeline(subset(nuclei_2dint_data$actin,nuclei_2dint_data$actin$trial=="T3"), cols,parameters_2d_int_data,paste("F-Actin ",nuc_labelling_2d_int_data,sep=""),"ACTIN_T3_")
    plot_white_bg_histogram_timeline(subset(nuclei_2dint_data$actin,nuclei_2dint_data$actin$trial=="T4"), cols,parameters_2d_int_data,paste("F-Actin ",nuc_labelling_2d_int_data,sep=""),"ACTIN_T4_")
    
    plot_white_bg_boxplot_timeline(subset(nuclei_2dint_data$actin,nuclei_2dint_data$actin$trial=="T1"), cols,parameters_2d_int_data,paste("F-Actin ",nuc_labelling_2d_int_data,sep=""),"ACTIN_T1_")
    plot_white_bg_boxplot_timeline(subset(nuclei_2dint_data$actin,nuclei_2dint_data$actin$trial=="T2"), cols,parameters_2d_int_data,paste("F-Actin ",nuc_labelling_2d_int_data,sep=""),"ACTIN_T2_")
    plot_white_bg_boxplot_timeline(subset(nuclei_2dint_data$actin,nuclei_2dint_data$actin$trial=="T3"), cols,parameters_2d_int_data,paste("F-Actin ",nuc_labelling_2d_int_data,sep=""),"ACTIN_T3_")
    plot_white_bg_boxplot_timeline(subset(nuclei_2dint_data$actin,nuclei_2dint_data$actin$trial=="T4"), cols,parameters_2d_int_data,paste("F-Actin ",nuc_labelling_2d_int_data,sep=""),"ACTIN_T4_")
    
    plot_white_bg_barplot_timeline_4_trials(nuclei_2dint_data$actin, cols,parameters_2d_int_data,paste("F-Actin ",nuc_labelling_2d_int_data,sep=""),"ACTIN_")
    
    
    plot_white_bg_histogram_timeline(subset(nuclei_2dint_data$LAMINAC,nuclei_2dint_data$LAMINAC$trial=="T1"), cols,parameters_2d_int_data,paste("LAMINAC ",nuc_labelling_2d_int_data,sep=""),"LAMINAC_T1_")
    plot_white_bg_histogram_timeline(subset(nuclei_2dint_data$LAMINAC,nuclei_2dint_data$LAMINAC$trial=="T2"), cols,parameters_2d_int_data,paste("LAMINAC ",nuc_labelling_2d_int_data,sep=""),"LAMINAC_T2_")
    plot_white_bg_histogram_timeline(subset(nuclei_2dint_data$LAMINAC,nuclei_2dint_data$LAMINAC$trial=="T3"), cols,parameters_2d_int_data,paste("LAMINAC ",nuc_labelling_2d_int_data,sep=""),"LAMINAC_T3_")
    plot_white_bg_histogram_timeline(subset(nuclei_2dint_data$LAMINAC,nuclei_2dint_data$LAMINAC$trial=="T4"), cols,parameters_2d_int_data,paste("LAMINAC ",nuc_labelling_2d_int_data,sep=""),"LAMINAC_T4_")
    
    plot_white_bg_boxplot_timeline(subset(nuclei_2dint_data$LAMINAC,nuclei_2dint_data$LAMINAC$trial=="T1"), cols,parameters_2d_int_data,paste("LAMINAC ",nuc_labelling_2d_int_data,sep=""),"LAMINAC_T1_")
    plot_white_bg_boxplot_timeline(subset(nuclei_2dint_data$LAMINAC,nuclei_2dint_data$LAMINAC$trial=="T2"), cols,parameters_2d_int_data,paste("LAMINAC ",nuc_labelling_2d_int_data,sep=""),"LAMINAC_T2_")
    plot_white_bg_boxplot_timeline(subset(nuclei_2dint_data$LAMINAC,nuclei_2dint_data$LAMINAC$trial=="T3"), cols,parameters_2d_int_data,paste("LAMINAC ",nuc_labelling_2d_int_data,sep=""),"LAMINAC_T3_")
    plot_white_bg_boxplot_timeline(subset(nuclei_2dint_data$LAMINAC,nuclei_2dint_data$LAMINAC$trial=="T4"), cols,parameters_2d_int_data,paste("LAMINAC ",nuc_labelling_2d_int_data,sep=""),"LAMINAC_T4_")
    
    plot_white_bg_barplot_timeline_4_trials(nuclei_2dint_data$LAMINAC, cols,parameters_2d_int_data,paste("LAMINAC ",nuc_labelling_2d_int_data,sep=""),"LAMINAC_")
    
    plot_white_bg_histogram_timeline(subset(nuclei_2dint_data$oct4,nuclei_2dint_data$oct4$trial=="T1"), cols,parameters_2d_int_data,paste("Oct4 ",nuc_labelling_2d_int_data,sep=""),"Oct4_T1_")
    plot_white_bg_histogram_timeline(subset(nuclei_2dint_data$oct4,nuclei_2dint_data$oct4$trial=="T2"), cols,parameters_2d_int_data,paste("Oct4 ",nuc_labelling_2d_int_data,sep=""),"Oct4_T2_")
    plot_white_bg_histogram_timeline(subset(nuclei_2dint_data$oct4,nuclei_2dint_data$oct4$trial=="T3"), cols,parameters_2d_int_data,paste("Oct4 ",nuc_labelling_2d_int_data,sep=""),"Oct4_T3_")
    plot_white_bg_histogram_timeline(subset(nuclei_2dint_data$oct4,nuclei_2dint_data$oct4$trial=="T4"), cols,parameters_2d_int_data,paste("Oct4 ",nuc_labelling_2d_int_data,sep=""),"Oct4_T4_")
    
    plot_white_bg_boxplot_timeline(subset(nuclei_2dint_data$oct4,nuclei_2dint_data$oct4$trial=="T1"), cols,parameters_2d_int_data,paste("Oct4 ",nuc_labelling_2d_int_data,sep=""),"Oct4_T1_")
    plot_white_bg_boxplot_timeline(subset(nuclei_2dint_data$oct4,nuclei_2dint_data$oct4$trial=="T2"), cols,parameters_2d_int_data,paste("Oct4 ",nuc_labelling_2d_int_data,sep=""),"Oct4_T2_")
    plot_white_bg_boxplot_timeline(subset(nuclei_2dint_data$oct4,nuclei_2dint_data$oct4$trial=="T3"), cols,parameters_2d_int_data,paste("Oct4 ",nuc_labelling_2d_int_data,sep=""),"Oct4_T3_")
    plot_white_bg_boxplot_timeline(subset(nuclei_2dint_data$oct4,nuclei_2dint_data$oct4$trial=="T4"), cols,parameters_2d_int_data,paste("Oct4 ",nuc_labelling_2d_int_data,sep=""),"Oct4_T4_")
    
    plot_white_bg_barplot_timeline_4_trials(nuclei_2dint_data$oct4, cols,parameters_2d_int_data,paste("Oct4 ",nuc_labelling_2d_int_data,sep=""),"Oct4_")
    
    
  }
  #Nuclear 3D_Int
  {
    dire<-paste(dird,"nuclear_3d_int/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_white_bg_histogram_timeline(subset(nuclei_3dint_data$dna,nuclei_3dint_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T1_")
    plot_white_bg_histogram_timeline(subset(nuclei_3dint_data$dna,nuclei_3dint_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T2_")
    plot_white_bg_histogram_timeline(subset(nuclei_3dint_data$dna,nuclei_3dint_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T3_")
    plot_white_bg_histogram_timeline(subset(nuclei_3dint_data$dna,nuclei_3dint_data$dna$trial=="T4"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T4_")
    
    plot_white_bg_boxplot_timeline(subset(nuclei_3dint_data$dna,nuclei_3dint_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T1_")
    plot_white_bg_boxplot_timeline(subset(nuclei_3dint_data$dna,nuclei_3dint_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T2_")
    plot_white_bg_boxplot_timeline(subset(nuclei_3dint_data$dna,nuclei_3dint_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T3_")
    plot_white_bg_boxplot_timeline(subset(nuclei_3dint_data$dna,nuclei_3dint_data$dna$trial=="T4"), cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_T4_")
    
    plot_white_bg_barplot_timeline_4_trials(nuclei_3dint_data$dna, cols,parameters_int_data,paste("DNA ",nuclabelling_int_data,sep=""),"DNA_")
    
    plot_white_bg_histogram_timeline(subset(nuclei_3dint_data$actin,nuclei_3dint_data$actin$trial=="T1"), cols,parameters_int_data,paste("F-Actin ",nuclabelling_int_data,sep=""),"ACTIN_T1_")
    plot_white_bg_histogram_timeline(subset(nuclei_3dint_data$actin,nuclei_3dint_data$actin$trial=="T2"), cols,parameters_int_data,paste("F-Actin ",nuclabelling_int_data,sep=""),"ACTIN_T2_")
    plot_white_bg_histogram_timeline(subset(nuclei_3dint_data$actin,nuclei_3dint_data$actin$trial=="T3"), cols,parameters_int_data,paste("F-Actin ",nuclabelling_int_data,sep=""),"ACTIN_T3_")
    plot_white_bg_histogram_timeline(subset(nuclei_3dint_data$actin,nuclei_3dint_data$actin$trial=="T4"), cols,parameters_int_data,paste("F-Actin ",nuclabelling_int_data,sep=""),"ACTIN_T4_")
    
    plot_white_bg_boxplot_timeline(subset(nuclei_3dint_data$actin,nuclei_3dint_data$actin$trial=="T1"), cols,parameters_int_data,paste("F-Actin ",nuclabelling_int_data,sep=""),"ACTIN_T1_")
    plot_white_bg_boxplot_timeline(subset(nuclei_3dint_data$actin,nuclei_3dint_data$actin$trial=="T2"), cols,parameters_int_data,paste("F-Actin ",nuclabelling_int_data,sep=""),"ACTIN_T2_")
    plot_white_bg_boxplot_timeline(subset(nuclei_3dint_data$actin,nuclei_3dint_data$actin$trial=="T3"), cols,parameters_int_data,paste("F-Actin ",nuclabelling_int_data,sep=""),"ACTIN_T3_")
    plot_white_bg_boxplot_timeline(subset(nuclei_3dint_data$actin,nuclei_3dint_data$actin$trial=="T4"), cols,parameters_int_data,paste("F-Actin ",nuclabelling_int_data,sep=""),"ACTIN_T4_")
    
    plot_white_bg_barplot_timeline_4_trials(nuclei_3dint_data$actin, cols,parameters_int_data,paste("F-Actin ",nuclabelling_int_data,sep=""),"ACTIN_")
    
    
    plot_white_bg_histogram_timeline(subset(nuclei_3dint_data$LAMINAC,nuclei_3dint_data$LAMINAC$trial=="T1"), cols,parameters_int_data,paste("LAMINAC ",nuclabelling_int_data,sep=""),"LAMINAC_T1_")
    plot_white_bg_histogram_timeline(subset(nuclei_3dint_data$LAMINAC,nuclei_3dint_data$LAMINAC$trial=="T2"), cols,parameters_int_data,paste("LAMINAC ",nuclabelling_int_data,sep=""),"LAMINAC_T2_")
    plot_white_bg_histogram_timeline(subset(nuclei_3dint_data$LAMINAC,nuclei_3dint_data$LAMINAC$trial=="T3"), cols,parameters_int_data,paste("LAMINAC ",nuclabelling_int_data,sep=""),"LAMINAC_T3_")
    plot_white_bg_histogram_timeline(subset(nuclei_3dint_data$LAMINAC,nuclei_3dint_data$LAMINAC$trial=="T4"), cols,parameters_int_data,paste("LAMINAC ",nuclabelling_int_data,sep=""),"LAMINAC_T4_")
    
    plot_white_bg_boxplot_timeline(subset(nuclei_3dint_data$LAMINAC,nuclei_3dint_data$LAMINAC$trial=="T1"), cols,parameters_int_data,paste("LAMINAC ",nuclabelling_int_data,sep=""),"LAMINAC_T1_")
    plot_white_bg_boxplot_timeline(subset(nuclei_3dint_data$LAMINAC,nuclei_3dint_data$LAMINAC$trial=="T2"), cols,parameters_int_data,paste("LAMINAC ",nuclabelling_int_data,sep=""),"LAMINAC_T2_")
    plot_white_bg_boxplot_timeline(subset(nuclei_3dint_data$LAMINAC,nuclei_3dint_data$LAMINAC$trial=="T3"), cols,parameters_int_data,paste("LAMINAC ",nuclabelling_int_data,sep=""),"LAMINAC_T3_")
    plot_white_bg_boxplot_timeline(subset(nuclei_3dint_data$LAMINAC,nuclei_3dint_data$LAMINAC$trial=="T4"), cols,parameters_int_data,paste("LAMINAC ",nuclabelling_int_data,sep=""),"LAMINAC_T4_")
    
    plot_white_bg_barplot_timeline_4_trials(nuclei_3dint_data$LAMINAC, cols,parameters_int_data,paste("LAMINAC ",nuclabelling_int_data,sep=""),"LAMINAC_")
    
    plot_white_bg_histogram_timeline(subset(nuclei_3dint_data$oct4,nuclei_3dint_data$oct4$trial=="T1"), cols,parameters_int_data,paste("Oct4 ",nuclabelling_int_data,sep=""),"Oct4_T1_")
    plot_white_bg_histogram_timeline(subset(nuclei_3dint_data$oct4,nuclei_3dint_data$oct4$trial=="T2"), cols,parameters_int_data,paste("Oct4 ",nuclabelling_int_data,sep=""),"Oct4_T2_")
    plot_white_bg_histogram_timeline(subset(nuclei_3dint_data$oct4,nuclei_3dint_data$oct4$trial=="T3"), cols,parameters_int_data,paste("Oct4 ",nuclabelling_int_data,sep=""),"Oct4_T3_")
    plot_white_bg_histogram_timeline(subset(nuclei_3dint_data$oct4,nuclei_3dint_data$oct4$trial=="T4"), cols,parameters_int_data,paste("Oct4 ",nuclabelling_int_data,sep=""),"Oct4_T4_")
    
    plot_white_bg_boxplot_timeline(subset(nuclei_3dint_data$oct4,nuclei_3dint_data$oct4$trial=="T1"), cols,parameters_int_data,paste("Oct4 ",nuclabelling_int_data,sep=""),"Oct4_T1_")
    plot_white_bg_boxplot_timeline(subset(nuclei_3dint_data$oct4,nuclei_3dint_data$oct4$trial=="T2"), cols,parameters_int_data,paste("Oct4 ",nuclabelling_int_data,sep=""),"Oct4_T2_")
    plot_white_bg_boxplot_timeline(subset(nuclei_3dint_data$oct4,nuclei_3dint_data$oct4$trial=="T3"), cols,parameters_int_data,paste("Oct4 ",nuclabelling_int_data,sep=""),"Oct4_T3_")
    plot_white_bg_boxplot_timeline(subset(nuclei_3dint_data$oct4,nuclei_3dint_data$oct4$trial=="T4"), cols,parameters_int_data,paste("Oct4 ",nuclabelling_int_data,sep=""),"Oct4_T4_")
    
    plot_white_bg_barplot_timeline_4_trials(nuclei_3dint_data$oct4, cols,parameters_int_data,paste("Oct4 ",nuclabelling_int_data,sep=""),"Oct4_")
    
    
  }
  #Cellular 3D_Int
  {
    dire<-paste(dird,"cellular_3d_int/",sep="")
    dir.create(dire)
    setwd(dire)
    plot_white_bg_histogram_timeline(subset(cell_3dint_data$dna,cell_3dint_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T1_")
    plot_white_bg_histogram_timeline(subset(cell_3dint_data$dna,cell_3dint_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T2_")
    plot_white_bg_histogram_timeline(subset(cell_3dint_data$dna,cell_3dint_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T3_")
    plot_white_bg_histogram_timeline(subset(cell_3dint_data$dna,cell_3dint_data$dna$trial=="T4"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T4_")
    
    plot_white_bg_boxplot_timeline(subset(cell_3dint_data$dna,cell_3dint_data$dna$trial=="T1"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T1_")
    plot_white_bg_boxplot_timeline(subset(cell_3dint_data$dna,cell_3dint_data$dna$trial=="T2"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T2_")
    plot_white_bg_boxplot_timeline(subset(cell_3dint_data$dna,cell_3dint_data$dna$trial=="T3"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T3_")
    plot_white_bg_boxplot_timeline(subset(cell_3dint_data$dna,cell_3dint_data$dna$trial=="T4"), cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_T4_")
    
    plot_white_bg_barplot_timeline_4_trials(cell_3dint_data$dna, cols,parameters_int_data,paste("DNA ",celllabelling_int_data,sep=""),"DNA_")
    
    plot_white_bg_histogram_timeline(subset(cell_3dint_data$actin,cell_3dint_data$actin$trial=="T1"), cols,parameters_int_data,paste("F-Actin ",celllabelling_int_data,sep=""),"ACTIN_T1_")
    plot_white_bg_histogram_timeline(subset(cell_3dint_data$actin,cell_3dint_data$actin$trial=="T2"), cols,parameters_int_data,paste("F-Actin ",celllabelling_int_data,sep=""),"ACTIN_T2_")
    plot_white_bg_histogram_timeline(subset(cell_3dint_data$actin,cell_3dint_data$actin$trial=="T3"), cols,parameters_int_data,paste("F-Actin ",celllabelling_int_data,sep=""),"ACTIN_T3_")
    plot_white_bg_histogram_timeline(subset(cell_3dint_data$actin,cell_3dint_data$actin$trial=="T4"), cols,parameters_int_data,paste("F-Actin ",celllabelling_int_data,sep=""),"ACTIN_T4_")
    
    plot_white_bg_boxplot_timeline(subset(cell_3dint_data$actin,cell_3dint_data$actin$trial=="T1"), cols,parameters_int_data,paste("F-Actin ",celllabelling_int_data,sep=""),"ACTIN_T1_")
    plot_white_bg_boxplot_timeline(subset(cell_3dint_data$actin,cell_3dint_data$actin$trial=="T2"), cols,parameters_int_data,paste("F-Actin ",celllabelling_int_data,sep=""),"ACTIN_T2_")
    plot_white_bg_boxplot_timeline(subset(cell_3dint_data$actin,cell_3dint_data$actin$trial=="T3"), cols,parameters_int_data,paste("F-Actin ",celllabelling_int_data,sep=""),"ACTIN_T3_")
    plot_white_bg_boxplot_timeline(subset(cell_3dint_data$actin,cell_3dint_data$actin$trial=="T4"), cols,parameters_int_data,paste("F-Actin ",celllabelling_int_data,sep=""),"ACTIN_T4_")
    
    plot_white_bg_barplot_timeline_4_trials(cell_3dint_data$actin, cols,parameters_int_data,paste("F-Actin ",celllabelling_int_data,sep=""),"ACTIN_")
    
    
    plot_white_bg_histogram_timeline(subset(cell_3dint_data$LAMINAC,cell_3dint_data$LAMINAC$trial=="T1"), cols,parameters_int_data,paste("LAMINAC ",celllabelling_int_data,sep=""),"LAMINAC_T1_")
    plot_white_bg_histogram_timeline(subset(cell_3dint_data$LAMINAC,cell_3dint_data$LAMINAC$trial=="T2"), cols,parameters_int_data,paste("LAMINAC ",celllabelling_int_data,sep=""),"LAMINAC_T2_")
    plot_white_bg_histogram_timeline(subset(cell_3dint_data$LAMINAC,cell_3dint_data$LAMINAC$trial=="T3"), cols,parameters_int_data,paste("LAMINAC ",celllabelling_int_data,sep=""),"LAMINAC_T3_")
    plot_white_bg_histogram_timeline(subset(cell_3dint_data$LAMINAC,cell_3dint_data$LAMINAC$trial=="T4"), cols,parameters_int_data,paste("LAMINAC ",celllabelling_int_data,sep=""),"LAMINAC_T4_")
    
    plot_white_bg_boxplot_timeline(subset(cell_3dint_data$LAMINAC,cell_3dint_data$LAMINAC$trial=="T1"), cols,parameters_int_data,paste("LAMINAC ",celllabelling_int_data,sep=""),"LAMINAC_T1_")
    plot_white_bg_boxplot_timeline(subset(cell_3dint_data$LAMINAC,cell_3dint_data$LAMINAC$trial=="T2"), cols,parameters_int_data,paste("LAMINAC ",celllabelling_int_data,sep=""),"LAMINAC_T2_")
    plot_white_bg_boxplot_timeline(subset(cell_3dint_data$LAMINAC,cell_3dint_data$LAMINAC$trial=="T3"), cols,parameters_int_data,paste("LAMINAC ",celllabelling_int_data,sep=""),"LAMINAC_T3_")
    plot_white_bg_boxplot_timeline(subset(cell_3dint_data$LAMINAC,cell_3dint_data$LAMINAC$trial=="T4"), cols,parameters_int_data,paste("LAMINAC ",celllabelling_int_data,sep=""),"LAMINAC_T4_")
    
    plot_white_bg_barplot_timeline_4_trials(cell_3dint_data$LAMINAC, cols,parameters_int_data,paste("LAMINAC ",celllabelling_int_data,sep=""),"LAMINAC_")
    
    plot_white_bg_histogram_timeline(subset(cell_3dint_data$oct4,cell_3dint_data$oct4$trial=="T1"), cols,parameters_int_data,paste("Oct4 ",celllabelling_int_data,sep=""),"Oct4_T1_")
    plot_white_bg_histogram_timeline(subset(cell_3dint_data$oct4,cell_3dint_data$oct4$trial=="T2"), cols,parameters_int_data,paste("Oct4 ",celllabelling_int_data,sep=""),"Oct4_T2_")
    plot_white_bg_histogram_timeline(subset(cell_3dint_data$oct4,cell_3dint_data$oct4$trial=="T3"), cols,parameters_int_data,paste("Oct4 ",celllabelling_int_data,sep=""),"Oct4_T3_")
    plot_white_bg_histogram_timeline(subset(cell_3dint_data$oct4,cell_3dint_data$oct4$trial=="T4"), cols,parameters_int_data,paste("Oct4 ",celllabelling_int_data,sep=""),"Oct4_T4_")
    
    plot_white_bg_boxplot_timeline(subset(cell_3dint_data$oct4,cell_3dint_data$oct4$trial=="T1"), cols,parameters_int_data,paste("Oct4 ",celllabelling_int_data,sep=""),"Oct4_T1_")
    plot_white_bg_boxplot_timeline(subset(cell_3dint_data$oct4,cell_3dint_data$oct4$trial=="T2"), cols,parameters_int_data,paste("Oct4 ",celllabelling_int_data,sep=""),"Oct4_T2_")
    plot_white_bg_boxplot_timeline(subset(cell_3dint_data$oct4,cell_3dint_data$oct4$trial=="T3"), cols,parameters_int_data,paste("Oct4 ",celllabelling_int_data,sep=""),"Oct4_T3_")
    plot_white_bg_boxplot_timeline(subset(cell_3dint_data$oct4,cell_3dint_data$oct4$trial=="T4"), cols,parameters_int_data,paste("Oct4 ",celllabelling_int_data,sep=""),"Oct4_T4_")
    
    plot_white_bg_barplot_timeline_4_trials(cell_3dint_data$oct4, cols,parameters_int_data,paste("Oct4 ",celllabelling_int_data,sep=""),"Oct4_")
    
    
  }
  
}

#write the collected data
{
  dir.create("E:/HMF3A_reprogramming/R_analysis/actin_LAMINAC_oct4/combined_data/")
  setwd("E:/HMF3A_reprogramming/R_analysis/actin_LAMINAC_oct4/combined_data/")
  
  write.csv(spheroid_combined, file="spheroid_combined.csv")
  write.csv(nucleus_combined, file="nucleus_combined.csv")
  
  write.csv(spheroid_2dint_data$actin, file="spheroid_2d_int_actin.csv")
  write.csv(spheroid_2dint_data$oct4, file="spheroid_2d_int_oct4.csv")
  write.csv(spheroid_2dint_data$LAMINAC, file="spheroid_2d_int_LAMINAC.csv")
  
  write.csv(spheroid_int_data$actin, file="spheroid_int_actin.csv")
  write.csv(spheroid_int_data$oct4, file="spheroid_int_oct4.csv")
  write.csv(spheroid_int_data$LAMINAC, file="spheroid_int_LAMINAC.csv")
  write.csv(spheroid_int_data$dna, file="spheroid_int_dna.csv")
  
  write.csv(nuclei_2dint_data$actin, file="nuclear_2d_int_actin.csv")
  write.csv(nuclei_2dint_data$oct4, file="nuclear_2d_int_oct4.csv")
  write.csv(nuclei_2dint_data$LAMINAC, file="nuclear_2d_int_LAMINAC.csv")
  
  write.csv(nuclei_3dint_data$actin, file="nuclear_int_actin.csv")
  write.csv(nuclei_3dint_data$oct4, file="nuclear_int_oct4.csv")
  write.csv(nuclei_3dint_data$LAMINAC, file="nuclear_int_LAMINAC.csv")
  write.csv(nuclei_3dint_data$dna, file="nuclear_int_dna.csv")
  
  write.csv(cell_3dint_data$actin, file="cellular_int_actin.csv")
  write.csv(cell_3dint_data$oct4, file="cellular_int_oct4.csv")
  write.csv(cell_3dint_data$LAMINAC, file="cellular_int_LAMINAC.csv")
  write.csv(cell_3dint_data$dna, file="cellular_int_dna.csv")
  
}



