#functions
source("E:/HMF3A_reprogramming/R_analysis/functions.R")

# initialise the experiment details and the subdirectories required
{
  path_to_ex<-c("E:/HMF3A_reprogramming/Timeline_H3k9ac/20191217_hmf_s13_repg_h3h9ac_568/",
                "E:/HMF3A_reprogramming/Timeline_H3k9ac/20191220_hmf_s12_repg_h3k9ac_568/",
                "E:/HMF3A_reprogramming/Timeline_H3k9ac/20191227_hmf_s17_oct4_488_ki67_568_h3k9ac/")
  
  exp<-c("20191217_hmf_s13_repg_h3h9ac_568","20191220_hmf_s12_repg_h3k9ac_568","20191227_hmf_s17_oct4_488_ki67_568_h3k9ac")
  
  data_types<-c("/3D geometrical data spheroid/",
                "/3D int_data spheroid/DNA/","/3D int_data spheroid/H3k9ac/")
  data_types_2d<-c("/2D_measures_spheroid/2D_spheroid_DNA.csv","/2D_measures_spheroid/2D_spheroid_H3k9ac.csv")
  nuc_data_types<-c("/3D geometrical data/","/3D ellipsoid/","/3D geometerical_simple/","/3D shape measure/",
                    "/3D int_data/DNA/","/3D int_data/H3K9AC/")
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
                       "/2D_measures_nuclei/2D_nucleusH3K9AC.csv")
  samples<-c("D0","D2","D4","D6","D8")
  trial<-c("T1","T2","T3")
  
  dird_apo<-"E:/HMF3A_reprogramming/R_analysis/H3k9ac/"
  dir.create(dird_apo)
  dird_H3k9ac<-paste(dird_apo,"spheroid_shells_H3k9ac", sep="")
  dird_dna<-paste(dird_apo,"spheroid_shells_dna", sep="")
  dir.create(dird_H3k9ac)
  dir.create(dird_dna)
}

# read in all the required features for spheroids
{
  plot(0:1,0:1)
  #T1
  {
    geometrical_data<-combine_sample_sets(path_to_ex[1],samples,data_types[1],"tsv",exp[1],trial[1])
    T1DNA_int_data<-combine_sample_sets(path_to_ex[1],samples,data_types[2],"tsv",exp[1],trial[1])
    T1H3k9ac_int_data<-combine_sample_sets(path_to_ex[1],samples,data_types[3],"tsv",exp[1],trial[1])
    geo_int_2D_data<-combine_sample_sets_2d(path_to_ex[1],samples,data_types_2d[1],exp[1],trial[1])
    
    T1_2d_data_H3k9ac<-combine_sample_sets_2d(path_to_ex[1],samples,data_types_2d[2],exp[1],trial[1])
    
    shell_H3k9ac_T1<-shell_analsysis(path_to_ex[1],samples,"/Shells/H3k9ac",exp[1],trial[1],dird_H3k9ac)
    shell_dna_T1<-shell_analsysis(path_to_ex[1],samples,"/Shells/DNA",exp[1],trial[1],dird_dna)
    
    T1_dna_axial<-axial_analysis(path_to_ex[1],samples,exp[1],trial[1],"/DNA/")
    T1_H3k9ac_axial<-axial_analysis(path_to_ex[1],samples,exp[1],trial[1],"/H3k9ac/")
    
    T1_H3k9ac_angle<-angle_analsysis(path_to_ex[1],samples,"/H3k9ac",exp[1],trial[1])
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
    T2H3k9ac_int_data<-combine_sample_sets(path_to_ex[2],samples,data_types[3],"tsv",exp[2],trial[2])
    geo_int_2D_data<-combine_sample_sets_2d(path_to_ex[2],samples,data_types_2d[1],exp[2],trial[2])
    
    T2_2d_data_H3k9ac<-combine_sample_sets_2d(path_to_ex[2],samples,data_types_2d[2],exp[2],trial[2])
    
    shell_H3k9ac_T2<-shell_analsysis(path_to_ex[2],samples,"/Shells/H3k9ac",exp[2],trial[2],dird_H3k9ac)
    shell_dna_T2<-shell_analsysis(path_to_ex[2],samples,"/Shells/DNA",exp[2],trial[2],dird_dna)
    
    T2_dna_axial<-axial_analysis(path_to_ex[2],samples,exp[2],trial[2],"/DNA/")
    T2_H3k9ac_axial<-axial_analysis(path_to_ex[2],samples,exp[2],trial[2],"/H3k9ac/")
    
    T2_H3k9ac_angle<-angle_analsysis(path_to_ex[2],samples,"/H3k9ac",exp[2],trial[2])
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
    T3H3k9ac_int_data<-combine_sample_sets(path_to_ex[3],samples,data_types[3],"tsv",exp[3],trial[3])
    geo_int_2D_data<-combine_sample_sets_2d(path_to_ex[3],samples,data_types_2d[1],exp[3],trial[3])
    
    T3_2d_data_H3k9ac<-combine_sample_sets_2d(path_to_ex[3],samples,data_types_2d[2],exp[3],trial[3])
    
    shell_H3k9ac_T3<-shell_analsysis(path_to_ex[3],samples,"/Shells/H3k9ac",exp[3],trial[3],dird_H3k9ac)
    shell_dna_T3<-shell_analsysis(path_to_ex[3],samples,"/Shells/DNA",exp[3],trial[3],dird_dna)
    
    T3_dna_axial<-axial_analysis(path_to_ex[3],samples,exp[3],trial[3],"/DNA/")
    T3_H3k9ac_axial<-axial_analysis(path_to_ex[3],samples,exp[3],trial[3],"/H3k9ac/")
    
    T3_H3k9ac_angle<-angle_analsysis(path_to_ex[3],samples,"/H3k9ac",exp[3],trial[3])
    T3_dna_angle<-angle_analsysis(path_to_ex[3],samples,"/DNA",exp[3],trial[3])
    
    geo_int_2D_data$Label<-substring(geo_int_2D_data$Label,5, nchar(geo_int_2D_data$Label)-4)
    geometrical_data$Label<-substring(geometrical_data$Label,0, nchar(geometrical_data$Label)-1)
    T3combined<-merge(geometrical_data,geo_int_2D_data, by="Label")
    
    rm(geometrical_data,geo_int_2D_data)
    
  }
  
  spheroid_2dint_data_H3k9ac<-rbind(T1_2d_data_H3k9ac,T2_2d_data_H3k9ac,T3_2d_data_H3k9ac)
  
  spheroid_int_data_H3k9ac<-rbind(T1H3k9ac_int_data,T2H3k9ac_int_data,T3H3k9ac_int_data)
  spheroid_int_data_dna<-rbind(T1DNA_int_data,T2DNA_int_data,T3DNA_int_data)
  
  spheroid_radial_H3k9ac<-rbind(shell_H3k9ac_T1,shell_H3k9ac_T2,shell_H3k9ac_T3)
  spheroid_radial_dna<-rbind(shell_dna_T1,shell_dna_T2,shell_dna_T3)
  
  spheroid_axial_dna<-rbind(T1_dna_axial,T2_dna_axial,T3_dna_axial)
  spheroid_axial_H3k9ac<-rbind(T1_H3k9ac_axial,T2_H3k9ac_axial,T3_H3k9ac_axial)
  
  spheroid_angle_H3k9ac<-rbind(T1_H3k9ac_angle,T2_H3k9ac_angle,T3_H3k9ac_angle)
  spheroid_angle_dna<-rbind(T1_dna_angle,T2_dna_angle,T3_dna_angle)
  
  spheroid_int_data<-list( H3k9ac=spheroid_int_data_H3k9ac,
                           dna=spheroid_int_data_dna)
  spheroid_radial<-list(H3k9ac=spheroid_radial_H3k9ac,
                        dna=spheroid_radial_dna)
  spheroid_axial<-list( H3k9ac=spheroid_axial_H3k9ac,
                        dna=spheroid_axial_dna)
  spheroid_angle<-list(H3k9ac=spheroid_angle_H3k9ac,
                       dna=spheroid_angle_dna)
  spheroid_2dint_data<-list( H3k9ac=spheroid_2dint_data_H3k9ac)
  
  
  spheroid_combined<-rbind(T1combined,T2combined,T3combined)
  colnames(spheroid_combined)[86:87]<-c("sample","trial")
  
  rm(T1combined,T2combined,T3combined,T1VIMENTIN_int_data,T2VIMENTIN_int_data,T3VIMENTIN_int_data,
     T1H3k9ac_int_data,T2H3k9ac_int_data,T3H3k9ac_int_data,T1ECAD_int_data,T2ECAD_int_data,T3ECAD_int_data,
     shell_H3k9ac_T1,shell_H3k9ac_T2,shell_H3k9ac_T3,shell_vimentin_T1,shell_vimentin_T2,shell_vimentin_T3,
     shell_ecad_T1,shell_ecad_T2,shell_ecad_T3,shell_dna_T1,shell_dna_T2,shell_dna_T3,
     T1_vimentin_axial,T2_vimentin_axial,T3_vimentin_axial,T1_dna_axial,T2_dna_axial,T3_dna_axial,
     T1_H3k9ac_axial,T2_H3k9ac_axial,T3_H3k9ac_axial,T1_ecad_axial,T2_ecad_axial,T3_ecad_axial,
     T1_vimentin_angle,T2_vimentin_angle,T3_vimentin_angle,T1_H3k9ac_angle,T2_H3k9ac_angle,T3_H3k9ac_angle,
     T1_ecad_angle,T2_ecad_angle,T3_ecad_angle,T1_dna_angle,T2_dna_angle,T3_dna_angle,
     T1DNA_int_data,T2DNA_int_data,T3DNA_int_data,
     spheroid_angle_vimentin, spheroid_angle_H3k9ac,spheroid_angle_ecad,spheroid_angle_dna,
     spheroid_axial_vimentin, spheroid_axial_H3k9ac,spheroid_axial_ecad,spheroid_axial_dna,
     spheroid_radial_vimentin, spheroid_radial_H3k9ac,spheroid_radial_ecad,spheroid_radial_dna,
     spheroid_int_data_vimentin, spheroid_int_data_H3k9ac,spheroid_int_data_ecad,spheroid_int_data_dna,
     T1_2d_data_VIMENTIN,T2_2d_data_VIMENTIN,T3_2d_data_VIMENTIN,T1_2d_data_H3k9ac,T2_2d_data_H3k9ac,T3_2d_data_H3k9ac,
     T1_2d_data_ECAD,T2_2d_data_ECAD,T3_2d_data_ECAD,spheroid_2dint_data_vimentin, spheroid_2dint_data_H3k9ac,
     spheroid_2dint_data_ecad)
  
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
    
    T1_H3K9AC_int_data<-combine_sample_sets(path_to_ex[1],samples,nuc_data_types[6],"tsv",exp[1],trial[1])
    
    T1_2d_H3K9AC_int_data<-combine_sample_sets_2d(path_to_ex[1],samples,nuc_data_types_2d[28],exp[1],trial[1])
    
    
    colnames(T1_compact_2D_data)[2]<-"Label"
    T1_compact_2D_data$Label<-substring(T1_compact_2D_data$Label,0, nchar(T1_compact_2D_data$Label)-4)
    T1_geo_int_2D_data$Label<-substring(T1_geo_int_2D_data$Label,5, nchar(T1_geo_int_2D_data$Label)-4)
    T1_geometrical_data$Label<-substring(T1_geometrical_data$Label,0, nchar(T1_geometrical_data$Label)-1)
    T1_DNA_int_data$Label<-substring(T1_DNA_int_data$Label,0, nchar(T1_DNA_int_data$Label)-1)
    T1_2d_H3K9AC_int_data$Label<-substring(T1_2d_H3K9AC_int_data$Label,5, nchar(T1_2d_H3K9AC_int_data$Label)-4)
    T1_H3K9AC_int_data$Label<-substring(T1_H3K9AC_int_data$Label,0, nchar(T1_H3K9AC_int_data$Label)-1)
    
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
    
    T2_H3K9AC_int_data<-combine_sample_sets(path_to_ex[2],samples,nuc_data_types[6],"tsv",exp[2],trial[2])
    
    T2_2d_H3K9AC_int_data<-combine_sample_sets_2d(path_to_ex[2],samples,nuc_data_types_2d[28],exp[2],trial[2])
    
    
    colnames(T2_compact_2D_data)[2]<-"Label"
    T2_compact_2D_data$Label<-substring(T2_compact_2D_data$Label,0, nchar(T2_compact_2D_data$Label)-4)
    T2_geo_int_2D_data$Label<-substring(T2_geo_int_2D_data$Label,5, nchar(T2_geo_int_2D_data$Label)-4)
    T2_geometrical_data$Label<-substring(T2_geometrical_data$Label,0, nchar(T2_geometrical_data$Label)-1)
    T2_2d_H3K9AC_int_data$Label<-substring(T2_2d_H3K9AC_int_data$Label,5, nchar(T2_2d_H3K9AC_int_data$Label)-4)
    T2_DNA_int_data$Label<-substring(T2_DNA_int_data$Label,0, nchar(T2_DNA_int_data$Label)-1)
    T2_H3K9AC_int_data$Label<-substring(T2_H3K9AC_int_data$Label,0, nchar(T2_H3K9AC_int_data$Label)-1)
    
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
    
    T3_H3K9AC_int_data<-combine_sample_sets(path_to_ex[3],samples,nuc_data_types[6],"tsv",exp[3],trial[3])
    
    T3_2d_H3K9AC_int_data<-combine_sample_sets_2d(path_to_ex[3],samples,nuc_data_types_2d[28],exp[3],trial[3])
    
    
    colnames(T3_compact_2D_data)[2]<-"Label"
    T3_compact_2D_data$Label<-substring(T3_compact_2D_data$Label,0, nchar(T3_compact_2D_data$Label)-4)
    T3_geo_int_2D_data$Label<-substring(T3_geo_int_2D_data$Label,5, nchar(T3_geo_int_2D_data$Label)-4)
    T3_geometrical_data$Label<-substring(T3_geometrical_data$Label,0, nchar(T3_geometrical_data$Label)-1)
    T3_2d_H3K9AC_int_data$Label<-substring(T3_2d_H3K9AC_int_data$Label,5, nchar(T3_2d_H3K9AC_int_data$Label)-4)
    T3_DNA_int_data$Label<-substring(T3_DNA_int_data$Label,0, nchar(T3_DNA_int_data$Label)-1)
    T3_H3K9AC_int_data$Label<-substring(T3_H3K9AC_int_data$Label,0, nchar(T3_H3K9AC_int_data$Label)-1)
    
    T3_combined<-merge(T3_geometrical_data,T3_geo_int_2D_data[,2:36], by="Label")
    T3_combined<-merge(T3_combined,T3_compact_2D_data[,2:8], by="Label")
    T3_combined<-merge(T3_combined,T3_DNA_int_data[,c(3:5,12:18)], by="Label")
    T3_combined<-merge(T3_combined,T3_edf_2D_data[,1:11], by="Label")
    T3_combined<-cbind(T3_combined,T3_gclm[,-c(1,128:131)])
    
    rm(T3_geometrical_data,T3_geo_int_2D_data,T3_compact_2D_data,T3_edf_2D_data,T3_gclm)
  }
  
  nuclei_2dint_data_H3K9AC<-rbind(T1_2d_H3K9AC_int_data,T2_2d_H3K9AC_int_data,T3_2d_H3K9AC_int_data)
  
  nuclei_int_data_H3K9AC<-rbind(T1_H3K9AC_int_data,T2_H3K9AC_int_data,T3_H3K9AC_int_data)
  nuclei_int_data_dna<-rbind(T1_DNA_int_data,T2_DNA_int_data,T3_DNA_int_data)
  
  
  nuclei_2dint_data<-list(H3K9AC=nuclei_2dint_data_H3K9AC)
  nuclei_3dint_data<-list(H3K9AC=nuclei_int_data_H3K9AC,dna=nuclei_int_data_dna)

  nucleus_combined<-rbind(T1_combined,T2_combined,T3_combined)
  
  rm(T1_2d_H3K9AC_int_data,T2_2d_H3K9AC_int_data,T3_2d_H3K9AC_int_data,T1_2d_OCT4_int_data,T2_2d_OCT4_int_data,T3_2d_OCT4_int_data,
     T1_2d_PMLC_int_data,T2_2d_PMLC_int_data,T3_2d_PMLC_int_data,T1_H3K9AC_int_data,T2_H3K9AC_int_data,T3_H3K9AC_int_data,
     T1_OCT4_int_data,T2_OCT4_int_data,T3_OCT4_int_data,T1_PMLC_int_data,T2_PMLC_int_data,T3_PMLC_int_data,
     T1_DNA_int_data,T2_DNA_int_data,T3_DNA_int_data,T1_cell_H3K9AC_int_data,T2_cell_H3K9AC_int_data,T3_cell_H3K9AC_int_data,
     T1_cell_OCT4_int_data,T2_cell_OCT4_int_data,T3_cell_OCT4_int_data,T1_cell_PMLC_int_data,T2_cell_PMLC_int_data,T3_cell_PMLC_int_data,
     T1_cell_DNA_int_data,T2_cell_DNA_int_data,T3_cell_DNA_int_data,T1_combined,T2_combined,T3_combined,
     nuclei_2dint_data_H3K9AC, nuclei_2dint_data_oct4,nuclei_2dint_data_pmlc,nuclei_int_data_H3K9AC, nuclei_int_data_oct4,
     nuclei_int_data_pmlc,nuclei_int_data_dna,cell_int_data_H3K9AC, cell_int_data_oct4,cell_int_data_pmlc,cell_int_data_dna)
  
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
  dird<-"E:/HMF3A_reprogramming/R_analysis/H3k9ac/black_bg_plots/"
  dir.create(dird)
  setwd(dird)
  cols<-colorRampPalette(c("red","yellow"))( 5) 
  
  #Spheoid Proj_Int
  {
    
    dire<-paste(dird,"spheroid_proj_int/",sep="")
    dir.create(dire)
    setwd(dire)
    
    
    plot_black_bg_histogram_timeline(subset(spheroid_2dint_data$H3k9ac,spheroid_2dint_data$H3k9ac$trial=="T1"), cols,parameters_2d_int_data,paste("H3k9ac ",labelling_2d_int_data,sep=""),"H3k9ac_T1_")
    plot_black_bg_histogram_timeline(subset(spheroid_2dint_data$H3k9ac,spheroid_2dint_data$H3k9ac$trial=="T2"), cols,parameters_2d_int_data,paste("H3k9ac ",labelling_2d_int_data,sep=""),"H3k9ac_T2_")
    plot_black_bg_histogram_timeline(subset(spheroid_2dint_data$H3k9ac,spheroid_2dint_data$H3k9ac$trial=="T3"), cols,parameters_2d_int_data,paste("H3k9ac ",labelling_2d_int_data,sep=""),"H3k9ac_T3_")
    
    plot_black_bg_boxplot_timeline(subset(spheroid_2dint_data$H3k9ac,spheroid_2dint_data$H3k9ac$trial=="T1"), cols,parameters_2d_int_data,paste("H3k9ac ",labelling_2d_int_data,sep=""),"H3k9ac_T1_")
    plot_black_bg_boxplot_timeline(subset(spheroid_2dint_data$H3k9ac,spheroid_2dint_data$H3k9ac$trial=="T2"), cols,parameters_2d_int_data,paste("H3k9ac ",labelling_2d_int_data,sep=""),"H3k9ac_T2_")
    plot_black_bg_boxplot_timeline(subset(spheroid_2dint_data$H3k9ac,spheroid_2dint_data$H3k9ac$trial=="T3"), cols,parameters_2d_int_data,paste("H3k9ac ",labelling_2d_int_data,sep=""),"H3k9ac_T3_")
    
    plot_black_bg_barplot_timeline(spheroid_2dint_data$H3k9ac, cols,parameters_2d_int_data,paste("H3k9ac ",labelling_2d_int_data,sep=""),"H3k9ac_")
    
    
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
    
    plot_black_bg_histogram_timeline(subset(spheroid_int_data$H3k9ac,spheroid_int_data$H3k9ac$trial=="T1"), cols,parameters_int_data,paste("H3k9ac ",labelling_int_data,sep=""),"H3k9ac_T1_")
    plot_black_bg_histogram_timeline(subset(spheroid_int_data$H3k9ac,spheroid_int_data$H3k9ac$trial=="T2"), cols,parameters_int_data,paste("H3k9ac ",labelling_int_data,sep=""),"H3k9ac_T2_")
    plot_black_bg_histogram_timeline(subset(spheroid_int_data$H3k9ac,spheroid_int_data$H3k9ac$trial=="T3"), cols,parameters_int_data,paste("H3k9ac ",labelling_int_data,sep=""),"H3k9ac_T3_")
    
    plot_black_bg_boxplot_timeline(subset(spheroid_int_data$H3k9ac,spheroid_int_data$H3k9ac$trial=="T1"), cols,parameters_int_data,paste("H3k9ac ",labelling_int_data,sep=""),"H3k9ac_T1_")
    plot_black_bg_boxplot_timeline(subset(spheroid_int_data$H3k9ac,spheroid_int_data$H3k9ac$trial=="T2"), cols,parameters_int_data,paste("H3k9ac ",labelling_int_data,sep=""),"H3k9ac_T2_")
    plot_black_bg_boxplot_timeline(subset(spheroid_int_data$H3k9ac,spheroid_int_data$H3k9ac$trial=="T3"), cols,parameters_int_data,paste("H3k9ac ",labelling_int_data,sep=""),"H3k9ac_T3_")
    
    plot_black_bg_barplot_timeline(spheroid_int_data$H3k9ac, cols,parameters_int_data,paste("H3k9ac ",labelling_int_data,sep=""),"H3k9ac_")
    
    
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
    
    
    plot_black_bg_histogram_timeline(subset(spheroid_angle$H3k9ac,spheroid_angle$H3k9ac$trial=="T1"), cols,parameters_angle_data,paste("H3k9ac ",labelling_angle_data,sep=""),"H3k9ac_T1_")
    plot_black_bg_histogram_timeline(subset(spheroid_angle$H3k9ac,spheroid_angle$H3k9ac$trial=="T2"), cols,parameters_angle_data,paste("H3k9ac ",labelling_angle_data,sep=""),"H3k9ac_T2_")
    plot_black_bg_histogram_timeline(subset(spheroid_angle$H3k9ac,spheroid_angle$H3k9ac$trial=="T3"), cols,parameters_angle_data,paste("H3k9ac ",labelling_angle_data,sep=""),"H3k9ac_T3_")
    
    plot_black_bg_boxplot_timeline(subset(spheroid_angle$H3k9ac,spheroid_angle$H3k9ac$trial=="T1"), cols,parameters_angle_data,paste("H3k9ac ",labelling_angle_data,sep=""),"H3k9ac_T1_")
    plot_black_bg_boxplot_timeline(subset(spheroid_angle$H3k9ac,spheroid_angle$H3k9ac$trial=="T2"), cols,parameters_angle_data,paste("H3k9ac ",labelling_angle_data,sep=""),"H3k9ac_T2_")
    plot_black_bg_boxplot_timeline(subset(spheroid_angle$H3k9ac,spheroid_angle$H3k9ac$trial=="T3"), cols,parameters_angle_data,paste("H3k9ac ",labelling_angle_data,sep=""),"H3k9ac_T3_")
    
    plot_black_bg_barplot_timeline(spheroid_angle$H3k9ac, cols,parameters_angle_data,paste("H3k9ac ",labelling_angle_data,sep=""),"H3k9ac_")
    
    
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
    
    plot_black_bg_histogram_timeline(subset(spheroid_axial$H3k9ac,spheroid_axial$H3k9ac$trial=="T1"), cols,parameters_axial_data,paste("H3k9ac ",labelling_axial_data,sep=""),"H3k9ac_T1_")
    plot_black_bg_histogram_timeline(subset(spheroid_axial$H3k9ac,spheroid_axial$H3k9ac$trial=="T2"), cols,parameters_axial_data,paste("H3k9ac ",labelling_axial_data,sep=""),"H3k9ac_T2_")
    plot_black_bg_histogram_timeline(subset(spheroid_axial$H3k9ac,spheroid_axial$H3k9ac$trial=="T3"), cols,parameters_axial_data,paste("H3k9ac ",labelling_axial_data,sep=""),"H3k9ac_T3_")
    
    plot_black_bg_boxplot_timeline(subset(spheroid_axial$H3k9ac,spheroid_axial$H3k9ac$trial=="T1"), cols,parameters_axial_data,paste("H3k9ac ",labelling_axial_data,sep=""),"H3k9ac_T1_")
    plot_black_bg_boxplot_timeline(subset(spheroid_axial$H3k9ac,spheroid_axial$H3k9ac$trial=="T2"), cols,parameters_axial_data,paste("H3k9ac ",labelling_axial_data,sep=""),"H3k9ac_T2_")
    plot_black_bg_boxplot_timeline(subset(spheroid_axial$H3k9ac,spheroid_axial$H3k9ac$trial=="T3"), cols,parameters_axial_data,paste("H3k9ac ",labelling_axial_data,sep=""),"H3k9ac_T3_")
    
    plot_black_bg_barplot_timeline(spheroid_axial$H3k9ac, cols,parameters_axial_data,paste("H3k9ac ",labelling_axial_data,sep=""),"H3k9ac_")
    
    
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
    
    plot_black_bg_histogram_timeline(subset(spheroid_radial$H3k9ac,spheroid_radial$H3k9ac$trial=="T1"), cols,parameters_radial_data,paste("H3k9ac ",labelling_radial_data,sep=""),"H3k9ac_T1_")
    plot_black_bg_histogram_timeline(subset(spheroid_radial$H3k9ac,spheroid_radial$H3k9ac$trial=="T2"), cols,parameters_radial_data,paste("H3k9ac ",labelling_radial_data,sep=""),"H3k9ac_T2_")
    plot_black_bg_histogram_timeline(subset(spheroid_radial$H3k9ac,spheroid_radial$H3k9ac$trial=="T3"), cols,parameters_radial_data,paste("H3k9ac ",labelling_radial_data,sep=""),"H3k9ac_T3_")
    
    plot_black_bg_boxplot_timeline(subset(spheroid_radial$H3k9ac,spheroid_radial$H3k9ac$trial=="T1"), cols,parameters_radial_data,paste("H3k9ac ",labelling_radial_data,sep=""),"H3k9ac_T1_")
    plot_black_bg_boxplot_timeline(subset(spheroid_radial$H3k9ac,spheroid_radial$H3k9ac$trial=="T2"), cols,parameters_radial_data,paste("H3k9ac ",labelling_radial_data,sep=""),"H3k9ac_T2_")
    plot_black_bg_boxplot_timeline(subset(spheroid_radial$H3k9ac,spheroid_radial$H3k9ac$trial=="T3"), cols,parameters_radial_data,paste("H3k9ac ",labelling_radial_data,sep=""),"H3k9ac_T3_")
    
    plot_black_bg_barplot_timeline(spheroid_radial$H3k9ac, cols,parameters_radial_data,paste("H3k9ac ",labelling_radial_data,sep=""),"H3k9ac_")
    
    
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
    
    plot_black_bg_histogram_timeline(subset(nuclei_2dint_data$H3K9AC,nuclei_2dint_data$H3K9AC$trial=="T1"), cols,parameters_2d_int_data,paste("H3K9AC ",nuc_labelling_2d_int_data,sep=""),"H3K9AC_T1_")
    plot_black_bg_histogram_timeline(subset(nuclei_2dint_data$H3K9AC,nuclei_2dint_data$H3K9AC$trial=="T2"), cols,parameters_2d_int_data,paste("H3K9AC ",nuc_labelling_2d_int_data,sep=""),"H3K9AC_T2_")
    plot_black_bg_histogram_timeline(subset(nuclei_2dint_data$H3K9AC,nuclei_2dint_data$H3K9AC$trial=="T3"), cols,parameters_2d_int_data,paste("H3K9AC ",nuc_labelling_2d_int_data,sep=""),"H3K9AC_T3_")
    
    plot_black_bg_boxplot_timeline(subset(nuclei_2dint_data$H3K9AC,nuclei_2dint_data$H3K9AC$trial=="T1"), cols,parameters_2d_int_data,paste("H3K9AC ",nuc_labelling_2d_int_data,sep=""),"H3K9AC_T1_")
    plot_black_bg_boxplot_timeline(subset(nuclei_2dint_data$H3K9AC,nuclei_2dint_data$H3K9AC$trial=="T2"), cols,parameters_2d_int_data,paste("H3K9AC ",nuc_labelling_2d_int_data,sep=""),"H3K9AC_T2_")
    plot_black_bg_boxplot_timeline(subset(nuclei_2dint_data$H3K9AC,nuclei_2dint_data$H3K9AC$trial=="T3"), cols,parameters_2d_int_data,paste("H3K9AC ",nuc_labelling_2d_int_data,sep=""),"H3K9AC_T3_")
    
    plot_black_bg_barplot_timeline(nuclei_2dint_data$H3K9AC, cols,parameters_2d_int_data,paste("H3K9AC ",nuc_labelling_2d_int_data,sep=""),"H3K9AC_")
    
    
    
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
    
    plot_black_bg_histogram_timeline(subset(nuclei_3dint_data$H3K9AC,nuclei_3dint_data$H3K9AC$trial=="T1"), cols,parameters_int_data,paste("H3K9AC ",nuclabelling_int_data,sep=""),"H3K9AC_T1_")
    plot_black_bg_histogram_timeline(subset(nuclei_3dint_data$H3K9AC,nuclei_3dint_data$H3K9AC$trial=="T2"), cols,parameters_int_data,paste("H3K9AC ",nuclabelling_int_data,sep=""),"H3K9AC_T2_")
    plot_black_bg_histogram_timeline(subset(nuclei_3dint_data$H3K9AC,nuclei_3dint_data$H3K9AC$trial=="T3"), cols,parameters_int_data,paste("H3K9AC ",nuclabelling_int_data,sep=""),"H3K9AC_T3_")
    
    plot_black_bg_boxplot_timeline(subset(nuclei_3dint_data$H3K9AC,nuclei_3dint_data$H3K9AC$trial=="T1"), cols,parameters_int_data,paste("H3K9AC ",nuclabelling_int_data,sep=""),"H3K9AC_T1_")
    plot_black_bg_boxplot_timeline(subset(nuclei_3dint_data$H3K9AC,nuclei_3dint_data$H3K9AC$trial=="T2"), cols,parameters_int_data,paste("H3K9AC ",nuclabelling_int_data,sep=""),"H3K9AC_T2_")
    plot_black_bg_boxplot_timeline(subset(nuclei_3dint_data$H3K9AC,nuclei_3dint_data$H3K9AC$trial=="T3"), cols,parameters_int_data,paste("H3K9AC ",nuclabelling_int_data,sep=""),"H3K9AC_T3_")
    
    plot_black_bg_barplot_timeline(nuclei_3dint_data$H3K9AC, cols,parameters_int_data,paste("H3K9AC ",nuclabelling_int_data,sep=""),"H3K9AC_")
    
    
    
  }
  
  dird<-"E:/HMF3A_reprogramming/R_analysis/H3k9ac/white_bg_plots/"
  dir.create(dird)
  setwd(dird)
  cols<-colorRampPalette(c("red","yellow"))( 5) 
  
  #Spheoid Proj_Int
  {
    
    dire<-paste(dird,"spheroid_proj_int/",sep="")
    dir.create(dire)
    setwd(dire)
    
    
    plot_white_bg_histogram_timeline(subset(spheroid_2dint_data$H3k9ac,spheroid_2dint_data$H3k9ac$trial=="T1"), cols,parameters_2d_int_data,paste("H3k9ac ",labelling_2d_int_data,sep=""),"H3k9ac_T1_")
    plot_white_bg_histogram_timeline(subset(spheroid_2dint_data$H3k9ac,spheroid_2dint_data$H3k9ac$trial=="T2"), cols,parameters_2d_int_data,paste("H3k9ac ",labelling_2d_int_data,sep=""),"H3k9ac_T2_")
    plot_white_bg_histogram_timeline(subset(spheroid_2dint_data$H3k9ac,spheroid_2dint_data$H3k9ac$trial=="T3"), cols,parameters_2d_int_data,paste("H3k9ac ",labelling_2d_int_data,sep=""),"H3k9ac_T3_")
    
    plot_white_bg_boxplot_timeline(subset(spheroid_2dint_data$H3k9ac,spheroid_2dint_data$H3k9ac$trial=="T1"), cols,parameters_2d_int_data,paste("H3k9ac ",labelling_2d_int_data,sep=""),"H3k9ac_T1_")
    plot_white_bg_boxplot_timeline(subset(spheroid_2dint_data$H3k9ac,spheroid_2dint_data$H3k9ac$trial=="T2"), cols,parameters_2d_int_data,paste("H3k9ac ",labelling_2d_int_data,sep=""),"H3k9ac_T2_")
    plot_white_bg_boxplot_timeline(subset(spheroid_2dint_data$H3k9ac,spheroid_2dint_data$H3k9ac$trial=="T3"), cols,parameters_2d_int_data,paste("H3k9ac ",labelling_2d_int_data,sep=""),"H3k9ac_T3_")
    
    plot_white_bg_barplot_timeline(spheroid_2dint_data$H3k9ac, cols,parameters_2d_int_data,paste("H3k9ac ",labelling_2d_int_data,sep=""),"H3k9ac_")
    
    
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
    
    plot_white_bg_histogram_timeline(subset(spheroid_int_data$H3k9ac,spheroid_int_data$H3k9ac$trial=="T1"), cols,parameters_int_data,paste("H3k9ac ",labelling_int_data,sep=""),"H3k9ac_T1_")
    plot_white_bg_histogram_timeline(subset(spheroid_int_data$H3k9ac,spheroid_int_data$H3k9ac$trial=="T2"), cols,parameters_int_data,paste("H3k9ac ",labelling_int_data,sep=""),"H3k9ac_T2_")
    plot_white_bg_histogram_timeline(subset(spheroid_int_data$H3k9ac,spheroid_int_data$H3k9ac$trial=="T3"), cols,parameters_int_data,paste("H3k9ac ",labelling_int_data,sep=""),"H3k9ac_T3_")
    
    plot_white_bg_boxplot_timeline(subset(spheroid_int_data$H3k9ac,spheroid_int_data$H3k9ac$trial=="T1"), cols,parameters_int_data,paste("H3k9ac ",labelling_int_data,sep=""),"H3k9ac_T1_")
    plot_white_bg_boxplot_timeline(subset(spheroid_int_data$H3k9ac,spheroid_int_data$H3k9ac$trial=="T2"), cols,parameters_int_data,paste("H3k9ac ",labelling_int_data,sep=""),"H3k9ac_T2_")
    plot_white_bg_boxplot_timeline(subset(spheroid_int_data$H3k9ac,spheroid_int_data$H3k9ac$trial=="T3"), cols,parameters_int_data,paste("H3k9ac ",labelling_int_data,sep=""),"H3k9ac_T3_")
    
    plot_white_bg_barplot_timeline(spheroid_int_data$H3k9ac, cols,parameters_int_data,paste("H3k9ac ",labelling_int_data,sep=""),"H3k9ac_")
    
    
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
    
    
    plot_white_bg_histogram_timeline(subset(spheroid_angle$H3k9ac,spheroid_angle$H3k9ac$trial=="T1"), cols,parameters_angle_data,paste("H3k9ac ",labelling_angle_data,sep=""),"H3k9ac_T1_")
    plot_white_bg_histogram_timeline(subset(spheroid_angle$H3k9ac,spheroid_angle$H3k9ac$trial=="T2"), cols,parameters_angle_data,paste("H3k9ac ",labelling_angle_data,sep=""),"H3k9ac_T2_")
    plot_white_bg_histogram_timeline(subset(spheroid_angle$H3k9ac,spheroid_angle$H3k9ac$trial=="T3"), cols,parameters_angle_data,paste("H3k9ac ",labelling_angle_data,sep=""),"H3k9ac_T3_")
    
    plot_white_bg_boxplot_timeline(subset(spheroid_angle$H3k9ac,spheroid_angle$H3k9ac$trial=="T1"), cols,parameters_angle_data,paste("H3k9ac ",labelling_angle_data,sep=""),"H3k9ac_T1_")
    plot_white_bg_boxplot_timeline(subset(spheroid_angle$H3k9ac,spheroid_angle$H3k9ac$trial=="T2"), cols,parameters_angle_data,paste("H3k9ac ",labelling_angle_data,sep=""),"H3k9ac_T2_")
    plot_white_bg_boxplot_timeline(subset(spheroid_angle$H3k9ac,spheroid_angle$H3k9ac$trial=="T3"), cols,parameters_angle_data,paste("H3k9ac ",labelling_angle_data,sep=""),"H3k9ac_T3_")
    
    plot_white_bg_barplot_timeline(spheroid_angle$H3k9ac, cols,parameters_angle_data,paste("H3k9ac ",labelling_angle_data,sep=""),"H3k9ac_")
    
    
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
    
    plot_white_bg_histogram_timeline(subset(spheroid_axial$H3k9ac,spheroid_axial$H3k9ac$trial=="T1"), cols,parameters_axial_data,paste("H3k9ac ",labelling_axial_data,sep=""),"H3k9ac_T1_")
    plot_white_bg_histogram_timeline(subset(spheroid_axial$H3k9ac,spheroid_axial$H3k9ac$trial=="T2"), cols,parameters_axial_data,paste("H3k9ac ",labelling_axial_data,sep=""),"H3k9ac_T2_")
    plot_white_bg_histogram_timeline(subset(spheroid_axial$H3k9ac,spheroid_axial$H3k9ac$trial=="T3"), cols,parameters_axial_data,paste("H3k9ac ",labelling_axial_data,sep=""),"H3k9ac_T3_")
    
    plot_white_bg_boxplot_timeline(subset(spheroid_axial$H3k9ac,spheroid_axial$H3k9ac$trial=="T1"), cols,parameters_axial_data,paste("H3k9ac ",labelling_axial_data,sep=""),"H3k9ac_T1_")
    plot_white_bg_boxplot_timeline(subset(spheroid_axial$H3k9ac,spheroid_axial$H3k9ac$trial=="T2"), cols,parameters_axial_data,paste("H3k9ac ",labelling_axial_data,sep=""),"H3k9ac_T2_")
    plot_white_bg_boxplot_timeline(subset(spheroid_axial$H3k9ac,spheroid_axial$H3k9ac$trial=="T3"), cols,parameters_axial_data,paste("H3k9ac ",labelling_axial_data,sep=""),"H3k9ac_T3_")
    
    plot_white_bg_barplot_timeline(spheroid_axial$H3k9ac, cols,parameters_axial_data,paste("H3k9ac ",labelling_axial_data,sep=""),"H3k9ac_")
    
    
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
    
    plot_white_bg_histogram_timeline(subset(spheroid_radial$H3k9ac,spheroid_radial$H3k9ac$trial=="T1"), cols,parameters_radial_data,paste("H3k9ac ",labelling_radial_data,sep=""),"H3k9ac_T1_")
    plot_white_bg_histogram_timeline(subset(spheroid_radial$H3k9ac,spheroid_radial$H3k9ac$trial=="T2"), cols,parameters_radial_data,paste("H3k9ac ",labelling_radial_data,sep=""),"H3k9ac_T2_")
    plot_white_bg_histogram_timeline(subset(spheroid_radial$H3k9ac,spheroid_radial$H3k9ac$trial=="T3"), cols,parameters_radial_data,paste("H3k9ac ",labelling_radial_data,sep=""),"H3k9ac_T3_")
    
    plot_white_bg_boxplot_timeline(subset(spheroid_radial$H3k9ac,spheroid_radial$H3k9ac$trial=="T1"), cols,parameters_radial_data,paste("H3k9ac ",labelling_radial_data,sep=""),"H3k9ac_T1_")
    plot_white_bg_boxplot_timeline(subset(spheroid_radial$H3k9ac,spheroid_radial$H3k9ac$trial=="T2"), cols,parameters_radial_data,paste("H3k9ac ",labelling_radial_data,sep=""),"H3k9ac_T2_")
    plot_white_bg_boxplot_timeline(subset(spheroid_radial$H3k9ac,spheroid_radial$H3k9ac$trial=="T3"), cols,parameters_radial_data,paste("H3k9ac ",labelling_radial_data,sep=""),"H3k9ac_T3_")
    
    plot_white_bg_barplot_timeline(spheroid_radial$H3k9ac, cols,parameters_radial_data,paste("H3k9ac ",labelling_radial_data,sep=""),"H3k9ac_")
    
    
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
    
    plot_white_bg_histogram_timeline(subset(nuclei_2dint_data$H3K9AC,nuclei_2dint_data$H3K9AC$trial=="T1"), cols,parameters_2d_int_data,paste("H3K9AC ",nuc_labelling_2d_int_data,sep=""),"H3K9AC_T1_")
    plot_white_bg_histogram_timeline(subset(nuclei_2dint_data$H3K9AC,nuclei_2dint_data$H3K9AC$trial=="T2"), cols,parameters_2d_int_data,paste("H3K9AC ",nuc_labelling_2d_int_data,sep=""),"H3K9AC_T2_")
    plot_white_bg_histogram_timeline(subset(nuclei_2dint_data$H3K9AC,nuclei_2dint_data$H3K9AC$trial=="T3"), cols,parameters_2d_int_data,paste("H3K9AC ",nuc_labelling_2d_int_data,sep=""),"H3K9AC_T3_")
    
    plot_white_bg_boxplot_timeline(subset(nuclei_2dint_data$H3K9AC,nuclei_2dint_data$H3K9AC$trial=="T1"), cols,parameters_2d_int_data,paste("H3K9AC ",nuc_labelling_2d_int_data,sep=""),"H3K9AC_T1_")
    plot_white_bg_boxplot_timeline(subset(nuclei_2dint_data$H3K9AC,nuclei_2dint_data$H3K9AC$trial=="T2"), cols,parameters_2d_int_data,paste("H3K9AC ",nuc_labelling_2d_int_data,sep=""),"H3K9AC_T2_")
    plot_white_bg_boxplot_timeline(subset(nuclei_2dint_data$H3K9AC,nuclei_2dint_data$H3K9AC$trial=="T3"), cols,parameters_2d_int_data,paste("H3K9AC ",nuc_labelling_2d_int_data,sep=""),"H3K9AC_T3_")
    
    plot_white_bg_barplot_timeline(nuclei_2dint_data$H3K9AC, cols,parameters_2d_int_data,paste("H3K9AC ",nuc_labelling_2d_int_data,sep=""),"H3K9AC_")
    
    
    
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
    
    plot_white_bg_histogram_timeline(subset(nuclei_3dint_data$H3K9AC,nuclei_3dint_data$H3K9AC$trial=="T1"), cols,parameters_int_data,paste("H3K9AC ",nuclabelling_int_data,sep=""),"H3K9AC_T1_")
    plot_white_bg_histogram_timeline(subset(nuclei_3dint_data$H3K9AC,nuclei_3dint_data$H3K9AC$trial=="T2"), cols,parameters_int_data,paste("H3K9AC ",nuclabelling_int_data,sep=""),"H3K9AC_T2_")
    plot_white_bg_histogram_timeline(subset(nuclei_3dint_data$H3K9AC,nuclei_3dint_data$H3K9AC$trial=="T3"), cols,parameters_int_data,paste("H3K9AC ",nuclabelling_int_data,sep=""),"H3K9AC_T3_")
    
    plot_white_bg_boxplot_timeline(subset(nuclei_3dint_data$H3K9AC,nuclei_3dint_data$H3K9AC$trial=="T1"), cols,parameters_int_data,paste("H3K9AC ",nuclabelling_int_data,sep=""),"H3K9AC_T1_")
    plot_white_bg_boxplot_timeline(subset(nuclei_3dint_data$H3K9AC,nuclei_3dint_data$H3K9AC$trial=="T2"), cols,parameters_int_data,paste("H3K9AC ",nuclabelling_int_data,sep=""),"H3K9AC_T2_")
    plot_white_bg_boxplot_timeline(subset(nuclei_3dint_data$H3K9AC,nuclei_3dint_data$H3K9AC$trial=="T3"), cols,parameters_int_data,paste("H3K9AC ",nuclabelling_int_data,sep=""),"H3K9AC_T3_")
    
    plot_white_bg_barplot_timeline(nuclei_3dint_data$H3K9AC, cols,parameters_int_data,paste("H3K9AC ",nuclabelling_int_data,sep=""),"H3K9AC_")
    
    
    
  }
}

#write the collected data
{
  dir.create("E:/HMF3A_reprogramming/R_analysis/H3K9AC/combined_data/")
  setwd("E:/HMF3A_reprogramming/R_analysis/H3K9AC/combined_data/")
  
  write.csv(spheroid_combined, file="spheroid_combined.csv")
  write.csv(nucleus_combined, file="nucleus_combined.csv")
  
  write.csv(spheroid_radial$H3k9ac, file="spheroid_radial_dist_H3K9AC.csv")
  
  write.csv(spheroid_2dint_data$H3k9ac, file="spheroid_2d_int_H3K9AC.csv")
  
  write.csv(spheroid_axial$H3k9ac, file="spheroid_axial_dist_H3K9AC.csv")
  write.csv(spheroid_axial$dna, file="spheroid_axial_dist_dna.csv")
  
  write.csv(spheroid_angle$H3k9ac, file="spheroid_angle_dist_H3K9AC.csv")
  write.csv(spheroid_angle$dna, file="spheroid_angle_dist_dna.csv")
  
  write.csv(spheroid_int_data$H3k9ac, file="spheroid_int_H3K9AC.csv")
  write.csv(spheroid_int_data$dna, file="spheroid_int_dna.csv")
  
  write.csv(nuclei_2dint_data$H3K9AC, file="nuclear_2d_int_H3K9AC.csv")
  
  write.csv(nuclei_3dint_data$H3K9AC, file="nuclear_int_H3K9AC.csv")
  write.csv(nuclei_3dint_data$dna, file="nuclear_int_dna.csv")
  
  
}