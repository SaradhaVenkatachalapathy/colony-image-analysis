
function run_nuclear_analysis(path_to_nuclear_analysis_macros_folder,path_to_experimental_folder, samples,channel_data,cytoplasm_channel,num_ch,sep) {
    //define path to the macros
    //move_raw=path_to_spheroid_analysis_macros_folder+"move_files_rawimages.ijm"
	binarize_nuc=path_to_nuclear_analysis_macros_folder+"binarize_the_nucleus.ijm"
	id_obj_3d=path_to_nuclear_analysis_macros_folder+"3dcrop_nuclei.ijm"
	
	meas_3d=path_to_nuclear_analysis_macros_folder+"measure_nucleus.ijm"
	meas_2d=path_to_nuclear_analysis_macros_folder+"projected_2dmeaure_nucleus.ijm"
	meas_ring=path_to_nuclear_analysis_macros_folder+"cell_ring.ijm"
	
	compaction=path_to_nuclear_analysis_macros_folder+"compaction_measures.ijm"
	edf=path_to_nuclear_analysis_macros_folder+"elliptic_fourier_descriptors.ijm.ijm"
	gclm=path_to_nuclear_analysis_macros_folder+"glcm_steps.ijm"
	
	viz_2d=path_to_nuclear_analysis_macros_folder+"zproject_and_padd_montage.ijm"

      for (i = 0; i < samples.length; i++) { // loop through the samples
			path_to_experimental_subfolder=path_to_experimental_folder+samples[i]+"\\"; 	//set path to subfolder.
			runMacro(move_raw,path_to_experimental_subfolder)								//move raw images to subfolder
			runMacro(binarize_nuc,path_to_experimental_subfolder)				//binarize nuclear
			runMacro(id_obj_3d,path_to_experimental_subfolder+sep+channel_data)					//identify 3d ojects
				
			runMacro(meas_2d,path_to_experimental_subfolder+sep+channel_data)						//measure 2d properties
			runMacro(meas_3d,path_to_experimental_subfolder+sep+channel_data)						//measure 3d properties
			
			runMacro(compaction,path_to_experimental_subfolder)			//compute the compaction
			runMacro(edf,path_to_experimental_subfolder)			//compute the axial distribbution
			runMacro(gclm,path_to_experimental_subfolder)			//compute the angular frequency distribution 
			runMacro(meas_ring,path_to_experimental_subfolder+sep+channel_data)						//measure 2d properties
			runMacro(viz_2d,path_to_experimental_subfolder+sep+channel_data+sep+path_to_nuclear_analysis_macros_folder)	//zproject and make a montage
   }
}

// set the arguments
samples=newArray("D0", "D2","D4","D6","D8"); 		//samples
sep=" ";								//separating charcter
channel_data="4 prot1 prot2 prot3";	//channel information
cytoplasm_channel="3";					//cytoplasm channel
nch=4;									// number of channels to be analysed

run_nuclear_analysis("E:\\nuclear_analysis\\","E:\\temp\\",samples,channel_data,cytoplasm_channel,nch,sep); 
run_nuclear_analysis("E:\\nuclear_analysis\\","E:\\temp\\",samples,channel_data,cytoplasm_channel,nch,sep); 
run_nuclear_analysis("E:\\nuclear_analysis\\","E:\\temp\\",samples,channel_data,cytoplasm_channel,nch,sep);


