/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////  SCRIPT WITH A FUNCTION TO RUN ALL THE SPHEREOID ANALYSIS PROGRAMS       											                   				   	
///////  WRITTEN BY: SARADHA VENKATACHALAPATHY                                                                                                                     
///////  ASSUMPTIONS: All the macros are saved in one folder
///////  DESCRIPTION: This script contains a function which accepts the path to the macro folder, the experimental folder, samples, channel information, cytoplasm channel,
///////					number of channels to be processed and the seperating character
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

function run_spheroid_analysis(path_to_spheroid_analysis_macros_folder,path_to_experimental_folder, samples,channel_data,cytoplasm_channel,num_ch,sep) {
    //define path to the macros
    move_raw=path_to_spheroid_analysis_macros_folder+"move_files_rawimages.ijm"
	vis_zproj=path_to_spheroid_analysis_macros_folder+"spheroid_zproject_raw.ijm"
	vis_3d=path_to_spheroid_analysis_macros_folder+"project_in3d.ijm"
	id_obj_3d=path_to_spheroid_analysis_macros_folder+"3d_objects_spheroid.ijm"
	meas_2d=path_to_spheroid_analysis_macros_folder+"spheroid_zproject_2d_measure.ijm"
	meas_3d=path_to_spheroid_analysis_macros_folder+"measures_spheroid.ijm"
	meas_shells=path_to_spheroid_analysis_macros_folder+"spheroid_shells.ijm"
	meas_zaxis=path_to_spheroid_analysis_macros_folder+"Spheroid_zaxis_profile.ijm"
	meas_angles=path_to_spheroid_analysis_macros_folder+"spheroid_angles.ijm"


      for (i = 0; i < samples.length; i++) { // loop through the samples
			path_to_experimental_subfolder=path_to_experimental_folder+samples[i]+"\\"; 	//set path to subfolder.
			runMacro(move_raw,path_to_experimental_subfolder)								//move raw images to subfolder
			runMacro(vis_zproj,path_to_experimental_subfolder+sep+channel_data)				//zproject and make a montage
			runMacro(vis_3d,path_to_experimental_subfolder+sep+num_ch)						//3d project the images
				
			runMacro(id_obj_3d,path_to_experimental_subfolder+sep+channel_data+sep+cytoplasm_channel)	//generate 3d objects of spheroid
			runMacro(meas_2d,path_to_experimental_subfolder+sep+channel_data)						//measure 2d properties
			runMacro(meas_3d,path_to_experimental_subfolder+sep+channel_data)						//measure 3d properties
			
			runMacro(meas_shells,path_to_experimental_subfolder+sep+channel_data)			//compute the radial distribution 
			runMacro(meas_zaxis,path_to_experimental_subfolder+sep+channel_data)			//compute the axial distribbution
			runMacro(meas_angles,path_to_experimental_subfolder+sep+channel_data)			//compute the angular frequency distribution 

   }


}

// set the arguments
samples=newArray("data", "dat1"); 		//samples
sep=" ";								//separating charcter
channel_data="4 OCT4 LAMINAC ACTIN";	//channel information
cytoplasm_channel="4";					//cytoplasm channel
nch=4;									// number of channels to be analysed

run_spheroid_analysis("E:\\spheroid_analysis\\","E:\\temp\\",samples,channel_data,cytoplasm_channel,nch,sep); // run the script