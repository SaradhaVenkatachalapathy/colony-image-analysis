/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////  SCRIPT TO IDENTIFY BINARIZE NUCLEI IN 3D      											                    
///////  WRITTEN BY: SARADHA VENKATACHALAPATHY                                                                                                                   
///////  ASSUMPTIONS: The input image is a confocal zstack and the channel 1 contains the nucleus 
///////  DESCRIPTION: The program opens nucleus channel, Enhance contrast, filter, threshold, binarizes the image. This is followed by computing the 3D distance map, threholding,
///////				  binarize and watershed. The saves the image in "after watershed" folder.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


dir = getArgument() // read in the argument

setBatchMode(true); // run in batch mode to save time and memory

// set paths
dirsa= dir + "rawimages"+ File.separator;
filenamesa=getFileList(dirsa);

dirw= dir + "after watershed"+ File.separator;
File.makeDirectory(dirw); 

//set the options for 3D objects counter and 3D ROI manager 

run("3D OC Options", "volume surface nb_of_obj._voxels nb_of_surf._voxels integrated_density mean_gray_value std_dev_gray_value median_gray_value minimum_gray_value maximum_gray_value centroid mean_distance_to_surface std_dev_distance_to_surface median_distance_to_surface centre_of_mass bounding_box dots_size=5 font_size=20 redirect_to=none");
run("3D Manager Options", "volume surface compactness fit_ellipse 3d_moments integrated_density mean_grey_value std_dev_grey_value mode_grey_value feret minimum_grey_value maximum_grey_value centroid_(pix) centroid_(unit) distance_to_surface centre_of_mass_(pix) centre_of_mass_(unit) bounding_box radial_distance surface_contact closest use distance_between_centers=10 distance_max_contact=1.80");


for(f=0;f<filenamesa.length;f++){ //loop through all images

	//open first channel of the raw image
	path=dirsa+filenamesa[f];
	run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT c_begin=1 c_end=1 c_step=1");
	imgName=getTitle(); 
	baseNameEnd=indexOf(imgName, ".nd2"); 
	baseName=substring(imgName, 0, baseNameEnd); 

	//Enhance contrast, filter, threshold, binarize and fill holes 
	run("Enhance Contrast...", "saturated=0.01 normalize process_all use");
	run("Kuwahara Filter", "sampling=5 stack");
	setAutoThreshold("Otsu dark stack");
	run("Make Binary", "method=Default background=Default");
	run("Invert LUT");
	//run("Erode", "stack");
	run("Erode", "stack");i
	run("Fill Holes", "stack");

	//Get the 3D distance map, then threhold, binarize and watershed
	run("3D Distance Map", "map=EDT image=baseName mask=Same threshold=1");
	setAutoThreshold("Otsu dark stack");
	run("Make Binary", "method=Default background=Default");
	run("Invert LUT");
	run("Erode", "stack");
	run("Watershed Irregular Features", "erosion=5 convexity_threshold=0.8 separator_size=0-15 stack");

	//Save the image
	saveAs("Tiff", dirw+baseName+".tiff"); 
	run("Close All");
}


