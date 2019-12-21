/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////  SCRIPT TO QUANTIFY GEOMETRICAL AND INTENSITY FEATURES
///////  WRITTEN BY: SARADHA VENKATACHALAPATHY                                                                                                                   
///////  ASSUMPTIONS: The input image is a confocal zstack and the channel 1 contains the nucleus. The raw files' extension is ".nd2". Spheroid volumes range from 5,0000-1,000,0000 kcu.microns. Tested on ImageJ v1.52p.
///////  DESCRIPTION: The script accepts the directory to a folder containing rawimages folder and the number of channels in the raw image that need to be analysed and the taget of the stains. 
///////				  It creates 3-6 subfolders: one for storing the geoemtrical features of the identified spheroids and others for storing intensity features of the selected image.
///////				  For all images in the identifed spheroids, it names and measures the geometrical and intensity features that are requied.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



args = split(getArgument()," "); 	// read in the array of argurments that are separated by a space
dirsa = args[0];	// The first argument is the directory
nchannels=args[1];	// The second argument is the number of channels
ch2_name=args[2];	// The third argument is the name of the protein in channel 2
ch3_name=args[3];	// The fourth argument is the name of the protein in channel 3
ch4_name=args[4];	// The fifth argument is the name of the protein in channel 4


setBatchMode(true); // run in batch mode to save time and memory

dirsa1=dirsa+"rawimages"+File.separator; // define the path to the folder containing the unprocessed or raw images
dirb= dirsa + "3d obj sphereoid"+ File.separator;  // define the path to the folder containing the image afer 3d segmentation to indentify spheroids

//set the options for 3D ROI manager 
run("3D Manager Options", "volume surface compactness fit_ellipse 3d_moments integrated_density mean_grey_value std_dev_grey_value mode_grey_value feret minimum_grey_value maximum_grey_value centroid_(pix) centroid_(unit) distance_to_surface centre_of_mass_(pix) centre_of_mass_(unit) bounding_box radial_distance surface_contact closest use exclude_objects_on_edges_xy distance_between_centers=10 distance_max_contact=1.80");

///creating new subdirectories for storing the processed data//
Dir_res_geo_m=dirsa + "3D geometrical data spheroid"+ File.separator;  // define the path to a new folder that will contain the 3d geoemetrical data for the spheroids
Dir_res_int=dirsa + "3D int_data spheroid"+ File.separator;   // define the path to a new folder that will contain subfolders below
Dir_res_ch1_int=Dir_res_int + "DNA" + File.separator;  // define the path to a new folder that will contain the 3d intensity data for the spheroids in channel 1/ DNA channel
Dir_res_ch2_int=Dir_res_int + ch2_name + File.separator;  // define the path to a new folder that will contain the 3d intensity data for the spheroids in channel 2
Dir_res_ch3_int=Dir_res_int + ch3_name + File.separator;  // define the path to a new folder that will contain the 3d intensity data for the spheroids in channel 3
Dir_res_ch4_int=Dir_res_int + ch4_name + File.separator;  // define the path to a new folder that will contain the 3d intensity data for the spheroids in channel 4
// create the paths defined above 
File.makeDirectory(Dir_res_geo_m); 
File.makeDirectory(Dir_res_int); 
File.makeDirectory(Dir_res_ch1_int); 
if(nchannels>1){ 
	File.makeDirectory(Dir_res_ch2_int); 
}
if(nchannels>2){ 
	File.makeDirectory(Dir_res_ch3_int); 
}
if(nchannels>3){ 
	File.makeDirectory(Dir_res_ch4_int); 
}

filenames=getFileList(dirb); // get the list of names of the images in the 3d objects folder

run("3D Manager"); // Open the 3D ROI manager
for(f=0;f<filenames.length;f++){
	
	// open the image containing 3d objects 
	open(dirb+filenames[f]);  //define the path
	// extract the the name of the file or remove file extension
	obj_image=getTitle(); //get the title and assign it. it will be a character string
	baseNameEnd=indexOf(obj_image, ".tiff");  //find the index of the string at with ".tiff" first appears NOTE: change this to correct file extension if working with other files 
	baseName=substring(obj_image, 0, baseNameEnd);  // get the substring such that the file extention gets removed from the character string 

	Ext.Manager3D_AddImage(); // add the segmented image
	Ext.Manager3D_Count(nb_obj);  // get the number of segmented rois

	name=newArray(nb_obj);//set an array for names

	for(p=0;p<nb_obj;p++){ // loop through all objects
		Ext.Manager3D_Bounding3D(p,x0,x1,y0,y1,z0,z1); // for the pth object, get the cordinates of the bounding rectangle
		Ext.Manager3D_MonoSelect(); // select only one object
		Ext.Manager3D_Select(p);	// select the pth object 
		name[p]=baseName+"_spheroid_"+x0+"_"+y0+"_"+z0+"_"+z1;  // set the name of the object and spheroid identifier
	}
	Ext.Manager3D_DeselectAll(); // deselct all
	Ext.Manager3D_Reset();
	run("Close All");
		
	// open the image containing 3d objects 
	open(dirb+filenames[f]);  //define the path
		run("Bin...", "x=4 y=4 z=1 bin=Average");
	
	Ext.Manager3D_AddImage(); // add the segmented image
	Ext.Manager3D_Count(nb_obj);  // get the number of segmented rois

	for(p=0;p<nb_obj;p++){ // loop through all objects
		Ext.Manager3D_MonoSelect(); // select only one object
		Ext.Manager3D_Select(p);	// select the pth object 
		Ext.Manager3D_Rename(name[p]); // rename the object
	}
	Ext.Manager3D_DeselectAll(); // deselct all
	
	Ext.Manager3D_SelectAll(); // select all objects
	Ext.Manager3D_Measure();   // measure geoemtric data
    Ext.Manager3D_SaveMeasure(Dir_res_geo_m+baseName+"_geometric.tsv");	//save the results
    Ext.Manager3D_CloseResult("M");	// close the results window
	run("Close All");	// close all open images

	name=baseName+".nd2"; // set the name with the proper extension
	path=dirsa1+name;	// set path
	run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT c_begin=1 c_end=1 c_step=1"); // open the first channel 
	run("Bin...", "x=4 y=4 z=1 bin=Average");
	Ext.Manager3D_Quantif(); // measure intensity 
	Ext.Manager3D_SaveQuantif(Dir_res_ch1_int+baseName+".tsv"); // save the results
	Ext.Manager3D_CloseResult("Q"); // close the results window
	run("Close All"); // close all open images

	if(nchannels>1){ // if channel 2 is needed
		run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT c_begin=2 c_end=2 c_step=1"); // open the second channel 
		run("Bin...", "x=4 y=4 z=1 bin=Average");
		Ext.Manager3D_Quantif(); // measure intensity 
		Ext.Manager3D_SaveQuantif(Dir_res_ch2_int+baseName+".tsv"); // save the results
		Ext.Manager3D_CloseResult("Q"); // close the results window
		run("Close All"); // close all open images
	}

	if(nchannels>2){ // if channel 3 is needed
		run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT c_begin=3 c_end=3 c_step=1"); // open the third channel 
		run("Bin...", "x=4 y=4 z=1 bin=Average");
		Ext.Manager3D_Quantif(); // measure intensity 
		Ext.Manager3D_SaveQuantif(Dir_res_ch3_int+baseName+".tsv"); // save the results
		Ext.Manager3D_CloseResult("Q"); // close the results window
		run("Close All"); // close all open images
	}
	
	if(nchannels>3){ // if channel 4 is needed
		run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT c_begin=4 c_end=4 c_step=1"); // open the fourth channel 
		run("Bin...", "x=4 y=4 z=1 bin=Average");
		Ext.Manager3D_Quantif(); // measure intensity 
		Ext.Manager3D_SaveQuantif(Dir_res_ch4_int+baseName+".tsv"); // save the results
		Ext.Manager3D_CloseResult("Q"); // close the results window
		run("Close All"); // close all open images
	}
	
	Ext.Manager3D_Reset();
	
	run("Collect Garbage"); // collect garbage
	run("Collect Garbage"); // collect garbage
	run("Collect Garbage"); // collect garbage

}
Ext.Manager3D_Close(); // close the 3D roi manager
