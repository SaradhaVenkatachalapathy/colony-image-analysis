/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////  SCRIPT TO OBTAIN CYTOPLASMIC INTENSITIES OF PROTEINS AS DEFINED BY A RING OF 2 MICRONS
///////  WRITTEN BY: SARADHA VENKATACHALAPATHY                                                                                                                   
///////  ASSUMPTIONS: The input image is a confocal zstack and the channel 1 contains the nucleus. 
///////  DESCRIPTION: Opens the image with 3D nuclear objects, dilate each nucleus upto 2 microns from its edge. Then it opens raw image, measure and save the results and images of each protein channel.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


args = split(getArgument()," "); 	// read in the array of argurments that are separated by a space
dirsa = args[0];	// The first argument is the directory
nchannels=args[1];	// The second argument is the number of channels
ch2_name=args[2];	// The third argument is the name of the protein in channel 2
ch3_name=args[3];	// The fourth argument is the name of the protein in channel 3
ch4_name=args[4];	// The fifth argument is the name of the protein in channel 4

//Set paths
dirb= dirsa + "3d obj nuclei"+ File.separator;
dirraw= dirsa + "rawimages"+ File.separator;

//make new directories for the images and results
dirbc= dirsa + "3d_objects_cell_2microns"+ File.separator;
dirc= dirsa + "cell_2microns_measure"+ File.separator;
dirc1=dirc + "DNA" +File.separator;  
dirc2=dirc + ch2_name +File.separator;  
dirc3=dirc + ch3_name +File.separator;  
dirc4=dirc + ch4_name +File.separator; 
File.makeDirectory(dirbc); 
File.makeDirectory(dirc);
File.makeDirectory(dirc1);
if(nchannels>1){ 
	File.makeDirectory(dirc2); 
	
}
if(nchannels>2){ 
	File.makeDirectory(dirc3); 
	
}
if(nchannels>3){ 
	File.makeDirectory(dirc4); 
	
}

//set the options for 3D objects counter and 3D ROI manager 
run("3D OC Options", "volume surface nb_of_obj._voxels nb_of_surf._voxels integrated_density mean_gray_value std_dev_gray_value median_gray_value minimum_gray_value maximum_gray_value centroid mean_distance_to_surface std_dev_distance_to_surface median_distance_to_surface centre_of_mass bounding_box dots_size=5 font_size=20 redirect_to=none");
run("3D Manager Options", "volume surface compactness fit_ellipse 3d_moments integrated_density mean_grey_value std_dev_grey_value mode_grey_value feret minimum_grey_value maximum_grey_value centroid_(pix) centroid_(unit) distance_to_surface centre_of_mass_(pix) centre_of_mass_(unit) bounding_box radial_distance surface_contact closest use distance_between_centers=10 distance_max_contact=1.80");


setBatchMode(true); // run in batch mode to save time and memory

filenames=getFileList(dirb);
baseName=newArray(filenames.length);
run("3D Manager"); // Open the 3D ROI manager

for(f=0;f<filenames.length;f++){ //loop through all images

	//Open segmented image
	open(dirb+filenames[f]);
	imgName=getTitle(); 
	baseNameEnd=indexOf(imgName, ".tiff"); 
	baseName[f]=substring(imgName, 0, baseNameEnd); 

	// Add nuclei to the 3D mananger and obtain the right names and close all
	Ext.Manager3D_AddImage();
	Ext.Manager3D_Count(nb_obj); 
	names=newArray(nb_obj);
	print("\\Clear"); 

	for(p=0;p<nb_obj;p++){
		t=p+1;
		//figureout the name of the nucleus 
		Ext.Manager3D_Bounding3D(p,x0,x1,y0,y1,z0,z1);
		zlows=z0+1;
		zhighs=z1+1;
		wids=x1-x0;
		heigs=y1-y0;
		Ext.Manager3D_MonoSelect();
		Ext.Manager3D_Select(p);
		names[p]=baseName[f]+"_nucleus_"+x0+"_"+y0+"_"+z0+"_"+z1;
	}
	run("Close All");
	Ext.Manager3D_Reset();

	
	//Open segmented image
	open(dirb+filenames[f]);
	imgName=getTitle(); 
	baseNameEnd=indexOf(imgName, ".tiff"); 
	baseName[f]=substring(imgName, 0, baseNameEnd); 


	//bin the image so that the pixel size is the similar across xy and z
	getVoxelSize(x, y, h, micron);
	fx=round(h/x);	// find the number of pixels to bin in x and y
	run("Bin...", "x=fx y=fx z=1 bin=Max"); // bin the image
	
	// calculate the number of dilations to get 2 microns
	getPixelSize(unit, pixelWidth, pixelHeight);
	num_dia=round(2/pixelHeight);

	// Add nuclei to the 3D mananger
	Ext.Manager3D_AddImage();
	Ext.Manager3D_Count(nb_obj); 
	
	// dilate the nucleus upto 2 microns from edge
	for(p=0;p<nb_obj;p++){
		t=p+1;
		//figureout the name of the nucleus 
		Ext.Manager3D_MonoSelect();
		Ext.Manager3D_Select(p);
		Ext.Manager3D_Rename(names[p]);
		for(k=1; k<=num_dia;k++){
			run("Dilate (3D)", "iso=t");
		}
	}
	saveAs("tiff",dirbc+baseName[f]); //save image
	run("Close All");

	//open raw image, measure and save the fields
	Ext.Manager3D_SelectAll();
	path=dirraw+baseName[f]+".nd2";
	//channel1
	run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT c_begin=1 c_end=1 c_step=1");
	getVoxelSize(x, y, h, micron);
	fx=round(h/x);	// find the number of pixels to bin in x and y
	run("Bin...", "x=fx y=fx z=1 bin=Max"); // bin the image
	Ext.Manager3D_Quantif(); //measure intensity
	Ext.Manager3D_SaveQuantif(dirc1 +baseName[f]+".tsv"); //save results
	Ext.Manager3D_CloseResult("Q");
	run("Close All");
	Ext.Manager3D_SelectAll();

	//channel2
	if(nchannels>1){
	run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT c_begin=2 c_end=2 c_step=1");
	getVoxelSize(x, y, h, micron);
	fx=round(h/x);	// find the number of pixels to bin in x and y
	run("Bin...", "x=fx y=fx z=1 bin=Max"); // bin the image
	Ext.Manager3D_Quantif();
	Ext.Manager3D_SaveQuantif(dirc2 +baseName[f]+".tsv");
	Ext.Manager3D_CloseResult("Q");
	run("Close All");
	Ext.Manager3D_SelectAll();
	}

	//channel3
	if(nchannels>2){
	run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT c_begin=3 c_end=3 c_step=1");
	getVoxelSize(x, y, h, micron);
	fx=round(h/x);	// find the number of pixels to bin in x and y
	run("Bin...", "x=fx y=fx z=1 bin=Max"); // bin the image
	Ext.Manager3D_Quantif();
	Ext.Manager3D_SaveQuantif(dirc3 +baseName[f]+".tsv");
	Ext.Manager3D_CloseResult("Q");
	run("Close All");
	Ext.Manager3D_SelectAll();
	}

	//channel4
	if(nchannels>3){
	run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT c_begin=4 c_end=4 c_step=1");
	getVoxelSize(x, y, h, micron);
	fx=round(h/x);	// find the number of pixels to bin in x and y
	run("Bin...", "x=fx y=fx z=1 bin=Max"); // bin the image
	Ext.Manager3D_Quantif();
	Ext.Manager3D_SaveQuantif(dirc4 +baseName[f]+".tsv");
	Ext.Manager3D_CloseResult("Q");
	run("Close All");
	Ext.Manager3D_SelectAll();
	}
	
	//Clean up
    run("Close All");
		Ext.Manager3D_Reset();
	
   	run("Clear Results");

   	run("Collect Garbage");
	
}

Ext.Manager3D_Close();
