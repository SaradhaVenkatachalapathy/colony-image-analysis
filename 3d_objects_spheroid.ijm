/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////  SCRIPT TO IDENTIFY, SEGMENT AND CROP 3D SPHEROIDS USING NUCLEAR STAIN 
///////  WRITTEN BY: SARADHA VENKATACHALAPATHY                                                                                                                   
///////  ASSUMPTIONS: The input image is a confocal zstack and the channel 1 contains the nucleus. The raw files' extension is ".nd2". Spheroid volumes range from 5,0000-1,000,0000 kcu.microns
///////  DESCRIPTION: The script accepts the directory to a folder containing rawimages folder and the number of channels in the raw image that need to be cropped. 
///////				  It creates 5 subfolders: one for storing the image afer 3d segmentation with identified spheroids and 4 for storing channels of the raw image cropped to represent each spheroid.
///////				  For all images in the rawimages folder, it smoothens, thresholds, binarizes and fill holes to obtain the binary image of a spheroid. Then identifies 3D objects of size ranging 
/////// 			  from 5,000-1,000,000 cu.microns. If there are spheroids in the image, then it saves the image containing the identified the objects, adds the images to the 3D roi manager. Then
///////				  for each identified object in the image, it gets the cordinates and crop the image along the 3D bounding box and thresholds to obtain the binary mask of the object. Then for each 
///////				  channel that is required, the raw image is cropped along the bounding box for each object and retains only sections of the image that are also present in the binary mask. This
///////				  ensures that in the cropped image there are not overlapping sections from adjacent spheroids. 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


args = split(getArgument()," "); 	// read in the array of argurments that are separated by a space
dirsa = args[0];	// The first argument is the directory
nchannels=args[1];	// The second argument is the number of channels
ch2_name=args[2];	// The third argument is the name of the protein in channel 2
ch3_name=args[3];	// The fourth argument is the name of the protein in channel 3
ch4_name=args[4];	// The fifth argument is the name of the protein in channel 4
actin_ch=args[5];

setBatchMode(true);	// run in batch mode to save time and memory

///creating new subdirectories for storing the processed data//
dirb= dirsa + "3d obj sphereoid"+ File.separator;	// define the path to a new folder that will contain the image afer 3d segmentation to indentify spheroids
dir= dirsa + "sphereoid_DNA"+ File.separator;		// define the path to a new folder that will contain the first channel of the raw image cropped to represent each spheroid
dirc2= dirsa + "sphereoid_"+ch2_name+ File.separator; 	// define the path to a new folder that will contain the second channel of the raw image cropped to represent each spheroid
dirc3= dirsa + "sphereoid_"+ch3_name+ File.separator; 	// define the path to a new folder that will contain the third channel of the raw image cropped to represent each spheroid
dirc4= dirsa + "sphereoid_"+ch4_name+ File.separator; 	// define the path to a new folder that will contain the fourth channel of the raw image cropped to represent each spheroid
// create the paths defined above 
File.makeDirectory(dirb);
File.makeDirectory(dir); 
if(nchannels>1){ 
	File.makeDirectory(dirc2); 
}
if(nchannels>2){ 
	File.makeDirectory(dirc3); 
}
if(nchannels>3){ 
	File.makeDirectory(dirc4); 
}

dirraw= dirsa + "rawimages"+ File.separator;	// define the path to the folder containing the unprocessed or raw images
filenamesa=getFileList(dirraw);		// get the list of names of the images in the rawimages folder


//set the options for 3D objects counter and 3D ROI manager 
run("3D OC Options", "volume surface nb_of_obj._voxels nb_of_surf._voxels integrated_density mean_gray_value std_dev_gray_value median_gray_value minimum_gray_value maximum_gray_value centroid mean_distance_to_surface std_dev_distance_to_surface median_distance_to_surface centre_of_mass bounding_box dots_size=5 font_size=20 redirect_to=none");
run("3D Manager Options", "volume surface compactness fit_ellipse 3d_moments integrated_density mean_grey_value std_dev_grey_value mode_grey_value feret minimum_grey_value maximum_grey_value centroid_(pix) centroid_(unit) distance_to_surface centre_of_mass_(pix) centre_of_mass_(unit) bounding_box radial_distance surface_contact closest use distance_between_centers=10 distance_max_contact=1.80");
run("3D Manager");	// Open the 3D ROI manager

for(f=0;f<filenamesa.length;f++){			//loop through all images

	// open the actin channel
	path=dirraw+filenamesa[f]; 	//define the path
	run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT c_begin=actin_ch c_end=actin_ch c_step=1");	// open the image

	// extract the the name of the file or remove file extension
	imgName=getTitle(); 	//get the title and assign it. it will be a character string
	baseNameEnd=indexOf(imgName, ".nd2"); 	//find the index of the string at with ".nd2" first appears NOTE: this only works for .nd2 files, change this to correct file extension if working with other files 
	baseName=substring(imgName, 0, baseNameEnd); // get the substring such that the file extention gets removed from the character string 

	//Smoothen, threshold, binarize and fill holes to obtain the binary image of a spheroid
	run("Gaussian Blur 3D...", "x=10 y=10 z=2"); // perform a 3D gaussian blur to smoothen over any inter-nuclear spaces 
	setAutoThreshold("Huang dark stack"); // Threshold the nucleus
	run("Make Binary", "method=Default background=Default");	// Binarize the nucleus
	run("Invert LUT");//invert LUTs
	run("Fill Holes", "stack"); //Fill any holes
	rename("actin");
	
	// open the nuclear channel (assummed to be channe1)  
	run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT c_begin=1 c_end=1 c_step=1");	// open the image
	//Smoothen, threshold, binarize and fill holes to obtain the binary image of a spheroid
	run("Gaussian Blur 3D...", "x=10 y=10 z=2"); // perform a 3D gaussian blur to smoothen over any inter-nuclear spaces 
	setAutoThreshold("Otsu dark stack"); // Threshold the nucleus
	run("Make Binary", "method=Default background=Default");	// Binarize the nucleus
	run("Invert LUT");//invert LUTs
	run("Fill Holes", "stack"); //Fill any holes
	rename("nuc");

	imageCalculator("OR create stack", "actin","nuc");
	run("Fill Holes", "stack");
	// Identify 3D objects of size ranging from 5,000-1,000,000 cu.microns
	getVoxelSize(width, height, depth, unit); 	// get the voxel dimentions
	a=1000000000/(width*height*depth);// obtain the number of pixels in a  volume of 1,000,000 cu.microns or  maximum value
	a_1=3000/(width*height*depth);// obtain the number of pixels in a  volume of 5,000 cu.microns or  minimum value
	run("3D Objects Counter", "threshold=128 min.=a_1 max.=a objects"); // Identify 3D objects using the 3D objects counter with a Size filter, display objects
	
	// Check if there are any 3d objects in the image 
	run("Z Project...", "projection=[Max Intensity]"); 	// project the image containg the objects
	run("Measure"); // measure
	close(); // close the projected image
	if(getResult("Mean",0)>0){    // if there is an object, then the mean intensity will be greater than 0
		
		saveAs("Tiff", dirb+baseName+".tiff"); // save the image containing the identified the objects
		obj_image=getTitle(); // get the image title
		
		Ext.Manager3D_AddImage();	// add the segmented image
		Ext.Manager3D_Count(nb_obj); // get the number of segmented rois
		run("Close All");	//close all images


		for(p=0;p<nb_obj;p++){ 		// loop through all objects 
			
			t=p+1; // temporary variable whos value is one more than the index

			//get the bounding box 
			Ext.Manager3D_Bounding3D(p,x0,x1,y0,y1,z0,z1);	// for the pth object, get the cordinates of the bounding rectangle
			zlows=z0+1; // the bottom z slice
			zhighs=z1+1; // the top z slice
			wids=x1-x0; // width along x axis
			heigs=y1-y0; //height along y axis

			// crop the image along the bounding box
			open(dirb+obj_image); // open the 3d objects image
			makeRectangle(x0, y0, wids, heigs); // make the 2D rectangle 
			run("Duplicate...", "duplicate range=zlows-zhighs"); // crop in z between the top and bottom slice of the spheroid

			// binarize to identify only the pth object: binary mask of the object
			setThreshold(t, t); // set threshold to identfy p-th object with gray level t
			run("Make Binary", "method=Default background=Default"); // binarize the object
			img1=getTitle(); // get the image title
			Ext.Manager3D_CloseResult("M"); // close the results window
		
			//For channel 1, crop the raw image along the bounding box of the pth object and retain only sections of the image that are also present in the binary mask of the pth object: to remove any overlapping sections
			run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT z_begin=zlows z_end=zhighs c_begin=1 c_end=1 c_step=1"); // open the first channel of the raw image between the top and bottom slice of the spheroid
			makeRectangle(x0, y0, wids, heigs);// make the 2D rectangle
			run("Duplicate...", "duplicate"); // crop
			img2=getTitle();// get the image title
			imageCalculator("AND create stack", img1,img2); // select region in the croped raw image that are also in the binary mask of the pth object
			run("Invert LUT"); // Invert LUT
			names=baseName+"_spheroid_"+x0+"_"+y0+"_"+z0+"_"+z1; // set the name of the image and spheroid identifier
			saveAs("tiff",dir+names);	// save the image
			Ext.Manager3D_CloseResult("M"); 	// close the results window 
		
			selectImage(img1); // select the binary mask image
			close("\\Others"); // close all others

			if(nchannels>1){ // if channel 2 is needed
				//For channel 2, crop the raw image along the bounding box of the pth object and retain only sections of the image that are also present in the binary mask of the pth object: to remove any overlapping sections
				run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT z_begin=zlows z_end=zhighs c_begin=2 c_end=2 c_step=1");  // open the second channel of the raw image between the top and bottom slice of the spheroid
				makeRectangle(x0, y0, wids, heigs);// make the 2D rectangle
				run("Duplicate...", "duplicate"); // crop 
				img2=getTitle();// get the image title
				imageCalculator("AND create stack", img1,img2);// select region in the croped raw image that are also in the binary mask of the pth object
				run("Invert LUT");// Invert LUT
				names=baseName+"_spheroid_"+x0+"_"+y0+"_"+z0+"_"+z1;// set the name of the image and spheroid identifier
				saveAs("tiff",dirc2+names);// save the image
				Ext.Manager3D_CloseResult("M");// close the results window 
		
				selectImage(img1); // select the binary mask image
				close("\\Others"); // close all others
			}
			if(nchannels>2){ // if channel 3 is needed
				//For channel 3, crop the raw image along the bounding box of the pth object and retain only sections of the image that are also present in the binary mask of the pth object: to remove any overlapping sections
				run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT z_begin=zlows z_end=zhighs c_begin=3 c_end=3 c_step=1"); 	// open the third channel of the raw image between the top and bottom slice of the spheroid
				makeRectangle(x0, y0, wids, heigs);// make the 2D rectangle
				run("Duplicate...", "duplicate");// crop
				img2=getTitle();// get the image title
				imageCalculator("AND create stack", img1,img2); // select region in the croped raw image that are also in the binary mask of the pth object
				run("Invert LUT");// Invert LUT
				names=baseName+"_spheroid_"+x0+"_"+y0+"_"+z0+"_"+z1; // set the name of the image and spheroid identifier
				saveAs("tiff",dirc3+names); // save the image
				Ext.Manager3D_CloseResult("M"); // close the results window 
		
				selectImage(img1); // select the binary mask image
				close("\\Others"); // close all others
			}
			if(nchannels>3){ // if channel 4 is needed
				//For channel 4, crop the raw image along the bounding box of the pth object and retain only sections of the image that are also present in the binary mask of the pth object: to remove any overlapping sections
				run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT z_begin=zlows z_end=zhighs c_begin=4 c_end=4 c_step=1");	// open the fourth channel of the raw image between the top and bottom slice of the spheroid
				makeRectangle(x0, y0, wids, heigs);// make the 2D rectangle
				run("Duplicate...", "duplicate");// crop
				img2=getTitle();// get the image title
				imageCalculator("AND create stack", img1,img2); // select region in the croped raw image that are also in the binary mask of the pth object
				run("Invert LUT");// Invert LUT
				names=baseName+"_spheroid_"+x0+"_"+y0+"_"+z0+"_"+z1; // set the name of the image and spheroid identifier
				saveAs("tiff",dirc4+names); // save the image
				Ext.Manager3D_CloseResult("M"); // close the results window 
		
		
				selectImage(img1); // select the binary mask image
				close("\\Others"); // close all others
			}
			run("Select None"); // get rid of all sections
     		run("Close All"); // close all open windows
		}
		run("Clear Results"); 	// clear the results table
	}
	run("Close All"); // close all open windows
  //  Ext.Manager3D_Close(); // close the 3D roi manager
  Ext.Manager3D_Reset();
   	run("Clear Results"); // clear the results table
}


 Ext.Manager3D_Close(); // close the 3D roi manager