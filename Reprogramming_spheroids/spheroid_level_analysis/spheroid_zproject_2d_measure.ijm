/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////  SCRIPT TO QUANTIFY THE PROJECTED GEOMETERICAL AND INTENSITY FEATURES
///////  WRITTEN BY: SARADHA VENKATACHALAPATHY                                                                                                                   
///////  ASSUMPTIONS: The input image is a confocal zstack and the channel 1 contains the nucleus. 
///////  DESCRIPTION: The script accepts the directory to a folder containing rawimages folder and the number of channels in the raw image that need to be analysed and the taget of the stains. 
///////				  It creates 1 subfolders: one for storing the features. It opens each image, does a sum projection and then measures geometrical and intensity features.  
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

args = split(getArgument()," "); 	// read in the array of argurments that are separated by a space
dirsa = args[0];	// The first argument is the directory
nchannels=args[1];	// The second argument is the number of channels
ch2_name=args[2];	// The third argument is the name of the protein in channel 2
ch3_name=args[3];	// The fourth argument is the name of the protein in channel 3
ch4_name=args[4];	// The fifth argument is the name of the protein in channel 4

// make a subdirectory to store the data
dirres= dirsa + "2D_measures_spheroid"+ File.separator;
File.makeDirectory(dirres); 

setBatchMode(true);	// run in batch mode to save time and memory

// set the measurements to be made//
run("Set Measurements...", "area mean standard modal min centroid center perimeter bounding fit shape feret's integrated median skewness kurtosis limit display redirect=None decimal=4");

dir1= dirsa + "sphereoid_DNA"+ File.separator;		// define the path to a new folder that will contain the first channel of the raw image cropped to represent each spheroid
list1 = getFileList(dir1); // get the list of the files
run("Clear Results");
for (i=0; i<list1.length; i++) {
	path = dir1+list1[i]; // set the path 
	open(path); //open image
	
	run("Z Project...", " projection=[Sum Slices]"); //Z project
	run("8-bit"); // convert to 8-bit image
	setThreshold(1, 255); 	// Set the threshold 
	
	run("Measure"); //Run the analysis
	run("Close All"); //close the images
}
run("Close All");
saveAs("Results",  dirres + "2D_spheroid_DNA.csv"); 
run("Clear Results");

if(nchannels>1){ 
	dir1= dirsa + "sphereoid_"+ch2_name+ File.separator; 	// define the path to a new folder that will contain the second channel of the raw image cropped to represent each spheroid
	list1 = getFileList(dir1); // get the list of the files

	for (i=0; i<list1.length; i++) {
		path = dir1+list1[i]; // set the path 
		open(path); //open image
	
		run("Z Project...", " projection=[Sum Slices]"); //Z project
		run("8-bit"); // convert to 8-bit image
		setThreshold(1, 255); 	// Set the threshold 
	
		run("Measure"); //Run the analysis
		run("Close All"); //close the images
	}
	run("Close All");
	saveAs("Results",  dirres + "2D_spheroid_"+ch2_name+".csv"); 
	run("Clear Results");
}

if(nchannels>2){ 
	dir1= dirsa + "sphereoid_"+ch3_name+ File.separator; 	// define the path to a new folder that will contain the third channel of the raw image cropped to represent each spheroid
	list1 = getFileList(dir1); // get the list of the files

	for (i=0; i<list1.length; i++) {
		path = dir1+list1[i]; // set the path 
		open(path); //open image
	
		run("Z Project...", " projection=[Sum Slices]"); //Z project
		run("8-bit"); // convert to 8-bit image
		setThreshold(1, 255); 	// Set the threshold 
	
		run("Measure"); //Run the analysis
		run("Close All"); //close the images
	}
	run("Close All");
	saveAs("Results",  dirres + "2D_spheroid_"+ch3_name+".csv"); 
	run("Clear Results");
}

if(nchannels>3){ 
	dir1= dirsa + "sphereoid_"+ch4_name+ File.separator; 	// define the path to a new folder that will contain the fourth channel of the raw image cropped to represent each spheroid
	list1 = getFileList(dir1); // get the list of the files

	for (i=0; i<list1.length; i++) {
		path = dir1+list1[i]; // set the path 
		open(path); //open image
	
		run("Z Project...", " projection=[Sum Slices]"); //Z project
		run("8-bit"); // convert to 8-bit image
		setThreshold(1, 255); 	// Set the threshold 
	
		run("Measure"); //Run the analysis
		run("Close All"); //close the images
	}
	run("Close All");
	saveAs("Results",  dirres + "2D_spheroid_"+ch4_name+".csv"); 
	run("Clear Results");
}