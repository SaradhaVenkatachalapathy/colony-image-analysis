/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////  SCRIPT TO COMPUTE INTENSITY OF THE SPHEROIDS ALONG THE Z AXIS
///////  WRITTEN BY: SARADHA VENKATACHALAPATHY                                                                                                                   
///////  ASSUMPTIONS: The input image is a confocal zstack and the channel 1 contains the nucleus.  Tested on ImageJ v1.52p.
///////  DESCRIPTION: The script accepts the directory to a folder containing segmented spheroids folder and the number of channels in the raw image that need to be analysed and the target of the stains. 
///////				  It creates 1 subfolders: one for storing the features. It opens each image, and measures geometrical and intensity features for each slice  
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


args = split(getArgument()," "); 	// read in the array of argurments that are separated by a space
dirsa = args[0];	// The first argument is the directory
nchannels=args[1];	// The second argument is the number of channels
ch2_name=args[2];	// The third argument is the name of the protein in channel 2
ch3_name=args[3];	// The fourth argument is the name of the protein in channel 3
ch4_name=args[4];	// The fifth argument is the name of the protein in channel 4

// make a subdirectory to store the data
dirres= dirsa + "stepwise_measures_spheroid"+ File.separator;
Dir_res_ch1_int=dirres + "DNA" + File.separator;  // define the path to a new folder that will contain the 3d intensity data for the spheroids in channel 1/ DNA channel
Dir_res_ch2_int=dirres + ch2_name + File.separator;  // define the path to a new folder that will contain the 3d intensity data for the spheroids in channel 2
Dir_res_ch3_int=dirres + ch3_name + File.separator;  // define the path to a new folder that will contain the 3d intensity data for the spheroids in channel 3
Dir_res_ch4_int=dirres + ch4_name + File.separator;  // define the path to a new folder that will contain the 3d intensity data for the spheroids in channel 4
// create the paths defined above 
File.makeDirectory(dirres); 
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
setBatchMode(true);	// run in batch mode to save time and memory
dir1= dirsa + "sphereoid_DNA"+ File.separator;		// define the path to a new folder that will contain the first channel of the raw image cropped to represent each spheroid
list1 = getFileList(dir1); // get the list of the files

for (i=0; i<list1.length; i++) {
	path = dir1+list1[i]; // set the path 
	open(path); //open image
	run("8-bit"); // convert to 8-bit image
	
	// extract the the name of the file or remove file extension
	imgName=getTitle(); 	//get the title and assign it. it will be a character string
	baseNameEnd=indexOf(imgName, ".tif"); 	//find the index of the string at with ".tif" first appears NOTE: this only works for .tif files, change this to correct file extension if working with other files 
	baseName=substring(imgName, 0, baseNameEnd); // get the substring such that the file extention gets removed from the character string 

	run("Clear Results");	//Clear Results
	setThreshold(1, 255); 	// Set the threshold 

	for (j=0; j<nSlices; j++){ //open all slices
		setSlice(j+1);	// set Slice
		run("Measure");		// Measure
	}
	run("Close All"); //close the images
	saveAs("Results",  Dir_res_ch1_int + baseName +".csv"); 
}

dir1= dirsa + "sphereoid_DNA"+ File.separator;		// define the path to a new folder that will contain the first channel of the raw image cropped to represent each spheroid
list1 = getFileList(dir1); // get the list of the files

for (i=0; i<list1.length; i++) {
	path = dir1+list1[i]; // set the path 
	open(path); //open image
	run("8-bit"); // convert to 8-bit image
	
	// extract the the name of the file or remove file extension
	imgName=getTitle(); 	//get the title and assign it. it will be a character string
	baseNameEnd=indexOf(imgName, ".tif"); 	//find the index of the string at with ".tif" first appears NOTE: this only works for .tif files, change this to correct file extension if working with other files 
	baseName=substring(imgName, 0, baseNameEnd); // get the substring such that the file extention gets removed from the character string 

	run("Clear Results");	//Clear Results
	setThreshold(1, 255); 	// Set the threshold 

	for (j=0; j<nSlices; j++){ //open all slices
		setSlice(j+1);	// set Slice
		run("Measure");		// Measure
	}
	run("Close All"); //close the images
	saveAs("Results",  Dir_res_ch1_int + baseName +".csv"); 
}

if(nchannels>1){ 
	dir1= dirsa + "sphereoid_"+ch2_name+ File.separator; 	// define the path to a new folder that will contain the second channel of the raw image cropped to represent each spheroid
	list1 = getFileList(dir1); // get the list of the files
	for (i=0; i<list1.length; i++) {
		path = dir1+list1[i]; // set the path 
		open(path); //open image
		run("8-bit"); // convert to 8-bit image
	
		// extract the the name of the file or remove file extension
		imgName=getTitle(); 	//get the title and assign it. it will be a character string
		baseNameEnd=indexOf(imgName, ".tif"); 	//find the index of the string at with ".tif" first appears NOTE: this only works for .tif files, change this to correct file extension if working with other files 
		baseName=substring(imgName, 0, baseNameEnd); // get the substring such that the file extention gets removed from the character string 

		run("Clear Results");	//Clear Results
		setThreshold(1, 255); 	// Set the threshold 

		for (j=0; j<nSlices; j++){ //open all slices
			setSlice(j+1);	// set Slice
			run("Measure");		// Measure
		}
		run("Close All"); //close the images
		saveAs("Results",  Dir_res_ch2_int + baseName +".csv"); 
	}
}

if(nchannels>2){ 
	dir1= dirsa + "sphereoid_"+ch3_name+ File.separator; 	// define the path to a new folder that will contain the second channel of the raw image cropped to represent each spheroid
	list1 = getFileList(dir1); // get the list of the files
	for (i=0; i<list1.length; i++) {
		path = dir1+list1[i]; // set the path 
		open(path); //open image
		run("8-bit"); // convert to 8-bit image
	
		// extract the the name of the file or remove file extension
		imgName=getTitle(); 	//get the title and assign it. it will be a character string
		baseNameEnd=indexOf(imgName, ".tif"); 	//find the index of the string at with ".nd2" first appears NOTE: this only works for .tif files, change this to correct file extension if working with other files 
		baseName=substring(imgName, 0, baseNameEnd); // get the substring such that the file extention gets removed from the character string 

		run("Clear Results");	//Clear Results
		setThreshold(1, 255); 	// Set the threshold 

		for (j=0; j<nSlices; j++){ //open all slices
			setSlice(j+1);	// set Slice
			run("Measure");		// Measure
		}
		run("Close All"); //close the images
		saveAs("Results",  Dir_res_ch3_int + baseName +".csv"); 
	}
}

if(nchannels>3){ 
	dir1= dirsa + "sphereoid_"+ch3_name+ File.separator; 	// define the path to a new folder that will contain the second channel of the raw image cropped to represent each spheroid
	list1 = getFileList(dir1); // get the list of the files
	for (i=0; i<list1.length; i++) {
		path = dir1+list1[i]; // set the path 
		open(path); //open image
		run("8-bit"); // convert to 8-bit image
	
		// extract the the name of the file or remove file extension
		imgName=getTitle(); 	//get the title and assign it. it will be a character string
		baseNameEnd=indexOf(imgName, ".tif"); 	//find the index of the string at with ".tif" first appears NOTE: this only works for .tif files, change this to correct file extension if working with other files 
		baseName=substring(imgName, 0, baseNameEnd); // get the substring such that the file extention gets removed from the character string 

		run("Clear Results");	//Clear Results
		setThreshold(1, 255); 	// Set the threshold 

		for (j=0; j<nSlices; j++){ //open all slices
			setSlice(j+1);	// set Slice
			run("Measure");		// Measure
		}
		run("Close All"); //close the images
		saveAs("Results",  Dir_res_ch4_int + baseName +".csv"); 
	}
}
