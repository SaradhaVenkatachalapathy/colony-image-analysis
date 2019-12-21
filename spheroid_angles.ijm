///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////  SCRIPT TO OOMPUTE ANGLULAR DISTRIBUTION FROM XZ PROJECTION OF A CONFOCAL STACK
///////  WRITTEN BY: SARADHA VENKATACHALAPATHY                                                                                                                   
///////  ASSUMPTIONS: The input image is a confocal zstack and the channel 1 contains the nucleus.  Tested on ImageJ v1.52p.
///////  DESCRIPTION: The script accepts the directory to a folder containing segmented spheroids folder and the number and names of channels in the raw image that need to be analysed and the target of the stains. 
///////				  It creates 1 subfolders: one for storing the features for each channel. It opens each image, obtains zprotection of the XY projection and computes the FFT of the image. Then it computes the 
///////				  frequency distribution of the angles in the image. 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


args = split(getArgument()," "); 	// read in the array of argurments that are separated by a space
dirsa = args[0];	// The first argument is the directory
nchannels=args[1];	// The second argument is the number of channels
ch2_name=args[2];	// The third argument is the name of the protein in channel 2
ch3_name=args[3];	// The fourth argument is the name of the protein in channel 3
ch4_name=args[4];	// The fifth argument is the name of the protein in channel 4

// make a subdirectory to store the data
dirres= dirsa + "angles_spheroid"+ File.separator;
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
dir1=dirsa + "rawimages"+ File.separator;
list1 = getFileList(dir1); // get the list of the files

run("Set Measurements...", "mean modal min integrated median redirect=None decimal=4"); // Set Measurements

for (i=0; i<list1.length; i++) {
	
	path = dir1+list1[i]; // set the path 
	run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT c_begin=1 c_end=1 c_step=1");	// open the image
	
	// extract the the name of the file or remove file extension
	imgName=getTitle(); 	//get the title and assign it. it will be a character string
	baseNameEnd=indexOf(imgName, ".nd2"); 	//find the index of the string at with ".nd2" first appears NOTE: this only works for .nd2 files, change this to correct file extension if working with other files 
	baseName=substring(imgName, 0, baseNameEnd); // get the substring such that the file extention gets removed from the character string 

	run("Reslice [/]...", "output=1.000 start=Top flip");	// obtain the XZ projection
	run("Z Project...", "projection=[Max Intensity]");		// z projection
	run("FFT");												// run FFT
	run("Rotate 90 Degrees Right");							// rotate the image

	getDimensions(width, height, channels, slices, frames);	//obtain the image dimensions
	for(index=0; index<90; index++){						//iterate over the angles 0-90
		makeLine(0, height/2+tan((PI/180)*index*1)*width, width, height/2-tan((PI/180)*index*1)*width,1); //draw a line selection through the middle of the image
		run("Measure");							// measure
		setResult("angle", index, 1*index); 	//set angle
		roiManager("add");						//add to roi manager
	}
	makeLine(width/2, height, width/2, 0,1);	// draw vertical line
	run("Measure");								//measure
	setResult("angle", index, 1*index);			//set angle
	roiManager("add");							//add to roi manager
	
	saveAs("Results",  Dir_res_ch1_int + baseName +".csv"); 	//save results
	run("Clear Results");	 									//clear results
	run("Close All");											//Close all iamges
	roiManager("reset");										//reset roi manager
}

if(nchannels>1){ 
for (i=0; i<list1.length; i++) {
	
	path = dir1+list1[i]; // set the path 
	run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT c_begin=2 c_end=2 c_step=1");	// open the image
	
	// extract the the name of the file or remove file extension
	imgName=getTitle(); 	//get the title and assign it. it will be a character string
	baseNameEnd=indexOf(imgName, ".nd2"); 	//find the index of the string at with ".nd2" first appears NOTE: this only works for .nd2 files, change this to correct file extension if working with other files 
	baseName=substring(imgName, 0, baseNameEnd); // get the substring such that the file extention gets removed from the character string 

	run("Reslice [/]...", "output=1.000 start=Top flip");	// obtain the XZ projection
	run("Z Project...", "projection=[Max Intensity]");		// z projection
	run("FFT");												// run FFT
	run("Rotate 90 Degrees Right");							// rotate the image

	getDimensions(width, height, channels, slices, frames);	//obtain the image dimensions
	for(index=0; index<90; index++){						//iterate over the angles 0-90
		makeLine(0, height/2+tan((PI/180)*index*1)*width, width, height/2-tan((PI/180)*index*1)*width,1); //draw a line selection through the middle of the image
		run("Measure");							// measure
		setResult("angle", index, 1*index); 	//set angle
		roiManager("add");						//add to roi manager
	}
	makeLine(width/2, height, width/2, 0,1);	// draw vertical line
	run("Measure");								//measure
	setResult("angle", index, 1*index);			//set angle
	roiManager("add");							//add to roi manager
	
	saveAs("Results",  Dir_res_ch1_int + baseName +".csv"); 	//save results
	run("Clear Results");	 									//clear results
	run("Close All");											//Close all iamges
	roiManager("reset");										//reset roi manager
}
}

if(nchannels>2){ 
for (i=0; i<list1.length; i++) {
	
	path = dir1+list1[i]; // set the path 
	run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT c_begin=3 c_end=3 c_step=1");	// open the image
	
	// extract the the name of the file or remove file extension
	imgName=getTitle(); 	//get the title and assign it. it will be a character string
	baseNameEnd=indexOf(imgName, ".nd2"); 	//find the index of the string at with ".nd2" first appears NOTE: this only works for .nd2 files, change this to correct file extension if working with other files 
	baseName=substring(imgName, 0, baseNameEnd); // get the substring such that the file extention gets removed from the character string 

	run("Reslice [/]...", "output=1.000 start=Top flip");	// obtain the XZ projection
	run("Z Project...", "projection=[Max Intensity]");		// z projection
	run("FFT");												// run FFT
	run("Rotate 90 Degrees Right");							// rotate the image

	getDimensions(width, height, channels, slices, frames);	//obtain the image dimensions
	for(index=0; index<90; index++){						//iterate over the angles 0-90
		makeLine(0, height/2+tan((PI/180)*index*1)*width, width, height/2-tan((PI/180)*index*1)*width,1); //draw a line selection through the middle of the image
		run("Measure");							// measure
		setResult("angle", index, 1*index); 	//set angle
		roiManager("add");						//add to roi manager
	}
	makeLine(width/2, height, width/2, 0,1);	// draw vertical line
	run("Measure");								//measure
	setResult("angle", index, 1*index);			//set angle
	roiManager("add");							//add to roi manager
	
	saveAs("Results",  Dir_res_ch1_int + baseName +".csv"); 	//save results
	run("Clear Results");	 									//clear results
	run("Close All");											//Close all iamges
	roiManager("reset");										//reset roi manager
}
}

if(nchannels>3){ 
for (i=0; i<list1.length; i++) {
	
	path = dir1+list1[i]; // set the path 
	run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT c_begin=4 c_end=4 c_step=1");	// open the image
	
	// extract the the name of the file or remove file extension
	imgName=getTitle(); 	//get the title and assign it. it will be a character string
	baseNameEnd=indexOf(imgName, ".nd2"); 	//find the index of the string at with ".nd2" first appears NOTE: this only works for .nd2 files, change this to correct file extension if working with other files 
	baseName=substring(imgName, 0, baseNameEnd); // get the substring such that the file extention gets removed from the character string 

	run("Reslice [/]...", "output=1.000 start=Top flip");	// obtain the XZ projection
	run("Z Project...", "projection=[Max Intensity]");		// z projection
	run("FFT");												// run FFT
	run("Rotate 90 Degrees Right");							// rotate the image

	getDimensions(width, height, channels, slices, frames);	//obtain the image dimensions
	for(index=0; index<90; index++){						//iterate over the angles 0-90
		makeLine(0, height/2+tan((PI/180)*index*1)*width, width, height/2-tan((PI/180)*index*1)*width,1); //draw a line selection through the middle of the image
		run("Measure");							// measure
		setResult("angle", index, 1*index); 	//set angle
		roiManager("add");						//add to roi manager
	}
	makeLine(width/2, height, width/2, 0,1);	// draw vertical line
	run("Measure");								//measure
	setResult("angle", index, 1*index);			//set angle
	roiManager("add");							//add to roi manager
	
	saveAs("Results",  Dir_res_ch1_int + baseName +".csv"); 	//save results
	run("Clear Results");	 									//clear results
	run("Close All");											//Close all iamges
	roiManager("reset");										//reset roi manager
}
}

