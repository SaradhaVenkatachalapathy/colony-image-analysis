/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////  SCRIPT TO 3D PROJECT ALL THE IMAGES 
///////  WRITTEN BY: SARADHA VENKATACHALAPATHY                                                                                                                   
///////  ASSUMPTIONS: The input image is a confocal zstack 
///////  DESCRIPTION: The script accepts the directory to a folder containing rawimages folder.
///////				  It creates 1 subfolder for storing the 3D-projected images.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

args = split(getArgument()," "); 	// read in the array of argurments that are separated by a space
dirsa = args[0];	// The first argument is the directory
nchannels=args[1];	// The second argument is the number of channels


// Make the new directory for projected image
newDir1 = dirsa + "project3d" + File.separator; 
File.makeDirectory(newDir1); 

dir= dirsa + "rawimages"+ File.separator; //set the path to the raw images
filename = getFileList(dir); // get the list of files

setBatchMode(true);	// run in batch mode to save time and memory

for (i=0; i<filename.length; i++) { 
	path=dir+filename[i]; // set path 
	run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT c_begin=1 c_end=nchannels c_step=1"); // open image

    // extract the the name of the file or remove file extension
	imgName=getTitle(); 	//get the title and assign it. it will be a character string
	baseNameEnd=indexOf(imgName, ".nd2"); 	//find the index of the string at with ".nd2" first appears NOTE: this only works for .nd2 files, change this to correct file extension if working with other files 
	baseName=substring(imgName, 0, baseNameEnd); // get the substring such that the file extention gets removed from the character string 

	// 3D project the image
	run("3D Project...", "projection=[Brightest Point] axis=X-Axis slice=1 initial=0 total=360 rotation=10 lower=1 upper=255 opacity=0 surface=100 interior=50 interpolate");

	saveAs("tiff", newDir1 + baseName); // save the image
    run("Close All"); //close the images
} 





	  