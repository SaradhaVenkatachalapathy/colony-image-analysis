/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////  SCRIPT TO PROJECT ALL THE IMAGES AND MAKE A MONTANGE FOR EACH IMAGE
///////  WRITTEN BY: SARADHA VENKATACHALAPATHY                                                                                                                   
///////  ASSUMPTIONS: The input image is a confocal zstack and the channel 1 contains the nucleus. 
///////  DESCRIPTION: The script accepts the directory to a folder containing rawimages folder and the number of channels in the raw image that need to be analysed and the taget of the stains. 
///////				  It creates 2 subfolders: one for storing the z-projected images. It opens each image, does a max intensity projection, stores the image. Then open all projected images 
///////				  and create a montage with a scale bar of 50 microns
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

args = split(getArgument()," "); 	// read in the array of argurments that are separated by a space
dirsa = args[0];	// The first argument is the directory
nchannels=args[1];	// The second argument is the number of channels
ch2_name=args[2];	// The third argument is the name of the protein in channel 2
ch3_name=args[3];	// The fourth argument is the name of the protein in channel 3
ch4_name=args[4];	// The fifth argument is the name of the protein in channel 4

dir= dirsa + "rawimages"+ File.separator; //set the path to the raw images
filename = getFileList(dir); // get the list of files

setBatchMode(true);	// run in batch mode to save time and memory

// Make the new directories for projected images and the montage
newDir1 = dirsa + "zproject" + File.separator; 
newDirm = dirsa + "Montage" + File.separator; 
File.makeDirectory(newDir1); 
File.makeDirectory(newDirm); 
 
/////// CHANNEL 1
// make sub directory for the first channel
newDir = newDir1 + "DNA" + File.separator; 
File.makeDirectory(newDir); 

for (i=0; i<filename.length; i++) { 
      path=dir+filename[i]; // set path 
	  run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT c_begin=1 c_end=1 c_step=1"); // open image

      // extract the the name of the file or remove file extension
	  imgName=getTitle(); 	//get the title and assign it. it will be a character string
	  baseNameEnd=indexOf(imgName, ".nd2"); 	//find the index of the string at with ".nd2" first appears NOTE: this only works for .nd2 files, change this to correct file extension if working with other files 
	  baseName=substring(imgName, 0, baseNameEnd); // get the substring such that the file extention gets removed from the character string 

      run("Z Project...", "projection=[Max Intensity]"); // project images
      saveAs("tiff", newDir + baseName); // save the image
      run("Close All"); //close the images
} 

run("Image Sequence...", "open=newDir sort"); // open the projected images as a stack
run("Scale Bar...", "width=50 font=14 color=White background=None location=[Lower Right] bold hide label"); // add scale bar of 50 microns
run("Make Montage...", "scale=1 border=3"); // Make a montage
saveAs("Jpeg", newDirm + "Montage_raw_DNA"); // save the image of the montage

/////// CHANNEL 2
if(nchannels>1){ 
	// make sub directory for the first channel
	newDir = newDir1 + ch2_name + File.separator; 
	File.makeDirectory(newDir); 

	for (i=0; i<filename.length; i++) { 
     	 path=dir+filename[i]; // set path 
	  	run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT c_begin=2 c_end=2 c_step=1"); // open image

      	// extract the the name of the file or remove file extension
	  	imgName=getTitle(); 	//get the title and assign it. it will be a character string
	  	baseNameEnd=indexOf(imgName, ".nd2"); 	//find the index of the string at with ".nd2" first appears NOTE: this only works for .nd2 files, change this to correct file extension if working with other files 
	  	baseName=substring(imgName, 0, baseNameEnd); // get the substring such that the file extention gets removed from the character string 

      	run("Z Project...", "projection=[Max Intensity]"); // project images
      	saveAs("tiff", newDir + baseName); // save the image
      	run("Close All"); //close the images
	} 

	run("Image Sequence...", "open=newDir sort"); // open the projected images as a stack
	run("Scale Bar...", "width=50 font=14 color=White background=None location=[Lower Right] bold hide label"); // add scale bar of 50 microns
	run("Make Montage...", "scale=1 border=3"); // Make a montage
	saveAs("Jpeg", newDirm + "Montage_raw_"+ch2_name); // save the image of the montage
}

/////// CHANNEL 3
if(nchannels>2){ 
	// make sub directory for the first channel
	newDir = newDir1 + ch3_name + File.separator; 
	File.makeDirectory(newDir); 

	for (i=0; i<filename.length; i++) { 
     	 path=dir+filename[i]; // set path 
	  	run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT c_begin=3 c_end=3 c_step=1"); // open image

      	// extract the the name of the file or remove file extension
	  	imgName=getTitle(); 	//get the title and assign it. it will be a character string
	  	baseNameEnd=indexOf(imgName, ".nd2"); 	//find the index of the string at with ".nd2" first appears NOTE: this only works for .nd2 files, change this to correct file extension if working with other files 
	  	baseName=substring(imgName, 0, baseNameEnd); // get the substring such that the file extention gets removed from the character string 

      	run("Z Project...", "projection=[Max Intensity]"); // project images
      	saveAs("tiff", newDir + baseName); // save the image
      	run("Close All"); //close the images
	} 

	run("Image Sequence...", "open=newDir sort"); // open the projected images as a stack
	run("Scale Bar...", "width=50 font=14 color=White background=None location=[Lower Right] bold hide label"); // add scale bar of 50 microns
	run("Make Montage...", "scale=1 border=3"); // Make a montage
	saveAs("Jpeg", newDirm + "Montage_raw_"+ch3_name); // save the image of the montage

}

/////// CHANNEL 4
if(nchannels>3){ 
	// make sub directory for the first channel
	newDir = newDir1 + ch4_name + File.separator; 
	File.makeDirectory(newDir); 

	for (i=0; i<filename.length; i++) { 
     	 path=dir+filename[i]; // set path 
	  	run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT c_begin=4 c_end=4 c_step=1"); // open image

      	// extract the the name of the file or remove file extension
	  	imgName=getTitle(); 	//get the title and assign it. it will be a character string
	  	baseNameEnd=indexOf(imgName, ".nd2"); 	//find the index of the string at with ".nd2" first appears NOTE: this only works for .nd2 files, change this to correct file extension if working with other files 
	  	baseName=substring(imgName, 0, baseNameEnd); // get the substring such that the file extention gets removed from the character string 

      	run("Z Project...", "projection=[Max Intensity]"); // project images
      	saveAs("tiff", newDir + baseName); // save the image
      	run("Close All"); //close the images
	} 

	run("Image Sequence...", "open=newDir sort"); // open the projected images as a stack
	run("Scale Bar...", "width=50 font=14 color=White background=None location=[Lower Right] bold hide label"); // add scale bar of 50 microns
	run("Make Montage...", "scale=1 border=3"); // Make a montage
	saveAs("Jpeg", newDirm + "Montage_raw_"+ch4_name); // save the image of the montage

}