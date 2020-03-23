/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////  SCRIPT COMPUTE THE ELLIPTIC FOURIER DESCRIPTORS OF AN IMAGE
///////  WRITTEN BY: SARADHA VENKATACHALAPATHY                                                                                                                   
///////  ASSUMPTIONS: The input image is a confocal zstack and the channel 1 contains the nucleus. 
///////  DESCRIPTION: The script opens images of individial nuclei, zprojects the confocal stack and computes the elliptic fourier descriptors of the nuclear boundary
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

dirsa = getArgument();	// The first argument is the directory

// make a subdirectory to store the data
dirres= dirsa + "2D_measures_nuclei"+ File.separator;
File.makeDirectory(dirres); 
setBatchMode(true);	// run in batch mode to save time and memory

dir1= dirsa + "indivisual_nuclei_DNA"+ File.separator;		// define the path to a new folder that will contain the first channel of the raw image cropped to represent each spheroid
list1 = getFileList(dir1); // get the list of the files

print("\\Clear"); //empties the Log
print("Label,EFD1,EFD2,EFD3,EFD4,EFD5,EFD6,EFD7,EFD8,EFD9,EFD10");

for (i=0; i<list1.length; i++) {
	path = dir1+list1[i]; // set the path 
	open(path); //open image

	//get the base name
	obj_image=getTitle();
	baseNameEnd=indexOf(obj_image, ".tif"); 
	baseName2=substring(obj_image, 0, baseNameEnd);
	
	//project the image and threshold
	if(nSlices>1){
	run("Z Project...", "projection=[Max Intensity]");
	setThreshold(1, 255);
	obj_image=getTitle();
	baseNameEnd=indexOf(obj_image, ".tif"); 
	baseName=substring(obj_image, 0, baseNameEnd);
	baseNameEnd=indexOf(baseName, " "); 
	if(baseNameEnd>0){
		baseName=substring(obj_image, 0, baseNameEnd);
	}

	//analyze particles to get roi of image
	run("Analyze Particles...", "add");
	roiManager("Select", 0);
	//compute elliptic fourer discriptors and save results
	run("Elliptic Fourier D.", "number=12 results");
	name="Results-EFD-"+baseName;
	IJ.renameResults(name,"Results");
	print(baseName2+","+getResult("efd", 2)+","+getResult("efd", 3)+","+getResult("efd", 4)+","+getResult("efd", 5)+","+getResult("efd", 6)+
	","+getResult("efd", 7)+","+getResult("efd", 8)+","+getResult("efd", 9)+","+getResult("efd", 10)+","+getResult("efd", 11));
	roiManager("reset");
	run("Clear Results");
	}
}

selectWindow("Log");  //select Log-window
saveAs("Text", dirres + "edf.csv");
print("\\Clear"); //empties the Log