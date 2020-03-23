/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////  SCRIPT TO COMPUTE THE GCLM FEATURES OF THE NUCLEUS
///////  WRITTEN BY: SARADHA VENKATACHALAPATHY                                                                                                                   
///////  ASSUMPTIONS: The input image is a confocal zstack and the channel 1 contains the nucleus. 
///////  DESCRIPTION: The script opens individual nuclei, zprojects, converts to an 8bit image and calculate the gclm texture features at 0, 90, 180 and 270 degrees with steps of 5,15,25,35 and 45 pixels.
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
print("Label,","Deg_0_Angular_Second_Moment,","Deg_0_Contrast,","Deg_0_Correlation,","Deg_0_Inverse_Difference_Moment,","Deg_0_Entropy,");

for (i=0; i<list1.length; i++) {
	path = dir1+list1[i]; // set the path 
	open(path); //open image
	obj_image=getTitle();
	baseNameEnd=indexOf(obj_image, ".tif"); 
	baseName=substring(obj_image, 0, baseNameEnd);
	if(nSlices>1){
		run("Z Project...", " projection=[Sum Slices]"); //Z project
		run("8-bit"); // convert to 8-bit image
		run("GLCM Texture", "enter="+5+ " select=[0 degrees] angular contrast correlation inverse entropy");
   		asm = getResult("Angular Second Moment",0); 
    	contrast = getResult("Contrast",0);
    	correlation = getResult("Correlation",0);
    	idm = getResult("Inverse Difference Moment   ",0); //Extra spaces needed due to source code error
    	entropy = getResult("Entropy",0);
    	print(baseName,",",asm,",",contrast,",",correlation,",",idm,",",entropy);
	}
	run("Close All"); //close the images
}
selectWindow("Log");  //select Log-window
saveAs("Text", dirres + "gclm_5_deg0.csv");

print("\\Clear"); //empties the Log
print("Label,","Deg_0_Angular_Second_Moment,","Deg_0_Contrast,","Deg_0_Correlation,","Deg_0_Inverse_Difference_Moment,","Deg_0_Entropy,");

for (i=0; i<list1.length; i++) {
	path = dir1+list1[i]; // set the path 
	open(path); //open image
	obj_image=getTitle();
	baseNameEnd=indexOf(obj_image, ".tif"); 
	baseName=substring(obj_image, 0, baseNameEnd);

	if(nSlices>1){
		run("Z Project...", " projection=[Sum Slices]"); //Z project
		run("8-bit"); // convert to 8-bit image
		run("GLCM Texture", "enter="+15+ " select=[0 degrees] angular contrast correlation inverse entropy");
    	asm = getResult("Angular Second Moment",0); 
    	contrast = getResult("Contrast",0);
    	correlation = getResult("Correlation",0);
    	idm = getResult("Inverse Difference Moment   ",0); //Extra spaces needed due to source code error
   		entropy = getResult("Entropy",0);
    	print(baseName,",",asm,",",contrast,",",correlation,",",idm,",",entropy);
	}
	run("Close All"); //close the images
}
selectWindow("Log");  //select Log-window
saveAs("Text", dirres + "gclm_15_deg0.csv")

print("\\Clear"); //empties the Log
print("Label,","Deg_0_Angular_Second_Moment,","Deg_0_Contrast,","Deg_0_Correlation,","Deg_0_Inverse_Difference_Moment,","Deg_0_Entropy,");

for (i=0; i<list1.length; i++) {
	path = dir1+list1[i]; // set the path 
	open(path); //open image
	obj_image=getTitle();
	baseNameEnd=indexOf(obj_image, ".tif"); 
	baseName=substring(obj_image, 0, baseNameEnd);
	if(nSlices>1){
	run("Z Project...", " projection=[Sum Slices]"); //Z project
	run("8-bit"); // convert to 8-bit image
	run("GLCM Texture", "enter="+25+ " select=[0 degrees] angular contrast correlation inverse entropy");
    asm = getResult("Angular Second Moment",0); 
    contrast = getResult("Contrast",0);
    correlation = getResult("Correlation",0);
    idm = getResult("Inverse Difference Moment   ",0); //Extra spaces needed due to source code error
    entropy = getResult("Entropy",0);
    print(baseName,",",asm,",",contrast,",",correlation,",",idm,",",entropy);
	}
	run("Close All"); //close the images
}
selectWindow("Log");  //select Log-window
saveAs("Text", dirres + "gclm_25_deg0.csv")

print("\\Clear"); //empties the Log
print("Label,","Deg_0_Angular_Second_Moment,","Deg_0_Contrast,","Deg_0_Correlation,","Deg_0_Inverse_Difference_Moment,","Deg_0_Entropy,");

for (i=0; i<list1.length; i++) {
	path = dir1+list1[i]; // set the path 
	open(path); //open image
	obj_image=getTitle();
	baseNameEnd=indexOf(obj_image, ".tif"); 
	baseName=substring(obj_image, 0, baseNameEnd);
	if(nSlices>1){
	run("Z Project...", " projection=[Sum Slices]"); //Z project
	run("8-bit"); // convert to 8-bit image
	run("GLCM Texture", "enter="+35+ " select=[0 degrees] angular contrast correlation inverse entropy");
    asm = getResult("Angular Second Moment",0); 
    contrast = getResult("Contrast",0);
    correlation = getResult("Correlation",0);
    idm = getResult("Inverse Difference Moment   ",0); //Extra spaces needed due to source code error
    entropy = getResult("Entropy",0);
    print(baseName,",",asm,",",contrast,",",correlation,",",idm,",",entropy);
	}
	run("Close All"); //close the images
}
selectWindow("Log");  //select Log-window
saveAs("Text", dirres + "gclm_35_deg0.csv")

print("\\Clear"); //empties the Log
print("Label,","Deg_0_Angular_Second_Moment,","Deg_0_Contrast,","Deg_0_Correlation,","Deg_0_Inverse_Difference_Moment,","Deg_0_Entropy,");

for (i=0; i<list1.length; i++) {
	path = dir1+list1[i]; // set the path 
	open(path); //open image
	obj_image=getTitle();
	baseNameEnd=indexOf(obj_image, ".tif"); 
	baseName=substring(obj_image, 0, baseNameEnd);
	if(nSlices>1){
	run("Z Project...", " projection=[Sum Slices]"); //Z project
	run("8-bit"); // convert to 8-bit image
	run("GLCM Texture", "enter="+45+ " select=[0 degrees] angular contrast correlation inverse entropy");
    asm = getResult("Angular Second Moment",0); 
    contrast = getResult("Contrast",0);
    correlation = getResult("Correlation",0);
    idm = getResult("Inverse Difference Moment   ",0); //Extra spaces needed due to source code error
    entropy = getResult("Entropy",0);
    print(baseName,",",asm,",",contrast,",",correlation,",",idm,",",entropy);
	}
	run("Close All"); //close the images
}
selectWindow("Log");  //select Log-window
saveAs("Text", dirres + "gclm_45_deg0.csv")

print("\\Clear"); //empties the Log
print("Label,","Deg_0_Angular_Second_Moment,","Deg_0_Contrast,","Deg_0_Correlation,","Deg_0_Inverse_Difference_Moment,","Deg_0_Entropy,");

for (i=0; i<list1.length; i++) {
	path = dir1+list1[i]; // set the path 
	open(path); //open image
	obj_image=getTitle();
	baseNameEnd=indexOf(obj_image, ".tif"); 
	baseName=substring(obj_image, 0, baseNameEnd);
	if(nSlices>1){
	run("Z Project...", " projection=[Sum Slices]"); //Z project
	run("8-bit"); // convert to 8-bit image
	run("GLCM Texture", "enter="+100+ " select=[0 degrees] angular contrast correlation inverse entropy");
    asm = getResult("Angular Second Moment",0); 
    contrast = getResult("Contrast",0);
    correlation = getResult("Correlation",0);
    idm = getResult("Inverse Difference Moment   ",0); //Extra spaces needed due to source code error
    entropy = getResult("Entropy",0);
    print(baseName,",",asm,",",contrast,",",correlation,",",idm,",",entropy);
	}
	run("Close All"); //close the images
}
selectWindow("Log");  //select Log-window
saveAs("Text", dirres + "gclm_100_deg0.csv")
////////////
print("\\Clear"); //empties the Log
print("Label,","Deg_90_Angular_Second_Moment,","Deg_90_Contrast,","Deg_90_Correlation,","Deg_90_Inverse_Difference_Moment,","Deg_90_Entropy,");

for (i=0; i<list1.length; i++) {
	path = dir1+list1[i]; // set the path 
	open(path); //open image
	obj_image=getTitle();
	baseNameEnd=indexOf(obj_image, ".tif"); 
	baseName=substring(obj_image, 0, baseNameEnd);
	if(nSlices>1){
	run("Z Project...", " projection=[Sum Slices]"); //Z project
	run("8-bit"); // convert to 8-bit image
	run("GLCM Texture", "enter="+5+ " select=[90 degrees] angular contrast correlation inverse entropy");
    asm = getResult("Angular Second Moment",0); 
    contrast = getResult("Contrast",0);
    correlation = getResult("Correlation",0);
    idm = getResult("Inverse Difference Moment   ",0); //Extra spaces needed due to source code error
    entropy = getResult("Entropy",0);
    print(baseName,",",asm,",",contrast,",",correlation,",",idm,",",entropy);
	}
	run("Close All"); //close the images
}
selectWindow("Log");  //select Log-window
saveAs("Text", dirres + "gclm_5_deg90.csv");

print("\\Clear"); //empties the Log
print("Label,","Deg_90_Angular_Second_Moment,","Deg_90_Contrast,","Deg_90_Correlation,","Deg_90_Inverse_Difference_Moment,","Deg_90_Entropy,");

for (i=0; i<list1.length; i++) {
	path = dir1+list1[i]; // set the path 
	open(path); //open image
	obj_image=getTitle();
	baseNameEnd=indexOf(obj_image, ".tif"); 
	baseName=substring(obj_image, 0, baseNameEnd);
	if(nSlices>1){
	run("Z Project...", " projection=[Sum Slices]"); //Z project
	run("8-bit"); // convert to 8-bit image
	run("GLCM Texture", "enter="+15+ " select=[90 degrees] angular contrast correlation inverse entropy");
    asm = getResult("Angular Second Moment",0); 
    contrast = getResult("Contrast",0);
    correlation = getResult("Correlation",0);
    idm = getResult("Inverse Difference Moment   ",0); //Extra spaces needed due to source code error
    entropy = getResult("Entropy",0);
    print(baseName,",",asm,",",contrast,",",correlation,",",idm,",",entropy);
	}
	run("Close All"); //close the images
}
selectWindow("Log");  //select Log-window
saveAs("Text", dirres + "gclm_15_deg90.csv");

print("\\Clear"); //empties the Log
print("Label,","Deg_90_Angular_Second_Moment,","Deg_90_Contrast,","Deg_90_Correlation,","Deg_90_Inverse_Difference_Moment,","Deg_90_Entropy,");

for (i=0; i<list1.length; i++) {
	path = dir1+list1[i]; // set the path 
	open(path); //open image
	obj_image=getTitle();
	baseNameEnd=indexOf(obj_image, ".tif"); 
	baseName=substring(obj_image, 0, baseNameEnd);
	if(nSlices>1){
	run("Z Project...", " projection=[Sum Slices]"); //Z project
	run("8-bit"); // convert to 8-bit image
	run("GLCM Texture", "enter="+25+ " select=[90 degrees] angular contrast correlation inverse entropy");
    asm = getResult("Angular Second Moment",0); 
    contrast = getResult("Contrast",0);
    correlation = getResult("Correlation",0);
    idm = getResult("Inverse Difference Moment   ",0); //Extra spaces needed due to source code error
    entropy = getResult("Entropy",0);
    print(baseName,",",asm,",",contrast,",",correlation,",",idm,",",entropy);
	}
	run("Close All"); //close the images
}
selectWindow("Log");  //select Log-window
saveAs("Text", dirres + "gclm_25_deg90.csv");

print("\\Clear"); //empties the Log
print("Label,","Deg_90_Angular_Second_Moment,","Deg_90_Contrast,","Deg_90_Correlation,","Deg_90_Inverse_Difference_Moment,","Deg_90_Entropy,");

for (i=0; i<list1.length; i++) {
	path = dir1+list1[i]; // set the path 
	open(path); //open image
	obj_image=getTitle();
	baseNameEnd=indexOf(obj_image, ".tif"); 
	baseName=substring(obj_image, 0, baseNameEnd);
	if(nSlices>1){
	run("Z Project...", " projection=[Sum Slices]"); //Z project
	run("8-bit"); // convert to 8-bit image
	run("GLCM Texture", "enter="+35+ " select=[90 degrees] angular contrast correlation inverse entropy");
    asm = getResult("Angular Second Moment",0); 
    contrast = getResult("Contrast",0);
    correlation = getResult("Correlation",0);
    idm = getResult("Inverse Difference Moment   ",0); //Extra spaces needed due to source code error
    entropy = getResult("Entropy",0);
    print(baseName,",",asm,",",contrast,",",correlation,",",idm,",",entropy);
	}
	run("Close All"); //close the images
}
selectWindow("Log");  //select Log-window
saveAs("Text", dirres + "gclm_35_deg90.csv");

print("\\Clear"); //empties the Log
print("Label,","Deg_90_Angular_Second_Moment,","Deg_90_Contrast,","Deg_90_Correlation,","Deg_90_Inverse_Difference_Moment,","Deg_90_Entropy,");

for (i=0; i<list1.length; i++) {
	path = dir1+list1[i]; // set the path 
	open(path); //open image
	obj_image=getTitle();
	baseNameEnd=indexOf(obj_image, ".tif"); 
	baseName=substring(obj_image, 0, baseNameEnd);
	if(nSlices>1){
	run("Z Project...", " projection=[Sum Slices]"); //Z project
	run("8-bit"); // convert to 8-bit image
	run("GLCM Texture", "enter="+45+ " select=[90 degrees] angular contrast correlation inverse entropy");
    asm = getResult("Angular Second Moment",0); 
    contrast = getResult("Contrast",0);
    correlation = getResult("Correlation",0);
    idm = getResult("Inverse Difference Moment   ",0); //Extra spaces needed due to source code error
    entropy = getResult("Entropy",0);
    print(baseName,",",asm,",",contrast,",",correlation,",",idm,",",entropy);
	}
	run("Close All"); //close the images
}
selectWindow("Log");  //select Log-window
saveAs("Text", dirres + "gclm_45_deg90.csv");

print("\\Clear"); //empties the Log
print("Label,","Deg_90_Angular_Second_Moment,","Deg_90_Contrast,","Deg_90_Correlation,","Deg_90_Inverse_Difference_Moment,","Deg_90_Entropy,");

for (i=0; i<list1.length; i++) {
	path = dir1+list1[i]; // set the path 
	open(path); //open image
	obj_image=getTitle();
	baseNameEnd=indexOf(obj_image, ".tif"); 
	baseName=substring(obj_image, 0, baseNameEnd);
	if(nSlices>1){
	run("Z Project...", " projection=[Sum Slices]"); //Z project
	run("8-bit"); // convert to 8-bit image
	run("GLCM Texture", "enter="+100+ " select=[90 degrees] angular contrast correlation inverse entropy");
    asm = getResult("Angular Second Moment",0); 
    contrast = getResult("Contrast",0);
    correlation = getResult("Correlation",0);
    idm = getResult("Inverse Difference Moment   ",0); //Extra spaces needed due to source code error
    entropy = getResult("Entropy",0);
    print(baseName,",",asm,",",contrast,",",correlation,",",idm,",",entropy);
	}
	run("Close All"); //close the images
}
selectWindow("Log");  //select Log-window
saveAs("Text", dirres + "gclm_100_deg90.csv");
////////////
print("\\Clear"); //empties the Log
print("Label,","Deg_180_Angular_Second_Moment,","Deg_180_Contrast,","Deg_180_Correlation,","Deg_180_Inverse_Difference_Moment,","Deg_180_Entropy,");
for (i=0; i<list1.length; i++) {
	path = dir1+list1[i]; // set the path 
	open(path); //open image
	obj_image=getTitle();
	baseNameEnd=indexOf(obj_image, ".tif"); 
	baseName=substring(obj_image, 0, baseNameEnd);
	if(nSlices>1){
	run("Z Project...", " projection=[Sum Slices]"); //Z project
	run("8-bit"); // convert to 8-bit image
	run("GLCM Texture", "enter="+5+ " select=[180 degrees] angular contrast correlation inverse entropy");
    asm = getResult("Angular Second Moment",0); 
    contrast = getResult("Contrast",0);
    correlation = getResult("Correlation",0);
    idm = getResult("Inverse Difference Moment   ",0); //Extra spaces needed due to source code error
    entropy = getResult("Entropy",0);
    print(baseName,",",asm,",",contrast,",",correlation,",",idm,",",entropy);
	}
	run("Close All"); //close the images
}
selectWindow("Log");  //select Log-window
saveAs("Text", dirres + "gclm_5_deg180.csv");

print("\\Clear"); //empties the Log
print("Label,","Deg_180_Angular_Second_Moment,","Deg_180_Contrast,","Deg_180_Correlation,","Deg_180_Inverse_Difference_Moment,","Deg_180_Entropy,");
for (i=0; i<list1.length; i++) {
	path = dir1+list1[i]; // set the path 
	open(path); //open image
	obj_image=getTitle();
	baseNameEnd=indexOf(obj_image, ".tif"); 
	baseName=substring(obj_image, 0, baseNameEnd);
	if(nSlices>1){
	run("Z Project...", " projection=[Sum Slices]"); //Z project
	run("8-bit"); // convert to 8-bit image
	run("GLCM Texture", "enter="+15+ " select=[180 degrees] angular contrast correlation inverse entropy");
    asm = getResult("Angular Second Moment",0); 
    contrast = getResult("Contrast",0);
    correlation = getResult("Correlation",0);
    idm = getResult("Inverse Difference Moment   ",0); //Extra spaces needed due to source code error
    entropy = getResult("Entropy",0);
    print(baseName,",",asm,",",contrast,",",correlation,",",idm,",",entropy);
	}
	run("Close All"); //close the images
}
selectWindow("Log");  //select Log-window
saveAs("Text", dirres + "gclm_15_deg180.csv");

print("\\Clear"); //empties the Log
print("Label,","Deg_180_Angular_Second_Moment,","Deg_180_Contrast,","Deg_180_Correlation,","Deg_180_Inverse_Difference_Moment,","Deg_180_Entropy,");
for (i=0; i<list1.length; i++) {
	path = dir1+list1[i]; // set the path 
	open(path); //open image
	obj_image=getTitle();
	baseNameEnd=indexOf(obj_image, ".tif"); 
	baseName=substring(obj_image, 0, baseNameEnd);
	if(nSlices>1){
	run("Z Project...", " projection=[Sum Slices]"); //Z project
	run("8-bit"); // convert to 8-bit image
	run("GLCM Texture", "enter="+25+ " select=[180 degrees] angular contrast correlation inverse entropy");
    asm = getResult("Angular Second Moment",0); 
    contrast = getResult("Contrast",0);
    correlation = getResult("Correlation",0);
    idm = getResult("Inverse Difference Moment   ",0); //Extra spaces needed due to source code error
    entropy = getResult("Entropy",0);
    print(baseName,",",asm,",",contrast,",",correlation,",",idm,",",entropy);
	}
	run("Close All"); //close the images
}
selectWindow("Log");  //select Log-window
saveAs("Text", dirres + "gclm_25_deg180.csv");

print("\\Clear"); //empties the Log
print("Label,","Deg_180_Angular_Second_Moment,","Deg_180_Contrast,","Deg_180_Correlation,","Deg_180_Inverse_Difference_Moment,","Deg_180_Entropy,");
for (i=0; i<list1.length; i++) {
	path = dir1+list1[i]; // set the path 
	open(path); //open image
	obj_image=getTitle();
	baseNameEnd=indexOf(obj_image, ".tif"); 
	baseName=substring(obj_image, 0, baseNameEnd);
	if(nSlices>1){
	run("Z Project...", " projection=[Sum Slices]"); //Z project
	run("8-bit"); // convert to 8-bit image
	run("GLCM Texture", "enter="+35+ " select=[180 degrees] angular contrast correlation inverse entropy");
    asm = getResult("Angular Second Moment",0); 
    contrast = getResult("Contrast",0);
    correlation = getResult("Correlation",0);
    idm = getResult("Inverse Difference Moment   ",0); //Extra spaces needed due to source code error
    entropy = getResult("Entropy",0);
    print(baseName,",",asm,",",contrast,",",correlation,",",idm,",",entropy);
	}
	run("Close All"); //close the images
}
selectWindow("Log");  //select Log-window
saveAs("Text", dirres + "gclm_35_deg180.csv");

print("\\Clear"); //empties the Log
print("Label,","Deg_180_Angular_Second_Moment,","Deg_180_Contrast,","Deg_180_Correlation,","Deg_180_Inverse_Difference_Moment,","Deg_180_Entropy,");
for (i=0; i<list1.length; i++) {
	path = dir1+list1[i]; // set the path 
	open(path); //open image
	obj_image=getTitle();
	baseNameEnd=indexOf(obj_image, ".tif"); 
	baseName=substring(obj_image, 0, baseNameEnd);
	if(nSlices>1){
	run("Z Project...", " projection=[Sum Slices]"); //Z project
	run("8-bit"); // convert to 8-bit image
	run("GLCM Texture", "enter="+45+ " select=[180 degrees] angular contrast correlation inverse entropy");
    asm = getResult("Angular Second Moment",0); 
    contrast = getResult("Contrast",0);
    correlation = getResult("Correlation",0);
    idm = getResult("Inverse Difference Moment   ",0); //Extra spaces needed due to source code error
    entropy = getResult("Entropy",0);
    print(baseName,",",asm,",",contrast,",",correlation,",",idm,",",entropy);
	}
	run("Close All"); //close the images
}
selectWindow("Log");  //select Log-window
saveAs("Text", dirres + "gclm_45_deg180.csv");

print("\\Clear"); //empties the Log
print("Label,","Deg_180_Angular_Second_Moment,","Deg_180_Contrast,","Deg_180_Correlation,","Deg_180_Inverse_Difference_Moment,","Deg_180_Entropy,");
for (i=0; i<list1.length; i++) {
	path = dir1+list1[i]; // set the path 
	open(path); //open image
	obj_image=getTitle();
	baseNameEnd=indexOf(obj_image, ".tif"); 
	baseName=substring(obj_image, 0, baseNameEnd);
	if(nSlices>1){
	run("Z Project...", " projection=[Sum Slices]"); //Z project
	run("8-bit"); // convert to 8-bit image
	run("GLCM Texture", "enter="+100+ " select=[180 degrees] angular contrast correlation inverse entropy");
    asm = getResult("Angular Second Moment",0); 
    contrast = getResult("Contrast",0);
    correlation = getResult("Correlation",0);
    idm = getResult("Inverse Difference Moment   ",0); //Extra spaces needed due to source code error
    entropy = getResult("Entropy",0);
    print(baseName,",",asm,",",contrast,",",correlation,",",idm,",",entropy);
	}
	run("Close All"); //close the images
}
selectWindow("Log");  //select Log-window
saveAs("Text", dirres + "gclm_100_deg180.csv");
////////////
print("\\Clear"); //empties the Log
print("Label,","Deg_270_Angular_Second_Moment,","Deg_270_Contrast,","Deg_270_Correlation,","Deg_270_Inverse_Difference_Moment,","Deg_270_Entropy,");
for (i=0; i<list1.length; i++) {
	path = dir1+list1[i]; // set the path 
	open(path); //open image
	obj_image=getTitle();
	baseNameEnd=indexOf(obj_image, ".tif"); 
	baseName=substring(obj_image, 0, baseNameEnd);
	if(nSlices>1){
	run("Z Project...", " projection=[Sum Slices]"); //Z project
	run("8-bit"); // convert to 8-bit image
	run("GLCM Texture", "enter="+5+ " select=[270 degrees] angular contrast correlation inverse entropy");
    asm = getResult("Angular Second Moment",0); 
    contrast = getResult("Contrast",0);
    correlation = getResult("Correlation",0);
    idm = getResult("Inverse Difference Moment   ",0); //Extra spaces needed due to source code error
    entropy = getResult("Entropy",0);
    print(baseName,",",asm,",",contrast,",",correlation,",",idm,",",entropy);
	}
	run("Close All"); //close the images
}
selectWindow("Log");  //select Log-window
saveAs("Text", dirres + "gclm_5_deg270.csv");

print("\\Clear"); //empties the Log
print("Label,","Deg_270_Angular_Second_Moment,","Deg_270_Contrast,","Deg_270_Correlation,","Deg_270_Inverse_Difference_Moment,","Deg_270_Entropy,");
for (i=0; i<list1.length; i++) {
	path = dir1+list1[i]; // set the path 
	open(path); //open image
	obj_image=getTitle();
	baseNameEnd=indexOf(obj_image, ".tif"); 
	baseName=substring(obj_image, 0, baseNameEnd);
	if(nSlices>1){
	run("Z Project...", " projection=[Sum Slices]"); //Z project
	run("8-bit"); // convert to 8-bit image
	run("GLCM Texture", "enter="+15+ " select=[270 degrees] angular contrast correlation inverse entropy");
    asm = getResult("Angular Second Moment",0); 
    contrast = getResult("Contrast",0);
    correlation = getResult("Correlation",0);
    idm = getResult("Inverse Difference Moment   ",0); //Extra spaces needed due to source code error
    entropy = getResult("Entropy",0);
    print(baseName,",",asm,",",contrast,",",correlation,",",idm,",",entropy);
	}
	run("Close All"); //close the images
}
selectWindow("Log");  //select Log-window
saveAs("Text", dirres + "gclm_15_deg270.csv");

print("\\Clear"); //empties the Log
print("Label,","Deg_270_Angular_Second_Moment,","Deg_270_Contrast,","Deg_270_Correlation,","Deg_270_Inverse_Difference_Moment,","Deg_270_Entropy,");
for (i=0; i<list1.length; i++) {
	path = dir1+list1[i]; // set the path 
	open(path); //open image
	obj_image=getTitle();
	baseNameEnd=indexOf(obj_image, ".tif"); 
	baseName=substring(obj_image, 0, baseNameEnd);
	if(nSlices>1){
	run("Z Project...", " projection=[Sum Slices]"); //Z project
	run("8-bit"); // convert to 8-bit image
	run("GLCM Texture", "enter="+25+ " select=[270 degrees] angular contrast correlation inverse entropy");
    asm = getResult("Angular Second Moment",0); 
    contrast = getResult("Contrast",0);
    correlation = getResult("Correlation",0);
    idm = getResult("Inverse Difference Moment   ",0); //Extra spaces needed due to source code error
    entropy = getResult("Entropy",0);
    print(baseName,",",asm,",",contrast,",",correlation,",",idm,",",entropy);
	}
	run("Close All"); //close the images
}
selectWindow("Log");  //select Log-window
saveAs("Text", dirres + "gclm_25_deg270.csv");

print("\\Clear"); //empties the Log
print("Label,","Deg_270_Angular_Second_Moment,","Deg_270_Contrast,","Deg_270_Correlation,","Deg_270_Inverse_Difference_Moment,","Deg_270_Entropy,");
for (i=0; i<list1.length; i++) {
	path = dir1+list1[i]; // set the path 
	open(path); //open image
	obj_image=getTitle();
	baseNameEnd=indexOf(obj_image, ".tif"); 
	baseName=substring(obj_image, 0, baseNameEnd);
	if(nSlices>1){
	run("Z Project...", " projection=[Sum Slices]"); //Z project
	run("8-bit"); // convert to 8-bit image
	run("GLCM Texture", "enter="+35+ " select=[270 degrees] angular contrast correlation inverse entropy");
    asm = getResult("Angular Second Moment",0); 
    contrast = getResult("Contrast",0);
    correlation = getResult("Correlation",0);
    idm = getResult("Inverse Difference Moment   ",0); //Extra spaces needed due to source code error
    entropy = getResult("Entropy",0);
    print(baseName,",",asm,",",contrast,",",correlation,",",idm,",",entropy);
	}
	run("Close All"); //close the images
}
selectWindow("Log");  //select Log-window
saveAs("Text", dirres + "gclm_35_deg270.csv");

print("\\Clear"); //empties the Log
print("Label,","Deg_270_Angular_Second_Moment,","Deg_270_Contrast,","Deg_270_Correlation,","Deg_270_Inverse_Difference_Moment,","Deg_270_Entropy,");
for (i=0; i<list1.length; i++) {
	path = dir1+list1[i]; // set the path 
	open(path); //open image
	obj_image=getTitle();
	baseNameEnd=indexOf(obj_image, ".tif"); 
	baseName=substring(obj_image, 0, baseNameEnd);
	if(nSlices>1){
	run("Z Project...", " projection=[Sum Slices]"); //Z project
	run("8-bit"); // convert to 8-bit image
	run("GLCM Texture", "enter="+45+ " select=[270 degrees] angular contrast correlation inverse entropy");
    asm = getResult("Angular Second Moment",0); 
    contrast = getResult("Contrast",0);
    correlation = getResult("Correlation",0);
    idm = getResult("Inverse Difference Moment   ",0); //Extra spaces needed due to source code error
    entropy = getResult("Entropy",0);
    print(baseName,",",asm,",",contrast,",",correlation,",",idm,",",entropy);
	}
	run("Close All"); //close the images
}
selectWindow("Log");  //select Log-window
saveAs("Text", dirres + "gclm_45_deg270.csv");

print("\\Clear"); //empties the Log
print("Label,","Deg_270_Angular_Second_Moment,","Deg_270_Contrast,","Deg_270_Correlation,","Deg_270_Inverse_Difference_Moment,","Deg_270_Entropy,");
for (i=0; i<list1.length; i++) {
	path = dir1+list1[i]; // set the path 
	open(path); //open image
	obj_image=getTitle();
	baseNameEnd=indexOf(obj_image, ".tif"); 
	baseName=substring(obj_image, 0, baseNameEnd);
	if(nSlices>1){
	run("Z Project...", " projection=[Sum Slices]"); //Z project
	run("8-bit"); // convert to 8-bit image
	run("GLCM Texture", "enter="+100+ " select=[270 degrees] angular contrast correlation inverse entropy");
    asm = getResult("Angular Second Moment",0); 
    contrast = getResult("Contrast",0);
    correlation = getResult("Correlation",0);
    idm = getResult("Inverse Difference Moment   ",0); //Extra spaces needed due to source code error
    entropy = getResult("Entropy",0);
    print(baseName,",",asm,",",contrast,",",correlation,",",idm,",",entropy);
	}
	run("Close All"); //close the images
}
selectWindow("Log");  //select Log-window
saveAs("Text", dirres + "gclm_100_deg270.csv");

