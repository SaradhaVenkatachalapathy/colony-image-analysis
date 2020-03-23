/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////  SCRIPT TO MAKE ALL IMAGES IN A FOLDER THE SAME SIZE BY PADDING
///////  WRITTEN BY: SARADHA VENKATACHALAPATHY                                                                                                                   
///////  ASSUMPTIONS: 
///////  DESCRIPTION: It gets the the maximum width and height of all the images in the input folder. Then the scripts resets the canvas of all images to 10+ maximum width/height
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

dirsa = getArgument(); //get directory of files

filename1 = getFileList(dirsa); 
setBatchMode(true);
//get the maximum width and height of all the images in the input folder
wid=0;
hei=0;
for (i=0; i<filename1.length; i++) { 
      path=dirsa+filename1[i];
      open(path);
	  if(getWidth() > wid){
	  	wid=getWidth();
	  }
	  if(getHeight() > hei){
	  	hei=getHeight();
	  }
      run("Close All"); 
       
	
}
//reset the canvas of all images to 10+ maximum width/height
wid=wid+10;
hei=hei+10;

for (i=0; i<filename1.length; i++) { 
      path=dirsa+filename1[i];
      open(path);
	  run("Canvas Size...", "width=wid height=hei position=Center zero");
	  saveAs("tiff", dirsa + getTitle()); 
      run("Close All"); 
	
}