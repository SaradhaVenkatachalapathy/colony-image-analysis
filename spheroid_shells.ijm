/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////  SCRIPT TO CREATE SHELLS BY ERODING A 3D OBJECT BY THE SAME DISTANCE IN X,Y AND Z AND MEASURING QUANTITITES IN THESE SHELLS
///////  WRITTEN BY: SARADHA VENKATACHALAPATHY                                                                                                                   
///////  ASSUMPTIONS: The input image is a confocal zstack and the channel 1 contains the nucleus. The raw files' extension is ".nd2". The erosion size is 5 microns in x,y and z
///////  DESCRIPTION: The script accepts the directory to a folder containing rawimages folder and the number of channels in the raw image that need to be analysed and the taget of the stains. 
///////				  It creates 3-6 subfolders: one for storing the geoemtrical features of the identified spheroids and others for storing intensity features of the selected image.
///////				  For all images in the identifed spheroids, it names and measures the geometrical and intensity features that are requiered. It bins the image so that the pixel size is 
///////				  the similar across xy and z and is equal to the erosion size. Then it erodes the object till the volume is is less than or equal to 100 cu.microns, saving the object after each 
//////				  iteration. Then it opens the requiured channels and computes the intensity features required.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

args = split(getArgument()," "); 	// read in the array of argurments that are separated by a space
dirsa = args[0];	// The first argument is the directory
nchannels=args[1];	// The second argument is the number of channels
ch2_name=args[2];	// The third argument is the name of the protein in channel 2
ch3_name=args[3];	// The fourth argument is the name of the protein in channel 3
ch4_name=args[4];	// The fifth argument is the name of the protein in channel 4
erode_size=5	// set step size to be 5 microns

dirsa= dirsa1 + "3d obj sphereoid"+ File.separator; // define the path to the folder containing the image afer 3d segmentation to indentify spheroids
filenames1=getFileList(dirsa); // get the list of names of the images in the 3d objects folder
dir_raw=dirsa1+"rawimages"+File.separator;  // define the path to the folder containing the unprocessed or raw images

///creating new subdirectories for storing the processed data//
dir2=dirsa1+"Shells"+File.separator;  // define the path to a new folder that will contain the following sub-folders for the 3d shells
newDir113 = dir2 + "geometric_measures"+ File.separator;  // define the path to a new folder that will contain the 3d geoemetrical data for the shells
newDir114= dir2 + "DNA" + File.separator; // define the path to a new folder that will contain the 3d intensity data for the spheroids in channel 1/ DNA channel
newDir124= dir2 + ch2_name+ File.separator; // define the path to a new folder that will contain the 3d intensity data for the spheroids in channel 2
newDir134= dir2 + ch3_name+ File.separator; // define the path to a new folder that will contain the 3d intensity data for the spheroids in channel 3 
newDir144= dir2 + ch4_name+ File.separator; // define the path to a new folder that will contain the 3d intensity data for the spheroids in channel 4 
File.makeDirectory(dir2); 
File.makeDirectory(newDir113); 
File.makeDirectory(newDir114); 
if(nchannels>1){ 
	File.makeDirectory(newDir124); 
}
if(nchannels>2){ 
	File.makeDirectory(newDir134); 
}
if(nchannels>3){ 
	File.makeDirectory(newDir144); 
}
setBatchMode(true); // process in batch mode

//set the options for 3D ROI manager, 3D objects counter and Measurements
run("3D OC Options", "volume surface nb_of_obj._voxels nb_of_surf._voxels centroid std_dev_distance_to_surface bounding_box show_masked_image_(redirection_requiered) dots_size=5 font_size=10 redirect_to=none");
run("3D Manager Options", "volume surface compactness fit_ellipse 3d_moments integrated_density mean_grey_value std_dev_grey_value mode_grey_value feret minimum_grey_value maximum_grey_value centroid_(pix) centroid_(unit) distance_to_surface centre_of_mass_(pix) centre_of_mass_(unit) bounding_box radial_distance surface_contact closest use distance_between_centers=10 distance_max_contact=1.80");
run("Set Measurements...", "area mean standard modal min centroid center perimeter bounding fit shape feret's integrated median skewness kurtosis area_fraction limit display redirect=None decimal=3");

for(i=0; i<filenames1.length; i++){
	
		// open the image containing 3d objects 
		path3=dirsa+filenames1[i];
		open(path3);
		run("8-bit"); // convert to 8 bit image
		
		// extract the the name of the file and remove file extension
		imgName1=getTitle(); //get the title and assign it. it will be a character string
		baseNameEnd=indexOf(imgName1, ".tiff");  //find the index of the string at with ".tiff" first appears NOTE: change this to correct file extension if working with other files 
		baseName=substring(imgName1, 0, baseNameEnd); // get the substring such that the file extention gets removed from the character string 

	
		//bin the image so that the pixel size is the similar across xy and z
		getVoxelSize(x, y, h, micron);
		fx=round(erode_size/x);	// find the number of pixels to bin in x and y
		fh=round(erode_size/h);	// find the number of pixels to bin in z
		run("Bin...", "x=fx y=fx z=fh bin=Average"); // bin the image
		//getVoxelSize(x, y, h, micron); // get the new voxel dimensions
		//print(imgName1,",",x,",",y,",",h); // display the new dimension
		
		//get filename, generate 3d object,  get the volume, 
		label=File.getName(path3); //get the file name
		selectWindow(imgName1); // select the 3d objects image
		setThreshold(1, 4095); // threshold
		run("Convert to Mask", "method=Default background=Default"); // conver to binary image
		run("3D Objects Counter", "threshold=1 slice=6 min.=0 max.=1000000000000 statistics objects"); // identify the object again
		selectWindow(imgName1); // select the image
		Vol=0;
		for(l=0;l<nResults;l++){
			Vol=Vol+getResult("Volume (micron^3)",0); // get the volume of the spheroid
		}

		//Add the objects to the 3dmanager, in case more than one spheroid exisits, merge and then name the object
		run("3D Manager"); // open 3d roi manager
		Ext.Manager3D_AddImage(); // add the image
		Ext.Manager3D_Count(nb_obj);  // get the nuumber of rois
		if(nb_obj>1){ // if more than one object exisit, merge them 
			Ext.Manager3D_SelectAll();
			Ext.Manager3D_Merge();

		}
		Ext.Manager3D_Select(0); //select the object
		Ext.Manager3D_Rename("spheroid_0"); // rename it nucleus_0
		
		//erode the object till the volume is is less than or equal to 100 cu.microns. 
		p=0;
		while(Vol>100){
			selectWindow(imgName1);
			run("Erode (3D)", "iso=255"); // erode the object in 3d by one unit 
			run("Z Project...", "projection=[Max Intensity]"); // project the eroded image
			//run("Threshold...");
			setThreshold(1, 4294967296); // set the threshold
			run("Measure");	// measure properties
			a=getResult("Area",0); // get the area
			
			if(a>0 ){ // if object exists, area will be greater than 0
				run("Close"); //close the projected image
				selectWindow(imgName1);// select eroded image
				run("3D Objects Counter", "threshold=1 min.=0 statistics objects"); // add the object
				Vol=0; 
				for(l=0;l<nResults;l++){
					Vol=Vol+getResult("Volume (micron^3)",0); // get the volme
				}
				Ext.Manager3D_AddImage(); // add the eroded object
				p=p+1; 
				print(getTitle(),",",Vol,",",p); // print the volume
				Ext.Manager3D_Count(nb_obj); 
				// if the number of objects is greater than expected, then merge the extra objects
				if(nb_obj>(p+1)){
					Ext.Manager3D_MultiSelect();
					for(k=p;k<nb_obj;k++){
						Ext.Manager3D_Select(k);
					}
					Ext.Manager3D_Merge();
					Ext.Manager3D_DeselectAll();		
				}
				Ext.Manager3D_Count(nb_obj); 
				Ext.Manager3D_Select((nb_obj-1));
				Ext.Manager3D_Rename("spheroid_"+p); // rename the object
				Ext.Manager3D_DeselectAll();
			}
			else{
				Vol=0;
			}
		}
		Ext.Manager3D_Count(p); 
		run("Close All");

		name=baseName+".nd2";
		path3=dir_raw+name; // set the path to the corresponding raw image
		//channel 1
		run("Bio-Formats", "open=path3 color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT c_begin=1 c_end=1 c_step=1");
		run("8-bit");
		//bin the image so that the pixel size is the similar across xy and z
		getVoxelSize(x, y, h, micron);
		fx=round(erode_size/x);
		fh=round(erode_size/h);
		run("Bin...", "x=fx y=fx z=fh bin=Average");
		
		// measure geometric properties
		Ext.Manager3D_DeselectAll();
		Ext.Manager3D_SelectAll();
    	Ext.Manager3D_Measure();
		Ext.Manager3D_SaveResult("M",newDir113+filenames1[i]+".tsv");
		Ext.Manager3D_CloseResult("M");

		// measuure DNA intensity
		Ext.Manager3D_DeselectAll();
		Ext.Manager3D_SelectAll();
		Ext.Manager3D_Quantif();
		Ext.Manager3D_SaveQuantif(newDir114+filenames1[i]+".tsv");
		Ext.Manager3D_CloseResult("Q");
		run("Close All");

		if(nchannels>1){ 
			//channel 2
			run("Bio-Formats", "open=path3 color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT c_begin=2 c_end=2 c_step=1");
			run("8-bit");
			//bin the image so that the pixel size is the similar across xy and z
			getVoxelSize(x, y, h, micron);
			fx=round(erode_size/x);
			fh=round(erode_size/h);
			run("Bin...", "x=fx y=fx z=fh bin=Average");
			
			Ext.Manager3D_DeselectAll();
			Ext.Manager3D_SelectAll();
			Ext.Manager3D_Quantif();
			Ext.Manager3D_SaveQuantif(newDir124+filenames1[i]+".tsv");
			Ext.Manager3D_CloseResult("Q");
			run("Close All");
		}

		if(nchannels>2){ 
			//channel 3
			run("Bio-Formats", "open=path3 color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT c_begin=3 c_end=3 c_step=1");
			run("8-bit");
			//bin the image so that the pixel size is the similar across xy and z
			getVoxelSize(x, y, h, micron);
			fx=round(erode_size/x);
			fh=round(erode_size/h);
			run("Bin...", "x=fx y=fx z=fh bin=Average");
		
			Ext.Manager3D_DeselectAll();
			Ext.Manager3D_SelectAll();
			Ext.Manager3D_Quantif();
			Ext.Manager3D_SaveQuantif(newDir134+filenames1[i]+".tsv");
			Ext.Manager3D_CloseResult("Q");
			run("Close All");
		}

		if(nchannels>3){ 
			//channel 4
			run("Bio-Formats", "open=path3 color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT c_begin=4 c_end=4 c_step=1");
			run("8-bit");
			//bin the image so that the pixel size is the similar across xy and z
			getVoxelSize(x, y, h, micron);
			fx=round(erode_size/x);
			fh=round(erode_size/h);
			run("Bin...", "x=fx y=fx z=fh bin=Average");
		
			Ext.Manager3D_DeselectAll();
			Ext.Manager3D_SelectAll();
			Ext.Manager3D_Quantif();
			Ext.Manager3D_SaveQuantif(newDir144+filenames1[i]+".tsv");
			Ext.Manager3D_CloseResult("Q");
			run("Close All");	
    		
		}

		Ext.Manager3D_Close();
    	run("Clear Results");
    	run("Collect Garbage");
		run("Collect Garbage");
	
	}






	