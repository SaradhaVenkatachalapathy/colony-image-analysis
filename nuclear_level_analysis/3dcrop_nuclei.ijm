/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////  SCRIPT TO IDENTIFY SEGMENT 3D NUCLEI FROM BINARY IMAGES      											                   
///////  WRITTEN BY: SARADHA VENKATACHALAPATHY                                                                                                                   
///////  ASSUMPTIONS: The input image is a confocal zstack and the channel 1 contains the nucleus 
///////  DESCRIPTION: The script accepts the directory to a folder containing rawimages folder and the number of channels in the raw image that need to be cropped. 
///////				  It creates 5 subfolders: one for storing the image afer 3d segmentation with identified nuclei and 4 for storing channels of the raw image cropped to represent each nuclei.
///////				  For all images in the after watersheld folder, it identifies 3D objects of size ranging from 50-2000 cu.microns. If there are nuclei in the image, then it saves 
/////// 			 the image containing the identified the objects, adds the images to the 3D roi manager. Then
///////				  for each identified object in the image, it gets the cordinates and crop the image along the 3D bounding box and thresholds to obtain the binary mask of the object. Then for each 
///////				  channel that is required, the raw image is cropped along the bounding box for each object and retains only sections of the image that are also present in the binary mask. This
///////				  ensures that in the cropped image there are not overlapping sections from adjacent nuclei. 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

args = split(getArgument()," "); 	// read in the array of argurments that are separated by a space
dirsa = args[0];	// The first argument is the directory
nchannels=args[1];	// The second argument is the number of channels
ch2_name=args[2];	// The third argument is the name of the protein in channel 2
ch3_name=args[3];	// The fourth argument is the name of the protein in channel 3
ch4_name=args[4];	// The fifth argument is the name of the protein in channel 4

setBatchMode(true);	// run in batch mode to save time and memory

///creating new subdirectories for storing the processed data//
dirb= dirsa + "3d obj nuclei"+ File.separator;	// define the path to a new folder that will contain the image afer 3d segmentation to indentify nuclei
dir= dirsa + "indivisual_nuclei_DNA"+ File.separator;		// define the path to a new folder that will contain the first channel of the raw image cropped to represent each nuclei
dirc2= dirsa + "indivisual_nuclei"+ch2_name+ File.separator; 	// define the path to a new folder that will contain the second channel of the raw image cropped to represent each nuclei
dirc3= dirsa + "indivisual_nuclei"+ch3_name+ File.separator; 	// define the path to a new folder that will contain the third channel of the raw image cropped to represent each nuclei
dirc4= dirsa + "indivisual_nuclei"+ch4_name+ File.separator; 	// define the path to a new folder that will contain the fourth channel of the raw image cropped to represent each nuclei
dir1= dirsa + "data"+ File.separator;
// create the paths defined above 
File.makeDirectory(dirb);
File.makeDirectory(dir); 
File.makeDirectory(dir1); 
if(nchannels>1){ 
	File.makeDirectory(dirc2); 
}
if(nchannels>2){ 
	File.makeDirectory(dirc3); 
}
if(nchannels>3){ 
	File.makeDirectory(dirc4); 
}

//set paths
//dirsa=getDirectory("Please choose the source directory");
dirw= dirsa + "after watershed"+ File.separator;
dirraw= dirsa + "rawimages"+ File.separator;

//set the options for 3D objects counter and 3D ROI manager 
run("3D OC Options", "volume surface nb_of_obj._voxels nb_of_surf._voxels integrated_density mean_gray_value std_dev_gray_value median_gray_value minimum_gray_value maximum_gray_value centroid mean_distance_to_surface std_dev_distance_to_surface median_distance_to_surface centre_of_mass bounding_box dots_size=5 font_size=20 redirect_to=none");
run("3D Manager Options", "volume surface compactness fit_ellipse 3d_moments integrated_density mean_grey_value std_dev_grey_value mode_grey_value feret minimum_grey_value maximum_grey_value centroid_(pix) centroid_(unit) distance_to_surface centre_of_mass_(pix) centre_of_mass_(unit) bounding_box radial_distance surface_contact closest use distance_between_centers=10 distance_max_contact=1.80");


setBatchMode(true);
filenames=getFileList(dirw);
filenamesraw=getFileList(dirraw);
baseName=newArray(filenames.length);

run("3D Manager");	// Open the 3D ROI manager
for(f=0;f<filenames.length;f++){	//loop through all images

	//open the watershed image and obtain the base name
	path=dirraw+filenamesraw[f];
	open(dirw+filenames[f]);
	imgName=getTitle(); 
	baseNameEnd=indexOf(imgName, ".tiff"); 
	baseName[f]=substring(imgName, 0, baseNameEnd); 
	
	//set the minimum and maximum volume of the nucleus
	getVoxelSize(width, height, depth, unit);
	a=2000/(width*height*depth);// maximum volume is 2000 cu.microns
	a_1=50/(width*height*depth);// minimum volume is 50 cu.microns
	
	run("3D Objects Counter", "threshold=128 slice=12 min.=a_1 max.=a objects"); //Identify 3d objctions with Size filter
	

	run("Z Project...", "projection=[Max Intensity]");
	run("Measure");
	close();
	if(getResult("Mean",0)>0){
		saveAs("Tiff", dirb+baseName[f]+".tiff"); 
		obj_image=getTitle();
		Ext.Manager3D_AddImage();
		Ext.Manager3D_Count(nb_obj); 
		run("Close All");

	
	zlow=newArray(nb_obj);
	zhigh=newArray(nb_obj);
	wid=newArray(nb_obj);
	heig=newArray(nb_obj);
	xv=newArray(nb_obj);
	yv=newArray(nb_obj);
	name=newArray(nb_obj);
	print("\\Clear"); 
	
	for(p=0;p<nb_obj;p++){
		t=p+1;
		Ext.Manager3D_Bounding3D(p,x0,x1,y0,y1,z0,z1);
		zlows=z0+1;
		zhighs=z1+1;
		wids=x1-x0;
		heigs=y1-y0;
		
		open(dirb+obj_image);
		makeRectangle(x0, y0, wids, heigs);
		run("Duplicate...", "duplicate range=zlows-zhighs");
		run("3D Convex Hull");//convex hull
		setThreshold(t, t);		
		run("Make Binary", "method=Default background=Default");
		img1=getTitle();
		Ext.Manager3D_CloseResult("M");
		
		//Channel1
		run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT z_begin=zlows z_end=zhighs c_begin=1 c_end=1 c_step=1");
		makeRectangle(x0, y0, wids, heigs);
		run("Duplicate...", "duplicate");
		img2=getTitle();
		imageCalculator("AND create stack", img1,img2);
		run("Invert LUT");
		names=baseName[f]+"_nucleus_"+x0+"_"+y0+"_"+z0+"_"+z1;
		saveAs("tiff",dir+names);
		imgName2=getTitle(); 
		Ext.Manager3D_CloseResult("M");
		
		selectImage(img1); 
		close("\\Others");

		if(nchannels>1){
		//Channel2
		run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT z_begin=zlows z_end=zhighs c_begin=2 c_end=2 c_step=1");
		makeRectangle(x0, y0, wids, heigs);
		run("Duplicate...", "duplicate");
		img2=getTitle();
		imageCalculator("AND create stack", img1,img2);
		run("Invert LUT");
		names=baseName[f]+"_nucleus_"+x0+"_"+y0+"_"+z0+"_"+z1;
		saveAs("tiff",dirc2+names);
		imgName2=getTitle(); 
		Ext.Manager3D_CloseResult("M");
		
		selectImage(img1); 
		close("\\Others");
		}

		if(nchannels>2){
		//Channel3
		run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT z_begin=zlows z_end=zhighs c_begin=3 c_end=3 c_step=1");
		makeRectangle(x0, y0, wids, heigs);
		run("Duplicate...", "duplicate");
		img2=getTitle();
		imageCalculator("AND create stack", img1,img2);
		run("Invert LUT");
		names=baseName[f]+"_nucleus_"+x0+"_"+y0+"_"+z0+"_"+z1;
		saveAs("tiff",dirc3+names);
		imgName2=getTitle(); 
		Ext.Manager3D_CloseResult("M");
		
		selectImage(img1); 
		close("\\Others");
		}
		
		if(nchannels>3){
		//Channel4
		run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT z_begin=zlows z_end=zhighs c_begin=4 c_end=4 c_step=1");
		makeRectangle(x0, y0, wids, heigs);
		run("Duplicate...", "duplicate");
		img2=getTitle();
		imageCalculator("AND create stack", img1,img2);
		run("Invert LUT");
		names=baseName[f]+"_nucleus_"+x0+"_"+y0+"_"+z0+"_"+z1;
		saveAs("tiff",dirc4+names);
		imgName2=getTitle(); 
		Ext.Manager3D_CloseResult("M");
		
		
		selectImage(img1); 
		close("\\Others");
		}
		
		
		
		xv[p]=x0;
		yv[p]=y0;
		wid[p]=wids;
		heig[p]=heigs;
		name[p]=names;
		zlow[p]=z0+1;
		zhigh[p]=z1+1;
		run("Select None");
     	run("Close All");


     	
	}


	saveAs("Results",dir1+baseName[f]+"_objects_statistics.csv" ); 
	run("Clear Results");
	for(ol=0; ol<nb_obj; ol++){
		setResult("Name",ol,name[ol]);
		setResult("BX0",ol,xv[ol]);
		setResult("BY0",ol,yv[ol]);
		setResult("Width",ol,wid[ol]);
		setResult("Height",ol,heig[ol]);
		setResult("zstart",ol,zlow[ol]);
		setResult("zend",ol,zhigh[ol]);
	}
	saveAs("Results",dir1+baseName[f]+"_log.csv" ); 	
	}
	
	
	run("Close All");
    Ext.Manager3D_Reset();
   	run("Clear Results");
}


Ext.Manager3D_Close();
