/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////  SCRIPT TO QUANTIFY GEOMETRICAL AND INTENSITY FEATURES
///////  WRITTEN BY: SARADHA VENKATACHALAPATHY                                                                                                                   
///////  ASSUMPTIONS: The input image is a confocal zstack and the channel 1 contains the nucleus. The raw files' extension is ".nd2". 
///////  DESCRIPTION: The script accepts the directory to a folder containing rawimages folder and the number of channels in the raw image that need to be analysed and the target of the stains. 
///////				  It creates 3-6 subfolders: one for storing the geoemtrical features of the identified nuclei and others for storing intensity features of the selected image.
///////				  For all images in the identifed nuclei, it names and measures the geometrical and intensity features that are requied.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


args = split(getArgument()," "); 	// read in the array of argurments that are separated by a space
dirsa = args[0];	// The first argument is the directory
nchannels=args[1];	// The second argument is the number of channels
ch2_name=args[2];	// The third argument is the name of the protein in channel 2
ch3_name=args[3];	// The fourth argument is the name of the protein in channel 3
ch4_name=args[4];	// The fifth argument is the name of the protein in channel 4


setBatchMode(true);
dirsa1=dirsa+"rawimages"+File.separator;
dirb= dirsa + "3d obj nuclei"+ File.separator;

run("3D Manager Options", "volume surface compactness fit_ellipse 3d_moments integrated_density mean_grey_value std_dev_grey_value mode_grey_value feret minimum_grey_value maximum_grey_value centroid_(pix) centroid_(unit) distance_to_surface centre_of_mass_(pix) centre_of_mass_(unit) bounding_box radial_distance surface_contact closest use distance_between_centers=10 distance_max_contact=1.80");

Dir_res_ellp=dirsa + "3D ellipsoid"+ File.separator;
Dir_res_geo=dirsa + "3D geometerical_simple"+ File.separator;
Dir_res_shape=dirsa + "3D shape measure"+ File.separator;
Dir_res_geo_m=dirsa + "3D geometrical data"+ File.separator;  
Dir_res_int=dirsa + "3D int_data"+ File.separator;  
Dir_res_ch1_int=Dir_res_int + "DNA" + File.separator;  
Dir_res_ch2_int=Dir_res_int + ch2_name +File.separator;  
Dir_res_ch3_int=Dir_res_int + ch3_name +File.separator;  
Dir_res_ch4_int=Dir_res_int + ch4_name +File.separator;  

File.makeDirectory(Dir_res_ellp); 
File.makeDirectory(Dir_res_geo); 
File.makeDirectory(Dir_res_shape); 
File.makeDirectory(Dir_res_geo_m); 
File.makeDirectory(Dir_res_int); 
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

filenames=getFileList(dirb);

for(f=0;f<filenames.length;f++){
	
	open(dirb+filenames[f]);
	obj_image=getTitle();
	baseNameEnd=indexOf(obj_image, ".tiff"); 
	baseName=substring(obj_image, 0, baseNameEnd); 
	
	run("3D Manager");
	Ext.Manager3D_AddImage();
	Ext.Manager3D_Count(nb_obj); 
	for(p=0;p<nb_obj;p++){
		Ext.Manager3D_Bounding3D(p,x0,x1,y0,y1,z0,z1);
		Ext.Manager3D_MonoSelect();
		Ext.Manager3D_Select(p);
		name=baseName+"_nucleus_"+x0+"_"+y0+"_"+z0+"_"+z1;
		Ext.Manager3D_Rename(name);
	}
	Ext.Manager3D_DeselectAll();

	
	Ext.Manager3D_DeselectAll();
	Ext.Manager3D_SelectAll();
	Ext.Manager3D_Measure();
    Ext.Manager3D_SaveMeasure(Dir_res_geo_m+baseName+"_geometric.tsv");
    Ext.Manager3D_CloseResult("M");
	run("Close All");

	name=baseName+".nd2";
	path=dirsa1+name;
	run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT c_begin=1 c_end=1 c_step=1");
	Ext.Manager3D_Quantif();
	Ext.Manager3D_SaveQuantif(Dir_res_ch1_int+baseName+".tsv");
	Ext.Manager3D_CloseResult("Q");
	run("Close All");

	if(nchannels>1){ 
	run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT c_begin=2 c_end=2 c_step=1");
	Ext.Manager3D_Quantif();
	Ext.Manager3D_SaveQuantif(Dir_res_ch2_int+baseName+".tsv");
	Ext.Manager3D_CloseResult("Q");
	run("Close All"); 
	}
	if(nchannels>2){ 
	run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT c_begin=3 c_end=3 c_step=1");
	Ext.Manager3D_Quantif();
	Ext.Manager3D_SaveQuantif(Dir_res_ch3_int+baseName+".tsv");
	Ext.Manager3D_CloseResult("Q");
	run("Close All");
	}
	if(nchannels>3){ 
	run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT c_begin=4 c_end=4 c_step=1");
	Ext.Manager3D_Quantif();
	Ext.Manager3D_SaveQuantif(Dir_res_ch4_int+baseName+".tsv");
	Ext.Manager3D_CloseResult("Q");
	run("Close All");
	}
	
	Ext.Manager3D_Reset();
	Ext.Manager3D_Close();
}
//3D geometric and shape measures

for(f=0;f<filenames.length;f++){
	open(dirb+filenames[f]);
	obj_image=getTitle();
	baseNameEnd=indexOf(obj_image, ".tiff"); 
	baseName=substring(obj_image, 0, baseNameEnd); 

	run("3D Ellipsoid Fitting", " ");
	selectWindow("Results");  //select Log-window 
	saveAs("Results", Dir_res_ellp+baseName+".csv");
	run("Clear Results");
	
	selectWindow(obj_image);   
	run("3D Geometrical Measure");
	selectWindow("Results");  //select Log-window 
	saveAs("Results", Dir_res_geo+baseName+".csv");
	run("Clear Results");
	
	selectWindow(obj_image);  
	run("3D Shape Measure");
	selectWindow("Results");  //select Log-window 
	saveAs("Results", Dir_res_shape+baseName+".csv");
	run("Clear Results");

	print("\\Clear");
	run("Close All");
}



