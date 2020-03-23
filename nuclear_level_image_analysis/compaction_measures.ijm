/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////  SCRIPT TO COMPUTE FEATURES THAT DESCRIBE THE COMPACTNESS OF THE DNA WITHIN THE NUCLEUS
///////  WRITTEN BY: SARADHA VENKATACHALAPATHY                                                                                                                   
///////  ASSUMPTIONS: The input image is a confocal zstack and the channel 1 contains the nucleus.
///////  DESCRIPTION: The script opens each image, thresholds the nucleus gets the nuclear feature. Then it performs a second threshold to obtain the Heterochromatin nodes (lower threshold= Mean +1.5* Stdev).
//////				  Then it computes the heterochromatin features and then saves the results in the folder "2D_measures_nuclei"
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

dirsa = getArgument(); //get the directory

//set paths
dir_dat= dirsa + "2D_measures_nuclei"+ File.separator;
setBatchMode(true);
dir1= dirsa + "indivisual_nuclei_DNA"+ File.separator;
list1 = getFileList(dir1);
//create variables
n=list1.length;
HCcontent=newArray(n);
ECcontent=newArray(n);
HCvolume=newArray(n);
ECvolume=newArray(n);
volume=newArray(n);
label=newArray(n);
//set options for measurements
run("Set Measurements...", "area mean standard modal min centroid center perimeter bounding fit shape feret's integrated median skewness kurtosis area_fraction stack limit display redirect=None decimal=3");
for (i=0; i<list1.length; i++) { //loop through all images

		//open image, convert to 8-bit and get he voxel characteristics
		path = dir1+list1[i];
		open(path);
		run("8-bit");
		getVoxelSize(width, height, depth, unit);

		//set threshold and analyse particles
		run("Clear Results");
		//run("Threshold...");
		setThreshold(1, 255);
		run("Analyze Particles...", "size=0-Infinity circularity=0.00-1.00 show=Nothing display stack");
		//obtain the variables of the nucleus
		num=nResults;
		volume[i]=0;
		nuc_area=newArray(num);
		low_thresh=newArray(num);
		slice_num=newArray(num);
		hcarea=newArray(num);
		tot_intensity=newArray(num);
		hc_intensity=newArray(num);
		
		for(k=0; k<num; k++){
			nuc_area[k]=getResult("Area",k);
			slice_num[k]=getResult("Slice",k);
			low_thresh[k]=getResult("Mean",k) +1.5*getResult("StdDev",k);
			volume[i]=volume[i]+getResult("Area",k)*depth;
			tot_intensity[k]=getResult("RawIntDen",k);
			
		}
		run("Clear Results");
		//threshold for HC node
		for(k=0; k<num; k++){
			setSlice(slice_num[k]);
			//run("Threshold...");
			setThreshold(low_thresh[k], 255);
			run("Measure");
		}
		for(k=0; k<num; k++){
			hcarea[k]=getResult("Area",k);
			hc_intensity[k]=getResult("RawIntDen",k);
		}

		//measure Heterochromatin features
		label[i]=File.getName(path);
		HCcontent[i]=0;
		ECcontent[i]=0;
		HCvolume[i]=0;
		ECvolume[i]=0;
		for (j=0; j<nResults; j++){
			HCvolume[i]=HCvolume[i]+hcarea[j];
			ECvolume[i]=ECvolume[i]+(nuc_area[j]-hcarea[j]);
			HCcontent[i]=HCcontent[i]+hc_intensity[j];
			ECcontent[i]=ECcontent[i]+(tot_intensity[j]-hc_intensity[j]);
			
		}
		close();
		HCvolume[i]=HCvolume[i]*depth;
		ECvolume[i]=ECvolume[i]*depth;
		run("Clear Results");
}
run("Clear Results");
//create the new results table
for (k=0; k<list1.length; k++) {
	setResult("Label1",k ,label[k]);
	setResult("HCcontent",k, HCcontent[k]);
	setResult("ECcontent",k, ECcontent[k]);
	setResult("HCvolume",k, HCvolume[k]);
	setResult("ECvolume",k, ECvolume[k]);
	setResult("HC_EC_content",k, HCcontent[k]/ECcontent[k]);
	setResult("HC_EC_volume",k, HCvolume[k]/ECvolume[k]);
	setResult("Volume",k,volume[k]);
}

saveAs("Results",dir_dat+"compaction.csv" ); 
run("Clear Results");