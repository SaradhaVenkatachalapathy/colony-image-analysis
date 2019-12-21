# Spheroid Analysis
These are programs writen to visualize and quantatively characterize spheroids. 

The first program to run is the _"move_files_rawimages.ijm"_. It moves all files in the input folder to a new folder it creates in the input directory called raw images.

For visualization of 3D spheroids, it is important to have 3D projections of the images to see the z profile along with a montage of zprojected images for visualizing the variability among the images aquired. Towards this end 2 programs are used. 
  1. The program _"project_in3d.ijm"_ perfroms a 3D projection of all the images. 
  2. The program _"spheroid_zproject_raw.ijm"_ projects the images along the Z axis and makes a montage for each channel. 
  
In order to then identify spheroids as 3D objects,we use _"3d_objects_spheroid.ijm"_. Here, the sheroids are identified using only the nuclear stain. After this step, a number of further analysis can be done. We have doen the following
  1. The program _"spheroid_zproject_2d_measure.ijm"_ quantifies the geometrical and intensity measures for a z-projected spheroid across channels 
  2. The program _"measures_spheroid.ijm"_ quantifies the 3D geometrical and intensity measures for spheroid across channels 
  3. The program _spheroid_shells.ijm"_ allows us to measure the radial distribution of intensities within a spheroid. It creates nesting shells by eroding the spheroid by the same distance in X, Y and Z thereby enabiling the quantification quantify of intensity distribution in these shells.
  4. The program _"Spheroid_zaxis_profile.ijm"_ allows us to measure the axial distribution of intensities within a spheroid. It measures the intensity along the z axis of the spheroid. 
  5. The program _"spheroid_zproject_2d_measure.ijm"_ quantifies the 2D geometrical and intensity measures for spheroid across channels from a Z projected image.
  6. The program _"spheroid_angles.ijm"_ quantifies the angles in a XZ projected image across defined channels

We can now run all the analysis using the function defined in _"run_spheroid_analysis.ijm"_.
