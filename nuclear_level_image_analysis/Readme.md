# Reprogramming Spheroid Nuclear level Analysis
These are programs written to visualize and qualitatively characterize nuclei of cells within a spheroid. 

The first program to run is the _"move_files_rawimages.ijm"_. It moves all files in the input folder to a new folder it creates in the input directory called raw images. 

In order to then identify nuclei as 3D objects, we use _"binarize_the_nucleus.ijm"_ to obtain binary images and _"3dcrop_nuclei.ijm"_ to segment and to obtain individual nuclear crops.

Below is a depiction of a typical nuclear segmentation of a 3D spheroid of cells.<br/>
<p align="center">
<img src='/nuclear_level_analysis/Sepgmentation.png' height='200' width='150'><br/>
</p>


After this step, a number of further analysis can be done. We have done the following<br/>
<p align="center">
<img src='/nuclear_level_analysis/Features.png' height='400' width='800'><br/>
</p>

Below are the programs that achieve the above. 

  1. The program _"projected_2dmeaure_nucleus.ijm"_ quantifies the geometrical and intensity measures for a z-projected nuclei across channels 
  2. The program _"measure_nucleus.ijm"_ quantifies the 3D geometrical and intensity measures for nuclei across channels 
  3. The program _"cell_ring.ijm"_ allows us to obtain cytoplasmic proteins as defined by a 2 µm ring from the nuclear edge.
  4. The program _"compaction_measures.ijm"_ allows us to measure chromatin compaction features.
  5. The program _"elliptic_fourier_descriptors.ijm"_ quantifies the elliptical fourier descriptors of 2D projected nuclear image. 
  6. The program _"glcm_steps.ijm"_ quantifies the GLCM matrix descriptors of 2D projected nuclear image. 

  The program _"zproject_and_padd_montage.ijm"_ projects the images along the Z axis and makes a montage for each channel. 

We can now run all the analysis using the function defined in _"run_nuclear_analysis.ijm"_.
