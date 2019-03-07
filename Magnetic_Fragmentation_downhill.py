
# coding: utf-8

# Import required modules and set up matplotlib display properties. Set up variables.
# 
# KL 25-2-19 Update to save out arrays of submap of HMI map. Update write_props to save out with 4 decimal points. 
# 
# KL 27-2-19 Update to save out arrays of labelled maps.
# 
# KL 5-3-19 Update to load in submaps made by make_submap_hmi.
# 
# KL 7-3-19 Updated and tested that results match previous ones. Slightly changed loop so that arrays are being split in pos and neg after the split_pol_thresh is called. Ready for Github.

# In[ ]:


# sunpy.map handles the coordinates of the data we are using
import sunpy.map
# numpy is a package that makes handling arrays much easier
import numpy as np
# matplotlib.pyplot is the standard plotting tool in python
import matplotlib.pyplot as plt
# skimage.measure performs the fragmentation
import skimage.measure as ms
# The next three imports are functions written just for this purpose
from split_pol_thresh import split_pol_thresh
from find_regions import find_regions
from write_props import write_props
# os is a library of operating system functions that we use to get filepaths
import os
import glob

# The next three lines are for plotting in a Jupyter notebook
# explanation of difference between: %matplotlib notebook and %matplotlib inline
# https://github.com/matplotlib/matplotlib/issues/4879
get_ipython().run_line_magic('matplotlib', 'inline')
plt.rcParams['figure.figsize'] = 11, 11
plt.rcParams.update({'font.size': 12})

fragmentation = 'downhill'


# In[ ]:


# Choose active region working on, data to be read, path to save the bulk properties and location of the submap
# * It does not accept '~/Dropbox/....' as input
# * The bulk_properties.txt contains the total number of fragments in each image, the total area in each image, 
# and the total flux in each image. We will need these to make plots like in Fraser's thesis.
# * The code allows the user to choose the bottom left and top right point of the observation window
# the submap will be created. So I choose around the active region I want to study. Then, inside find_regions.py,
# it will read in the date of each HMI folder and differencially rotate where the bl and tr should be for each
# new map. Choose the bottom left and top right of the initial observation window using solar coordinates 
#(in arcsecs from disk centre)

act_region ='SPoCA21717'


# In[ ]:


br_thresholds = [50, 100] # Gauss 

for i in range(len(br_thresholds)):
    threshold = br_thresholds[i]

    #pc
    main = '/scratch/mshtarii/dina/output_'+ act_region +'/'
    subarraysdir = '/scratch/mshtarii/dina/hmi_submap_arrays/'+ act_region +'_hmi_submap_arrays/'

    #mac
    #datapath = '/Users/dina/phd-mac/large_folders/mag_frag_project/'+ act_region +'_HMI_data/'
    #main = '/Users/dina/Desktop/from_pc/output_'+ act_region +'/'
    #subarraysdir = '//Users/dina/phd-mac/large_folders/mag_frag_project/SPoCA21717_hmi_submap_arrays/'

    propsdir = main +'frag_props_'+ fragmentation +'_'+ str(threshold) +'G_a0_4dp/'    
    bulkdir = main +'bulk_props_'+ str(threshold) +'G/'
    labelarraysdir = main + 'labelled_maps_arrays_' +fragmentation +'_'+ str(threshold) +'G_a0_4dp/'

    if not os.path.exists(propsdir):
        os.makedirs(propsdir)

    if not os.path.exists(bulkdir):
        os.makedirs(bulkdir)

    if not os.path.exists(labelarraysdir):
        os.makedirs(labelarraysdir)    
        
        
    # Load the HMI submaps       
    files_to_load=sorted(glob.glob(subarraysdir +"*"))
    
    # Read in the first submap and keep its date
    filename = files_to_load[0]
    submap_area = sunpy.map.Map(filename)
    initial_obs_time = submap_area.date
    
    # For each submap, start the fragmentation process: find fragments and create text document that represents the properties of fragments in each image file.
    for image in range(0,1):#len(files_to_load)):
        filename =files_to_load[image]
        submap = sunpy.map.Map(filename)

        # Calculate the time difference in the form [days.fraction_of_day] from the original image in days
        # The brackets at total_seconds() indicate that this is a function attached to the variable
        new_obs_time = submap.date
        timediff = (new_obs_time - initial_obs_time).total_seconds()/86400
        
        ###############################CONTINUE ON WEDNESDAY#######################################
        
        # Add a padding around the submap by converting the last pixel of each side of the submap to zero.
        # Set all of the edge pixels in the submaps to zero. If not and there is a label on the edge, when the code
        # will say "look at col+1 (or row+1)" I will be asking it to look out of the image and give back an error.
        submap_area_data = np.pad(submap.data, ((1, 1), (1, 1)), 'constant', constant_values=(0,0))
        
        # Split the submap into two: positive and negative.
        # Strip out all negative data in the positive submap, and vice versa
        pos_submap_data = split_pol_thresh(submap_area_data, threshold, 'pos')
        neg_submap_data = split_pol_thresh(submap_area_data, threshold, 'neg')
        
        ################
        #Optionally: plot out the HMI maps    
#         downsubmapsdir = main +'hmi_submaps_'+ str(threshold) +'G/'
#         if not os.path.exists(ldownsubmapsdir):
#         os.makedirs(downsubmapsdir)    
#         pos_submap_area.plot()
#         imagefilename =  submapsdir + 'HMI_'+str(threshold)+'_G_' + str(image).zfill(4)
#         plt.savefig(imagefilename, dpi = 300)
#         plt.clf()
        ###################

        # Find the regions in the positive and negative magnetogram data
        pos_region_frame, num_pos_regions = find_regions(pos_submap_data)
        neg_region_frame, num_neg_regions = find_regions(neg_submap_data)

        # Convert the data into label frames that scikit-image can use
        pos_labeled_frame, pos_num_labels = ms.label(pos_region_frame.astype(int),
                                                     return_num = True,
                                                     connectivity = 2)
        neg_labeled_frame, neg_num_labels = ms.label(neg_region_frame.astype(int),
                                                     return_num = True,
                                                     connectivity = 2)

        # Save them out to use later
        np.save(labelarraysdir  +'labelled_map_arrays_'+ new_obs_time.strftime('%Y%m%d_%H%M%S').zfill(4) +'_p', pos_labeled_frame)
        np.save(labelarraysdir  +'labelled_map_arrays_'+ new_obs_time.strftime('%Y%m%d_%H%M%S').zfill(4) +'_n', neg_labeled_frame)

        #Use the scikit-image method 'regionprops' to get various properties on
        #the fragments in the image data
        #Cache=False means that no properties will be computed until they are specifically called in write_props
        pos_properties = ms.regionprops(pos_labeled_frame,
                                        intensity_image= submap.data, cache=False)
        neg_properties = ms.regionprops(neg_labeled_frame,
                                        intensity_image= submap.data, cache=False)

        # Write the properties to files. Note that to add new properties,
        # you need to do so in the 'write_props' function
        write_props(pos_properties, 'p', image, new_obs_time, submap, act_region, threshold, main, fragmentation, propsdir, bulkdir)
        write_props(neg_properties, 'n', image, new_obs_time, submap, act_region, threshold, main, fragmentation, propsdir, bulkdir)

