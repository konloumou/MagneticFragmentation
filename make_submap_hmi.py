
# coding: utf-8

# Iain's code of loading in HMI data, selecting  submap and saving it out as a fits file. Start the pipeline with this one. Run on the cluster as a python code.

# In[ ]:


import glob
from datetime import datetime
import os

import astropy.units as u
import sunpy.map
from sunpy.coordinates import frames
from astropy.coordinates import SkyCoord

from window_corner_rotation import window_corner_rotation

import warnings
warnings.filterwarnings("ignore")


# In[ ]:


act_region = 'SPoCA21717'

datapath = '/scratch/mshtarii/dina/magfrag_data/SPoCA21717_HMI_data/'
hmiimg = sorted(glob.glob(datapath+"/hmi.m*fits"))

main = '/scratch/mshtarii/dina/output_'+ act_region +'/'

subarraysdir = main +'hmi_submap_arrays/'
outname ='hmisub_'

if not os.path.exists(subarraysdir):
        os.makedirs(subarraysdir)

print(len(hmiimg))


# In[ ]:


# box at initial time
bottom_left = [-221,-430] 
top_right = [-170,-387]

data0=sunpy.map.Map(hmiimg[0]).rotate()

initial_obs_time = data0.date

bl = SkyCoord(bottom_left[0] * u.arcsec,
              bottom_left[1] * u.arcsec, frame = frames.Helioprojective, \
              obstime = initial_obs_time)  
tr = SkyCoord(top_right[0] * u.arcsec,
              top_right[1] * u.arcsec, frame = frames.Helioprojective, \
              obstime = initial_obs_time) 
bl_stonyhurst = bl.transform_to(frames.HeliographicStonyhurst)
tr_stonyhurst = tr.transform_to(frames.HeliographicStonyhurst)

# Create a submap using those coordinates
submap = data0.submap(bl, tr)

submap.save(subarraysdir+outname+submap.date.strftime("%Y%m%d_%H%M%S")+'.fits',overwrite=True)


# In[ ]:


for i in range(1,len(hmiimg)):
    print(str(i+1)+' of '+str(len(hmiimg))+' '+hmiimg[i])
    data=sunpy.map.Map(hmiimg[i]).rotate()
    timediff = (data.date - initial_obs_time).total_seconds()/86400
    [bl, tr] = window_corner_rotation(bl_stonyhurst, tr_stonyhurst, timediff, data.date)


    # Create a submap using those coordinates
    submap = data.submap(bl, tr)

    submap.save(subarraysdir+outname+submap.date.strftime("%Y%m%d_%H%M%S")+'.fits',overwrite=True)

