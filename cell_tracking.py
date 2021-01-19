# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 19:52:30 2017

@author: sarahkeegan and keriabermudez

This is a modification of the cell_traking.py written by Sarah

"""
import os
#os.chdir('/Users/keriabermudez/Dropbox/GitHub/cell-segmentation-analysis/')

import cell_segmentation as cellseg
from skimage import io as sio
import numpy as np
from skimage import img_as_ubyte, exposure
import pandas as pd
import os
#tracks cells from input file given initial tracking windows.

all_measurements = pd.DataFrame()
group = "WT" #"KO"
dir_ = "./PCNA_CYD1_Videos/"+group+ "/"
intensity_to_meas = 'cd1' #'pcna'
#tw_df = pd.read_table(dir_ + "tracking_windows_jack.csv",sep=',')

xlsx = pd.ExcelFile(dir_ + group+"_table.xlsx")
tw_df = pd.read_excel(xlsx, 'Blad1')


sample_nums = tw_df['Sample_Number'].unique()
results_folder = "results_"+intensity_to_meas+"_Blad1"

#%%
rescale_cyd1 = False

for sample_num in sample_nums:
    #for each sample, load movie files
    #then track each cell that has been identified by a tracking window
    cur_tws = tw_df.query("Sample_Number == "+str(sample_num))
    
    image_stack_pcna = sio.imread("%ssample_%.2d_PCNA.tif" % (dir_,sample_num))
    image_stack_cyd1 = sio.imread("%ssample_%.2d_CYD1.tif" % (dir_,sample_num))
    
    #color stack, in 8-bit                    
    image_stack_pcna_8bit = np.zeros(image_stack_pcna.shape, dtype='uint8')
    for i in range(image_stack_pcna.shape[0]):
        img1 = exposure.rescale_intensity(image_stack_pcna[i])
        img1 = img_as_ubyte(img1)
        image_stack_pcna_8bit[i] = img1
    color_stack_pcna_8bit = np.stack((image_stack_pcna_8bit,
                        np.zeros(shape=image_stack_pcna_8bit.shape, dtype=image_stack_pcna_8bit.dtype),
                        np.zeros(shape=image_stack_pcna_8bit.shape, dtype=image_stack_pcna_8bit.dtype)),axis=-1)
    
    #color stact cd1
    image_stack_cyd1_8bit = np.zeros(image_stack_cyd1.shape, dtype='uint8')
    for i in range(image_stack_cyd1.shape[0]):
        img1 = exposure.rescale_intensity(image_stack_cyd1[i])
        img1 = img_as_ubyte(img1)
        image_stack_cyd1_8bit[i] = img1
    color_stack_cyd1_8bit = np.stack((np.zeros(shape=image_stack_cyd1_8bit.shape, dtype=image_stack_cyd1_8bit.dtype),
                        image_stack_cyd1_8bit,
                        np.zeros(shape=image_stack_cyd1_8bit.shape, dtype=image_stack_cyd1_8bit.dtype)),axis=-1)                   
    if(rescale_cyd1):
        image_stack_cyd1_rescl = np.zeros(image_stack_cyd1.shape, dtype='uint16')
        for i in range(image_stack_cyd1.shape[0]):
            image_stack_cyd1_rescl[i] = exposure.rescale_intensity(image_stack_cyd1[i])
        intensity_stack = image_stack_cyd1_rescl
    else:
        intensity_stack = image_stack_cyd1

    print ("Tracking for sample " + str(sample_num))
    for index, row in cur_tws.iterrows():
        
        track_window = (int(row['X']),int(row['Y']),int(row['Width']),int(row['Height']))
        start_s = int(row['StartS'])-1
        end_s = int(row['EndS'])-1
        
        if start_s == 0:
            continue
        z = int(row['Start_Frame']-1)
        label = int(row['Label'])
        print ("label " + str(label) + "...")
        zstack = color_stack_pcna_8bit[:,:,:,0].copy()

        if intensity_to_meas == 'cd1':
            ct2 = cellseg.cell_tracking(intensity_stack, image_stack_pcna, color_stack_cyd1_8bit) #intensity, tracking, display
        else:
            ct2 = cellseg.cell_tracking(image_stack_pcna, image_stack_pcna, color_stack_pcna_8bit) #intensity, tracking, display

        ct2.set_watershed_param(smooth_distance = True, kernel = 10, min_distance=15) # kernel = 10, min_distance=10 or 15 !
        (stack_graph,measurements) = ct2.track_window_graph(label, z, track_window,
                                    use_camshift=False, enhance_bool = True, blur_bool = True, kernel_size = 11, size = 2.5,vline_x =start_s )
        #Save video
        sio.imsave("%s%s/tracked_window_sample_%.2d_%s.tif" % (dir_,results_folder,sample_num,label),stack_graph)

        #save measurements
        measurements['Sample_Number'] = sample_num
        measurements['Label'] = label

        max_inten = measurements.mean_intensity.max()
        min_inten = measurements.mean_intensity.min()
        measurements['norm_mean_inten'] =  (measurements['mean_intensity']- min_inten)/(max_inten- min_inten)
        measurements['offset_z'] = measurements.z - start_s
        all_measurements = all_measurements.append(measurements)

#output results
all_measurements.index = range(len(all_measurements.index))
if intensity_to_meas == 'cd1':   
    all_measurements.rename(columns={'mean_intensity': 'cyd1_mean_intensity',}, inplace=True)
else:   
    all_measurements.rename(columns={'mean_intensity': 'pcna_mean_intensity',}, inplace=True)
all_measurements.to_csv(dir_+results_folder+"/tracked_cells.csv")

