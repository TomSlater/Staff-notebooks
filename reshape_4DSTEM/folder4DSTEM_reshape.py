# -*- coding: utf-8 -*-
"""
Created on Tue Oct 30 13:50:17 2018

@author: eha56862

This is to auto-convert a time-stamped 4DSTEM folder into a similar folder structure
in the processing folder of the same visit. The following files get saved for each stack:
    - reshaped 4DSTEM HDF5 file
    - The above file binned by 4 in the diffraction plane
    - TIFF file and PNG file of incoherent BF reconstruction
    - TIFF file and PNG file of sparsed sum of the diffraction patterns
User can choose via the parameters set below to use the fly-back over-exposed pixel 
to reshape the data or by providing the scan dimensions. If the fly-back reshape fails, 
it reshapes using the scan size. So scan_X and scan_Y values should be provided.

This file can be run by:
   > runfile('.../folder4DSTEM_reshape.py')
or from python:
    import os
    os.chdir(r’the path containing the attached .py files’)
    exec(open('folder4DSTEM_reshape.py').read())
"""

# User parameters:
# =============================================================================
fly_back = True
frame_size = False
# If using frame_size to reshape:
scan_X = 256  # number of lines
scan_Y = 256  # number of probe positions per line
# =============================================================================
import reshape_4DSTEM_funcs as reshape
import hyperspy.api as hs
import os
import time

[acquisition_folder, processing_folder] = reshape.folders_setup()

# For converting a series of mib datasets to hdf5
os.chdir(acquisition_folder)
for dirName in os.listdir(acquisition_folder):
    time0 = time.time()
    print('********************************************************')
    print('Currently active in this directory: %s' % dirName)
    try:
        for root, dirs, files in os.walk(acquisition_folder + '\\' + dirName, topdown = False):
            if os.path.exists(processing_folder + r'\\' + dirName):
                # If the folder already exists assumes the hdf5 is already saved and moves on
                print(dirName + ' Folder already exists!') 
            else:
                for file_name in files:
                    if file_name.endswith('.hdr'):
                        img_flag = 0
                        # checks for correct file type
                        print(file_name)
                        os.chdir(acquisition_folder + '\\' + dirName)
                        # load the data lazily into hyperspy 
                        dp = hs.load(file_name, lazy = True)
                        print('loaded to hyperspy')
                        # checks to see if it is a multi-frame data before reshaping
                        if any(dp.axes_manager.navigation_axes):
                                # attampt to reshape the data 
                                if fly_back:
                                    try:
                                        dp = reshape.reshape_4DSTEM_FlyBack(dp, scan_X)
                                        print('Data reshaped to: '+ str(dp.axes_manager.navigation_shape))
                                    except :
                                        print('Data reshape using flyback failed! - Reshaping using scan size instead.')
                                        # if reshape fales bin data dowm 
                                        dp = reshape.reshape_4DSTEM_FrameSize(dp, scan_X, scan_Y)
                                        dp_bin = dp.rebin(scale = (1,1,4,4))
                                        img_flag = 1
                                elif frame_size:
                                    dp = reshape.reshape_4DSTEM_FrameSize(dp, scan_X, scan_Y)
                                    print('Data reshaped to: '+ str(dp.axes_manager.navigation_shape))
                                else:
                                    print('=====================================================')
                                    print('One of the reshaping algorithms need to be selected!')
                                    print('----------Data not reshaped ---------- ')
                                    print('=====================================================')
                                    
                                dp_bin = dp.rebin(scale = (1,1,4,4))
                                dp_bin.compute()
                                
                                #calculate incoherent bright-field image
                                ibf = dp_bin.sum(axis = dp_bin.axes_manager.signal_axes)
                                ibf = hs.signals.Signal2D(ibf)
                                
                                # sum dp image of a subset dataset
                                dp_subset = dp.inav[0::int(dp.axes_manager[0].size / 50), 0::int(dp.axes_manager[0].size / 50)]
                                sum_dp_subset = dp_subset.sum()
                                
                                img_flag = 1
                # save data
                os.chdir(processing_folder)
                       
                os.mkdir(dirName)
                os.chdir(processing_folder + r'\\' + dirName)
                print('saving hdf5 : ' + file_name +'.hdf5')
                dp.save(file_name, extension = 'hdf5')
                print('saved hdf5 : ' + file_name +'.hdf5')
                
                print('saving binned data: ' + file_name + '_binned.hdf5')
                dp_bin.save('binned_' + file_name, extension = 'hdf5')
                print('saved binned data: binned_' + file_name + '.hdf5')

                if img_flag == 1:
                    print('saving average diffraction pattern')
                    file_dp = file_name.rpartition('.')[0]+ '_subset_dp'
                    sum_dp_subset.save(file_dp, extension = 'tiff')
                    sum_dp_subset.save(file_dp, extension = 'png')
                    print('saving ibf image')
                    file_ibf =  file_name.rpartition('.')[0]+ '_ibf'
                    ibf.save(file_ibf, extension = 'tiff')
                    ibf.save(file_ibf, extension = 'png')
                    
                print('Total time elapsed (seconds): ', int(time.time() - time0))
    except :
        print('============================================================')
        print('Something went wrong with ' + file_name)
        print('Skipping to next file')
        print('============================================================')
        continue

