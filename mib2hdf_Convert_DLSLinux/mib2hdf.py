import argparse
from IdentifyPotentialConversions import check_differences

def convert(beamline, year, visit, mib_files, folder = None):
    """    
    This is to convert a set of time-stamped 4DSTEM folders (mib_files) into a
    similar folder structure in the processing folder of the same visit. The
    following files get saved for each stack:
        - reshaped 4DSTEM HDF5 file
        - The above file binned by 4 in the diffraction plane
        - TIFF file and PNG file of incoherent BF reconstruction
        - TIFF file and PNG file of sparsed sum of the diffraction patterns
    The code gets the floor of square_root of th total number of frames and
    uses that as scan_x / scan_y values to initially try to reshape using the
    fly_back pixel. If that fails, it just reshapes using those values.
    """
    import reshape_4DSTEM_funcs as reshape
    import hyperspy.api as hs
    import os
    import time
    import numpy as np
    
    if folder:
        raw_location = os.path.join('/dls',beamline,'data', year, visit, 'Merlin', os.path.relpath(folder))
    else:
        raw_location = os.path.join('/dls',beamline,'data', year, visit, 'Merlin')  
        
    proc_location = os.path.join('/dls',beamline,'data', year, visit, 'processing', 'Merlin')
    if not os.path.exists(proc_location):
        os.mkdir(proc_location)
    
    os.chdir(raw_location)
    
    # get together the full paths of the datasets
    mib_files_locations = []
    for i, file in enumerate(mib_files):
        if folder:
            mib_files_locations.append(os.path.join(raw_location, os.path.relpath(mib_files[i][0].split('/')[-1]))) 
        else:  
            mib_files_locations.append(os.path.join(raw_location, os.path.relpath(mib_files[i][0]))) 
    
    # main loop 
    for j, dirName in enumerate(mib_files_locations):
        # time0 = time.time()
        print('********************************************************')
        print('Currently active in this directory: %s' % dirName)
        try:
            if os.path.exists(os.path.join(proc_location, os.path.relpath(mib_files[j][0]))):
                # If the folder already exists assumes the hdf5 is already saved and moves on
                print(os.path.join(proc_location, os.path.relpath(mib_files[j][0])) + ' Folder already exists!')
            else:
                os.makedirs(os.path.join(proc_location, os.path.relpath(mib_files[j][0])))    
                for root, dirs, files in os.walk(dirName):
                    for file_name in files:
                        if file_name.endswith('.hdr'):
                            img_flag = 0
                            # checks for correct file type
                            print('Converting the file:  ', file_name)
                            os.chdir(dirName)
                            # load the data lazily into hyperspy 
                            dp = hs.load(file_name, lazy = True)
                            # Calculate scan_X
                            scan_X = int(np.sqrt(dp.axes_manager[0].size))
                            print('Data loaded to hyperspy')
                            # checks to see if it is a multi-frame data before reshaping
                            if any(dp.axes_manager.navigation_axes):
                            # attampt to reshape the data 
                                try:
                                    dp = reshape.reshape_4DSTEM_FlyBack(dp, scan_X)
                                    print('Data reshaped using flyback pixel to: '+ str(dp.axes_manager.navigation_shape))
                                    dp_bin = dp.rebin(scale = (1,1,4,4))
                                except:
                                    print('Data reshape using flyback pixel failed! - Reshaping using scan size instead.')
                                    # if reshape fales bin data dowm 
                                    dp = reshape.reshape_4DSTEM_FrameSize(dp, scan_X, scan_X)
                                    dp_bin = dp.rebin(scale = (1,1,4,4))
                                    img_flag = 1
                                
                                # out of lazy to compute images
                                dp_bin.compute()
                                
                                #calculate incoherent bright-field image
                                ibf = dp_bin.sum(axis = dp_bin.axes_manager.signal_axes)
                                ibf = hs.signals.Signal2D(ibf)
                                
                                # sum dp image of a subset dataset
                                dp_subset = dp.inav[0::int(dp.axes_manager[0].size / 50), 0::int(dp.axes_manager[0].size / 50)]
                                sum_dp_subset = dp_subset.sum()
                                img_flag = 1
                    # save data
                    os.chdir(os.path.join(proc_location, os.path.relpath(mib_files[j][0])))
                               
                    print('Saving hdf5 : ' + file_name.rpartition('.')[0] +'.hdf5')
                    dp.save(file_name, extension = 'hdf5')
                    print('Saved hdf5 : ' + file_name.rpartition('.')[0] +'.hdf5')
                        
                    print('Saving binned data: ' + file_name.rpartition('.')[0] + '_binned.hdf5')
                    dp_bin.save('binned_' + file_name, extension = 'hdf5')
                    print('Saved binned data: binned_' + file_name.rpartition('.')[0] + '.hdf5')
            
                    if img_flag == 1:
                        print('Saving average diffraction pattern')
                        file_dp = file_name.rpartition('.')[0]+ '_subset_dp'
                        sum_dp_subset.save(file_dp, extension = 'tiff')
                        sum_dp_subset.save(file_dp, extension = 'png')
                        print('Saving ibf image')
                        file_ibf =  file_name.rpartition('.')[0]+ '_ibf'
                        ibf.save(file_ibf, extension = 'tiff')
                        ibf.save(file_ibf, extension = 'png')
                            
               # print('Total time elapsed (seconds): ', int(time.time() - time0))
        except:
            print('============================================================')
            print('Something went wrong with ' + file_name)
            print('Skipping to next file')
            print('============================================================')
            continue

def main(beamline, year, visit, folder = None):
    
    [to_convert, mib_files] = check_differences(beamline, year, visit, folder)
    if bool(to_convert):
        convert(beamline, year, visit, mib_files, folder)
    else:
        print('Nothing to convert here!')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('beamline', help='Beamline name')
    parser.add_argument('year', help='Year')
    parser.add_argument('visit', help='Session visit code')
    parser.add_argument('folder', nargs= '?',help='OPTION to add a specific folder within a visit \
                        Merlin folder to look for data, e.g. sample1/dataset1/')
    v_help = "Display all debug log messages"
    parser.add_argument("-v", "--verbose", help=v_help, action="store_true",
                        default=False)

    args = parser.parse_args()
    main(args.beamline, args.year, args.visit, args.folder)
