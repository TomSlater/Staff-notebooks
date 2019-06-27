import argparse
from IdentifyPotentialConversions import check_differences
import gc
from mib_np_import import *

def max_contrast8(d):
    data = d.data
    data = data - data.min()
    data = data * (255 / data.max())
    d.data = data
    return d


def convert(beamline, year, visit, mib_to_convert, STEM_flag, scan_X, folder):
    """    
    This is to convert a set of time-stamped 4DSTEM folders (mib_files) into a
    similar folder structure in the processing folder of the same visit. The
    following files get saved for each stack:
        - reshaped 4DSTEM HDF5 file
        - The above file binned by 4 in the diffraction plane
        - TIFF file and JPG file of incoherent BF reconstruction
        - TIFF file and JPG file of sparsed sum of the diffraction patterns
    The code gets scan_X value to initially try to reshape using the
    fly_back pixel. If that fails, it just reshapes using those values.
    If STEM_flag is passed as 0 , i.e. TEM data, it only saves a sum image
    and the HDF5
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
    
    # getting the raw data folders as a set
    data_folders = []
    for path in mib_to_convert:
        data_folders.append(os.path.join(*path.split('/')[:-1]))
    data_folders_set = set(data_folders)
    
    # main loop 
    for mib_path in list(data_folders_set):
        #time0 = time.time()
        print('********************************************************')
        print('Currently active in this directory: %s' % mib_path.split('/')[-1])
        
        # get number of mib files in this folder
        mib_num = 0
        mib_list = []
        for file in os.listdir('/'+ mib_path):
            if file.endswith('mib'):
                mib_num += 1
                mib_list.append(file)
        #print(mib_num)
        #print(mib_list)
        #print('STEM_flag: ', STEM_flag)
        
        if mib_num == 1:
            try:
                #print(mib_path)
                #print(mib_list[0])
                dp = mib_np_reader('/' +mib_path + '/'+ mib_list[0])
                #print('loaded the mib file')
                
            except ValueError:
                print('file could not be read into an array!')
            if (STEM_flag == 0 or STEM_flag == '0'): # if it is TEM data
                process_path = proc_location +'/'+ os.path.join(*mib_path.split('/')[6:])
                if not os.path.exists(process_path):
                    os.makedirs(process_path)
                dp.save(process_path + '/' +mib_list[0], extension = 'hdf5')
                dp_sum = max_contrast8(dp.sum())
                dp_sum.save(process_path + '/' +mib_list[0]+'_sum', extension = 'jpg')
            else:
                process_path = proc_location +'/'+ os.path.join(*mib_path.split('/')[6:])
                if not os.path.exists(process_path):
                    os.makedirs(process_path)
                img_flag = 0
                
                print('Data loaded to hyperspy')
                # checks to see if it is a multi-frame data before reshaping
                if dp.axes_manager[0].size > 1:
                # attampt to reshape the data 
                    try:
                        dp = reshape.reshape_4DSTEM_FlyBack(dp, scan_X)
                        print('Data reshaped using flyback pixel to: '+ str(dp.axes_manager.navigation_shape))
                        dp_bin = dp.rebin(scale = (1,1,4,4))
                    except:
                        print('Data reshape using flyback pixel failed! - Reshaping using scan size instead.')
                        num_frames = dp.axes_manager[0].size
                        dp = reshape.reshape_4DSTEM_FrameSize(dp, scan_X, int(num_frames / scan_X))
                        dp_bin = dp.rebin(scale = (1,1,4,4))
                        img_flag = 1
                            
                    # out of lazy to compute images
                    dp_bin.compute(progressbar = False)                                
                    #calculate incoherent bright-field image
                    ibf = dp_bin.sum(axis = dp_bin.axes_manager.signal_axes)
                    ibf = hs.signals.Signal2D(ibf)
                    ibf = max_contrast8(ibf)
                 
                    # sum dp image of a subset dataset
                    dp_subset = dp.inav[0::int(dp.axes_manager[0].size / 50), 0::int(dp.axes_manager[0].size / 50)]
                    sum_dp_subset = dp_subset.sum()
                    sum_dp_subset = max_contrast8(sum_dp_subset)
                    img_flag = 1
                    
                    # save data
                
                    if img_flag == 1:
                        print('Saving average diffraction pattern')
                        file_dp = mib_list[0].rpartition('.')[0]+ '_subset_dp'
                        sum_dp_subset.save(process_path+'/'+file_dp, extension = 'tiff')
                        sum_dp_subset.save(process_path+'/'+file_dp, extension = 'jpg')
                        print('Saving ibf image')
                        file_ibf =  mib_list[0].rpartition('.')[0]+ '_ibf'
                        ibf.save(process_path+'/'+file_ibf, extension = 'tiff')
                        ibf.save(process_path+'/'+file_ibf, extension = 'jpg')
                                  
                print('Saving hdf5 : ' + mib_list[0].rpartition('.')[0] +'.hdf5')
                dp.save(process_path+'/'+mib_list[0], extension = 'hdf5')
                print('Saved hdf5 : ' + mib_list[0].rpartition('.')[0] +'.hdf5')
                    
                if dp.axes_manager[0].size > 1:
                    print('Saving binned data: ' + mib_list[0].rpartition('.')[0] + '_binned.hdf5')
                    dp_bin.save(process_path+ '/'+'binned_' + mib_list[0], extension = 'hdf5')
                    print('Saved binned data: binned_' + mib_list[0].rpartition('.')[0] + '.hdf5')
                    del dp_bin
                
                del dp
                gc.collect()
                
            
        elif mib_num > 1:
            if (STEM_flag == 0 or STEM_flag == '0'):
                process_path = proc_location +'/'+ os.path.join(*mib_path.split('/')[6:])
                if not os.path.exists(process_path):
                    os.makedirs(process_path)
                for k, file in enumerate(mib_list):
                    print(mib_path)
                    print(file)
                    dp = mib_np_reader('/' +mib_path + '/'+ file)
                    
                    dp.save(process_path + '/' +file, extension = 'hdf5')
                    dp_sum = max_contrast8(dp.sum())
                    dp_sum.save(process_path + '/' +file+'_sum', extension = 'jpg')
                    

def watch_convert(beamline, year, visit, STEM_flag, scan_X, folder):
    
    [to_convert, mib_files, mib_to_convert] = check_differences(beamline, year, visit, folder)
    #holder for raw data path
    import os, time
    if folder:
        raw_location = os.path.join('/dls',beamline,'data', year, visit, 'Merlin', os.path.relpath(folder))
    else:
        raw_location = os.path.join('/dls',beamline,'data', year, visit, 'Merlin')
    if bool(to_convert):
        convert(beamline, year, visit, mib_to_convert, STEM_flag, scan_X, folder)
    else:
        #watch_check = input('Do you want to keep watching this folder? (Y/N)')
        watch_check = 'Y'
        if (watch_check == 'Y' or watch_check == 'y'):
            print(raw_location)
            path_to_watch = raw_location
            [to_convert, mib_files, mib_to_convert] = check_differences(beamline, year, visit, folder)
            before = dict ([(f, None) for f in os.listdir (path_to_watch)])
            while True:
                time.sleep (60)
                after = dict ([(f, None) for f in os.listdir (path_to_watch)])
                added = [f for f in after if not f in before]
                if added: 
                    print("Added dataset: ", ", ".join (added))
                    # print(added[-1])
                    new_data_folder = os.listdir (path = os.path.join(path_to_watch, added[-1]))
                    for f in new_data_folder:
                        if f.endswith('mib'):
                            wait_flag = 1
                            while wait_flag == 1:
                                try:
                                    print('file name: ', f)
                                    # f_size = os.stat(os.path.join(path_to_watch, added[-1], f)).st_size
                                    # f_size = os.path.getsize(os.path.join(path_to_watch, added[-1], f))
                                    # above give the file size from source
                                    # but this below throws an error while copy is not complete:
                                    f_size = os.path.getsize(f)
                                    # print(f_size)
                                    print('file size: ', f_size)
                                    wait_flag = 0
                                except FileNotFoundError:
                                    time.sleep(20)
                                    print('waiting for mib data to copy!!')
                                    print(os.path.isfile(os.path.join(path_to_watch, added[-1], f)))
                                    if os.path.isfile(os.path.join(path_to_watch, added[-1], f)):
                                        wait_flag = 0
                                    else:
                                        pass
           
                    [to_convert, mib_files, mib_to_convert] = check_differences(beamline, year, visit, folder)
                    convert(beamline, year, visit, mib_files, STEM_flag, scan_X, folder)
                    # print('here!')
                before = after

def main(beamline, year, visit, STEM_flag = 1, scan_X = 256, folder = None):
    #print('Active in main!')
    #print(STEM_flag)
    
    watch_convert(beamline, year, visit, STEM_flag , scan_X, folder)
            

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('beamline', help='Beamline name')
    parser.add_argument('year', help='Year')
    parser.add_argument('visit', help='Session visit code')
    parser.add_argument('STEM_flag', nargs= '?',default=1, help='OPTION to specify non-STEM Medipix data \
                        - default is STEM')
    parser.add_argument('scan_X', nargs= '?',default=256,help='OPTION to specify number of probe positions in \
                        x direction- default is 256')
    parser.add_argument('folder', nargs= '?',default=None, help='OPTION to add a specific folder within a visit \
                        Merlin folder to look for data, e.g. sample1/dataset1/')
    v_help = "Display all debug log messages"
    parser.add_argument("-v", "--verbose", help=v_help, action="store_true",
                        default=False)

    args = parser.parse_args()
    
    main(args.beamline, args.year, args.visit, args.STEM_flag, args.scan_X, args.folder)
    #watch_convert(args.beamline, args.year, args.visit, args.folder, args.STEM_flag, args.scan_X)
