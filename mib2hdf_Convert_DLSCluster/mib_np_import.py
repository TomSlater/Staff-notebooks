#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 15:02:01 2019
This script imports Medipix mib files into numpy and passes as a lazy
signal to hyperspy.
It works on Raw and not-RAW files and can accommodate single chip and
quad choip data.
It only needs the mib file for import, i.e. ignores the .HDR file.

@author: eha56862
"""
import os
import numpy as np
import h5py
import dask.array as da
import dask
import hyperspy.api as hs

def _manageHeader(fname):
    ''' Getting all the necessary information from the header  of the mib file '''
    
    Header = str()
    with open(fname,'rb') as input:
        aByte= input.read(1)
        Header += str(aByte.decode('ascii'))
   
        # This gets rid of the header 
        while aByte and ord(aByte) != 0: 
       
            aByte= input.read(1)
            Header += str(aByte.decode('ascii'))
    
    ####################################################
    #print(' Header str ',header_str)
    #header = ImageHeader(header_str)
    elements_in_header = Header.split(',')
    
    #print(elements_in_header)
    
    #Nframes = int(elements_in_header[1]) always returns 00001
    #print(elements_in_header[2])
    
    DataOffset = int(elements_in_header[2])
    
    NChips = int(elements_in_header[3])
    #if elements_in_header[4] == '1024':
        #self.PixelXdim = int(elements_in_header[4])//2
        #self.PixelYdim = int(elements_in_header[5])*2
    #else:
    #self.PixelXdim = int(elements_in_header[4])
    #self.PixelYdim = int(elements_in_header[5])
    #print(self.PixelXdim)
    PixelDepthInFile= elements_in_header[6]
    sensorLayout = elements_in_header[7].strip()    
    Timestamp = elements_in_header[9]
    shuttertime = float(elements_in_header[10])
    #GainMode= elements_in_header[13]
    
    if PixelDepthInFile == 'R64':
        bitdepth =int(elements_in_header[18]) # RAW
    elif PixelDepthInFile =='U16':
        bitdepth =12
    elif PixelDepthInFile =='U08':
        bitdepth =6

        #example output for 6bit 256*256 4DSTEM data:
        #(768, 4, 'R64', '2x2', '2019-06-14 11:46:12.607836', 0.0002, 6)
        #example output for 12bit single frame nor RAW:
        #(768, 4, 'U16', '2x2', '2019-06-06 11:12:42.001309', 0.001, 12)
        
    hdr = (DataOffset,NChips,PixelDepthInFile,sensorLayout,Timestamp,shuttertime,bitdepth)
    return hdr

def parse_hdr(fp):
    """Parse information from mib file header info from _manageHeader function.
    Accepts file object 'fp' a mib file. Returns dictionary hdr_info.
    """
    hdr_info = {}
    
    read_hdr = _manageHeader(fp)

    #set the array size of the chip 

    if read_hdr[3] == '1x1':
        hdr_info['width'] = 256
        hdr_info['height'] = 256
    elif read_hdr[3] == '2x2':
        hdr_info['width'] = 512
        hdr_info['height'] = 512
    
    hdr_info['Assembly Size'] = read_hdr[3]
    #convert frames to depth
    #hdr_info['depth'] = int(hdr_info['Frames in Acquisition (Number)'])
    #set mib offset
    hdr_info['offset'] = read_hdr[0]
    #set data-type
    hdr_info['data-type'] = 'unsigned'
    #set data-length
    if read_hdr[6] == '1':
		#binary data recorded as 8 bit numbers
        hdr_info['data-length'] = '8'
    else:
		#changes 6 to 8 , 12 to 16 and 24 to 32 bit
        cd_int = int(read_hdr[6])
        hdr_info['data-length'] = str(int((cd_int + cd_int/3) ))
        
    hdr_info['Counter Depth (number)'] = int(read_hdr[6])
    if read_hdr[2] =='R64':
        hdr_info['raw'] = 'R64'
    else:
        hdr_info['raw'] = 'MIB'
    #set byte order
    hdr_info['byte-order'] = 'dont-care'
    #set record by to stack of images
    hdr_info['record-by'] = 'image'
    
    
    #set title to file name
    hdr_info['title'] = fp.split('.')[0]
    #set time and date
    #Adding the try argument to accommodate the new hdr formatting as of April 2018
    try:
        year, month, day_time = read_hdr[4].split('-')
        day , time = day_time.split(' ')
        hdr_info['date'] = year + month + day
        hdr_info['time'] = time
    except:
        day, month, year_time = read_hdr[4].split('/')
        year , time = year_time.split(' ')
        hdr_info['date'] = year + month + day
        hdr_info['time'] = time
        
    hdr_info['data offset'] = read_hdr[0]
    #hdr_info['raw'] = _manageHeader(fp.name[:-3]+'mib')[:2]
		
    # print(hdr_info)
    return hdr_info



def read_mib(hdr_info, fp, mmap_mode='r', save_hdf = False, path = None):
    """Read the raw file object 'fp' based on the information given in the
    'hdr_info' dictionary.

    Parameters
    ----------
    hdr_info: dict
        A dictionary containing the keywords as parsed by read_hdr
    fp:
    mmap_mode: {None, 'r+', 'r', 'w+', 'c'}, optional
    If not None, then memory-map the file, using the given mode
    (see `numpy.memmap`).  The mode has no effect for pickled or
    zipped files.
    A memory-mapped array is stored on disk, and not directly loaded
    into memory.  However, it can be accessed and sliced like any
    ndarray.  Memory mapping is especially useful for accessing
    small fragments of large files without reading the entire file
    into memory.


    """
    #raw mib files sizes in bytes for single frames 
    
    reader_offset = 0
    if hdr_info['Assembly Size'] == '2x2':
        mib_file_size_dict = {
        '1': 33536,
        '6': 262912,
        '12': 525056,
        }
    if hdr_info['Assembly Size'] == '1x1':
        mib_file_size_dict = {
        '1': 8576,
        '6': 65920,
        '12': 131456,
        }
    
    width = hdr_info['width']
    height = hdr_info['height']
    #depth = hdr_info['depth']
    offset = hdr_info['offset']
    data_length = hdr_info['data-length']
    data_type = hdr_info['data-type']
    endian = hdr_info['byte-order']
    record_by = hdr_info['record-by']
    
    file_size = os.path.getsize(fp[:-3]+'mib')
    if hdr_info['raw'] == 'R64':
        #file_size = os.path.getsize(fp[:-3]+'mib')
        #print(file_size)
        single_frame = mib_file_size_dict.get(str(hdr_info['Counter Depth (number)']))
        depth = int(file_size / single_frame)
    elif hdr_info['raw'] == 'MIB':
        if hdr_info['Counter Depth (number)'] =='1':
            
            single_frame = mib_file_size_dict.get('6')  # 1 bit and 6 bit non-raw frames have the same size
            depth = int(file_size / single_frame)
        else:
            single_frame = mib_file_size_dict.get(str(hdr_info['Counter Depth (number)']))
            depth = int(file_size / single_frame)
            #print(depth)
            

    if data_type == 'signed':
        data_type = 'int'
    elif data_type == 'unsigned':
        data_type = 'uint'
    elif data_type == 'float':
        pass
    else:
        raise TypeError('Unknown "data-type" string.')

    #if endian == 'big-endian':
    #    endian = '>'
    #elif endian == 'little-endian':
    #    endian = '<'
    #else:
    #    endian = '='
    # mib data always big-endian 
    endian = '>'
    #print(data_length)
    data_type += str(int(data_length))
    data_type = np.dtype(data_type)
    data_type = data_type.newbyteorder(endian)
    
    #set header number of bits
#    if hdr_info['Assembly Size'] == '1x1':
#        if data_length == '16':
#            hdr_multiplier = (int(data_length)/8)**-1
#            hdr_bits = int(hdr_info['data offset'] * hdr_multiplier)
#        else:
#            hdr_multiplier = (int(data_length)/8)**-1
#            hdr_bits = int(hdr_info['data offset'] * hdr_multiplier)
#    if hdr_info['Assembly Size'] == '2x2':
#        hdr_multiplier = (int(data_length)/8)**-1
#        hdr_bits = int(hdr_info['data offset'] * hdr_multiplier)
    hdr_multiplier = (int(data_length)/8)**-1
    hdr_bits = int(hdr_info['data offset'] * hdr_multiplier)
    #print(hdr_bits)
	
    #print(hdr_bits, data_length, hdr_multiplier)
	
    data = np.memmap(fp,
                     offset=reader_offset,
                     dtype=data_type,
                     mode=mmap_mode)
    #print(data.shape)

    if record_by == 'vector':   # spectral image
        size = (height, width, depth)
        try:
            data = data.reshape(size)
        # in case of incomplete frame:    
        except ValueError:
            if hdr_info['raw'] == 'R64':
                file_size = os.getsize(fp[:-3]+'mib')
                #print(file_size)
                single_frame = mib_file_size_dict.get(str(hdr_info['Counter Depth (number)']))
                correct_size = int(file_size / single_frame)
                #print(correct_size)
                data = data.reshape(correct_size)
            
        
    elif record_by == 'image':  # stack of images
        width_height = width * height
        #print('width_height :',width_height)
        #a_width, a_height = round(np.sqrt(depth)), depth/ (round(np.sqrt(depth)))
        size = (depth, height, width)
        #print('size: ', size)
        #remove headers at the beginning of each frame and reshape

        if hdr_info['raw'] == 'R64':
            file_size = os.path.getsize(fp[:-3]+'mib')
            #print(file_size)
            single_frame = mib_file_size_dict.get(str(hdr_info['Counter Depth (number)']))
            correct_size = int(file_size / single_frame)
            data = data.reshape(-1, width_height + hdr_bits)[:,-width_height:].reshape(correct_size, width, height)
                
        
        
        if hdr_info['raw'] == 'R64':
            if hdr_info['Counter Depth (number)'] == 24 or  hdr_info['Counter Depth (number)'] == 12:
                COLS = 4

            if hdr_info['Counter Depth (number)'] == 1:
                COLS = 64

            if hdr_info['Counter Depth (number)'] == 6:
                COLS = 8
                
            #print(hdr_info['depth'])
            # print(correct_size)
            
            data = data.reshape(depth,height * (height//COLS) , COLS )
                
        
            data = np.flip(data,2)
            
            if hdr_info['Assembly Size'] == '2x2':
                
                try:
                    data = data.reshape(depth,512 // 2, 512 * 2 )
                except ValueError:
                    data = data.reshape(correct_size,512 // 2, 512 * 2 )
                
                det1 = data[:,:,0:256]
                det2 = data[:,:,256:512]
                det3 = data[:,:,512:512+256]
                det4 = data[:,:,512+256:]
                
                det3 = np.flip(det3,2)
                det3 = np.flip(det3,1)
                
                det4 = np.flip(det4,2)
                det4 = np.flip(det4,1)
                
                data = np.concatenate((np.concatenate((det1,det3),1),np.concatenate((det2,det4),1)),2)
                
        if hdr_info['Assembly Size'] == '1x1':
            #print(depth)
            correct_size = int(file_size / single_frame)
            #print(correct_size, data.shape)
            data = data.reshape(-1, width_height + hdr_bits)[:,-width_height:].reshape(correct_size, width, height)
            data = data.reshape(depth,256, 256 )
        
    elif record_by == 'dont-care':  # stack of images
        size = (height, width)
        data = data.reshape(size)
        
    da_data = da.from_array(data)
        
    if save_hdf:
        os.chdir(path)
        f = h5py.File('raw_data', 'w')
        f.create_dataset('dataset1', data = da_data)
        f.close()
      
    #print(da_data.shape)    
    data_hs = hs.signals.Signal2D(da_data).as_lazy()
    return data_hs

def mib_np_reader(mib_filename):
    hdr_stuff = parse_hdr(mib_filename)
    data = read_mib(hdr_stuff, mib_filename)
    
    return data
##Testing
#import numpy as np
#import os
#import hyperspy.api as hs
#os.chdir(r'Y:\2019\mg22317-14\Merlin\raw\20190614 114607')
##os.chdir(r'/dls/e02/data/2019/mg22317-14/Merlin/raw/20190614 125820')
#fp = 'ceria_CL230um_CL11cm_256_2-6bit_800us_01.mib'
##fp = 'ceria_CL230um_CL11cm_256_2-6bit_600us_56576_03.mib'
#test = parse_hdr(fp)
#d = read_mib(test, fp)
#s = hs.signals.Signal2D(d).as_lazy()
