import numpy as np
from IPython.core.debugger import Tracer; debug_here = Tracer()
import shutil
import os

def change_param_file(param_file,change_params):
    """Change the param file of an Xnmr file directory.
    A helper function for save_xnmr_bin
    """
    f=open(param_file,'r')
    lines = f.readlines()
    f.close()

    #There are actually two npts in the param file
    #that need to be changed.. this special case is taken care of here
    if 'npts' in change_params:
        change_params['acq_npts'] = change_params['npts']

    try:
        for j in range(len(lines)):
            x = lines[j].find('=')
            if x > 0:
                param = lines[j][0:x-1]
                if param in change_params:
                    newval = change_params[param]
                    lines[j] = param+' = '+newval+'\n'
    except:
        debug_here()
        raise Exception('Error in change_param_file function')

    fr=open(param_file,'w')
    fr.writelines(lines)
    fr.close()
    return



def save_xnmr_bin_fid(fname,data):
    """Save an FID to Xnmr binary format.
    This is designed to save to an existing Xnmr file directory, not create a new one.
    Use when you've done some processing of Xnmr data in python and you want to
    save it to a new Xnmr file.

    It will create a backup directory in the Xnmr file directory and copy over everything
    there prior to overwriting.

    Note that this overwrites 'data' and 'data.fid' with the same data. Therefore, loading
    'raw' or 'processed' data uxing read_xnmr_bin() will return the same thing.
    """

    change_params = {
        'npts':'',
        'na2':''}
    
    #Figure out the dimensions
    if data.ndim==2:
        if data.shape[0] > data.shape[1]:
            raise Exception('FIDs should be in rows, not in columns. First dimension should be smaller than the 2nd')
        change_params['na2'] = str(data.shape[0])
    else:
        change_params['na2'] = 1
    change_params['npts']=str(data.shape[1])
    
    #Transform to a linear array
    real = (data.real).flatten()
    imag = (data.imag).flatten()
    real = real.astype('float32')
    imag = imag.astype('float32')
    raw_out = np.empty(len(real)+len(imag),dtype=np.float32)

    raw_out[0::2] = real
    raw_out[1::2] = imag

    #Back up the parameter and data files (if not already done so)
    if os.path.exists(fname+'/backups'):
        print("Param backups exists")
    else:
        os.makedirs(fname+'/backups')
        shutil.copy(fname+'/params',fname+'/backups/params')
        print('Param file backed up')
        shutil.copy(fname+'/data.fid',fname+'/backups/data.fid')
        print('raw data (data.fid) backed up.')
        if os.path.exists(fname+'/data'):
            shutil.move(fname+'/data',fname+'/backups/data')
            print('Processed data files removed (moved to backup folder)')

    #Write out the data (will overwrite)
    try:
        raw_out.tofile(fname+'/data.fid')
        raw_out.tofile(fname+'/data')
    except:
        raise Exception('Error writing out the array as data.fid')
    print('Wrote data.fid')
    
    change_param_file(fname+'/params',change_params)
    print('Changed param file')

    print('save_xnmr_bin.py completed')
    return
    
