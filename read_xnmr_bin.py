def read_xnmr_bin(fname,what='processed',silence=False):
    """Read binary Xnmr files directly into Numpy arrays.
    Nnmr data folders have two binary data files: 'data' and 'data.fid'. 'data' is the raw data, as acquired by the spectrometer. This will be the FID. 'data.fid' is the processed data. If you open up an Xnmr file, do some processing, then save it,  'data.fid' is created. Depending on what processing you do, this could be a spectrum or an FID.

    To return the raw data, pass what='raw'
    To return the processed data (could be a spectrum or FID), pass what='processed'

    There is one caveat: you could end up with an Xnmr file which only has 'data.fid' and no 'data' file. This occurs if you do some processing (maybe something like add and subtract), then save the result into a new Xnmr file. In this case, the function will return an error if you ask for what='raw'.

    Set silence=True to turn off all of the information messages

    fname is the path to the Xnmr file (a directory)
    """
    import numpy as np
    if not silence:
        print('Reading ', fname)
    f = open(fname+'/params','r')
    a = f.readlines()
    f.close()
    
    npts= 0
    na2 = 0
    dwell = 0
    na = 0
    for i in range(len(a)):
        ind = a[i].strip().find('npts =')
        if(ind > -0.5):
            npts = int(a[i].strip()[ind+6:])

        ind = a[i].strip().find('na =')
        if(ind > -0.5):
            na = int(a[i].strip()[ind+4:])

        ind = a[i].strip().find('na2 =')
        if(ind > -0.5):
            na2 = int(a[i].strip()[ind+5:])

        ind = a[i].strip().find('dwell =')
        if(ind > -0.5):
            dwell = float(a[i].strip()[ind+7:])*1e-6
    
    print('npts=%i, na=%i, na2=%i, dwell=%.0e' % (npts,na,na2,dwell))
    
    if( (npts==0) or (na2==0) or (dwell<1e-12)):
        raise Exception('ERROR: read_Xnmr. Got 0 na, na2, or dwell')
    
    
    if(what=='processed'):
        a=np.fromfile(fname+'/data',np.float32)
    elif(what=='raw'):
        try:
            a=np.fromfile(fname+'/data.fid',np.float32)
        except:
            raise Exception('Could not find data.fid. Xnmr file may have been saved as a processed spectrum.')
    else:
        raise Exception('ERROR: read_xnmr_bin. Need to specify "processed" or "raw"')

    real = a[0::2].reshape((na2,-1))
    imag = a[1::2].reshape((na2,-1))
    
    if(((real.shape[0] != na2) or (imag.shape[0] != na2)) and not what=='processed'):
        raise Exception('ERROR: readXnmr. reshaped real/imag matrices dimension mismatch na2/npts')
    
    #now, calculate the frequencies
    if(what=='processed'):
        if not silence:
            print('Returning processed data')
        f = open(fname+'/proc_params','r')
        l = f.readlines()
        f.close()

        spect = True
        for line in l:
            if line.strip().find('FT:') > -0.5:
                if int(line[3:].strip()) == 0:
                    spect = False

        if spect:
            if not silence:
                print('Processed data is frequency domain')
            freqs = -1/dwell*(np.arange(npts)-npts/2)/npts
            return freqs,real,imag,'freq'
        else:
            if not silence:
                print('Processed data is time domain')
            time = np.arange(npts)*dwell
            return time,real,imag,'time'


    else:
        if not silence:
            print('Returning unprocessed time domain data')
        time = np.arange(npts)*dwell
        return time,real,imag,'time'

