def read_xnmr_bin(fname,what='processed',silence=False):
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
        print('ERROR: read_Xnmr. Got 0 na, na2, or dwell')
        return -1
    
    
    if(what=='processed'):
        a=np.fromfile(fname+'/data',np.float32)
    elif(what=='raw'):
        a=np.fromfile(fname+'/data.fid',np.float32)
    else:
        print('ERROR: read_xnmr_bin. Need to specify "processed" or "raw"')
        return -1

    real = a[0::2].reshape((na2,-1))
    imag = a[1::2].reshape((na2,-1))
    
    if(((real.shape[0] != na2) or (imag.shape[0] != na2)) and not what=='processed'):
        print('ERROR: readXnmr. reshaped real/imag matrices dimension mismatch na2/npts')
        return -1
    
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

