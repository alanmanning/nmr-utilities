# -*- coding: utf-8 -*-



#######################################################
##### USEFUL FUNCTIONS
#######################################################

def get_ppm(nu_unknown,nu_known,delta_known):
    """Calculate the ppm using an external reference.
    delta_known = the shift of the known reference (eg. if this is 10 ppm, delta_known=10e-6)
    nu_known = the frequency corresponding to the shift of the known reference (in Hz)
    nu_known = the frequencies to convert to ppm
    """
    assert abs(delta_known) < 1, "delta_known is not in ppm, its in absolute units (eg. if this is 10 ppm, delta_known=10e-6)"
    
    nu_ref = nu_known/(delta_known+1)
    return (nu_unknown - nu_ref)/nu_ref

def convert_tensor(iso=None,aniso=None,asym=None,xx=None,yy=None,zz=None):
    """Calculate between tensor conventions
    See: http://anorganik.uni-tuebingen.de/klaus/nmr/index.php?p=conventions/csa/csa
    Use with convert_tensor(iso=arg1,aniso=arg2,asym=arg3) to get xx, yy, zz
    Use with convert_tensor(xx=arg1,yy=arg2,zz=arg3) to get iso, aniso, asym
    """
    import numpy as np

    if (iso is not None) and (aniso is not None) and (asym is not None) and (xx is None) and (yy is None) and (zz is None):
        aniso_red = 2*aniso/3. #reduced anisotropy
        shift = np.zeros(3)
        shift[0] = iso + aniso_red
        shift[1] = iso - aniso_red*(1-asym)/2
        shift[2] = iso - aniso_red*(1+asym)/2

        #Convert to weird Haeberlen convention ordering
        inds = np.argsort(np.abs(shift - iso))
        yy = shift[inds[0]]
        xx = shift[inds[1]] 
        zz = shift[inds[2]]

        #test:
        t_iso = np.average((xx,yy,zz))
        t_aniso1 = (zz - t_iso)*3./2.
        t_aniso2 = zz - (xx+yy)/2.
        t_asym = (yy-xx)/(t_aniso1*2/3)

        assert np.abs(t_iso-iso)<1e-10, 'iso messed up'
        assert np.abs(t_aniso1-t_aniso2)<1e-10, 'two ways of calculating aniso messed up'
        assert np.abs(t_aniso1-aniso)<1e-10, 'aniso messed up'
        assert np.abs(t_asym - asym)<1e-10, 'asym messed up'

        shift = {
            'yy' : shift[inds[0]],
            'xx' : shift[inds[1]], 
            'zz' : shift[inds[2]],
        }

        return shift

    elif (iso is None) and (aniso is None) and (asym is None) and (xx is not None) and (yy is not None) and (zz is not None):

        iso = (xx+yy+zz)/3.0

        vals = np.array([xx,yy,zz])
        inds = np.argsort(np.abs(vals-iso))
        yy = vals[inds[0]]
        xx = vals[inds[1]]
        zz = vals[inds[2]]

        aniso_red = zz-iso
        aniso = 3*aniso_red/2
        asym = (yy-xx)/aniso_red

        shift = {
            'iso' : iso,
            'aniso' : aniso,
            'asym' : asym,
        }

        return shift

    else:
        raise Exception('You passed the arguments incorrectly')

#######################################################
##### FUNCTIONS WHICH ARE PROBABLY NOT USEFUL ANYMORE
#######################################################
# These were written to deal with the older, Vnmr data. The 400 MHz
# spectrometer now has VnmrJ software.

def plot_spectra_varian(filenames,labels,xlimits,ref_trans_freq,ref_shift_hz,ref_shift_ppm,samp_trans_freq=-1,exportfilenames=[]):
    import matplotlib.pyplot as plt
    import numpy as np
    
    print(exportfilenames)
   
    thefig = plt.figure();
    thefig.set_facecolor("#FFFFFF");
    for i in range(len(filenames)):
        f = filenames[i]
        if len(filenames) == len(exportfilenames):
            expf = exportfilenames[i]
            vals = read_1d_xnmr('/home/manning/nmr/nmr_data/' + expf,False,False)
        else:
            vals = read_1d_xnmr('/home/manning/nmr/nmr_data/' + f + '.fidexport.txt',False,False)
        samp_freqs = vals[0];
        samp_sigs = vals[1]/max(vals[1])
        
        if(samp_trans_freq < 0):
            samp_trans_freq = get_varian_param('/home/manning/nmr/nmr_data/' + f + '.fid','sfrq')*1e6
        samp_ppms = get_ppm(samp_trans_freq, samp_freqs, ref_trans_freq, ref_shift_hz, ref_shift_ppm)
        plt.subplot(len(filenames),1,i+1)
        plt.plot(samp_ppms,samp_sigs,label=labels[i])
        if len(xlimits) > 1:
            plt.xlim(xlimits)
        plt.grid()
        if i < len(filenames) -1:
            plt.gca().axes.set_xticklabels([])
        plt.gca().axes.set_yticklabels([])
        plt.gca().invert_xaxis()
        plt.legend()
    
    plt.gca().axes.get_xaxis().set_visible(True)
    plt.show(block=False)
    plt.ylabel('Relative intensity')
    plt.xlabel('PPM')


def get_varian_param(filename,paramname,is_string=False):
    from os import system

    grep_string = 'grep -A 1 "^' + paramname + ' " ' + filename + '/procpar'
    system(grep_string + '> .gtmp')
    grep_output = open('.gtmp') #stupid way of doing it but necessary on a mac.
    lines = grep_output.readlines()
    grep_output.close()
    l=lines[1].replace("\n","").strip()
    
    
    
    if(is_string==False):
        vals = [float(x) for x in l.split(" ")]
    else:
        vals = l.replace('"','').strip();
        
    if len(vals) > 2:
        return vals[1:]
    else:
        return vals[-1]

def get_xnmr_param(filename,paramname,is_string=False):
    from os import system
    import numpy as np

    #print filename
    grep_string = 'grep -w ' + paramname + ' ' + filename + '/params'
    system(grep_string + '> .gtmp')
    grep_output = open('.gtmp') #stupid way of doing it but necessary on a mac.
    lines = grep_output.readlines()
    grep_output.close()
    out = np.zeros(len(lines))
    for i in range(len(lines)):
        l = lines[i].replace('\n','')
        l = l[l.index('=')+2:]
        out[i] = float(l)
    return out

def get_varian_time_start(filename):
    import datetime as dt
    import numpy as np
    tc = get_varian_param(filename,'time_complete',is_string=True)
    time_completed = dt.datetime(int(tc[1:5]),int(tc[5:7]),int(tc[7:9]),int(tc[10:12]),int(tc[12:14]),int(tc[14:16]))
    nt = get_varian_param(filename,'nt')
    pad = get_varian_param(filename,'pad')
    d1 = get_varian_param(filename,'d1')
    
    l=1
    nt_ar = False
    pad_ar = False
    if isinstance(nt,list):
        nt_ar = True
        l = len(nt)
    if isinstance(pad,list):
        pad_ar = True
        if nt_ar and (len(pad) != l):
            input('ERROR--pad and nt lengths aren''t the same -- (nmrtools.get_varian_time_run)')
        l = len(pad)
    
    if l==1:
        time_start = time_completed - dt.timedelta(seconds=pad + nt*d1)
        return time_start
    else:
        if nt_ar and not pad_ar:
            nt = np.asarray(nt)
            pad = pad*np.ones(len(nt))
        elif pad_ar and not nt_ar:
            pad = np.asarray(pad)
            nt = nt*np.ones(len(pad))
        elif pad_ar and nt_ar:
            pad = np.asarray(pad)
            nt = np.asarray(nt)
        else:
            input('ERROR, arrays mixed up -- nmrtools.get_varian_time_run')
        
        time_start = time_completed - dt.timedelta(seconds=sum(pad) + sum(nt*d1))
        times = [time_start]
        for i in range(1,len(pad)):
            the_pad = pad[i]
            the_nt = nt[i]
            old_time = times[-1]
            times.append(old_time + dt.timedelta(seconds=the_pad + the_nt*d1))
        return times
            
            
    

def get_varian_time_completed(filename):
    import datetime as dt
    import numpy as np
    tc = get_varian_param(filename,'time_complete',is_string=True)
    time_completed = dt.datetime(int(tc[1:5]),int(tc[5:7]),int(tc[7:9]),int(tc[10:12]),int(tc[12:14]),int(tc[14:16]))
    return time_completed

def get_varian_time_midscan(filename):
    import datetime as datetime
    ts = get_varian_time_start(filename)
    tc = get_varian_time_completed(filename)
    
    
    if isinstance(ts,list):
        tmid = []
        for i in range(len(ts)):
            tmid.append(datetime.timedelta(seconds=0.5*(tc[i] - ts[i]).total_seconds()))
    else:
        tmid = datetime.timedelta(seconds=0.5*(tc - ts).total_seconds())
    
    return tmid
        
# def calc_ppm(delta_ref, nu_ref, nu_samp):
#     nu_0 = nu_ref/(delta_ref+1)
#     delta = (nu_samp - nu_0)/nu_0
#     return delta



## Old get_ppm: this has a bug that introduces a small (~1 ppm error)        
# def get_ppm(samp_trans_freq, samp_freqs, ref_trans_freq, ref_shift_hz, ref_shift_ppm):
#     print('Depreciated and slightly incorrect. Use calc_ppm instead!')
#     print('FIX GET_PPM!!!')
#     f=samp_freqs+samp_trans_freq
#     ref=ref_shift_hz+ref_trans_freq
#     ref_shift_ppm_scaled = ref_shift_ppm/10**6
#     f0 = ref/(1+ref_shift_ppm_scaled)
#     shift = (f-f0)/f0
#     ppm = shift*10**6
#     return ppm
    
def read_1d_xnmr(filename,plotflag,loadall):
    import numpy
    print("read_1d_xnmr: loading...")
    
    m=numpy.genfromtxt(filename,skip_header=1)
    hz = m[:,1]
    real = m[:,2]
    
    if(plotflag):
        from pylab import plot, show
        plot(hz,real,'g')
        show(block=False)
        
    if(loadall):
        points = m[:,0]
        imag = m[:,3]
        return [points,hz,real,imag]
    else:
        return [hz,real]
    
    
def read_2d_xnmr(filename):
    import numpy as np
    #from IPython.core.debugger import Tracer; debug_here = Tracer()

    print("read_2d_xnmr: loading...", end=' ')
    
    m = np.transpose(np.loadtxt(filename,skiprows=1,delimiter=' '))
  
    #if 2d time domain data is being read in
    if m.shape[0] == 6:
    #x point, time, y point, time, real, imag 
        x_point_num = m[0]
        x_t = m[1]
        y_point_num = m[2]
        y_t = m[3]
        real = m[4]
        imag = m[5]
       
        real = np.reshape(real,(len(np.unique(y_point_num)),-1)) 
        imag = np.reshape(imag,(len(np.unique(y_point_num)),-1)) 
        x_t = np.reshape(x_t,(len(np.unique(y_point_num)),-1))

        return x_t,real,imag
    
    #if no 2d time domain data is being read in (for experiments
    #aquired in Xnmr, generally)
    elif m.shape[0] == 5: 
    #x point, time, y point, real, imag
        x_point_num = m[0]
        x_t = m[1]
        y_point_num = m[2]
        real = m[3]
        imag = m[4]
        real = np.reshape(real, (len(np.unique(y_point_num)),-1))
        imag = np.reshape(imag, (len(np.unique(y_point_num)),-1))
        x_t = np.reshape(x_t,(len(np.unique(y_point_num)),-1))
        return x_t,real,imag
        
    else:
        print('UNKNOWN FORMAT')
        debug_here()
   # m=numpy.genfromtxt(filename,skip_header=1)
   # print m.shape 
   # if len(m) == 5: #x point, time, y point, real, imag
   #    x = m[:,0]
   #     t = m[:,1]
   #     y = m[:,2]
   #     real = m[:,3]
   #     imag = m[:,4]
   #     print '# unique x=',numpy.unique(x)
   #     print '# unique y=',numpy.unique(y)
   # else: #x point, freq, y point, time, real, imag
   #    x = m[:,0]
   #     f = m[:,1]
   #     y = m[:,2]
   #     t = m[:,3]
   #     real = m[:,4]
   #     imag = m[:,5]
   # pts = m[:,0]
   # xlen = numpy.diff(pts,n=1).argmin() + 1
   # ylen = pts.size/xlen
   # if(abs(ylen-round(ylen)) > 0.0001):
   #     import warnings
   #     warnings.warn("Xnmr data not read correctly")
   # ylen=int(round(ylen))
   # print "loaded. size = " + str(xlen) + " x " + str(ylen)
   # 
   # #x point, freq, y point, time, real, imag

   # hz=[]
   # for i in range(1,ylen+1):
   #     hz.append(m[(i-1)*xlen:i*xlen,1])  
   # real=[]
   # for i in range(1,ylen+1):
   #     real.append(m[(i-1)*xlen:i*xlen,4])

   # if(loadall):
   #     x_points=[]
   #     for i in range(1,ylen+1):
   #         x_points.append(m[(i-1)*xlen:i*xlen,0])
   #     y_points=[]
   #     for i in range(1,ylen+1):
   #         y_points.append(m[(i-1)*xlen:i*xlen,2])
   #     time=[]
   #     for i in range(1,ylen+1):
   #         time.append(m[(i-1)*xlen:i*xlen,3])        
   #     imag=[]
   #     for i in range(1,ylen+1):
   #         imag.append(m[(i-1)*xlen:i*xlen,5])
   #         
   #     return [x_points, hz, y_points, time, real, imag]
   # else:
   #     return [hz,real]   
    
    
def find_main_peak_limits(sig,scale_down,scale_up,plotflag):
    import numpy
    
    max_ind = sig.argmax();
    ulim = int(max_ind + scale_down*sig.size)
    llim = int(max_ind - scale_up*sig.size)
    
    peak_ok = (ulim < (sig.size - 1)) and (llim > 0)
    
    if(peak_ok):
        peak_height = numpy.mean(sig[round(max_ind - 0.01*sig.size):round(max_ind + 0.01*sig.size)])
        baseline = (numpy.mean(sig[0:llim]) + numpy.mean(sig[ulim:-1]))/2
        st_dev = (numpy.std(sig[0:llim]) + numpy.std(sig[ulim:-1]))/2
        peak_ok = peak_ok and ((baseline + st_dev)/peak_height < 1.0)
    
    if(not(peak_ok)):
        llim = ulim = 0
        if(plotflag):
            from pylab import plot, show
            plot(sig,'b',[llim,llim],[sig.min(),sig.max()],'r',[ulim,ulim],[sig.min(),sig.max()],'r',)
            show()
        
    if(plotflag):
        from pylab import plot, show
        plot(sig,'g',[llim,llim],[min(sig),max(sig)],'r',[ulim,ulim],[min(sig),max(sig)],'r')
        show(block=False)
        
    return [llim,ulim]
    

def integrate_peaks(filename,plotflag):
    import numpy
    print("integrate_peaks with" + filename)
    
    [hz,sig] = read_2d_xnmr(filename,False,False)

    sig_ints = []
    for i in range(0,len(hz)):
        [llim,ulim] = find_main_peak_limits(sig[i],0.1,0.1,False)
        sig_ints.append(sig[i][llim:ulim].sum())
        if (llim == ulim) or (sig_ints[-1] < 0):
            sig_ints[-1] = -1
            
    if(plotflag):
        from pylab import plot, show
        plot(sig_ints)
        show()
        
    sig_ints = numpy.array(sig_ints)
    fid=open(filename + ".peakint","w")
    sig_ints.tofile(fid,sep=",",format="%e")
    fid.close()

    return sig_ints
    

def get_diffusion_values(procpars,peakints,outname):
    import numpy
    from pylab import plot, show
    

    gscale =  0.163698e-4
    gamma   = 267.522e+6
    
    
    old_q = numpy.array([0])
    allDs = []
    allDELTAs = []
    
    for i in range(0,len(peakints)):
        pfile = procpars[i]
        intfile = peakints[i]
        print("\n\n=================" + intfile)
        
        DELTA  = get_varian_param(pfile,"BigT")
        g      = gscale*numpy.asarray(get_varian_param(pfile,"gzlvl1"))
        delta  = get_varian_param(pfile,"gt1")
        
        sig = numpy.genfromtxt(intfile,delimiter=",")
        
        ok = False
        good_vals = sig > 0
        if numpy.count_nonzero(good_vals) > 3:
            q      = gamma*g[good_vals]*delta/(2*numpy.pi)
            ylin  = numpy.log(sig[good_vals])
            xlin  = q*q
            [m,c] = numpy.polyfit(xlin,ylin,1)
            D     = -m/((2*numpy.pi)**2*(DELTA-delta/3))*10**4
            #plot(xlin, ylin,'.k')
            #plot(xlin, xlin*m + c, 'b')
            #show();
            
            if D<0:
                print("BAD FIT, NOT USING")
            elif (len(old_q) > 1) and (old_q != q).any():
                print("BIG PROBLEM, q'S ARE DIFFERENT!")
            else:
                ok = True
                print("D = " + str(float(D)))
                allDELTAs.append(DELTA)
                allDs.append(D)
                
        if not ok:
            print("NO GOOD DATA IN FILE")
            allDELTAs.append(DELTA)
            allDs.append(-1)

    fid=open(outname,"w")
    numpy.asarray(allDELTAs).tofile(fid,sep=",",format="%e")
    fid.close()
    fid=open(outname,"a")
    fid.writelines("\n")
    numpy.asarray(allDs).tofile(fid,sep=",",format="%e")
    fid.close()
    
     
def get_model_echo(gamma,Delta,g,delta,Dbulk,R):
    import numpy
    exparg = -2*gamma*gamma*g*g*s
    R = numpy.exp(exparg)
    return R
    
def get_model_sum(Delta,delta,Dbulk,R):
    import numpy
    eqn_roots = [2.08157597781810,5.94036999057271, 9.20584014293666, 12.4044450219020,15.5792364103872, 18.7426455847748, 21.8996964794928, 25.0528252809930,28.2033610039524, 31.3520917265645,  34.4995149213670, 37.6459603230864]
    alpha_roots = numpy.asarray(eqn_roots)/R    
    s=0
    for j in range(1):
        alpha2 = alpha_roots[j]*alpha_roots[j]
        alpha2_R2  = alpha2*R*R
        alpha2_D   = alpha2*Dbulk
        term1 = 1/(alpha2*(alpha2_R2-2))
        term2 = 2*delta/alpha2_D
        term3 = numpy.exp(-alpha2_D*(Delta - delta))
        term4 = 2*numpy.exp(-alpha2_D*Delta)
        term5 = 2*numpy.exp(alpha2_D*delta)
        term6 = numpy.exp(-alpha2_D*(Delta+delta))
        topfrac = 2 + term3 - term4 - term5 + term6
        brackets = term2 - topfrac/alpha2_D**2
        s =  s + term1*brackets
    return s

def get_model_diffusion(Delta,delta,R,Dbulk):
    import numpy as np
    Deff = get_model_sum(Delta,delta,Dbulk,R)/(delta**2*(Delta-delta/3))
    np.nan_to_num(Deff)
    return Deff   
     
        

def plot_matrix(mat,xar=[],yar=[],xlabel=[],ylabel=[],title=[]):
    import numpy as np
    import pylab as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    fig,ax = plt.subplots(1,1)
    vmin = np.floor(np.min(mat))
    vmax = np.ceil(np.max(mat))
    print(title,'vmin=',vmin,'vmax',vmax)
    if len(xar)<1e-6:
        xar = list(range(mat.shape[1]))
    if len(yar)<1e-6:
        yar = list(range(mat.shape[0]))
    im = ax.imshow(mat, cmap='jet',
        vmin=vmin,vmax=vmax,
        extent=[xar[0],xar[-1],yar[-1],yar[0]],aspect='auto')

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    fig.suptitle(title,size=16)

    #divider = make_axes_locatable(ax)
    #cax = divider.append_axes("right", size="20%", pad=0.05)
    #cbar = plt.colorbar(im, cax=cax, format="%.2f")
    plt.show()
