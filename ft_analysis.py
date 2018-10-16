import numpy as np
from IPython.core.debugger import Tracer; debug_here = Tracer()
import pyfftw
import pickle
from scipy.optimize import minimize_scalar
import pylab as plt

#def do_baseline_correct1(data):
#    return data - np.average(data)

def do_baseline_correct2(data,points):
    """Do a baseline correction on a spectrum.
    You probably don't want to use this function, it never worked too well.
    """
    if data.ndim==2:
        if data.shape[0] > data.shape[1]:
            raise Exception('ft_analysis.do_ft: Data expected in rows, not columns. data.shape[0] > data.shape[1]!')

    if data.ndim==2:
        if len(points)==0:
            offsets = 0.5*np.average(data[:,0:200],axis=-1)+np.average(data[:,0:200],axis=-1)
        else:
            offsets = np.average(data[:,points[0]:points[1]+1],axis=-1)
        return offsets,data-offsets[:,np.newaxis]
    else:
        if len(points)==0:
            offset = 0.5*(np.average(data[0:200])+np.average(data[-200:])) 
        else:
            offset = np.average(data[points[0]:points[1]+1])
        return offset,data - offset 

def do_broadening(t,data,gf):
    """Apply Gaussian broadening
    """
    if data.ndim==2:
        if data.shape[0] > data.shape[1]:
            raise Exception('ft_analysis.do_ft: Data expected in rows, not columns. data.shape[0] > data.shape[1]!')

    if data.ndim==2:
        t=t[np.newaxis,:]
    return data*np.exp(-t*np.pi*gf/1.6651*t*np.pi*gf/1.6651)

def get_freqs(n,dwell):
    """Get the frequencies associated with do_ft_fftw or do_ft.
    n=number of datapoints acquired in the FID.
    dwell=dwell of the FID
    """
    return np.fft.fftshift(np.fft.fftfreq(n,dwell))

def do_ft(t,data,gf=None,bc=None):
    """Calculate the spectrum from an FID using Numpy's FT libraries.
    Don't use this, it's very slow. Use do_ft_fftw below instead.
    """
    if data.ndim==2:
        if data.shape[0] > data.shape[1]:
            raise Exception('ft_analysis.do_ft: Data expected in rows, not columns. data.shape[0] > data.shape[1]!')


    if gf is not None:
        data = do_broadening(t,data,gf)
    ft = np.fft.ifft(data,axis=-1)
    ft = np.fft.fftshift(ft,axes=-1)
    if bc is not None:
        bcor,ft = do_baseline_correct2(ft,bc)
        return bcor,ft
    else:
        return ft

def do_ft_fftw(data,kind='fft',wisdom_path='/home/alan/Dropbox/nmr-analysis/fftw_wisdom.bin'):
    """Calculate the spectrum from an FID using FFTW.
    Use this for calculating FT's. Can handle one and two-dimensional data.
    
    FFTW uses a "wisdom" file, which stores FFTW's knowledge of the fastest way to 
    perform the FT for a specific array. This will vary depending on the datatype, the size, and
    the CPU. You'll need to change the location of the wisdom file to somewhere on your machine.

    kind='fft' (forward transform) or 'ifft' (inverse transform) NOTE: These are actually
    mapped to the opposite transforms in FFTW to keep it consistent with the way Xnmr
    does the transforms.
    """
    if data.ndim==2:
        if data.shape[0] > data.shape[1]:
            raise Exception('ft_analysis.do_ft: Data expected in rows, not columns. data.shape[0] > data.shape[1]!')

    try:
        wisdom = pickle.load(open(wisdom_path,'rb'))
    except:
        print('no stored wisdom')
    else:
        ret = pyfftw.import_wisdom(wisdom)

    dtype = data.dtype
    if data.dtype == np.complex128:
        dtype = 'complex128'
    elif data.dtype == np.complex256:
        dtype = 'complex256'
    elif data.dtype == np.complex64:
        dtype = 'complex64'
    elif data.dtype == np.float128:
        dtype = 'float128'
    elif data.dtype == np.float64:
        dtype = 'float64'
    elif data.dtype == np.float32:
        dtype = 'float32'
    else:
        raise Exception('Add support for more dtypes!')

    a = pyfftw.empty_aligned(data.shape, dtype=dtype)
    a[:] = data

    # The kind and the FT performed seems backwards. That's just how Xnmr does it
    # so I'm following that convention
    if kind=='fft':
        b = pyfftw.interfaces.numpy_fft.ifft(a,planner_effort='FFTW_MEASURE',threads=1)
    elif kind=='ifft':
        b = pyfftw.interfaces.numpy_fft.fft(a,planner_effort='FFTW_MEASURE',threads=1)
    else:
        raise Exception('do_ft_fftw: Invalid kind of ft. got kind='+str(kind))
    

    pickle.dump(pyfftw.export_wisdom(),open(wisdom_path,'wb'))

    return b

def baseline_correct(spectra,freqs):
    """Do a baseline correction on a spectrum.
    You probably don't want to use this function, it never worked too well.
    """
    sw = np.max(freqs) - np.min(freqs)
    window = 0.05*sw
    i1 = np.logical_and(freqs<np.max(freqs),freqs>np.max(freqs)-window)
    i2 = np.logical_and(freqs>np.min(freqs),freqs<np.min(freqs)+window)
    if spectra.ndim==2:
        if spectra.shape[0] > spectra.shape[1]:
            raise Exception('ft_analysis.real_spectra_to_fid: Data expected in rows, not columns. data.shape[0] > data.shape[1]!')
        bcors = np.median( np.concatenate((spectra[:,i1],spectra[:,i2]), axis=1), axis=1)[:,np.newaxis]
    elif spectra.ndim==1:
        try:
            bcors = np.median(np.concatenate((spectra[i1],spectra[i2])))
        except:
            debug_here()
    return bcors

def real_spectra_to_fid(spectra,freqs):
    """Convert a spectrum which only contains the real part to a complex FID.
    This is possible because of the Kramers-Kronig relations / Sokhotski–Plemelj theorem / Hilbert transform
    (the same thing with different name).
    """
    if spectra.ndim==2:
        if spectra.shape[0] > spectra.shape[1]:
            raise Exception('ft_analysis.real_spectra_to_fid: Data expected in rows, not columns. data.shape[0] > data.shape[1]!')
    if (freqs[3:8] < 0).any():
        raise Exception('Freqs and data should not be shifted')
    if (np.abs(spectra.imag) > 1e-9).any():
        raise Exception('Spectrum should not have any imaginary parts. Use regular ifft.')

    assert np.sum(spectra.imag)==0, 'Spectra must have zero imaginary component'

    fid = do_ft_fftw(np.real(spectra),kind='ifft')
    fid += np.sign(freqs)*fid
    #fid = np.fft.fft(np.fft.ifftshift(spectra.real+hilbert(spectra.real,axis=-1).imag*1j,axes=-1),axis=-1)
    return fid



def auto_phase_imag_to_zero(t,fid,power=4,full_output=False,plot=False):
    """Auto-phase the FID.
    This function autophases the FID by phasing until the imaginary part
    at t=0 is 0. It uses polynomial fitting (with a degree determined by power).
    """

    def resid(phase0,x,y):
        p = np.polyfit(x,(y*np.exp(-1j*phase0/180*np.pi)).imag,power)
        return p[power]**2

    t_start = 1*t
    fid_start = 1*fid
    result = minimize_scalar(resid, args=(t_start,fid_start), bounds=(-3,3),method='bounded', tol=None, options=None)

    fid_start *= np.exp(-1j*result['x']/180*np.pi)
    p = np.polyfit(t_start,fid_start.imag,power)
    t_new = np.arange(0,t_start[-1],t_start[1]-t_start[0])
    fit = np.zeros_like(t_new)
    for i in range(len(p)):
        fit += p[i]*np.power(t_new,power-i)

    if plot:
        plt.plot(t_start,fid_start.imag,'.b')
        plt.plot()
        plt.plot(t_new,fit,'-g')
        plt.show()

    # print('Auto-phase by backwards prediction: phase0=%f' % result['x'])
    if full_output:
        return result['x'], p
    else:
        return result['x']


def auto_phase_real_to_mag(fid_start,start_sign):
    """Auto-phase the FID.
    This function autophases the FID by setting the real part equal to the
    imaginary part. If the real part should be negative (eg. after an inversion
    pulse), then start_sign should be 1.
    """
    if start_sign not in [-1,1]:
        raise Exception('auto_phase_real_to_mag: start_sign must be 1 or -1')
    
    a=np.average(fid_start)
    phase0 = np.arctan(a.imag/a.real*start_sign)*180/np.pi
    print('Auto-phase by matching real with magnitude: phase0=%f' % phase0)
    return phase0