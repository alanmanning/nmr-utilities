from IPython.core.debugger import Tracer; debug_here = Tracer()
import numpy as np
from scipy.special import wofz
#File locations
drop_dir = "/home/alan/Dropbox/"
import imp
ft_analysis = imp.load_source("ft_analysis",drop_dir+"nmr-analysis/"+"ft_analysis.py")

#Load the shared c library and tell python what
#types to pass it.
import ctypes
from numpy.ctypeslib import ndpointer
lib = ctypes.cdll.LoadLibrary(drop_dir+'nmr-analysis/lineshapes_cfuncs.so')
sl_t_cfunc = lib.sl_t
sl_t_cfunc.restype = ctypes.c_int
sl_t_cfunc.argtypes = [ctypes.c_double,ctypes.c_double, ctypes.c_double, ctypes.c_double,
    ndpointer(ctypes.c_double,flags="C_CONTIGUOUS"),
    ndpointer(ctypes.c_double,flags="C_CONTIGUOUS"),
    ndpointer(ctypes.c_double,flags="C_CONTIGUOUS"),
    ctypes.c_int,ctypes.c_int,
    ndpointer(ctypes.c_double,flags="C_CONTIGUOUS"),
    ndpointer(ctypes.c_double,flags="C_CONTIGUOUS")]
sl_t_lorentz_cfunc = lib.sl_t_lorentz
sl_t_lorentz_cfunc.restype = ctypes.c_int
sl_t_lorentz_cfunc.argtypes = [ctypes.c_double,ctypes.c_double, ctypes.c_double, ctypes.c_double,
    ndpointer(ctypes.c_double,flags="C_CONTIGUOUS"),
    ndpointer(ctypes.c_double,flags="C_CONTIGUOUS"),
    ndpointer(ctypes.c_double,flags="C_CONTIGUOUS"),
    ctypes.c_int,ctypes.c_int,
    ndpointer(ctypes.c_double,flags="C_CONTIGUOUS"),
    ndpointer(ctypes.c_double,flags="C_CONTIGUOUS")]
sl_f_cfunc = lib.sl_f
sl_f_cfunc.restype = ctypes.c_int
sl_f_cfunc.argtypes = [ctypes.c_double,ctypes.c_double, ctypes.c_double, ctypes.c_double,
    ndpointer(ctypes.c_double,flags="C_CONTIGUOUS"),
    ndpointer(ctypes.c_double,flags="C_CONTIGUOUS"),
    ndpointer(ctypes.c_double,flags="C_CONTIGUOUS"),
    ctypes.c_int,ctypes.c_int,
    ndpointer(ctypes.c_double,flags="C_CONTIGUOUS")]

########
######## Python functions
########

def py_gauss_lor_t(amp,lorentz_wid, gauss_wid, center, t):
    lorentz = py_lor_t(1,lorentz_wid,0,t)
    gauss   = py_gauss_t(1,gauss_wid,0,t)
    return amp*lorentz*gauss*np.exp(-1j*2*np.pi*center*t)

# def py_gauss_lor_f(amp,lorentz_wid, gauss_wid, center, t):
#     if t[0]<0:
#         raise Exception('Need to pass time, not frequency (FT is used, not direct calculation)')
#     sig1 = py_gauss_lor_t(amp,lorentz_wid,gauss_wid,center,t)
#     sig2 = ft_analysis.do_ft_fftw(sig1,kind='fft')
#     return np.fft.fftshift(sig2)

def py_gauss_lor_f(amp,lorentz_wid, gauss_wid, center, f):
    #Computation of the Voigt profile using the Faddeeva function
    z = ((f-center)+1j*lorentz_wid)/gauss_wid/np.sqrt(2)
    return amp*wofz(z).real/gauss_wid/np.sqrt(2*np.pi)

def py_gauss_t(amp,wid,center,t):
    # #This is normalized for integral from t=0..inf
    # alpha = (2*np.pi*wid)/4/np.sqrt(np.log(2))
    # gauss = np.exp(-alpha**2*t**2)/np.sqrt(np.pi)/2/alpha
    #Not normalized
    alpha = (2*np.pi*wid)/4/np.sqrt(np.log(2))
    gauss = np.exp(-alpha**2*t**2)
    return amp*gauss*np.exp(-1j*2*np.pi*center*t)

def py_lor_t(amp,wid,center,t):
    # #This is normalized for integral from t=0..inf
    # lorentz = np.exp(-t/t2)/t2 
    #Not normalized
    lorentz = np.exp(-t*wid*np.pi)
    return amp*lorentz*np.exp(-1j*2*np.pi*center*t)

# def py_lor_f(amp,t2,center,f):
#     lorentz = amp/np.pi*( 1/t2 )/ ((1/t2**2) + (2*np.pi*(center-f))**2)
#     return amp*lorentz

# def py_gauss_f(amp,wid,center,f):
#     alpha = (2*np.pi*wid)/4/np.sqrt(np.log(2))
#     gauss = 1/np.sqrt(2)/alpha*np.exp(-((f-center)*2*np.pi/2/alpha)**2)
#     return amp*gauss

# def py_sl_f(amp,wid, wid0, center, p2, weights, freq):
#     wids = np.sqrt(wid0**2 + (wid*p2[np.newaxis,:])**2) 
#     alpha = (2*np.pi*wids)/4/np.sqrt(np.log(2))
#     gauss = NUMPY_IFFT_SCALE/np.sqrt(2)/alpha*np.exp(-((freq[:,np.newaxis]-center)*2*np.pi/2/alpha)**2)*weights[np.newaxis,:]
#     return amp*np.sum(gauss,axis=1)

# def py_sl_t(amp,wid, wid0, center,p2, weights, t):
#     wids = np.sqrt(wid0**2 + (wid*p2)**2) 
#     alpha = (2*np.pi*wids)/4/np.sqrt(np.log(2))
#     gauss = np.exp(-alpha**2*t[:,np.newaxis]**2)*weights[np.newaxis,:]
#     return amp*np.sum(gauss,axis=1)*np.exp(-1j*2*np.pi*center*t)

########
######## C helper functions (to do easy multithreading have to pass lists)
########

def c_sl_t(args):
    amp = args[0]
    wid = args[1]
    wid0 = args[2]
    center = args[3]
    p2 = args[4]
    weights = args[5]
    t = args[6]
    fid_re = np.empty(len(t))
    fid_im = np.empty(len(t))
    sl_t_cfunc(amp, wid, wid0,    center,  p2,      weights, t,    len(t),len(p2), fid_re, fid_im)
    return fid_re +1j*fid_im

def c_sl_t_lorentz(args):
    amp = args[0]
    wid = args[1]
    wid0 = args[2]
    center = args[3]
    p2 = args[4]
    weights = args[5]
    t = args[6]
    fid_re = np.empty(len(t))
    fid_im = np.empty(len(t))
    sl_t_lorentz_cfunc(amp, wid, wid0,    center,  p2,      weights, t,    len(t),len(p2), fid_re, fid_im)
    return fid_re +1j*fid_im
    
def c_sl_f(args):
    amp = args[0]
    wid = args[1]
    wid0 = args[2]
    center = args[3]
    p2 = args[4]
    weights = args[5]
    freq = args[6]
    spect_re = np.empty(len(freq))
    sl_f_cfunc(amp, wid, wid0, center,     p2,       weights, freq, len(freq),len(p2), spect_re)
    return spect_re

########
######## Testing functions
########
def test():
    import pylab as plt
    import time
    drop_dir = "/home/alan/Dropbox/"
    import imp
    ft_analysis = imp.load_source("ft_analysis",drop_dir+"nmr-analysis/"+"ft_analysis.py")
    

    t = 1e-6*np.arange(130e3)

    #Lorentzian
    t2  = 1e-3
    amp = 10
    center = -512





    # t = 1e-6*np.arange(130e3)
    # freq = np.linspace(-250e3,250e3,len(t))
    # wid = 20e3
    # wid0 = 1.0
    # center=-300
    # p2 = 0.5*(3*np.cos(np.linspace(0,np.pi/2,600))**2 - 1)
    # weights = np.ones(len(p2))
    
    # t1 = time.time()
    # fid1 = py_sl_t(1.0,wid, wid0, center, p2, weights, t)
    # t2 = time.time()
    # fid2 = c_sl_t([1.0,wid, wid0, center, p2, weights, t])
    # t3 = time.time()
    
    # print('Python took %g ms' % (1000*(t2-t1)))
    # print('C took %g ms' % (1000*(t3-t2)))
    
    # if np.allclose(fid1,fid2):
    #     print('FID test OK!')
    # else:
    #     print('FID test failed.')
    #     plt.plot(fid1.real,'-b',linewidth=2,alpha=0.5,label='Python')
    #     plt.plot(fid1.imag,'--b',linewidth=2,alpha=0.5)
    #     plt.plot(fid2.real,'-g',linewidth=2,alpha=0.5,label='C')
    #     plt.plot(fid2.imag,'--g',linewidth=2,alpha=0.5)
    #     plt.legend(frameon=False)
    #     plt.show(block=False)
    #     debug_here()
    
    # t1 = time.time()
    # spect1 = py_sl_f(1.0,wid, wid0, center, p2, weights, freq)
    # t2 = time.time()
    # spect2 = c_sl_f([1.0,wid, wid0, center, p2, weights, freq])
    # t3 = time.time()
    
    # print('Python took %g ms' % (1000*(t2-t1)))
    # print('C took %g ms' % (1000*(t3-t2)))
    
    # if np.allclose(spect1,spect2):
    #     print('Spect test OK!')
    # else:
    #     print('Spect test failed.')
    #     plt.plot(freq,spect1,'-b',linewidth=2,alpha=0.5, label='Python')
    #     plt.plot(freq,spect1,'-g',linewidth=2,alpha=0.5, label='C')
    #     plt.legend(frameon=False)
    #     plt.show(block=False)
    #     debug_here()
    







