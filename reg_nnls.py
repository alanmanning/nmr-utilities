import scipy.optimize as opt
import numpy as np
import lmfit
import pylab as plt

def fit_reg_nnls(fname=None,data=None,t2=np.logspace(-3,0.3,200),
    init_lambda=1,doplot=False,writefiles=False,chisq_misfit=[1.02,1.025]):
    """ Perform regularized NNLS fitting to CPMG curves to extract T2 distributions.
    This is a pretty clunky function. It works, but there's a lot of extra junk.

    You can either load data from a txt file (fname) or passed to data directly.
    If passing to data, it expects data[0] to be CPMG decay time, data[1] to be CPMG decay amplitude.

    t2s are the distribution x-values.

    init_lambda is the initial regularization parameter guess.

    doplot says whether to plot the output (useful for debugging)

    writefiles: an option to write out txt files with the distribution. Probably not useful.

    chisqr_misfit = regularized NNLS always fits worse than non-regularized NNLS. Lambda
    (the regularization parameter) is adjusted until the regularized NNLS chi-squared is
    a certain percentage worse than the non-regularized NNLS. This target percentage falls
    within the range defined by chisqr_misfit. The default value of [1.02,1.025] means the
    regularized NNLS value will be 2.0 to 2.5% worse than the non-regularized value.

    RETURNS:
    fit_nreg_y, fit_reg_y, t2, x_nreg, x_reg, final_lambda, chisq_reg, chisq_nreg

    fit_nreg_y = the fit to CPMG curve from non-regularized NNLS
    fit_reg_y = the fit to CPMG curve from regularized NNLS
    t2 = the t2 distribution x-values
    x_nreg = the t2 distribution from non-regularized NNLS (confusingly called "x", when it's more of a "y".)
    x_reg = the same for regularized NNLS
    final_lambda = the final lambda value used in the regularized NNLS
    chisqr_reg = the chi-squared value from the regularized NNLS
    chisqr_nreg = the chis-squared value from the non-regularized NNLS

    If you just want the regularized NNLS distribution, plot x_reg vs. t2!
    """


    def try_lambda(params,C,d,chisq_nreg=-1,chisq_misfit=-1,return_x=False):
        trial_lambda = params['lambda_val']*1
        C_local = 1*C
        C_local[len(data_y):,:] = C_local[len(data_y):,:]*trial_lambda #create C = col[A,lambda*I]
        x,chisq = opt.nnls(C_local,d)
    
        if(return_x): #if we aren't using this function for fitting
            return x,chisq
        else: 
            #we are using this function for fitting with COBYLA routine
            #compare chisq, return non-zero with linear function if outside intended range
    
            rel =  chisq/chisq_nreg
            # print('(%f,%f) ' %(trial_lambda,rel),end='')
            #if (rel >= chisq_misfit[0]) and (rel <= chisq_misfit[1]):
            #    ret = 0
            #elif rel < chisq_misfit[0]:
            #    ret = (chisq_misfit[0]-rel)**2
            #else: # rel>misfit[1]
            #    ret = (rel-chisq_misfit[1])**2
            ret = (rel-0.5*chisq_misfit[0]-0.5*chisq_misfit[1])**2
            #ret = (rel - chisq_misfit[0])**2
            # print('trial_lambda=',trial_lambda,' ret=',ret)
            return ret

    if fname is not None:
        # print('Loading data from file')
        time, data_y = np.loadtxt(fname,unpack=True)
    elif data is not None:
        # print('Using data passed to function')
        time, data_y = data
    else:
        raise Exception('Need to provide an fname to load data from or pass [time,cpmg_amp] directly to data')

    
    #NNLS is ||Ax-data_y||^2 subject to x>=0
    t2_2d,time_2d = np.meshgrid(t2,time)
    A = np.exp(-time_2d/t2_2d)
#    A[:,-1] = 0*A[:,-1] + 1 #last row is for constant offset
    
    #do the NNLS (non-regularized) to get a chisq value to compare to
    # print('DOING NON-REGULARIZED NNLS') 
    x_nreg,chisq_nreg = opt.nnls(A,data_y)
    fit_nreg_y = np.dot(A,x_nreg)
    # print('done!')
    
    #Now construct matrices for regularized NNLS
    C = np.concatenate((A,np.eye(len(t2)))) #C = col[A,lambda*I]... but we will multiply the bottom part by lambda later
    d = np.concatenate((data_y,np.zeros(len(t2)))) #create d = col[data_y, 0]
    
    param0 = lmfit.Parameters()
    param0.add('lambda_val',value=init_lambda,min=0)
    # print('DOING REGULARIZED NNLS') 
    # print('(trial lambda,chisq misfit ratio)=',end='')
    result = lmfit.minimize(try_lambda,param0,args=(C,d,chisq_nreg,chisq_misfit,False),method='nelder')
    final_lambda = result.params['lambda_val']*1
    # print('done!')
    
    #OK, we know what lambda we want. Now calculate the result (we've already done this, but it shouldn't take to long to do again)
    x_reg,chisq_reg = try_lambda(result.params,C,d,return_x=True)
    if(np.logical_or(chisq_reg/chisq_nreg < chisq_misfit[0], chisq_reg/chisq_nreg > chisq_misfit[1])):
        print('Fitting failed: final chisq not within misfit range!')
    # print('final lambda=%f, offset/max_amp=%f' %(final_lambda,x_reg[-1]/np.max(data_y)) )
    fit_reg_y = np.dot(A,x_reg)
    
    #Output some files if wanted
    if writefiles:
        print('Saving output text files')
        np.savetxt('nonreg-nnls-fit-'+fname,np.transpose([time,fit_nreg_y]),header='Echo time (s), non-regularized NNLS fit')
        np.savetxt('nonreg-nnls-t2-distribution-'+fname,np.transpose([t2,x_nreg]),header='T2 time (s), amplitude (exponential multiplying factor)')
        np.savetxt('reg-nnls-t2-distribution-'+fname,np.transpose([t2,x_reg]),header='T2 time (s), amplitude (exponential multiplying factor)')

    #make plots if wanted
    if doplot:
        plt.close('all')
        fig = plt.figure()
        ax = fig.add_subplot(111)
#        ax.plot(time,data_y,'.',label='data')
        ax.plot(time,fit_reg_y-data_y,'--k',linewidth=2,label=('regularized NNLS\nlambda=%f' % final_lambda))
        ax.plot(time,fit_nreg_y-data_y,':g',linewidth=2,label='non-regularized NNLS')
        ax.set_xlabel('Echo time (s)')
        ax.set_ylabel('Echo amplitude residuals')
        ax.legend(frameon=False)
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(t2,x_reg,'-b',linewidth=2,label='T2 distribution')
        ax.plot(t2,x_reg,'.k',linewidth=2)
        plt.legend(frameon=False,loc='upper left')
        ax.set_xscale('log')
        ax.set_xlabel('T2 times')
        ax.set_ylabel('Amplitudes')
        ax.set_title('Regularized NNLS T2 distribution')

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(t2,x_nreg,'-b',linewidth=2,label='T2 distribution')
        ax.plot(t2,x_nreg,'.k',linewidth=2)
        plt.legend(frameon=False,loc='upper left')
        ax.set_xscale('log')
        ax.set_xlabel('T2 times')
        ax.set_ylabel('Amplitudes')
        ax.set_title('Non-Regularized NNLS T2 distribution')
        plt.show()

    return fit_nreg_y, fit_reg_y, t2, x_nreg, x_reg, final_lambda, chisq_reg, chisq_nreg

