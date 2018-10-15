def fit_reg_nnls(fname,data=[],t2s=[],init_lambda=1,doplot=False,writefiles=False,chisq_misfit=None):
    import scipy.optimize as opt
    import numpy as np
    from IPython.core.debugger import Tracer; debug_here = Tracer()
    import lmfit
    import pylab as plt

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
    
    if len(t2s)==0:
        # print('No t2s given. Using default: t2 = np.logspace(-3,0,500)')
        t2 = np.logspace(-3,0,200) #1ms to 1s, 180 t2s as in Barta et al. Last one is for constant offset term
    else:
        # print('Using the t2s given')
        t2 = t2s 
        
    if chisq_misfit is None:
        chisq_misfit = [1.02,1.025] #limits for chisq misfit of regularized compared to non-regularized. Usually 2-2.5%
    
    if len(data) == 0:
        # print('Loading data from file')
        time, data_y = np.loadtxt(fname,unpack=True)
    else:
        # print('Using data passed to function')
        time, data_y = data

    
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

