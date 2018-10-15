import numpy as np
import lmfit
import pylab as plt
from IPython.core.debugger import Tracer; debug_here = Tracer()

def do_mc_error(
                resid_func,
                resid_func_kws,
                fit_params,
                n_iter,
                fit_data,
                noise_sd=None,
                fitter_kws={},
                randomize_init_params=False,
                noise_func=None,
                noise_func_args=None,
                ):

    if 'exp_data' not in resid_func_kws:
        raise Exception('One of the arguments of the residual function has to be "exp_data"!')

    param_vals = {}
    for pn in fit_params:
        param_vals[pn] = np.zeros(n_iter)

    #Sanity checks if randomizing initial params each iteration
    if randomize_init_params:
        if fitter_kws.get('method','') == 'differential_evolution':
            Print('do_mc_error: randomize_init_params=True but using differential_evolution. Starting params will be randomized by the solver anyway')
        else:
            for pn in fit_params:
                if (fit_params[pn].vary) and (fit_params[pn].expr is None) and (np.isinf(fit_params[pn].min) or np.isinf(fit_params[pn].max)):
                    raise Exception('do_mc_error: randomize_init_params=True but param {0:s} has infinity limits'.format(pn))
    # fit_data = resid_func_kws['exp_data'] - resid_func(fit_params,**resid_func_kws)

    #If no noise function is given, use the noise standard deviation
    if noise_func is None:
        assert noise_sd is not None, 'If no noise_func is passed, the stdev of the residual noise must be (noise_sd).'
        noise_func = lambda size,stdev : stdev*np.random.standard_normal(size=size)
        noise_func_args = {'size' : fit_data.shape, 'stdev' : noise_sd}

    num_at_lims = 0
    for i in range(n_iter):
        print('{0:d}'.format(i),end=',',flush=True)


        if randomize_init_params:
            for pn in fit_params:
                if fit_params[pn].vary and fit_params[pn].expr is None:
                    span = fit_params[pn].max - fit_params[pn].min
                    fit_params[pn].value = fit_params[pn].min + 0.05*span + np.random.rand()*0.9*span

        resid_func_kws['exp_data'] = fit_data + noise_func(**noise_func_args)

        # plt.plot(resid_func_kws['exp_data'].flatten())
        # plt.plot(fit_data.flatten())
        # plt.show()
        results = lmfit.minimize(resid_func,fit_params,kws=resid_func_kws,**fitter_kws)

        #Check to make sure results aren't pegged at one of the parameter limits. If this happens a lot,
        #the parameters probably aren't very well constrained. Increase the limits or change the model
        for pn in param_vals:
            param_vals[pn][i] = results.params[pn].value
            if not results.params[pn].vary:
                continue
            llim,ulim = results.params[pn].min,results.params[pn].max
            if np.isfinite(llim) and (results.params[pn].value-llim)<1e-3:
                num_at_lims+=1
                # print('\ndo_mc_error: Lower limit warning! parameter {0:s} has value of {1:f} and min {2:f}. Consider decreasing the minimum.\n'.format(pn,results.params[pn].value,llim))
            if np.isfinite(ulim) and (ulim-results.params[pn].value)<1e-3:
                num_at_lims+=1
                # print('\ndo_mc_error: Upper limit warning! parameter {0:s} has value of {1:f} and max {2:f}. Consider increasing the maximum.\n'.format(pn,results.params[pn].value,ulim))

    print('\nAt lims {0:d} times out of {1:d} iterations'.format(num_at_lims,n_iter))
    for pn in param_vals:
        fit_params.pop(pn)
        fit_params.add(pn,value=np.average(param_vals[pn]),vary=False)
        fit_params[pn].stderr = np.std(param_vals[pn])
        # # fit_params[pn].min = -np.inf
        # # fit_params[pn].max = np.inf
        # fit_params[pn].value = np.average(param_vals[pn])
        # fit_params[pn].stderr = np.std(param_vals[pn])
        # # fit_params[pn].vary = False


    return fit_params