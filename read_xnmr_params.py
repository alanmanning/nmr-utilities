"""
This script reads the params from the Xnmr param file and returns
them as a dictionary. Things that can be converted to numbers will be
returned as a float, otherwise they're left as strings. Arrayed parameters
will be returned as lists
-APM 23 Feb 2017

"""

from IPython.core.debugger import Tracer; debug_here = Tracer()

def read_xnmr_params(fname,debug=False):
    f=open(fname+'/params','r')
    lines = f.readlines()
    f.close()

    params = dict()
    params['seqfil'] = lines[0].strip()
    try:
        for i in range(1,len(lines)):
            x = lines[i].find('=')
            if x > 0:
                pname = lines[i][0:x-1]
                pstring = lines[i][x+1:].strip()
                isfloat = True
                try:
                    pval = float(pstring)
                except:
                    isfloat = False
                    pval = pstring
                if not (pname in params): #new param, so lets add it to the dictionary
                    params[pname] = pval
                else:
                    if not isinstance(params[pname],list):
                        params[pname] = [params[pname]]
                    params[pname] += [pval]
    except:
        print('Error reading param file. Issue with line: ' + lines[i])
        debug_here()
        return 0
    if debug:
        return params,lines
    else:
        return params
