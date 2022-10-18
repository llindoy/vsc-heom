import numpy as np
from matplotlib import pyplot as plt
import scipy.optimize
import glob

def monoExp(x, m, t):
    return -np.abs(t) * x + m

def get_rate(fname, L, skip = 1, color=None):
    try:
        data = np.genfromtxt(fname, skip_header=skip)
        d2 = data[:, 1]
        t = data[:, 0]

        p0 = (0.5, 0.1) # start with values near those we expect
        N = d2.shape[0]
        params, cv = scipy.optimize.curve_fit(monoExp, t[(N*4)//5:], np.log(d2[(N*4)//5:]-0.5), p0)
        m, k = params

        beta = 1052.8
        Eb = 0.3/ 27.2114

        assert(d2.shape[0] > 100)

        diff = np.diff(d2)/(t[1]-t[0])
        kappa = -diff/(d2[1:] + (d2[1:]-1))
        N = d2.shape[0]//5
        kappas = np.convolve(kappa, np.ones(N)/N, mode='valid')
        print(L, np.abs(k), kappa[-1], kappas[-1])
        if not color is None:
            plt.semilogx(t[1:], kappa, label="K_"+str(L))
            plt.semilogx(t[:-N], kappas, '--', label="aveK_"+str(L))
        else:
            plt.semilogx(t[1:], kappa, label="K_"+str(L), color=color)
            plt.semilogx(t[:-N], kappas, '--', label="aveK_"+str(L), color=color)    
    except: 
        print(fname)
        return None

file = "examples/rate_theory_calculations/heom_cavity/run_1193.18181818_15.out"

print("filelabel exp_fit instantaneous averaged")
get_rate(file, "1193", color='k')#, color='k')#, color='k')

plt.legend()
plt.show()
