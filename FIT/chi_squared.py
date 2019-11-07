import matplotlib.pyplot as plt
import numpy as np
import math
from scipy import integrate, stats, optimize
import time
import iminuit


class Chi_Squared():
    # no __init__

    def chi_squared_stat(params,f,data,bins,sigma=''): # params are parameters of f
        
        x_fit = f(bins,*params)

        if sigma == '':
            sigma_data = np.sqrt(data) # poisson error
        else:
            sigma_data = sigma

        np.place(sigma_data,sigma_data<1.,1.) 

        appo = (x_fit - data)/ sigma_data
        appo = np.power(appo,2)
        
        return appo.sum()
