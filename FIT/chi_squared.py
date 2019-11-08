import matplotlib.pyplot as plt
import numpy as np
import math
from scipy import integrate, stats, optimize
import time
import iminuit


class Chi_Squared():
    # no __init__
    def __init__ (self,f,data,bins,sigma=''):

        self.f = f
        self.data = data
        self.bins = bins

        if sigma == '':
            self.sigma_data = np.sqrt(data) # poisson error
        else:
            self.sigma_data = sigma

        np.place(self.sigma_data,self.sigma_data<1.,1.) 
        

    '''def chi_squared_stat(self,params,f,data,bins,sigma=''): # params are parameters of f
        
        x_fit = f(bins,*params)

        if sigma == '':
            sigma_data = np.sqrt(data) # poisson error
        else:
            sigma_data = sigma

        np.place(sigma_data,sigma_data<1.,1.) 

        appo = (x_fit - data)/ sigma_data
        appo = np.power(appo,2)
        
        return appo.sum() '''

    def chi_squared_stat(self,params): 

        x_fit = self.f(self.bins,*params)

        appo = (x_fit - self.data)/ self.sigma_data
        appo = np.power(appo,2)
        
        return appo.sum() 

    def chi_squared_stat_minuit(self,mu,sigma): 

        x_fit = self.f(self.bins,mu,sigma)

        appo = (x_fit - self.data)/ self.sigma_data
        appo = np.power(appo,2)
        
        return appo.sum() 






