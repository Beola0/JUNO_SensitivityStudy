import matplotlib.pyplot as plt
import numpy as np
import math, scipy
from scipy import integrate, stats, optimize
import time
import iminuit
from test_statistic import Test_Statistic


class Chi_Squared(Test_Statistic):
    # no __init__
    def __init__ (self,f,x_data,bins,sigma=''):

        self.f = f
        self.x_data = x_data
        self.bins = bins

        if sigma == '':
            self.sigma_data = np.sqrt(x_data) # poisson error
        else:
            self.sigma_data = sigma
        
        # check if there are sigma = 0 
        np.place(self.sigma_data,self.sigma_data<1.,3.) 

        self.V = np.ones((len(bins),len(bins)))
        self.V_inv = np.ones((len(bins),len(bins)))

    def set_invert_mat(self,V):
        self.V = V
        self.V_inv = np.linalg.inv(self.V)

        M = len(V[0,:])        
        U,s,Vh = scipy.linalg.svd(self.V)
        S = np.zeros((M,M))
        S_inv = np.zeros((M,M))
        for i in np.arange(0,M): 
            S[i,i] = s[i] 
            S_inv[i,i] = 1./s[i]

        appo = np.dot(Vh,S_inv)
        V_inv = np.dot(appo,U)

        return U,s,Vh, V_inv
                  

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

    def chi_squared_V_(self,params): # to use with scipy.optimize.minimize()

        x_fit = self.f(self.bins,*params)
        X = self.x_data - x_fit

        M = len(self.bins)
        appo = np.zeros(M)
        for j in np.arange(0,M):
            appo[j] = (self.V_inv[j,:] * X).sum()
        appo1 = X * appo

        #appo = (x_fit - self.data)/ self.sigma_data
        #appo = np.power(appo,2)
        
        return appo1.sum() 

    def chi_squared_stat_minuit(self,mu,sigma): # to use with Minuit

        x_fit = self.f(self.bins,mu,sigma)

        appo = (x_fit - self.x_data)/ self.sigma_data
        appo = np.power(appo,2)
        
        return appo.sum() 

    def chi_squared_V(self,mu,sigma):

        x_fit = self.f(self.bins,mu,sigma)
        X = self.x_data - x_fit

        M = len(self.bins)
        appo = np.zeros(M)
        for j in np.arange(0,M):
            appo[j] = (self.V_inv[j,:] * X).sum()
        appo1 = X * appo

        #appo = np.dot(self.V_inv,X)
        #appo1 = np.dot(X,appo)

        return appo1.sum()






