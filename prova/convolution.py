import numpy as np
import matplotlib.pyplot as plt
import math


class Convolution ():

    #def __init__(self): # not necessary


    def Gaussian (self,x,sigma):
        const = math.sqrt(2*math.pi)
        appo = 1/const/sigma * np.exp(- np.power(x,2)/2./(sigma**2))
        return appo

    # convolution of f with a gaussian in the range E

    # using np.convolve()
    def np_conv(self,f,E,sigma,plot_this=False): 
        g = self.Gaussian(E-np.mean(E),sigma)
        self.conv_np = np.convolve(f,g,mode='same')

        if plot_this:

            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(E,f,'b',label='spectrum')
            ax.plot(E,g,'r',label='gaussian')
            ax.grid()
            ax.legend()

            fig1 = plt.figure()
            ax1 = fig1.add_subplot(111)
            ax1.plot(E,self.conv_np,'b',label='convolution (numpy)')
            ax1.grid()
            ax1.legend()

        return self.conv_np

    # using a numerical method
    def numerical_conv(self,f,E,sigma,plot_this=False):
        g = self.Gaussian(E-np.mean(E),sigma)
        self.conv_num = np.zeros(len(E))
        n=0
        for E0 in E:
            appo = self.Gaussian(E-E0,sigma)
            prod = appo * f
            self.conv_num[n] = prod.sum() # da sostituire con integrale simpson
            n += 1

        if plot_this:

            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(E,f,'b',label='spectrum')
            ax.plot(E,g,'r',label='gaussian')
            ax.grid()
            ax.legend()

            fig1 = plt.figure()
            ax1 = fig1.add_subplot(111)
            ax1.plot(E,self.conv_num,'b',label='convolution (numerical)')
            ax1.grid()
            ax1.legend()

        return self.conv_num        












