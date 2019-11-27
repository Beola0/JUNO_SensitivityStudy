import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import integrate
from scipy import signal


class Convolution ():

    def Gaussian (self,x,sigma):
        const = math.sqrt(2*math.pi)
        appo = 1/const/sigma * np.exp(- np.power(x,2)/2./(sigma**2))
        return appo

    ### convolution of a function f with a gaussian in a given range E

    ### using np.convolve()
    def np_conv(self,f,E,sigma='',a='',b='',plot_this=False,plot_start=False): 

        Evis = E - 0.8

        if sigma != '':
            
            g = self.Gaussian(Evis-np.mean(Evis),sigma)
            self.conv_np = np.convolve(f,g,mode='same')
            #self.conv_np = signal.fftconvolve(f,g,mode='same') # same results as np.convolve()

        if a != '' or b != '':
        
            rad = a**2 / Evis + b**2
            sigma_Evis = np.sqrt(rad) * Evis

            g = self.Gaussian(Evis-np.mean(Evis),sigma_Evis)
            self.conv_np = np.convolve(f,g,mode='same')
            #self.conv_np = signal.fftconvolve(f,g,mode='same')

        if plot_this:

            fig1 = plt.figure()
            ax1 = fig1.add_subplot(111)
            ax1.plot(Evis,self.conv_np,'b',linewidth=1,label='Convolved spectrum')
            ax1.set_xlabel(r'$\text{E}_{\text{vis}}$ [\si{MeV}]')
            ax1.set_ylabel(r'N($\bar{\nu}$) [arb. unit]')
            ax1.set_title(r'Convolution (numpy) with a Gaussian' + '\nwith variable width')
            ax1.grid()
            ax1.legend()

        if plot_start:

            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(Evis,f,'b',linewidth=1,label='Unconvolved spectrum')
            ax.set_xlabel(r'$\text{E}_{\text{dep}}$ [\si{MeV}]')
            ax.set_ylabel(r'N($\bar{\nu}$) [arb. unit]')
            ax.set_title(r'Starting spectrum')
            ax.grid()
            ax.legend()

        return self.conv_np


    ### using a numerical method: numerical convolution (over Evis)
    ### need (sigma) OR (a or b)
    def numerical_conv(self,f,E,sigma='',a='',b='',plot_this=False,plot_start=False):

        Evis = E - 0.8

        if sigma != '' and (a != '' or b != ''):
            sigma = ''
            print("Both 'sigma' and 'a-b' have been given. 'sigma' will be ignored, 'a-b' will be used.")

        # gaussian with fixed or given width
        if sigma != '':

            self.conv_num = np.zeros(len(Evis))
            n=0
            for E0 in Evis:
                appo = self.Gaussian(Evis-E0,sigma)
                prod = appo * f
                self.conv_num[n] = integrate.simps(prod,Evis) 
                n += 1

        # gaussian with variable width, set by a (stochastic term) and b (constant term)
        if a != '' or b != '':
        
            rad = a**2 / Evis + b**2
            sigma_Evis = np.sqrt(rad) * Evis

            self.conv_num = np.zeros(len(Evis))
            n=0
            for E0 in Evis:
                appo = self.Gaussian(Evis-E0,sigma_Evis)
                prod = appo * f
                self.conv_num[n] = integrate.simps(prod,Evis) 
                n += 1

        if plot_this:

            fig1 = plt.figure()
            ax1 = fig1.add_subplot(111)
            ax1.plot(Evis,self.conv_num,'b',linewidth=1,label='Convolved spectrum')
            ax1.set_xlabel(r'$\text{E}_{\text{vis}}$ [\si{MeV}]')
            ax1.set_ylabel(r'N($\bar{\nu}$) [arb. unit]')
            ax1.set_ylim(-0.005,0.095)
            ax1.set_title(r'Numerical convolution with a Gaussian' + '\nwith variable width')
            ax1.grid()
            ax1.legend()

        if plot_start:

            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(Evis,f,'b',linewidth=1,label='Unconvolved spectrum')
            ax.set_xlabel(r'$\text{E}_{\text{dep}}$ [\si{MeV}]')
            ax.set_ylabel(r'N($\bar{\nu}$) [arb. unit]')
            ax.set_ylim(-0.005,0.095)
            ax.set_title(r'Starting spectrum')
            ax.grid()
            ax.legend()

        return self.conv_num        








