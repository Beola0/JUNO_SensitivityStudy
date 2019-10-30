import sys
sys.path.insert(0,'/Users/beatricejelmini/Desktop/juno_scripts/Reactor')
sys.path.insert(0,'/Users/beatricejelmini/Desktop/juno_scripts/Oscillation')
sys.path.insert(0,'/Users/beatricejelmini/Desktop/juno_scripts/Spectrum')

import latex # to use latex font and siunitx package
from reactor import Reactor_Spectrum # import class
from oscillation import Oscillation_Prob 
from spectrum import Oscillated_Spectrum
from convolution import Convolution

import math
import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
#import time


def Gaussian (x,sigma):
    const = math.sqrt(2*math.pi)
    appo = 1/const/sigma * np.exp(- np.power(x,2)/2./(sigma**2))
    return appo



#time_start = time.process_time_ns()

E = np.arange(1.806,10.01,0.01) # in MeV
Evis = E - 0.8

spectrum = Oscillated_Spectrum(t12=0.310, m21=7.39*10**(-5), t13_N=0.02240, m3l_N=2.525*10**(-3), t13_I=0.02263, m3l_I=-2.512*10**(-3))

spectrum.osc_spectrum(E,0,plot_this=False)
zero = np.zeros(len(E))

fig = plt.figure()
ax = fig.add_subplot(111,projection = '3d')
ax.plot(E,zero,spectrum.norm_osc_spect_N,'b',linewidth=1.,label='NO')


a = 0.029
b = 0.008

Gs = np.empty((len(E),len(E)))
n_=0
for E0 in Evis:
    appo = math.pow(a,2) / E0 + math.pow(b,2)
    sigma = math.sqrt(appo) * E0
    g = Gaussian(Evis - E0,sigma)
    g_appo = g * spectrum.norm_osc_spect_N[n_]
    Gs[n_,] = g_appo
    n_ += 1

g_sum = np.empty(len(E))
g_sum_simps = np.empty(len(E))
for n0 in np.arange(0,len(E)):
    g_sum[n0] = Gs[:,n0].sum()
    g_sum_simps[n0] = integrate.simps(Gs[:,n0],E)

const = (Evis[-1] - Evis[0]) / len(Evis)
#ax.plot(zero,Evis,g_sum*const,'g',linewidth=1.,label='convolution')
ax.plot(zero,Evis,g_sum_simps,'r',linewidth=1.,label='convolution')

n = np.array([20,70,120,170,220,270,320,370,420,470,520,570,620,670,720]) # al MeV e al mezzo MeV
for n0 in n:

    x = E[n0] # fixed Enu
    xvis = x - 0.8
    appo = math.pow(a,2) / xvis + math.pow(b,2)
    sigma = math.sqrt(appo) * xvis
    #sigma = 1.
    g = Gaussian (Evis-xvis,sigma)
    g_appo = g * spectrum.norm_osc_spect_N[n0]
    x_appo = np.empty(len(E))
    x_appo.fill(x)
    ax.plot(x_appo,Evis,g_appo,linewidth=1.,label=r'sigma = %.4f $10^{-2}$' % (sigma*100))

#time_el = time.process_time_ns()
#print('elapsed time: ' + str(time_el*10**(-6)) + ' ms')

fig.suptitle('Convolution via a graphic method')
ax.set_xlabel(r'$\text{E}_{\nu}$ [\si{MeV}]')
ax.set_xlim(0.,11.)
ax.set_ylabel(r'$\text{E}_{\text{vis}}^{\text{obs}}$ [\si{MeV}]')
ax.set_ylim(0.,11.)
ax.set_zlabel(r'N [arb. unit]')
#ax.set_zlim(-0.005,0.095)
ax.legend(prop={'size': 8})
fig.savefig('3d_convolution.pdf',format='pdf')

fig1 = plt.figure()
fig1.suptitle(r'3d convolution')
ax1 = fig1.add_subplot(111)
#ax1.plot(Evis,spectrum.norm_osc_spect_N,'b',linewidth=1.,label='NO')
#ax1.plot(Evis,g_sum*const,'g',linewidth=1.,label='convolution - num')
ax1.plot(Evis,g_sum_simps,'b',linewidth=1.,label='convolution - simps')
ax1.set_xlabel(r'$\text{E}_{\text{vis}}^{\text{obs}}$ [\si{MeV}]')
ax1.set_ylabel(r'N [arb. unit]')
ax1.set_ylim(-0.005,0.095)
ax1.grid()
ax1.legend()
fig1.savefig('conv_3d.pdf',format='pdf')

plt.ion()
plt.show()
