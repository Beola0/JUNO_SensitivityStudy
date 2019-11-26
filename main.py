import os
cwd = os.getcwd()

import sys
sys.path.insert(0,cwd+'/Reactor')
sys.path.insert(0,cwd+'/Oscillation')
sys.path.insert(0,cwd+'/Spectrum')
sys.path.insert(0,cwd+'/Exp_Resolution')

import latex 
from reactor import Reactor_Spectrum 
from oscillation import Oscillation_Prob 
from spectrum import Oscillated_Spectrum
from convolution import Convolution

import numpy as np
import matplotlib.pyplot as plt
import time

##### MAIN PROGRAM #####

#time_start = time.process_time_ns()

E = np.arange(1.806,10.01,0.001) # in MeV



### creation of reactor's spectrum, class Reactor_Spectrum
### set plot_this=True to plot the reactor's flux, the IBD cross-section and the reactor's spectrum, respectively

react = Reactor_Spectrum()
flux = react.flux(E,plot_this=False) 
xsec = react.cross_section(E,plot_this=False)
reactor_spectrum = react.unosc_spectrum(E,plot_this=False)



### evaluation of the survival probability, class Oscillation_Prob
### set plot_this=True to plot the survibal probability
### the survival probability is plotted both as function of L/E and of E
### NO and N refer to the Normal Ordering; IO and I refer to the Inverted Ordering

### input: sin^2(theta_12), deltam^2_21, NO: sin^2(theta_13), deltam^2_3l, IO: sin^2(theta_13), deltam^2_3l
### values from JHEP01 (2019) 106, table 1 (http://www.nu-fit.org/?q=node/8) (see also for the notation)
prob = Oscillation_Prob(t12=0.310, m21=7.39*10**(-5), t13_N=0.02240, m3l_N=2.525*10**(-3), t13_I=0.02263, m3l_I=-2.512*10**(-3))
prob_N, prob_I = prob.eval_prob(E,0,plot_this=False) # 1 for NO, -1 for IO, 0 for both (plotting)



### creation of the oscillated spectrum, class Oscillated_Spectrum
### set plot_this=True to plot the oscillated spectrum
### set plot_un=True to plot also the unoscillated spectrum

### parameters from nu-fit.org as above (JHEP01 (2019) 106, table 1, with SK)
### input: sin^2(theta_12), deltam^2_21, NO: sin^2(theta_13), deltam^2_3l, IO: sin^2(theta_13), deltam^2_3l
spectrum = Oscillated_Spectrum(t12=0.310, m21=7.39*10**(-5), t13_N=0.02240, m3l_N=2.525*10**(-3), t13_I=0.02263, m3l_I=-2.512*10**(-3))
spectrum_N, spectrum_I = spectrum.osc_spectrum(E,-1,plot_this=False,plot_un=False) # 1 for NO, -1 for IO, 0 for both (plotting)



### oscillated spectrum with experimental resolution via numerical convolution, class Convolution used in the class Spectrum
### set plot_this=True to plot the numerical convolution
### set plot_start=True to plot the starting spectrum
### for further reference see https://arxiv.org/abs/1210.8141

a = 0.029 # stochastic term
b = 0.008 # constant term 
resol_spect_N, resol_spect_I = spectrum.resol_spectrum (E,a,b,0,plot_this=False) # 1 for NO, -1 for IO, 0 for both (plotting)

#elapsed_time = time.process_time_ns() - time_start
#print('elapsed time: '+str(elapsed_time*10**(-6))+' ms')

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(E-0.8,resol_spect_N,'b',linewidth=1,label='NO')
ax.plot(E-0.8,resol_spect_I,'r--',linewidth=1,label='IO')
ax.set_xlabel(r'$\text{E}_{\text{vis}}$ [\si{MeV}]')
ax.set_ylabel(r'N($\bar{\nu}$) [arb. unit]')
ax.set_ylim(-0.005,0.095)
ax.set_title(r'Antineutrino spectrum' + '\nwith finite energy resolution')
ax.text(8.05,0.05,r'a = \SI{%.1f}{\percent}' % (a*100) + '\nb = \SI{%.1f}{\percent}' % (b*100))
ax.legend()
ax.grid()
fig.savefig('oscillated_spectrum.pdf',format='pdf',transparent=True) 
print('Saved plot in oscillated_spectrum.pdf')

plt.ion()
plt.show()







