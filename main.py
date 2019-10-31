import sys
sys.path.insert(0,'/Users/beatricejelmini/Desktop/JUNO/JUNO_codes/Reactor')
sys.path.insert(0,'/Users/beatricejelmini/Desktop/JUNO/JUNO_codes/Oscillation')
sys.path.insert(0,'/Users/beatricejelmini/Desktop/JUNO/JUNO_codes/Spectrum')
sys.path.insert(0,'/Users/beatricejelmini/Desktop/JUNO/JUNO_codes/Exp_Resolution')

import latex 
from reactor import Reactor_Spectrum 
from oscillation import Oscillation_Prob 
from spectrum import Oscillated_Spectrum
from convolution import Convolution

import numpy as np
import matplotlib.pyplot as plt
import time

## MAIN PROGRAM ##

time_start = time.process_time_ns()

E = np.arange(1.806,10.01,0.001) # in MeV


# creation of reactor's spectrum

react = Reactor_Spectrum()
flux = react.flux(E,plot_this=False)
xsec = react.cross_section(E,plot_this=False)
reactor_spectrum = react.unosc_spectrum(E,plot_this=False)


# evaluation of the survival probability

# input: sin2(theta12), deltam2_21, NO: sin2(theta13), deltam2_3l, IO: sin2(theta13), deltam2_3l
# values from nu-fit.org: JHEP01 (2019) 106, table 1, (Spaniards)
prob = Oscillation_Prob(t12=0.310, m21=7.39*10**(-5), t13_N=0.02240, m3l_N=2.525*10**(-3), t13_I=0.02263, m3l_I=-2.512*10**(-3))
prob_N, prob_I = prob.eval_prob(E,0,plot_this=False)


# oscillated spectrum

# parameters from nu-fit.org (JHEP01 (2019) 106, table 1, with SK)
# input: sin2(theta12), deltam2_21, NO: sin2(theta13), deltam2_3l, IO: sin2(theta13), deltam2_3l
spectrum = Oscillated_Spectrum(t12=0.310, m21=7.39*10**(-5), t13_N=0.02240, m3l_N=2.525*10**(-3), t13_I=0.02263, m3l_I=-2.512*10**(-3))
spectrum_N, spectrum_I = spectrum.osc_spectrum(E,0,plot_this=True)

# experimental resolution
a = 0.029
b = 0.008
convol = Convolution()
conv_num_N = convol.numerical_conv(spectrum_N,E,a=a,b=b,plot_this=False)
conv_num_I = convol.numerical_conv(spectrum_I,E,a=a,b=b,plot_this=False)

#elapsed_time = time.process_time_ns() - time_start
#print('elapsed time: '+str(elapsed_time*10**(-6))+' ms')

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(E-0.8,conv_num_N,'b',linewidth=1,label='NO')
ax.plot(E-0.8,conv_num_I,'r--',linewidth=1,label='IO')
ax.set_xlabel(r'$\text{E}_{\text{vis}}^{\text{obs}}$ [\si{MeV}]')
ax.set_ylabel(r'N($\bar{\nu}$) [arb. unit]')
ax.set_ylim(-0.005,0.095)
ax.set_title(r'Antineutrino spectrum with finite energy resolution' + '\n(numerical convolution)')
ax.text(8.05,0.05,r'a = \SI{%.1f}{\percent}' % (a*100) + '\nb = \SI{%.1f}{\percent}' % (b*100))
ax.legend()
ax.grid()
fig.savefig('osc_spectrum_w_resolution_num.pdf',format='pdf',transparent=True) 

plt.ion()
plt.show()







