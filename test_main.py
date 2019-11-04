# to import something from a file.py in a different directory
import sys
sys.path.insert(0,'/Users/beatricejelmini/Desktop/JUNO/JUNO_codes/Reactor')
sys.path.insert(0,'/Users/beatricejelmini/Desktop/JUNO/JUNO_codes/Oscillation')
sys.path.insert(0,'/Users/beatricejelmini/Desktop/JUNO/JUNO_codes/Spectrum')
sys.path.insert(0,'/Users/beatricejelmini/Desktop/JUNO/JUNO_codes/Exp_Resolution')

import latex # to use latex font and siunitx package
from reactor import Reactor_Spectrum # import class
from oscillation import Oscillation_Prob 
from spectrum import Oscillated_Spectrum
from convolution import Convolution

import numpy as np
import matplotlib.pyplot as plt
import time


## MAIN PROGRAM ##

time_start = time.process_time_ns()

E = np.arange(1.806,10.01,0.01) # in MeV

# Creation of the reactor's spectrum
#react = Reactor_Spectrum()
# IBD threshold is 1.8 MeV  
#flux = react.flux(E,plot_this=True)
#xsec = react.cross_section(E,plot_this=True)
#react.unosc_spectrum(E,plot_this=True)

# input: sin2(theta12), deltam2_21, NO: sin2(theta13), deltam2_3l, IO: sin2(theta13), deltam2_3l
# values from nu-fit.org: JHEP01 (2019) 106, table 1, (Spaniards)
prob = Oscillation_Prob(t12=0.310, m21=7.39*10**(-5), t13_N=0.02240, m3l_N=2.525*10**(-3), t13_I=0.02263, m3l_I=-2.512*10**(-3))

#prob.eval_prob(E,0,plot_this=True)
#prob.eval_prob_jhep(E,0,plot_this=True)


# parameters from nu-fit.org (JHEP01 (2019) 106, table 1, with SK)
# input: sin2(theta12), deltam2_21, NO: sin2(theta13), deltam2_3l, IO: sin2(theta13), deltam2_3l
spectrum = Oscillated_Spectrum(t12=0.310, m21=7.39*10**(-5), t13_N=0.02240, m3l_N=2.525*10**(-3), t13_I=0.02263, m3l_I=-2.512*10**(-3))
#print(str(time.clock()-time_start))

#spectrum.unosc_spectrum(E,plot_this=True)
#spectrum.eval_prob(E,0,plot_this=True)

spectrum.osc_spectrum(E,0,plot_this=False)

a = 0.029
b = 0.008
#Evis = E - 0.8
#appo = a**2 / Evis + b**2
#sigma_Evis = np.sqrt(appo) * Evis

convol = Convolution()
conv_np_N = convol.np_conv(spectrum.norm_osc_spect_N,E,a=a,b=b,plot_this=False)
conv_np_I = convol.np_conv(spectrum.norm_osc_spect_I,E,a=a,b=b,plot_this=False)
conv_num_N = convol.numerical_conv(spectrum.norm_osc_spect_N,E,a=a,b=b,plot_this=False)
conv_num_I = convol.numerical_conv(spectrum.norm_osc_spect_I,E,a=a,b=b,plot_this=False)

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

'''sig = np.array([0.04,0.05,0.06,0.07,0.08,0.09,0.1])
for num in sig:
    ax.plot(E-0.8,convol.numerical_conv(spectrum.norm_osc_spect_N,E,sigma=num), linewidth=0.7,label='sigma = %.2f' % (num)) 
ax.legend() '''


const = (E[-1] - E[0])/len(E)
fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.plot(E-0.8,conv_np_N*const,'b',linewidth=1,label='NO')
ax1.plot(E-0.8,conv_np_I*const,'r--',linewidth=1,label='IO')
ax1.set_xlabel(r'$\text{E}_{\text{vis}}^{\text{obs}}$ [\si{MeV}]')
ax1.set_ylabel(r'N($\bar{\nu}$) [arb. unit]')
ax1.set_ylim(-0.005,0.095)
ax1.set_title(r'Antineutrino spectrum with finite energy resolution' + '\n(numpy.convolve())')
ax1.text(8.05,0.05,r'a = \SI{%.1f}{\percent}' % (a*100) + '\nb = \SI{%.1f}{\percent}' % (b*100))
ax1.legend()
ax1.grid()
fig1.savefig('osc_spectrum_w_resolution_np.pdf',format='pdf',transparent=True)

'''sig = np.array([0.04,0.05,0.06,0.07,0.08,0.09,0.1])
for num in sig:
    ax1.plot(E-0.8,convol.np_conv(spectrum.norm_osc_spect_N,E,sigma=num) * const, linewidth=0.7,label='sigma = %.2f' % (num)) 
ax1.legend() '''




'''fig = plt.figure()
#fig.suptitle(r'Antineutrino spectrum')
ax = fig.add_subplot(111)
ax.legend()
ax.grid()
ax.set_xlabel(r'$\text{E}_{\nu}$ [\si{MeV}]')
ax.set_ylabel(r'N($\bar{\nu}$) [arb. unit]')
ax.set_ylim(-0.005,0.095)
ax.set_title(r'Antineutrino spectrum')
ax.plot(E,spectrum.norm_osc_spect_N,'b',linewidth=1,label='NO')
ax.plot(E,spectrum.norm_osc_spect_I,'r--',linewidth=1,label='IO')
ax.legend()
fig.savefig('Spectrum/osc_spectrum_N_I.pdf',format='pdf',transparent=True)  '''

plt.ion()
plt.show()




