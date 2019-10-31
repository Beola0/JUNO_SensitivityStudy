import sys
sys.path.insert(0,'/Users/beatricejelmini/Desktop/JUNO/JUNO_codes')

import math
import numpy as np
import matplotlib.pyplot as plt
import latex


def Gaussian (x,sigma):
    const = math.sqrt(2*math.pi)
    appo = 1/const/sigma * np.exp(- np.power(x,2)/2./(sigma**2))
    return appo


a = 0.029
b = 0.008

E_nu = np.arange(1.806,10.01,0.01) # in MeV
Evis = E_nu - 0.8

appo = math.pow(a,2) / Evis + math.pow(b,2)
sigma_Evis = np.sqrt(appo) * Evis

#g = Gaussian(Evis-5,0.05)
#g_var = Gaussian(Evis-5,sigma_Evis)

fig = plt.figure()
ax = fig.add_subplot(111)

ax.set_xlabel(r'$\text{E}_{\text{vis}}$ [\si{MeV}]')
ax.set_xlim(0,10)
ax.set_ylabel(r'$\sigma_\text{E}$ [\si{MeV}]',color='b')
ax.plot(Evis,sigma_Evis,'b',linewidth=1)
ax.tick_params(axis='y', labelcolor='b')

ax.text(4.05,0.105,r'a = \SI{%.1f}{\percent}' %(a*100) + '\nb = \SI{%.1f}{\percent}' % (b*100))
ax.set_title('Energy resolution')
ax.grid()

ax1 = ax.twinx() # new axes with same x
ax1.set_ylabel(r'$\sigma_\text{E}/E~(\si{\percent})$',color='r')
ax1.plot(Evis,sigma_Evis/Evis*100,'r',linewidth=1)
ax1.tick_params(axis='y', labelcolor='r')

fig.tight_layout()

fig.savefig('resolution.pdf',format='pdf',transparent=True)

plt.ion()
plt.show()


