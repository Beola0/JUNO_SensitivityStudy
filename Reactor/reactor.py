#import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import simps


def reactor_exp(x,a,b,c): # input x is an array
    res = np.exp(a-b*x-c*(x**2))
    return res


class Reactor_Spectrum:

    def __init__(self):
        # fixed (for now) fission fractions
        self.f235u = 0.58
        self.f239pu = 0.30
        self.f238u = 0.07
        self.f241pu = 0.05

    ### this method evaluates the antineutrino flux
    ### for further reference: http://inspirehep.net/record/25814/ and https://arxiv.org/abs/0807.3203, eq. (2)
    def flux(self,E,plot_this=False): 
        u235 = reactor_exp(E,0.870,0.160,0.091)
        pu239 = reactor_exp(E,0.896,0.239,0.0981)
        u238 = reactor_exp(E,0.976,0.162,0.0790)
        pu241 = reactor_exp(E,0.793,0.080,0.1085)

        self.tot_flux = self.f235u * u235 + self.f239pu * pu239 + self.f238u * u238 + self.f241pu * pu241

        if plot_this:

            fig = plt.figure()
            #fig.suptitle('Reactor antineutrino flux')
            ax = fig.add_subplot(111)

            ax.plot(E,self.f235u * u235,'b',linewidth=1,label=r'$^{235}$U')
            ax.plot(E,self.f239pu * pu239,'r',linewidth=1,label=r'$^{239}$Pu')
            ax.plot(E,self.f238u * u238,'g',linewidth=1,label=r'$^{238}$U')
            ax.plot(E,self.f241pu * pu241,'y',linewidth=1,label=r'$^{241}$Pu')
            ax.plot(E,self.tot_flux,'k',linewidth=1,label='total')

            ax.legend()
            ax.grid()
            ax.set_xlabel(r'$\text{E}_{\nu}$ [$\si{MeV}$]')
            ax.set_ylabel(r'$\Phi_{\nu}$')
            ax.set_title(r'Reactor antineutrino flux')

            plt.savefig('Reactor/flux.pdf',format='pdf',transparent=True)
            print('\nThe plot has been saved in Reactor/flux.pdf\n')

        return self.tot_flux

    ### this method evaluates the IBD cross section
    ### for further reference: Strumia, Vissani, https://arxiv.org/abs/astro-ph/0302055, eq. (25)
    def cross_section(self,E,plot_this=False): # E is the neutrino's energy
        alpha = -0.07056
        beta = 0.02018
        gamma = -0.001953
        Delta = 1.293 # MeV, mass(n)-mass(p)
        m_e = 0.511 # MeV
        const = 10**(-43) # cm^2

        E_e = np.subtract(E,Delta)  # electron's energy

        appo = np.power(E_e,2) - m_e**2
        p_e = np.sqrt(appo) # electron's momentum

        appo_exp = alpha + beta * np.log(E) + gamma * np.power(np.log(E),3)
        E_exp = np.power(E,appo_exp)

        self.x_sec = const * p_e * E_e * E_exp

        if plot_this:

            fig = plt.figure()
            #fig.suptitle('IBD cross section')
            ax = fig.add_subplot(111)
            
            ax.plot(E,self.x_sec,'k',linewidth=1,label='IBD cross section')
            ax.grid()
            ax.set_xlabel(r'$\text{E}_{\nu}$ [\si{MeV}]')
            ax.set_ylabel(r'$\sigma_{\text{IBD}}$ [\si{\centi\meter\squared}]')
            ax.set_title(r'IBD cross section')

            plt.savefig('Reactor/cross_section.pdf',format='pdf',transparent=True)
            print('\nThe plot has been saved in Reactor/cross_section.pdf\n')

        return self.x_sec

    # this method combines the flux and the cross section to obtain the spectrum
    def unosc_spectrum(self,E,plot_this=False):
        self.flux(E,plot_this=False)
        self.cross_section(E,plot_this=False)

        self.spectrum_un = self.tot_flux * self.x_sec
        integral = simps(self.spectrum_un,E)
        self.norm_spectrum_un = self.spectrum_un/integral

        if plot_this:

            fig = plt.figure()
            #fig.suptitle('Unoscillated reactor spectrum')
            ax = fig.add_subplot(111)

            #ax.plot(E,self.spectrum_un,'b',linewidth=1,label='spectrum')
            ax.plot(E,self.norm_spectrum_un,'k',linewidth=1,label='spectrum')
            ax.grid()
            ax.set_xlabel(r'$\text{E}_{\nu}$ [\si{MeV}]')
            ax.set_ylabel(r'N($\nu$) [arb. unit]')
            ax.set_title(r'Unoscillated reactor spectrum')

            plt.savefig('Reactor/unoscillated_spectrum.pdf',format='pdf',transparent=True)
            print('\nThe plot has been saved in Reactor/unoscillated_spectrum.pdf\n')

        return self.norm_spectrum_un

        





