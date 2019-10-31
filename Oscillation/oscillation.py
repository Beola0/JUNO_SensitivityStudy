import matplotlib.pyplot as plt
import numpy as np
import math


class Oscillation_Prob: # survival probability

    def __init__(self,t12,m21,t13_N,m3l_N,t13_I,m3l_I):
        self.baseline = 52.5 # [km], fixed for now

        self.sin2_12 = t12
        self.deltam_21 = m21 # [eV^2]

        # Normal Ordering
        self.sin2_13_N = t13_N
        self.deltam_3l_N = m3l_N # [eV^2]

        # Inverted Ordering
        self.sin2_13_I = t13_I
        self.deltam_3l_I = m3l_I # [eV^2]


    # this function evaluates the oscillatory term
    def sin2 (self,x,dm2): 
        appo = 1.27 * dm2 * x  # x in [m/MeV]
        return np.power(np.sin(appo),2)


    # this function evalutates the survival prob as function of L/E and E 
    def eval_prob(self,E,ordering,plot_this=False): # formula from PRD 78, 2008, Chinese

        if (ordering < -1) or (ordering > 1):
            print('Error: use 1 for NO, -1 for IO, 0 for both')
            return -1

        x = np.arange(0.01,500000,1.) # [m/MeV]
        x_E = self.baseline*1000 / E # [m/MeV]

        # Normal Ordering
        A_N = math.pow(1-self.sin2_13_N,2) * 4. * self.sin2_12 * (1 - self.sin2_12)
        B_N = (1 - self.sin2_12) * 4 * self.sin2_13_N * (1 - self.sin2_13_N)
        C_N = self.sin2_12 * 4 * self.sin2_13_N * (1 - self.sin2_13_N)

        self.prob_N = 1. - A_N * self.sin2(x,self.deltam_21) - B_N * self.sin2(x,self.deltam_3l_N) - C_N * self.sin2(x,self.deltam_3l_N-self.deltam_21)
        self.prob_E_N = 1. - A_N * self.sin2(x_E,self.deltam_21) - B_N * self.sin2(x_E,self.deltam_3l_N) - C_N * self.sin2(x_E,self.deltam_3l_N-self.deltam_21)

        # Inverted Ordering
        A_I = math.pow(1-self.sin2_13_I,2) * 4. * self.sin2_12 * (1 - self.sin2_12)
        B_I = (1 - self.sin2_12) * 4 * self.sin2_13_I * (1 - self.sin2_13_I)
        C_I = self.sin2_12 * 4 * self.sin2_13_I * (1 - self.sin2_13_I)

        self.prob_I = 1. - A_I * self.sin2(x,self.deltam_21) - B_I * self.sin2(x,self.deltam_3l_I+self.deltam_21) - C_I * self.sin2(x,self.deltam_3l_I)
        self.prob_E_I = 1. - A_I * self.sin2(x_E,self.deltam_21) - B_I * self.sin2(x_E,self.deltam_3l_I+self.deltam_21) - C_I * self.sin2(x_E,self.deltam_3l_I)

        if plot_this:
            
            fig = plt.figure()
            #fig.suptitle(r'Survival probability') # as function of L/E
            ax = fig.add_subplot(111)
            ax.grid()
            ax.set_xlim(left=0.02,right=80)
            ax.set_xlabel(r'L/$\text{E}_{\nu}$ [\si[per-mode=symbol]{\kilo\meter\per\MeV}]')
            ax.set_ylabel(r'P ($\bar{\nu}_{e} \rightarrow \bar{\nu}_{e}$)')
            ax.set_title(r'Survival probability')

            fig1 = plt.figure()
            #fig1.suptitle(r'Survival probability') # as function of E
            ax1 = fig1.add_subplot(111)
            ax1.grid()
            ax1.set_xlabel(r'$\text{E}_{\nu}$ [\si[per-mode=symbol]{\MeV}]')
            ax1.set_ylabel(r'P ($\bar{\nu}_{e} \rightarrow \bar{\nu}_{e}$)')
            ax1.set_title(r'Survival probability')

            if ordering == 1 or ordering == 0: # NO

                ax.semilogx(x/1000,self.prob_N,'b',linewidth=1,label='NO')
                ax.legend()
                fig.savefig('Oscillation/prob_N_LE.pdf',format='pdf',transparent=True)
                ax1.plot(E,self.prob_E_N,'b',linewidth=1,label='NO')
                ax1.legend()
                fig1.savefig('Oscillation/prob_N_E.pdf',format='pdf',transparent=True)

            if ordering == -1: # IO
                ax.semilogx(x/1000,self.prob_I,'r',linewidth=1,label='IO')
                ax.legend()
                fig.savefig('Oscillation/prob_I_LE.pdf',format='pdf',transparent=True)
                ax1.plot(E,self.prob_E_I,'r',linewidth=1,label='IO')
                ax1.legend()
                fig1.savefig('Oscillation/prob_I_E.pdf',format='pdf',transparent=True)

            if ordering == 0: # IO
                ax.semilogx(x/1000,self.prob_I,'r--',linewidth=1,label='IO')
                ax.legend()
                fig.savefig('Oscillation/prob_LE.pdf',format='pdf',transparent=True)
                ax1.plot(E,self.prob_E_I,'r--',linewidth=1,label='IO')
                ax1.legend()
                fig1.savefig('Oscillation/prob_E.pdf',format='pdf',transparent=True)

        return self.prob_E_N, self.prob_E_I


