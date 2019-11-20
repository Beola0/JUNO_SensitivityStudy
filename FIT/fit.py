import sys
sys.path.insert(0,'/Users/beatricejelmini/Desktop/JUNO/JUNO_codes')

import latex
import matplotlib.pyplot as plt
import numpy as np
import math, scipy
import time
from chi_squared import Chi_Squared
from iminuit import Minuit
from scipy import stats
#from scipy.linalg import cholesky
from mpl_toolkits.axes_grid1 import make_axes_locatable
from numpy import linalg as LA


def Gaussian (x,mu,sigma):
    const = math.sqrt(2*math.pi) * sigma
    return 1/const * np.exp(- np.power(x-mu,2)/2./(sigma**2)) * 1000.


# MAIN PROGRAM

time_start_c = time.perf_counter_ns()
bins = np.arange(-2.5,2.6,0.1)



### without poisson fluctuations

x_asimov = Gaussian(bins,mu=0.,sigma=1.)
params0 = np.ndarray((2,),buffer=np.array([0.,1.])) # mu and sigma

chi_sq_g = Chi_Squared(Gaussian,x_asimov,bins)
m_g = Minuit(chi_sq_g.chi_squared_stat_minuit,mu=params0[0],sigma=params0[1],pedantic=False)
m_g.migrad()
m_g.minos()
x_fit = Gaussian(bins,m_g.values['mu'],m_g.values['sigma'])
cov_mat = m_g.np_matrix()
rho = cov_mat[0,1] / math.sqrt(cov_mat[0,0]) / math.sqrt(cov_mat[1,1])

print('\nMinimization results:'+'\nparameters:')
print(m_g.values)
print('covariance matrix:')
print(m_g.np_matrix())
print('correlation: ' + str(rho)) 

'''fig = plt.figure(figsize=[12.5, 7.5])
ax = fig.add_subplot(121)
ax.set_title(r'Asimov dataset')
entries, edges, _ = ax.hist(bins,len(bins),weights=x_asimov,color='b',histtype='step',label=r'histo')
#bin_centers = np.zeros(len(bins)) 
#for n in np.arange(len(bins)): 
#    bin_centers[n] = (edges[n+1] + edges[n])/2.      
ax.errorbar(bins,entries,yerr=np.sqrt(entries),fmt='r.',elinewidth=1.,capsize=1.5,markersize=3.,label=r'error bars')
ax.plot(bins,x_fit,'g',linewidth=1.5,label=r'fit')
ax.set_ylabel(r'N')
ax.legend()
ax.grid()

# contour plot
mu = np.arange(-0.1,0.101,0.001) # M
sigma = np.arange(0.9,1.101,0.001) # N
chi = np.ndarray((len(sigma),len(mu))) # N * M

for n in np.arange(len(sigma)):
    for m in np.arange(len(mu)):
        chi[n,m] = chi_sq_g.chi_squared_stat_minuit(mu[m],sigma[n])

ax2 = fig.add_subplot(122)
ax2.set_title(r'Contour plot')
#m_g.draw_mncontour('mu','sigma')
ax2.plot(m_g.values['mu'],m_g.values['sigma'],'k.',markersize=3.)
level = np.array([1,4,9,16,25]) 
chi_appo = chi - chi.min()
cs = ax2.contourf(mu,sigma,chi_appo,level) 
ax2.set_xlabel(r'$\mu$')
ax2.set_ylabel(r'$\sigma$')
cbar = fig.colorbar(cs)
cbar.ax.set_ylabel(r'$\Delta \chi^2$')
ax2.grid() '''



'''### with poisson fluctuations

x_poisson = np.ndarray(len(x_asimov))
n=0
for x in x_asimov:
    x_poisson[n] = np.random.poisson(x,1)
    n += 1

chi_sq_gp = Chi_Squared(Gaussian,x_poisson,bins)
params0_p = np.ndarray((2,),buffer=np.array([0.,1.])) # mu and sigma

m_p = Minuit(chi_sq_gp.chi_squared_stat_minuit,mu=params0_p[0],sigma=params0_p[1],pedantic=False) #fix_mu = True 
m_p.migrad()
m_p.minos()
x_fit_p = Gaussian(bins,m_p.values['mu'],m_p.values['sigma'])
cov_mat_p = m_p.np_matrix()
rho_p = cov_mat_p[0,1] / math.sqrt(cov_mat_p[0,0]) / math.sqrt(cov_mat_p[1,1])

print('\nMinimization results (with poisson fluctuations):'+'\nparameters:')
print(m_p.values)
print('covariance matrix:')
print(m_p.np_matrix())
print('correlation: ' + str(rho_p))

fig3 = plt.figure(figsize=[12.5, 7.5])
ax3 = fig3.add_subplot(121)
ax3.set_title(r'Dataset with poisson fluctuations')
entries_p, edges, _ = ax3.hist(bins,len(bins),weights=x_poisson,color='b',histtype='step',label=r'histo')
ax3.errorbar(bins,entries_p,yerr=np.sqrt(entries_p),fmt='r.',elinewidth=1.,capsize=1.5,markersize=3.,label=r'error bars')
ax3.plot(bins,x_fit_p,'g',linewidth=1.,label=r'fit')
ax3.legend()
ax3.grid()

mu_p = np.arange(-0.1,0.101,0.001) # M
sigma_p = np.arange(0.9,1.101,0.001) # N
chi_p = np.ndarray((len(sigma_p),len(mu_p))) # N * M

for n in np.arange(len(sigma_p)):
    for m in np.arange(len(mu_p)):
        chi_p[n,m] = chi_sq_gp.chi_squared_stat_minuit(mu_p[m],sigma_p[n])

ax4 = fig3.add_subplot(122)
ax4.set_title(r'Contour plot (poisson fluctuations)')
#m_p.draw_mncontour('mu','sigma',nsigma=5)  
ax4.plot(m_p.values['mu'],m_p.values['sigma'],'k.',markersize=3.)
level = np.array([1,4,9,16,25]) 
chi_p_appo = chi_p - chi_p.min()
cs = ax4.contourf(mu_p,sigma_p,chi_p_appo,level) 
ax4.set_xlabel(r'$\mu$')
ax4.set_ylabel(r'$\sigma$')
cbar = fig3.colorbar(cs)
cbar.ax.set_ylabel(r'$\Delta \chi^2$')
ax4.grid() '''

#time_el_c = time.perf_counter_ns()
#print('\nelapsed time: ' + str(time_el_c*10**(-6)) + ' ms') 



'''### test of mu and sigma distribution obtained from the minimization of 1000 samples with Poisson fluctuations

dim = 1000
mu_array = np.zeros(dim)
sigma_array = np.zeros(dim)

for n0 in np.arange(0,dim):
    x_poisson = np.ndarray(len(x_asimov))
    n=0
    for x in x_asimov:
        x_poisson[n] = np.random.poisson(x,1)
        n += 1
    chi_sq_ = Chi_Squared(Gaussian,x_poisson,bins)
    m_ = Minuit(chi_sq_.chi_squared_stat_minuit,mu=params0_p[0],sigma=params0_p[1],pedantic=False) 
    m_.migrad()
    mu_array[n0] = m_.values['mu']
    sigma_array[n0] = m_.values['sigma']

fig5 = plt.figure(figsize=[12.5, 7.5])

N_bins = 25
ax5 = fig5.add_subplot(121)
ax5.set_title(r'$\mu$ distribution')
ax5.hist(mu_array,N_bins,color='b',histtype='step')
ax5.set_xlabel(r'$\mu$')
ax5.set_ylabel(r'N')
ax5.grid()

ax6 = fig5.add_subplot(122)
ax6.set_title(r'$\sigma$ distribution')
ax6.hist(sigma_array,N_bins,color='r',histtype='step')
ax6.set_xlabel(r'$\sigma$')
ax6.set_ylabel(r'N')
ax6.grid()'''



'''### construction of covariance matrix (without correlation) for poisson fluctuations, only statistical uncertainty

M = len(bins)
V = np.full((M,M),0.)
V_2 = np.full((M,M),0.)
N_samples = 100.
x_data_samples = np.empty((int(N_samples),M))
R = np.full((M,M),0.)

sigma_j = np.full((M),0.)
sigma_k = np.full((M),0.)

N_ = 0
for N0 in np.arange(0,N_samples):

    x_poisson = np.ndarray(len(x_asimov))
    n=0
    for x in x_asimov:
        x_poisson[n] = np.random.poisson(x,1)
        n += 1

    for j in np.arange(0,M):
        for k in np.arange(0,M):   
            V[j,k] += (x_poisson[j] - x_asimov[j]) * (x_poisson[k] - x_asimov[k])

    x_data_samples[N_,] = x_poisson
    N_ += 1


for j in np.arange(0,M):
    appo = ((x_data_samples[:,j] -  x_asimov[j])**2).sum()
    sigma_j[j] = math.sqrt(appo)

for k in np.arange(0,M):
    appo = ((x_data_samples[:,k] -  x_asimov[k])**2).sum()
    sigma_k[k] = math.sqrt(appo)

for j in np.arange(0,M):
    for k in np.arange(0,M):   
        s_j = (x_data_samples[:,j] - x_asimov[j])**2
        s_j = s_j.sum()
        s_k = (x_data_samples[:,k] - x_asimov[k])**2
        s_k = s_k.sum()
        V_2[j,k] = ((x_data_samples[:,j] - x_asimov[j]) * (x_data_samples[:,k] - x_asimov[k])).sum()
        R[j,k] = V_2[j,k] / (np.sqrt(s_j) * np.sqrt(s_k))

V = V/N_samples
V_2 = V_2/N_samples
np.savetxt('cov_mat.txt',V)
V_inv = np.linalg.inv(V) '''

# poisson correlation matrix
M = len(bins)
V_p = np.full((M,M),0.)
n=0
for x0 in x_asimov:
    V_p[n,n] = x0
    n += 1

#np.place(V_p,np.diag(V_p)<1.,1.) 
V_p_inv = np.linalg.inv(V_p)

### correlation matrix loaded from file

b = np.loadtxt('corr_mu_sigma.txt',delimiter=',',unpack=False) 
V = np.loadtxt('corr_matrix.txt',delimiter=',',unpack=False) 
V_inv = np.linalg.inv(V)

x_data = x_asimov
x_data_bis = Gaussian(bins,b[0,4],b[1,4])
params0 = np.ndarray((2,),buffer=np.array([0.,1.])) 

'''def fun (mu,sigma):
    x_fit = Gaussian(bins,mu,sigma)
    X = x_data - x_fit
    M = len(bins)
    appo = np.zeros(M)
    for j in np.arange(0,M):
        appo[j] = (V_inv[j,:] * X).sum()
    appo1 = X * appo

    return appo1.sum() '''

chi_sq = Chi_Squared(Gaussian,x_data,bins,sigma='')
chi_sq.set_invert_mat(V)
m_c = Minuit(chi_sq.chi_squared_V,mu=params0[0],sigma=params0[1],pedantic=False)
m_c.migrad()
#m_c.minos()
x_fit_c = Gaussian(bins,m_c.values['mu'],m_c.values['sigma'])
cov_mat = m_c.np_matrix()
rho = cov_mat[0,1] / math.sqrt(cov_mat[0,0]) / math.sqrt(cov_mat[1,1])

print('\nMinimization results:'+'\nparameters:')
print(m_c.values)
print('covariance matrix:')
print(m_c.np_matrix())
print('correlation: ' + str(rho)) 

fig_c = plt.figure(figsize=[12.5, 7.5])
fig_c.suptitle(r'Dataset with mu-sigma correlation')
ax_c1 = fig_c.add_subplot(121)
ax_c1.set_title(r'Dataset')
entries, edges, _ = ax_c1.hist(bins,len(bins),weights=x_data,color='b',histtype='step',label=r'histo')    
#x_c1.errorbar(bins,entries,yerr=np.sqrt(entries),fmt='r.',elinewidth=1.,capsize=1.5,markersize=3.,label=r'error bars')
ax_c1.plot(bins,x_fit_c,'g',linewidth=1.5,label=r'fit')
ax_c1.set_ylabel(r'N')
ax_c1.legend()
ax_c1.grid() 

# contour plot
mu = np.arange(-0.1,0.101,0.001) # M
sigma = np.arange(0.9,1.101,0.001) # N
mu = np.arange(-0.13,0.132,0.002) # M
sigma = np.arange(0.87,1.132,0.002) # N
#mu = np.arange(-0.06,0.0605,0.0005) # M
#sigma = np.arange(0.94,1.0605,0.0005) # N
#mu = np.arange(-0.01,0.0115,0.0005) # M
#sigma = np.arange(0.993,1.0075,0.0005) # N
chi_c = np.ndarray((len(sigma),len(mu))) # N * M

for n in np.arange(len(sigma)):
    for m in np.arange(len(mu)):
        #chi_c[n,m] = fun(mu[m],sigma[n])
        chi_c[n,m] = chi_sq.chi_squared_V(mu[m],sigma[n])

ax_c2 = fig_c.add_subplot(122)
ax_c2.set_title(r'Contour plot')
#m_c.draw_mncontour('mu','sigma')
ax_c2.plot(m_c.values['mu'],m_c.values['sigma'],'k.',markersize=3.)
level = np.array([1,4,9,16,25]) 
chi_appo = chi_c - chi_c.min()
cs = ax_c2.contourf(mu,sigma,chi_appo,level) 
#ax_c2.contour(mu, sigma, chi_appo, cs.levels, colors='k')
ax_c2.set_xlabel(r'$\mu$')
ax_c2.set_ylabel(r'$\sigma$')
cbar = fig_c.colorbar(cs)
cbar.ax.set_ylabel(r'$\Delta \chi^2$')
ax_c2.grid()


plt.ion()
plt.show()



