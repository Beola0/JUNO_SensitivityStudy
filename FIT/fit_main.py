import matplotlib.pyplot as plt
import numpy as np
import math
from scipy import integrate, stats, optimize
import time
import iminuit
from chi_squared import Chi_Squared


def Gaussian (x,mu,sigma):
    const = math.sqrt(2*math.pi) * sigma
    return 1/const * np.exp(- np.power(x-mu,2)/2./(sigma**2)) * 1000.

def chi_squared (params,f,x_data,bins):
    #mu = params[0]
    #sigma = params[1]
    x_fit = f(bins,*params)
    sigma_data = np.sqrt(x_data)
    np.place(sigma_data,sigma_data<1.,1.) 
    appo = (x_fit - x_data)/ sigma_data
    appo = np.power(appo,2)
    return appo.sum()

time_start = time.process_time_ns()
# main program
bins = np.arange(-5.,5.2,0.1)

chi_sq = Chi_Squared()

# without poisson fluctuations

x_asimov = Gaussian(bins,mu=0.,sigma=1.)

params0 = np.ndarray((2,),buffer=np.array([0.1,1.])) # mu and sigma
res = optimize.minimize(chi_squared,x0=params0,args=(Gaussian,x_asimov,bins),method='BFGS')
x_fit = Gaussian(bins,res.x[0],res.x[1])

print('Minimization results:')
print(res)

mu = np.arange(-0.06,0.061,0.001) # M
sigma = np.arange(0.94,1.061,0.001) # N
chi = np.ndarray((len(sigma),len(mu))) # N * M

for n in np.arange(len(sigma)):
    for m in np.arange(len(mu)):
        chi[n,m] = chi_squared((mu[m],sigma[n]),Gaussian,x_asimov,bins)

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
ax2.set_title('Asimov dataset')
entries, edges, _ = ax2.hist(bins,len(bins),weights=x_asimov,color='b',histtype='step',label='histo')
ax2.errorbar(bins,entries,yerr=np.sqrt(entries),fmt='r.',elinewidth=1.,capsize=1.5,markersize=3.,label='error bars')
ax2.plot(bins,x_fit,'g',linewidth=1.5,label='fit')
ax2.legend()
ax2.grid()

fig3 = plt.figure()
ax3 = fig3.add_subplot(111)
ax3.set_title('Contour plot')
ax3.plot(res.x[0],res.x[1],'k.',markersize=3.)
level = np.array([1,4,9,16,25]) 
cs = ax3.contour(mu,sigma,chi,50) 
ax3.set_xlabel('mu')
ax3.set_ylabel('sigma')
cbar = fig3.colorbar(cs)
cbar.ax.set_ylabel('chi squared')
ax3.grid()



# with poisson fluctuations

x_poisson = np.ndarray(len(x_asimov))
n=0
for x in x_asimov:
    x_poisson[n] = np.random.poisson(x,1)
    n += 1

params0_p = np.ndarray((2,),buffer=np.array([0.,1.])) # mu and sigma
res_p = optimize.minimize(chi_squared,x0=params0_p,args=(Gaussian,x_poisson,bins),method='BFGS')
res_iminuit = iminuit.minimize(chi_squared,x0=params0_p,args=(Gaussian,x_poisson,bins))
m_ = res_iminuit.minuit
x_fit_p = Gaussian(bins,res_p.x[0],res_p.x[1])

print('\nMinimization results (with poisson fluctuations):')
print(res_p)

print('\nMinimization results (with poisson fluctuations) - iminuit:')
print(res_iminuit)

mu_p = np.arange(-0.06,0.061,0.001) # M
sigma_p = np.arange(0.94,1.061,0.001) # N
chi_p = np.ndarray((len(sigma),len(mu))) # N * M

for n in np.arange(len(sigma_p)):
    for m in np.arange(len(mu_p)):
        chi_p[n,m] = chi_squared((mu_p[m],sigma_p[n]),Gaussian,x_poisson,bins)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_title('Dataset with poisson fluctuations')
entries_p, edges, _ = ax.hist(bins,len(bins),weights=x_poisson,color='b',histtype='step',label='histo')
ax.errorbar(bins,entries_p,yerr=np.sqrt(entries_p),fmt='r.',elinewidth=1.,capsize=1.5,markersize=3.,label='error bars')
ax.plot(bins,x_fit_p,'g',linewidth=1.,label='fit')
ax.legend()
ax.grid()

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.set_title('Contour plot with poisson fluctuations')
ax1.plot(res_p.x[0],res_p.x[1],'k.',markersize=3.)
level = np.array([1,4,9,16,25]) 
cs = ax1.contour(mu_p,sigma_p,chi_p,50) 
ax1.set_xlabel('mu')
ax1.set_ylabel('sigma')
cbar = fig1.colorbar(cs)
cbar.ax.set_ylabel('chi squared')
ax1.grid()

fig4 = plt.figure()
ax4 = fig4.add_subplot(111)
ax4.set_title('Contour plot - iminuit')
m_.draw_mncontour('x0','x1',nsigma=5)  
ax4.set_xlabel('mu') 
ax4.set_ylabel('sigma')  
ax4.grid() 

time_el = time.process_time_ns()
print('elapsed time: ' + str(time_el*10**(-6)) + ' ms')


plt.ion()
plt.show()



'''popt,pcov = optimize.curve_fit(Gaussian,bins,x_asimov)
rho = pcov[0,1] / math.sqrt(pcov[0,0]) / math.sqrt(pcov[1,1])
x_fit = Gaussian(bins,popt[0],popt[1])
print('\nmu: ' + str(popt[0]) + '\nsigma: ' + str(popt[1]) )
print('\ncovariance matrix:')
print(pcov)
print('correlation: ' + str(rho))

# chi square
chisq, p = stats.chisquare(x_asimov,x_fit)
print('\nscipy.stats.chisquare')
print('chisquare: ' + str(chisq))
print('p-value: ' + str(p))

# log likelihood
stat, p_stat = stats.power_divergence(x_asimov,x_fit,lambda_=0)
print('\nscipy.stats.power_divergence')
print('stat: ' + str(stat))
print('p-value: ' + str(p_stat)) '''

