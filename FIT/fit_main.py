import matplotlib.pyplot as plt
import numpy as np
import math
from scipy import integrate, stats, optimize


def Gaussian (x,mu,sigma):
    const = math.sqrt(2*math.pi) * sigma
    return 1/const * np.exp(- np.power(x-mu,2)/2./(sigma**2)) * 1000.


x = np.arange(-5.,5.2,0.2)
g = Gaussian(x,mu=0.,sigma=1.)

popt,pcov = optimize.curve_fit(Gaussian,x,g)
rho = pcov[0,1] / math.sqrt(pcov[0,0]) / math.sqrt(pcov[1,1])
g_exp = Gaussian(x,popt[0],popt[1])
print('\nmu: ' + str(popt[0]) + '\nsigma: ' + str(popt[1]) )
print('\ncovariance matrix:')
print(pcov)
print('correlation: ' + str(rho))

# chi square
chisq, p = stats.chisquare(g,g_exp)
print('\nscipy.stats.chisquare')
print('chisquare: ' + str(chisq))
print('p-value: ' + str(p))


# log likelihood
stat, p_stat = stats.power_divergence(g,g_exp,lambda_=0)
print('\nscipy.stats.power_divergence')
print('stat: ' + str(stat))
print('p-value: ' + str(p_stat))





# plot as a function
'''fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(x,g,'r',marker='.',linewidth=1.,label='gaussian') '''


# plot as a histogram
'''fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.bar(x,g,0.1,color='b')
#ax1.hist(x,len(x),weights=g,color='r') '''


#Nbins = (x[-1] - x[0])/ 0.1
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
ax2.set_title('Asimov dataset')
#ax2.plot(x,g,'r',linewidth=1.)
entries, edges, _ = ax2.hist(x,len(x),weights=g,color='b',histtype='step',label='histo')
ax2.errorbar(x,entries,yerr=np.sqrt(entries),fmt='r.',elinewidth=1.,capsize=1.5,markersize=3.,label='error bars')
ax2.plot(x,g_exp,'g',linewidth=1.5,label='fit')
#ax2.hist(x,len(x),weights=x/10.,color='r')
ax2.legend()
ax2.grid()

plt.ion()
plt.show()
