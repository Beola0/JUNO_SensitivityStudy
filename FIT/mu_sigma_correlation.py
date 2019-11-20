import sys
sys.path.insert(0,'/Users/beatricejelmini/Desktop/JUNO/JUNO_codes')

import latex
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy import stats
from mpl_toolkits.axes_grid1 import make_axes_locatable
from numpy import linalg as LA


def Gaussian (x,mu,sigma):
    const = math.sqrt(2*math.pi) * sigma
    return 1/const * np.exp(- np.power(x-mu,2)/2./(sigma**2)) * 1000.


bins = np.arange(-2.5,2.6,0.1)
x_asimov = Gaussian(bins,mu=0.,sigma=1.)

### construction of covariance matrix (with mu-sigma correlation)

N_samples = 5000
a = stats.norm.rvs(size=(2,N_samples))
corr = 0.5
b =  np.ndarray((2,N_samples)) # mu and sigma correlated
b[0,:] = a[0,:]
b[1,:] = corr * a[0,:] + np.sqrt(1-corr**2) * a[1,:] 
mu_err = 0.01
sigma_err = 0.007

# rescaling the normal random variables
mu_l = b[0,:].min()
mu_r = b[0,:].max()
b[0,:] = (b[0,:]-mu_l) * 2*(10*mu_err) / (mu_r - mu_l) - (10*mu_err)
sigma_l = b[1,:].min()
sigma_r = b[1,:].max()
b[1,:] = (b[1,:]-sigma_l) * 2*(10*sigma_err) / (sigma_r - sigma_l) - (10*sigma_err) + 1.

np.savetxt('corr_mu_sigma.txt',b,delimiter=',')

fig_ = plt.figure()
ax_ = fig_.add_subplot(111)
# generation of N samples
M = len(bins)
X_samples_asimov = np.empty((M,N_samples))
X_bins = np.empty((M,N_samples))
for N0 in np.arange(0,N_samples):
    X_samples_asimov[:,N0] = Gaussian(bins,mu=b[0,N0],sigma=b[1,N0])
    #X_bins[:,N0] = bins
    ax_.hist(bins,M,weights=X_samples_asimov[:,N0],histtype='step')
    
#ax_.hist(X_bins,M,weights=X_samples_asimov,histtype='step',stacked=True, fill=False)

# construction of the covariance matric
V = np.zeros((M,M),dtype='float32')
for j in np.arange(0,M):
    for k in np.arange(0,M):   
        # sum over samples
        V[j,k] = ((X_samples_asimov[j,:] - x_asimov[j]) * (X_samples_asimov[k,:] - x_asimov[k])).sum()

V = V/N_samples
# inversion of the covariance matrix
V_inv = np.linalg.inv(V)

np.savetxt('corr_matrix.txt',V,delimiter=',')

# plot of correlated parameter
fig, axScatter = plt.subplots(figsize=(8., 6.5))
fig.suptitle(r'$\mu - \sigma$ correlation, $\rho$ = %.2f' % (corr))
axScatter.scatter(b[0], b[1], s=10., c='b', marker='.')
axScatter.set_aspect(1.)
axScatter.set_xlabel(r'$\mu$')
axScatter.set_ylabel(r'$\sigma$')
divider = make_axes_locatable(axScatter)
axHistx = divider.append_axes("top", 1.7, pad=0.1, sharex=axScatter)
axHisty = divider.append_axes("right", 1.7, pad=0.1, sharey=axScatter)
plt.setp(axHistx.get_xticklabels() + axHisty.get_yticklabels(), visible=False)
axHistx.hist(b[0],bins=10,color='b',histtype='step')
axHisty.hist(b[1],bins=10,color='b',orientation='horizontal',histtype='step')
fig.savefig('mu_sigma_corr.pdf', format='pdf')
plt.draw()

# plot of the covariance matrix
fig_mat = plt.figure()
ax_mat = fig_mat.add_subplot(111)
ax_mat.set_title(r'Covariance matrix - $\rho$ = %.2f' % (corr))
im = ax_mat.imshow(V,cmap='brg',interpolation='none',origin={'lower','lower'}) 
bar = fig_mat.colorbar(im)
fig_mat.savefig('cov_mat.pdf', format='pdf')


plt.ion()
plt.show()
