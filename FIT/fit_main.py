import matplotlib.pyplot as plt
import numpy as np
import math
from scipy import integrate

def Gaussian (x,mu,sigma):
    const = math.sqrt(2*math.pi) * sigma
    return 1/const * np.exp(- np.power(x-mu,2)/2./(sigma**2))


x = np.arange(-5.,5.,0.1)
g = Gaussian(x,mu=0.,sigma=1.)



# plot as a function
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(x,g,'r',marker='.',linewidth=1.,label='gaussian')


# plot as a histogram
fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.bar(x,g,0.1)

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)

'''for x0 in x:
    appo = Gaussian(x0,mu=0.,sigma=1.)
    ar = np.full(int(str(appo*100)[:2]),x0) '''
ax2.hist(g)

plt.ion()
plt.show()
