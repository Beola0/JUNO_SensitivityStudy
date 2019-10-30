from convolution import Convolution
import numpy as np
from scipy import integrate
from scipy import signal
import matplotlib.pyplot as plt
import math
import time

def Gaussiana (x):
    const = math.sqrt(2*math.pi)
    sigma = 0.4
    appo = 1/const/sigma * np.exp(- np.power(x,2)/2./(sigma**2))
    #norm = scipy.integrate.simp(appo,x)
    return appo

time_start = time.process_time_ns()

imp = signal.unit_impulse(101,[20,80])
    
x = np.arange(0.,10.1,0.1)
g = Gaussiana (x-5)
#print(integrate.simps(g,x))

convol = Convolution()
res_np = convol.np_conv(imp,x,0.4)
res_num = convol.numerical_conv(imp,x,0.4)

#g2 = Gaussiana (x)
conv1 = np.convolve(imp,g,mode='same')
conv2 = signal.convolve(imp,g,mode='same',method='direct')

# convolution step by step
t = np.arange(0,10.,0.1)
conv = np.zeros(len(t))
n=0
for t0 in t:
    appo = Gaussiana(x-t0)
    prod = appo * imp
    conv[n] = prod.sum()
    n += 1

el_time = time.process_time_ns() - time_start
print('elapsed time: ',str(el_time*10**(-6)),' ms')


#fig = plt.figure()
#ax = fig.add_subplot(111)
#ax.plot(x,imp,'b')
#ax.plot(x,g,'r')
#ax.plot(x,conv1,'g')
#ax.grid()

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.plot(x,conv1,'b',label='np')
ax1.plot(x,res_np,'r',label='class - np')
#ax1.plot(x,conv2,'r',label='scipy - same')
ax1.plot(t,conv,'g',label='num conv')
ax1.plot(x,res_num,'k',label='class - num conv')
ax1.grid()
ax1.legend()

plt.ion()
plt.show()
