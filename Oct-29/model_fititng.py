import h5py
import numpy as np
import pylab as plt
from scipy.optimize import curve_fit
f1 = h5py.File('data.hdf','r')

x = f1['data/xpos'][:]
y = f1['data/ypos'][:]


plt.plot(x,y, 'o', label = 'raw data')
#plt.show()

def fit1(x,a,b,c, d, e):
    return a * np.cos(b * x) - d*np.exp(c*x) + e/(x)

params, pcov = curve_fit(fit1,x,y)
print(params)
fit_func = fit1(x,*params)

chi2 = np.sum((fit_func - y)**2)
print(chi2)
a = chi2/(len(fit_func) - len(params) - 1)
b = np.sum((y - np.mean(y))**2)/(len(y) - 1)
R2 = 1 - a/b

print(R2)
x = np.sort(x)
plt.plot(x,fit1(x,*params), 'r-', label = 'fitted')
plt.legend()
plt.show()