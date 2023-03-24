# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
from scipy import special
import numpy as np 
from scipy.integrate import dblquad
from scipy.integrate import quad
import scipy.constants as cst

lamda = 1.55 #en microns
k = (2*np.pi/lamda)
def e_guide_moyen(kc, d, phi):
    return 1/(d*kc)*(np.sin((kc*d/2)+phi)-np.sin((-d*kc/2)+phi))

def kc(nc, neff): 
    return k*(nc**2-neff**2)

def asub(neff, nsub): 
    return k*(neff**2-nsub**2)

def asup(neff, nsup):
    return k*(neff**2-nsup**2)

def phi(kc, asup, asub): 
    return (1/2)*(np.arctan(asup/kc)-np.arctan(asub/kc)) 

def w(w0,z, zr):
    return w0*(1+(z/zr)**2)**(1/2)

#%% Dans le guide 10 microns
nc = 1.80
neff = 1.796
nsub = 1.46
nsup = 1 # ou 1
d= 10 #dimension du guide en microns 

kc = kc(nc, neff)
asup = asup(neff, nsup)
asub = asub(neff, nsub)
phi = phi(kc, asup, asub)

print(e_guide_moyen(kc, d, phi))


#%% Dans le guide 1 microns
nc = 1.80
neff = 1.62
nsub = 1.46
nsup = 1 # ou 1
d= 1 #dimension du guide 

kc = kc(nc, neff)
asup = asup(neff, nsup)
asub = asub(neff, nsub)
phi = phi(kc, asup, asub)

print(e_guide_moyen(kc, d, phi))
#%%Volume Bulk
L = 10000 #1cm
w0 = 5


zr = np.pi*w0**2*1.80/lamda

volume = np.pi*quad(lambda z:w(w0,z,zr)**2,-L/2, L/2)[0]
print(volume)

#%%Dans le bulk 
L = 10000 
w0 = 5




zr = np.pi*(w0**2)*1.80/lamda

def gaussian(r,z):  
    return (r/w(w0,z, zr))*np.exp(-r**2/w(w0,z, zr)**2)

average =((2*np.pi*w0)/volume)*dblquad(gaussian,-L/2, L/2,0,lambda z:w(w0,z,zr))[0]

print(average)

#%% taux de pompage et différence de population
#On se remet en SI
T1 = 11.4E-3
T2= 1E-6


def W(rabi): 
    return rabi**2*T2/2
    

rabi_bulk =  2*np.pi*690E3
rabi_guide =rabi_bulk*6 #Corrielli 2016 sur PrYSO


Wguide = W(rabi_guide)
Wbulk = W(rabi_bulk)


def n(t,pompe): 
    return (1/(1+2*T1*pompe))+(0.5-(1/(1+2*T1*pompe)))*np.exp(-(2*pompe+(1/T1))*t)

nguide = n(500E-6, Wguide)
nbulk = n(500E-6, Wbulk)