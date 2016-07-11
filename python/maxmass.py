##  this code does the following:
##  read in table of max BH mass vs metals
##  interpolate this data
##  makes a uniform probability function of mass for a given metallicity
##  calculates the average chirp mass for this metallicity

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import scipy.interpolate as interpol

metals,maxbh = np.loadtxt("mass_max_of_z.dat",unpack=True)

#print metals
#print maxbh

metalgrid = np.linspace(min(metals),max(metals),1000)


f = interpol.interp1d(metals,maxbh)
maxmassgrid = f(metalgrid)
#print f(0.1)
#print "******"

"""
plt.plot(metals,maxbh)
plt.xlabel("metals")
plt.ylabel("max BH mass")
plt.show()
### ok we are good with the interp. ###
"""

input_metal = 0.1
#sigma = 5.0
"""thismetal = np.abs(metalgrid - input_metal).argmin()
print thismetal
place = thismetal #np.where(metalgrid == thismetal)
print place
maxmass = maxmassgrid[place]
print maxmass"""

"""
def probability(input_metal,sigma):
    place = np.abs(metalgrid - input_metal).argmin()
    maxmass = maxmassgrid[place]
    #print place,maxmass
    prob = np.exp(-np.power(maxmassgrid-maxmass,2.0) / (2*np.power(sigma,2.)))
    return prob

dist = probability(input_metal,sigma)

plt.plot(maxmassgrid,dist)
plt.xlabel("BH mass")
plt.ylabel("probability")
plt.title(" input metallicity of "+str(input_metal))
plt.show()                
"""

def make_uniform_func(metal_now):
    MmaxNow = f(metal_now)
    #print MmaxNow
    nm = 1./(MmaxNow - 3)  # normalization
    return lambda x, Nm = nm, A = MmaxNow: np.where( np.logical_and(x<A , x > 3), nm* np.ones(len(x)), np.zeros(len(x)))


prob_function = make_uniform_func(input_metal) # input is a metallicity
newmassgrid = np.linspace(0,100,50)
#print np.logical_and(newmassgrid < 80 , newmassgrid > 3)
p = prob_function(newmassgrid)
plt.plot(newmassgrid,p)
#print p
plt.show()


import scipy
def average_mass(input_probability_distribution):
    masses = np.linspace(0,100,50)
    prob = input_probability_distribution(masses)
    result = scipy.integrate.simps(prob*np.power(masses, 15./6), masses)
    return np.power(result, 6/15.)

the_answer = average_mass(prob_function)
print "the mean mass for the metallicity",input_metal, " is" ,the_answer
