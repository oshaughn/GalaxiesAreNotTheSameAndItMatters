
import numpy as np
#import astropy
from astropy.cosmology import WMAP9 as cosmo   # http://astropy.readthedocs.org/en/latest/cosmology/

rosDebug =True

jbTimeScale = 38.7811  #  = = 1.223*10^18 Convert[Second/(Giga Year), 1]
jbStarMass =  26676*0.3
jbGasMass = 26676.0  # Msun?

class StarCatalog:
    """
    Class containing all information loaded from star output files.
    """
    def __init__(self, fname, z):
        if rosDebug:
            print "Loading ", fname,
        # Data format: 
        self.dat = np.loadtxt(fname)
        self.z = float(z)
        self.mass = len(self.dat)*jbStarMass  # I suspect this is wrong because I don't remember how to distinguish between gas        
        if rosDebug:
            print " ....finished, saved as z= ", z, " with log M = ", np.log10(self.mass)
        # links, so no cost
        self.Zs = self.dat[:,0]  # Star particle birth 
        self.Zg = self.dat[:,1]  # Mean metallicity of nearest 32 particles
        self.tb  = self.dat[:,2]  # time in Gyr.  Beware cosmology
        # if rosDebug:
        #     print  "    cosmo time ", cosmo.age(self.z), " versus max t in file ", np.max(self.tb)    # check time conversion and cosmology
        self.weights = np.ones(len(self.dat))*1.0/len(self.dat)  # treat all particles equally

    def time(self):
        return cosmo.age(self.z)

    def integrate(self,fn):   # fn takes 3 arguments: Zs, Zg, and t in Gyr
        return np.sum(fn(self.Zs,self.Zg, self.tb)*self.weights)   # integrate an arbitrary expression of these quantities. Use for adding in time and Z dependent formation


def cumulative_P(vals):
    vals_sorted = np.sort(np.copy(vals))
    P = np.arange(len(vals))*1.0/len(vals)
    return vals_sorted,P

def cumulative_P_weighted(vals, wts):
    idx_sorted_util  = np.lexsort((np.arange(len(vals)), vals))
    vals_sorted  = np.array([vals[k] for k in idx_sorted_util])
    wts  = np.array([wts[k] for k in idx_sorted_util])
    cum_wts = np.cumsum(wts)
    cum_wts = cum_wts/cum_wts[-1]

    return vals_sorted, cum_wts
    

