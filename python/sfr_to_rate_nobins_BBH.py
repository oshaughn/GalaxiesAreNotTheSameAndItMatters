#! /usr/bin/env python
#  sfr_to_rate_nobins.py 
#  GOAL
#      Take SFR history and spit out an event rate, using different response functions
#      --time-exponent 0   : use uniform in time, starting after 10 Myr and running to 
#      --time-exponent 1   : use 1/t delay time, start (DEFAULT)
#      --time-max    100    : time max in Gyr for response function (DEFAULT).  If time-exponent<1, time-max=10Gyr
# 
#      --starmetal-file
#      --Z-exponent 0
#
#  USAGE
#    python sfr_to_rate_nobins.py --starmetal-file ~/PersonalJBArchive/boring-galaxy/z0.016.starmetal.dat --Z-exponent 2
#    python sfr_to_rate_nobins.py --starmetal-file ~/PersonalJBArchive/boring-galaxy/z0.016.starmetal.dat
#
# PLOTS
#    Third plot is money plot: shows you where the mergers are coming from, for different exponent choices.

import matplotlib.pyplot as plt
import sys
import numpy as np
from scipy.integrate import simps

import bisect
from scipy.interpolate import interp1d
from scipy.interpolate import griddata
from scipy.integrate import quad

import ParseStars

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--starmetal-file", default=None)
parser.add_argument("--type-bbh",action='store_true',help="Volume function selection bias (and Z exponent?) appropriate to BBH")
parser.add_argument("--type-bhns",action='store_true',help="Volume function selection bias (and Z exponent?) appropriate to BH-NS")
parser.add_argument("--time-exponent",default=1,type=float)
parser.add_argument("--Z-exponent",default=None,type=float)
parser.add_argument("--verbose", action='store_true')
opts = parser.parse_args()
if not opts.starmetal_file:
    print "  Please specify an SFR file "


# Maximum mass information (single BH)

metals,maxbh = np.loadtxt("mass_max_of_z.dat",unpack=True)
metals = metals[::-1]  # write in order
maxbh = maxbh[::-1]
metals = np.concatenate( (np.array([1e-18]), metals,[1e9]))
maxbh = np.concatenate( ([maxbh[0]],maxbh,[maxbh[-1]] ))

#plt.plot(metals, maxbh);plt.legend(); plt.xlabel("metals") ; plt.ylabel("max BH mass")
#plt.savefig('inputs.png'); plt.clf
#plt.plot(np.log10(metals), maxbh);plt.legend(); plt.xlabel("metals") ; plt.ylabel("max BH mass")
#plt.savefig('inputs-scaled.png')

#maxbh_of_logZ = interp1d(np.log10(metals),maxbh,bounds_error=True,fill_value = np.max(maxbh),kind='linear')
#maxbh_of_Z = lambda x: maxbh_of_logZ(np.log10(x))
maxbh_of_Z = interp1d(metals,maxbh,bounds_error=False)#,bounds_error=True,fill_value =np.max(maxbh),kind='linear',copy=True)

#print np.min(np.log10(metals)), np.max(np.log10(metals)), maxbh_of_logZ(0.0)


#print np.min(metals),np.max(metals), maxbh_of_Z(1)
#plt.plot(metals,maxbh,label='raw')
# Zvals = np.linspace(0.001, 2,100)
# plt.plot(Zvals, maxbh_of_Z(Zvals),label='fit')
# print maxbh_of_Z(Zvals)
# plt.legend()
# plt.xlabel("metals")
# plt.ylabel("max BH mass")
# plt.savefig('tmp.png')
# plt.show()
# sys.exit(0)

def mc(m1,m2):
    return np.power(m1*m2,3./5.)/np.power(m1+m2, 1./5.)

def average_func(func,input_probability_distribution):
    masses = np.linspace(0.1,100,500)  # zero total mass causes problems with mchirp definition
    prob = input_probability_distribution(masses)
    result = simps(prob*func(masses), masses)
    return result

def make_uniform_func(metal_now):
    MmaxNow = maxbh_of_Z(metal_now)
#    print MmaxNow
    nm = 1./(MmaxNow - 3)  # normalization
    return lambda x, Nm = nm, A = MmaxNow: np.where( np.logical_and(x<A , x > 3), nm* np.ones(len(x)), np.zeros(len(x)))


def make_power_func(metal_now,mass_exponent=2):
    MmaxNow = maxbh_of_Z(metal_now)
    ret, err = quad( lambda x: np.power(x,-1*mass_exponent),3,MmaxNow)
    nm  =ret
#    print metal_now, MmaxNow, nm
    return lambda x, Nm = nm, A = MmaxNow,p=mass_exponent: np.where( np.logical_and(x<A , x > 3), nm* np.power(x,-1.*p), np.zeros(len(x)))

def Vaveraged_bns(Z):
    mc_av =mc(1.4, 1.4)
    return np.power(mc_av, 15./6.)

def Vaveraged_bhns(Z):
    prob_Mbh = make_uniform_func(Z)
    mc_av = np.power(average_func(lambda x: np.power(mc(1.4,x),15./6), prob_Mbh),5./16.)
    nm_check = average_func(lambda x: 1, prob_Mbh)
#    print " internal ", Z, mc_av, nm_check
    return np.power(mc_av, 15./6.)/nm_check

def Vaveraged_bhbh(Z):
    prob_Mbh = make_power_func(Z)
    mc_av = np.power(average_func(lambda x: np.power(mc(x,x),15./6), prob_Mbh),5./16.)
    nm_check = average_func(lambda x: 1, prob_Mbh)
#    print " internal ", Z, mc_av, nm_check
    return np.power(mc_av, 15./6.)/nm_check


# val = Vaveraged_bhns(1)
# prob_Mbh = make_uniform_func(1)
# mc_av = np.power(average_func(lambda x: np.power(mc(1.4,x),15./6), prob_Mbh),5./16.)
# print 1, maxbh_of_Z(1), mc(1.4,maxbh_of_Z(1)),    mc_av
# sys.exit(0)

#prob_Mbh = make_power_func(1)
#print prob_Mbh(np.linspace(1,30,20))
#mc_av = np.power(average_func(lambda x: np.power(mc(x,x),15./6), prob_Mbh),5./16.)
#print 1, maxbh_of_Z(1), mc_av
#sys.exit(0)

# Get data for three key functions, and interpolate it
if opts.verbose:
    print " Building 'V' factor used later...here is the specific table "
    print "  log10 Z/Zsun     V_bns  V_bhns  V_bhbh "
dat_logZoverZsun = np.linspace(-9,7,400)
dat_nsns = np.zeros(len(dat_logZoverZsun))
dat_bhns=np.zeros(len(dat_logZoverZsun))
dat_bhbh=np.zeros(len(dat_logZoverZsun))
for indx in np.arange(len(dat_nsns)):
    lZ =dat_logZoverZsun[indx]
    v1=dat_nsns[indx] = Vaveraged_bns(np.power(10,lZ))
    v2=dat_bhns[indx] = Vaveraged_bhns(np.power(10,lZ))
    v3=dat_bhbh[indx] = Vaveraged_bhbh(np.power(10,lZ))
    if opts.verbose:
        print lZ, v1,v2,v3

fn_nsns= interp1d(dat_logZoverZsun,dat_nsns)
fn_bhns= interp1d(dat_logZoverZsun,dat_bhns)
fn_bhbh= interp1d(dat_logZoverZsun,dat_bhbh)
#sys.exit(0)

# SFR files:
#    - not guaranteed to have same # per Z bin ! 
if opts.verbose:
    print " File: Loading ", opts.starmetal_file
dat = np.loadtxt(opts.starmetal_file,skiprows=1)
if opts.verbose:
    print " File: Range check on Z ", np.min(dat[:,0]), np.max(dat[:,0])
# format of this file: 3 columns, star metallicity, nearby gas metallicity, formation time of star
lambda0 = 1e-3

Zsun_val_stellar = 0.02    # Problem: galaxy plot has zsun = 0.01, but KB has Zsun=0.02
Zsun_val_JB =0.01

fac_convert = Zsun_val_stellar/Zsun_val_JB

def volume_function():
    """
    volume_function: compute the metallicity-dependent factor proportional to V \propto \mc^(15/6)
    """
    # we will add something real here
    if not opts.type_bbh and not opts.type_bhns:
        print " Volume function: using NSNS"
        return fn_nsns(fac_convert*dat[:,0]/Zsun_val_stellar)
    elif opts.type_bhns:
        print " Volume function: using BHNS"
        return fn_bhns(fac_convert*dat[:,0]/Zsun_val_stellar)
    else:
        print " Volume function: using BHBH"
        return fn_bhbh(fac_convert*dat[:,0]/Zsun_val_stellar)

def metal_function(z_exp):
    # number of merging compact binaries in a starburst depends on metallicity
    return np.minimum(np.power((1e-8+fac_convert*dat[:,0])/Zsun_val_stellar,-1.*z_exp), 1e3*np.ones(len(dat)))

def time_function(time_exp):
    """
    time_function: dP/dt \propto 1/t for t > 0.01 Gyr
    EXCEPT for BH-BH, suppress formation at early times
    """
    # where does this come from?
    standard_delay =  np.where(dat[:,2]>0.01,1./np.power(dat[:,2]+1e-4, time_exp),np.zeros(len(dat)))
    if opts.type_bbh:
        standard_delay =  np.where(dat[:,2]>0.01,1./np.power(dat[:,2]+1e-4, 1),np.zeros(len(dat)))   # start with 1/t, to avodi confusion
        standard_delay *= np.where(fac_convert*dat[:,0]>0.25*Zsun_val_stellar, 0.1*(dat[:,0]/1e4),np.ones(len(dat)))  # suppression factor: convert to delay time which is uniform delay time over 10^4 Gyr, and includes 0.1 of all BHBH
        return standard_delay
    else:
        return standard_delay


if opts.type_bbh:
    opts.Z_exponent = 2
if not opts.type_bbh and not opts.type_bhns:
    opts.Z_exponent =0   # checkme
if opts.type_bhns:
    opts.Z_exponent =0   # checkme


if opts.verbose:
    print " ------- "
if not (opts.Z_exponent is None):
    # print a specific value
    weights = metal_function(opts.Z_exponent) * time_function(opts.time_exponent)*volume_function()
    print opts.starmetal_file, opts.Z_exponent,  len(dat)*ParseStars.jbStarMass/1e10, np.sum(weights)*lambda0*ParseStars.jbStarMass/1e10 #, np.sum(weights)/len(dat)
    hist, bin_edges = np.histogram(dat[:,2], weights=weights,density=True,bins=40)
    bin_centers = (bin_edges[:-1]+bin_edges[1:])/2
    np.savetxt("hist_data.dat", np.array([bin_centers,hist]).T)
    plt.plot(bin_centers, hist)
    plt.xlabel("t (Gyr)")
    plt.ylabel("dN/dt")
    plt.savefig("hist_result.png"); plt.clf()
else:
    # Loop over all from 0 to 3 in steps of 0.1
    for Zexp in np.linspace(0,3, 30):
        weights = metal_function(Zexp)*time_function( opts.time_exponent)*volume_function()
        print opts.starmetal_file, Zexp,  len(dat)*ParseStars.jbStarMass/1e10, np.sum(weights)*lambda0*ParseStars.jbStarMass/1e10, np.sum(weights)/len(dat)
        hist, bin_edges = np.histogram(dat[:,0], weights=weights)
        bin_centers = (bin_edges[:-1]+bin_edges[1:])/2
        plt.plot(bin_centers, hist)
        plt.savefig("hist_result_" + str(Zexp) + ".png")

