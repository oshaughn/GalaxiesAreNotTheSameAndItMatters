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
#    python sfr_to_rate.py --starmetal-file ~/PersonalJBArchive/boring-galaxy/sfrhist_metalbin_boring.txt --plot
#
# PLOTS
#    Third plot is money plot: shows you where the mergers are coming from, for different exponent choices.

import numpy as np

import bisect
from scipy.interpolate import griddata

import ParseStars

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--starmetal-file", default=None)
parser.add_argument("--time-exponent",default=1,type=float)
parser.add_argument("--Z-exponent",default=0,type=float)
opts = parser.parse_args()
if not opts.starmetal_file:
    print "  Please specify an SFR file "


# SFR files:
#    - not guaranteed to have same # per Z bin ! 
dat = np.loadtxt(opts.starmetal_file,skiprows=1)
lambda0 = 1e-3

if opts.Z_exponent:
    # print a specific value
    weights = np.minimum(np.power((1e-8+dat[:,0])/0.02,-1.*opts.Z_exponent), 1e3*np.ones(len(dat))) *1./np.power(0.01+dat[:,2], opts.time_exponent)
    print opts.starmetal_file, opts.Z_exponent,  len(dat)*ParseStars.jbStarMass/1e10, np.sum(weights)*lambda0*ParseStars.jbStarMass/1e10, np.sum(weights)/len(dat)
else:
    # Loop over all from 0 to 3 in steps of 0.1
    for Zexp in np.linspace(0,3, 30):
        weights = np.minimum(np.power((1e-8+dat[:,0])/0.02,-1.*Zexp), 1e3*np.ones(len(dat))) *1./np.power(0.01+dat[:,2], opts.time_exponent)
        print opts.starmetal_file, Zexp,  len(dat)*ParseStars.jbStarMass/1e10, np.sum(weights)*lambda0*ParseStars.jbStarMass/1e10, np.sum(weights)/len(dat)
