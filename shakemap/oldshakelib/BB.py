#!/usr/bin/env python

############################################################################
# Module to apply amplitude and sigma modifications when converting 
# from one peak ground motion definition to another, per Beyer and
# Bommer, 2006, BSSA, 96(4A).
#
# Functions:
#
# param = one of 'pga', 'pgv', 'psa03', 'psa10', 'psa30'
#
# val = getAM2Max_pgm(param)
# val = getRand2Max_pgm(param)
# val = getGM2Max_pgm(param)
#	These functions return a multiplier to convert Arithmetic Mean,
#	Random Component, or Geometric Mean (respectively) ground motion
#	amplitudes to Maximum Component amplitude.
#
# sigma = getAM2Max_sigma(param, sigma)
# sigma = getRand2Max_sigma(param, sigma)
# sigma = getGM2Max_sigma(param, sigma)
#	These functions take a sigma (in linear space) and convert it from
#	Arithmetic Mean, Random Component, or Geometric Mean (respectively)
#	to Maximum Component sigma.
#
############################################################################

import numpy as np

#define some global variables
# Geometric Mean to Arithmetic Mean (& vice versa)
columns = ['c12','c34','R']
gm2amvals = [[1.00,1.00,0.01],
             [1.00,1.00,0.01],
             [1.00,1.00,0.01414],
             [1.00,1.00,0.02],
             [1.00,1.00,0.02]]
GM2AM = {'pga':dict(zip(columns,gm2amvals[0][:])),
         'pgv':dict(zip(columns,gm2amvals[1][:])),
         'psa03':dict(zip(columns,gm2amvals[2][:])),
         'psa10':dict(zip(columns,gm2amvals[3][:])),
         'psa30':dict(zip(columns,gm2amvals[4][:]))}

# Geometric Mean to Random component (& vice versa)
gm2randvals = [[1.00,1.03,0.07],
               [1.00,1.03,0.09],
               [1.00,1.05,0.08656],
               [1.00,1.05,0.11],
               [1.00,1.05,0.11]]
GM2Rand = {'pga':dict(zip(columns,gm2randvals[0][:])),
           'pgv':dict(zip(columns,gm2randvals[1][:])),
           'psa03':dict(zip(columns,gm2randvals[2][:])),
           'psa10':dict(zip(columns,gm2randvals[3][:])),
           'psa30':dict(zip(columns,gm2randvals[4][:]))}

# Geometric Mean to Envelope (i.e., maximum component) (& vice versa)
gm2envvals = [[1.10,    1.02, 0.05],
              [1.15,    1.03, 0.06],
              [1.14140, 1.02, 0.05242],
              [1.20,    1.02, 0.07],
              [1.20,    1.02, 0.07]]
GM2Env = {'pga':dict(zip(columns,gm2envvals[0][:])),
          'pgv':dict(zip(columns,gm2envvals[1][:])),
          'psa03':dict(zip(columns,gm2envvals[2][:])),
          'psa10':dict(zip(columns,gm2envvals[3][:])),
          'psa30':dict(zip(columns,gm2envvals[4][:]))}


def getAM2Max_pgm(param):
    return getBBpgm(param,GM2Env) / getBBpgm(param,GM2AM)

def getRand2Max_pgm(param):
    return getBBpgm(param,GM2Env) / getBBpgm(param,GM2Rand)

def getGM2Max_pgm(param):
    return getBBpgm(param,GM2Env);

def getAM2Max_sigma(param, sigma):
    sigma = np.log10(sigma)
    sigma = inverseBBsigma(param, GM2AM, sigma)
    sigma = applyBBsigma(param, GM2Env, sigma)
    return np.power(10,sigma)

def getRand2Max_sigma(param, sigma):
    sigma = np.log10(sigma)
    sigma = inverseBBsigma(param, GM2Rand, sigma)
    sigma = applyBBsigma(param, GM2Env, sigma)
    return np.power(10,sigma)

def getGM2Max_sigma(param, sigma):
    sigma = np.log10(sigma)
    sigma = applyBBsigma(param, GM2Env, sigma)
    return np.power(10,sigma)


####################################################################
# Internal routines
####################################################################

# use vars qw ( $GM2AM $GM2Rand $GM2Env );
# use enum qw( c12 R c34 );

def getBBpgm(param, type):
  return type[param]['c12']

def applyBBsigma(param, type, value):
    R = type[param]['R']
    c34 = type[param]['c34']
    return np.sqrt(value * value * R * R + c34 * c34)

def inverseBBsigma(param, type, value):
    R = type[param]['R']
    c34 = type[param]['c34']
    return np.sqrt((value * value - c34 * c34) / (R * R))
