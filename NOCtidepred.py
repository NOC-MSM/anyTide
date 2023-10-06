#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
NOCtidepred.py

Harmonic reconstruction.
Reconstructs for a vector of times.
Reconstructs for either an array of locations or a point.

Based on matlab code (NOCtidepred.m) from Simon Williams: sdwil@noc.ac.uk

Author: jelt@noc.ac.uk
Date: 13 Jul 2018
"""
# Clean
# For point reconstruction
import numpy as np
import math # atan2
import pandas as pd
import datetime
from os import path # fix the path
import matplotlib.pyplot as plt

#import anyTide_Cwrapper

# Additional for map reconstruction
from netCDF4 import Dataset
import sys
sys.path.insert(0, "../DEV_jelt/NEMO_diag/IT/")
#import internaltideharmonics_NEMO as ITh
from AMM60_tools import findJI
from AMM60_tools import harmonictable

##########################################################################
## Functions
##########################################################################

def set_names_phases():

    # The following phase speeds and names are indexed by what I'm called the
    # Doodson number, which start at one: Doodson(SA)=1.
    # Note that for python indexing these are called with indexing starting from
    # zero so that names[ doodson(SA)-1 ] = SA
    sig0 = [0.0410686,   0.0821373,   0.5443747,   1.0158958,
            1.0980330,  12.8542862,  12.9271398,  13.3986609,
           13.4715145,  13.9430356,  14.0251729,  14.4920521,
           14.5695476,  14.9178647,  14.9589314,  15.0000000,
           15.0410686,  15.0821353,  15.1232059,  15.5125897,
           15.5854433,  16.0569644,  16.1391017,  27.3416965,
           27.4238338,  27.8953548,  27.9682085,  28.4397295,
           28.5125832,  28.9019670,  28.9841042,  29.0662415,
           29.4556253,  29.5284789,  29.9589333,  30.0000000,
           30.0410667,  30.0821373,  30.5443747,  30.6265120,
           31.0158958,  42.9271398,  43.4761564,  43.9430356,
           44.0251729,  45.0410686,  57.4238338,  57.9682085,
           58.4397295,  58.9841042,  59.0662415,  60.0000000,
           60.0821373,  86.4079380,  86.9523127,  87.4238338,
           87.9682085,  88.0503458,  88.9841042,  89.0662415,
           26.4079380,  26.8701754,  26.9523127,  27.5059711,
           28.3575923,  29.9178627,  31.0887494,  42.3827651,
           43.0092771,  44.5695476,  56.8701754,  56.9523127,
           57.8860712,  71.9112441,  72.4602606,  73.0092771,
           84.8476675,  85.3920423,  85.8542797,  85.9364170,
           86.3258007,  86.4807917,  86.8701754,  87.4966874,
           88.5125832,  88.5947205, 114.8476675, 115.3920423,
          115.9364170, 116.4079380, 116.9523127, 117.0344500,
          117.5059711, 117.9682085, 118.0503458, 145.9364170,
          146.9523127, 174.3761465, 174.9205212, 175.9364170,
           27.4966874,  27.8860712,  28.9430356,  29.0251729,
           30.4715211,  31.0980330,  56.4079380,  57.4966874,
           58.5125832,  59.5284789,  28.3986629,  28.4807962,
           72.9271398,  74.0251729,  29.5284789,   0.0000000,
            0.0000000,   0.0000000,   0.0000000,   0.0000000]

    names = ['SA',   'SSA', 'MM', 'MSF', 'MF', '2Q1', 'SIG1', 'Q1', 'RO1',
             'O1',   'MP1', 'M1', 'CHI1', 'PI1', 'P1', 'S1', 'K1', 'PSI1', 'PHI1',
             'TH1',  'J1', 'SO1', 'OO1', 'OQ2', 'MNS2', '2N2', 'MU2', 'N2', 'NU2',
             'OP2',  'M2', 'MKS2', 'LAM2', 'L2', 'T2', 'S2', 'R2', 'K2', 'MSN2',
             'KJ2',  '2SM2', 'MO3', 'M3', 'SO3', 'MK3', 'SK3', 'MN4', 'M4', 'SN4',
             'MS4',  'MK4', 'S4', 'SK4', '2MN6', 'M6', 'MSN6', '2MS6', '2MK6', '2SM6',
             'MSK6', '2MN2S2', '3MSK2', '3M2S2', 'MNK2S2', 'SNK2', '2SK2', '2MS2N2',
             'MQ3',  '2MP3', '2MQ3', '3MK4', '3MS4', '2MSK4', '3MK5', 'M5', '3MO5',
             '2MNS6','3MNS6', '4MK6', '4MS6', '2MSNK6', '2MV6', '3MSK6', '4MN6', '3MSN6',
             'MKL6', '2MN8', '3MN8', 'M8', '2MSN8', '3MS8', '3MK8', 'MSNK8', '2MS8',
             '2MSK8','4MS10', '3M2S10', '4MSN12', '5MS12', '4M2S12', 'MVS2', '2MK2',
             'MA2',  'MB2', 'MSV2', 'SKM2', '2MNS4', 'MV4', '3MN4', '2MSN4', 'NA2',
             'NB2',  'MSO5', 'MSK5', '2MN2']

    return names, np.array(sig0) # transform sig0 from a 'list' into an 'array'

##########################################################################
def phamp0fast(mjdn):
#
# Calculates nodal amplitude factors f, and phases (including nodal
# corrections) for the standard list of constituents at time 00:00
# on day mjdn
#
# Input: mjdn (integer] = modified julian day number
#
# Output:
#
# f (double array length ncmax] = nodal amplitude factors
# v (double array length ncmax] = phases (degrees) including nodal corrections
# optionally return u, nodal correction to phase
#
# calculate ecliptic mean longitudes of moon (s), sun(h), lunar perigee (p),
# lunar ascending node (en) and perihelion (p1) at 00:00 on mjdn

      [s,h,p,en,p1] = longfindfast(mjdn)

# calculate nodal amplitude factors and phase corrections at 00:00 on mjdn

      [u,f] =  ufsetfast(p,en)

# calculate mean phases at 00:00 on mjdn

      v = vsetfast(s,h,p,p1)

# sum phases

      v = v + u

      return f,v

#########################################################################
def longfindfast(mjdn):
# =================================================================
#
#     Calculate astronomical arguments for tides
#     This version from Chris Hughes 28 June 2002
#
# =================================================================
#
# Input: mjdn = modified julian day number
#        iout = integer unit for output (0 produces no output)
# Output: s = mean longitude of Moon
#         h = mean longitude of Sun
#         p = mean longitude of lunar perigee
#        EN = mean longitude of lunar node (point where plane of
#             lunar orbit crosses plane of equator)
#        p1 = mean longitude of perihelion
#
#        All in degrees,  calculated at 00:00 on that day.
#
#
# Time is calculated internally from zero at 12:00 1-Jan-2000 TDT
# Longitudes at t=0 and formulae for other times taken from
# Cartwright, D.E., 1985: Tidal Prediction and Modern Time Scales,
# International Hydrographic Review, Monaco, 62(1), 127-138.
#
# The formulae are for time (TT) in (atomic clock) seconds, converted to
# Julian centuries (86400*365.25*100 s = 1 Julian century). Due to
# variations in earth rotation rate, this differs from civil time by
# an amount which increases with distance from the reference time - by
# less than a minute within a century, increasing to about 6 hours at
# 750 BC. 6 hours in the position of the moon translates to an M2 phase
# shift of about 6 degrees. They are part of a Taylor series and will
# drift from reality at long time scales (centuries).
#
# To account for this, two time variables are used:
#
# It is assumed that the measurement time of day is t = civil time (UTC, less than
#    a second from UT1 which is actually a measure of earth rotation)
#
# tt = Terrestrial Time (the scale used to derive the position of the sun
#    and moon) is the actual (atomic) time elapsed. Units of tt are Julian centuries
#    (36525*86400 seconds) measured relative to 12:00 1 Jan 2000 (J2000)
#
# A quadratic correction is applied to convert from t to tt: deltat = tt - t,
#   where deltat = a + b*tt + c*tt^2 (in practice using t instead of tt)
#
# c is 31.0 (taken from Ch 14 of Stephenson: Historical eclipses and earth rotation)
#   this corresponds to a lengthening of length of day at a rate 1.7 ms per century
#   (2*36525*31 = 0.001697)
#
# b is 90.0, chosen to ensure a good fit to measured deltat since 1900 and
#   reasonable fit to earlier astronomical observations and eclipse data going
#   back to 750 BC. This gives an excellent fit for 1940-2000 (2-3 seconds)
#   but the error rises to ~10 s earlier in 20th century and by late 2015.
#   Errors remain within about 30 minutes for times back to 750 BC.
#
# a is chosen by insisting that deltat = 32.184 s at t0 = 00:00 1 Jan 1958,
#   since that was the actual value of deltat when TAI was set equal to UT2
#   (an approximation of UT1 with seasonal terms removed). Since
#   mjd of this time is 36204.0, and mjd of J2000 is 51544.5d0, this gives
#   t0 = (36204.0d0-51544.5d0)/36525.d0 in Julian centuries, and hence
#   a = 32.184d0 - b*t0 - c*t0**2
#
# The size of errors here (seconds to minutes) should only be interpreted
# in terms of the distance moved by the sun and moon relative to the stars
# over that period. The biggest issue is the moon, which moves at about 0.55
# degrees per hour.
#
#
    cycle=360.0
    c = 32.0
    b = 90.0
    t0 = (36204.0-51544.5)/36525.0
    a = 32.184 - b * t0 - c * t0**2
#
# calculate t (julian centuries UTC after 12:00, 1 Jan 2000)
# 51544 = mjdn of 1 Jan 2000
#

    t = (mjdn-51544-0.5)/36525

#
# calculate correction to t for earth rotation variations and
# difference between TDT and TAI (TDT-TAI=32.184s)
#     (dt in seconds, set dt=32.284s at 00:00 1-Jan-1958)
#
#     mjdn for that date is 36204
#
#      t0 = (36204.0d0-51544.5d0)/36525.d0
#      c = 31.0d0
#      b = 90.0d0
#      a = 32.184d0 - b*t0 - c*t0**2

    dt = a +  b*t  +  c*t**2

    tt = t + dt/(86400.0*36525.0)

    s   = 218.3166 + 481267.8811*tt - 0.0019*tt**2
    h   = 280.4661 +  36000.7698*tt + 0.0003*tt**2
    p   =  83.3532 +   4069.0136*tt - 0.0106*tt**2
    EN  = 125.0445 -   1934.1364*tt + 0.0018*tt**2
    p1  = 282.9384 +      1.7194*tt + 0.0002*tt**2

    s  = np.mod( s,cycle)
    h  = np.mod( h,cycle)
    p  = np.mod( p,cycle)
    EN = np.mod(EN,cycle)
    p1 = np.mod(p1,cycle)

    s[ s < 0] += cycle
    h[ h < 0] += cycle
    p[ p < 0] += cycle
    EN[EN < 0] += cycle
    p1[p1 < 0] += cycle

    return s,h,p,EN,p1

##########################################################################
def ufsetfast(p,en):
#
# Computes nodal adjustment factors for amplitude (f) and phase
# (u, degrees) given ecliptic longitudes of lunar perigee (p) and
# lunar ascending node (en) in degrees.
#
# Input:
#
#  p (double) - ecliptic longitude of lunar perigee (degrees)
# en (double) - ecliptic longitude of lunar ascending node (degrees)
#
# Output:
#
# u (double array length ncmax) - list of nodal phase adustments (degrees)
#                                 for the standard list of constituents
# f (double array length ncmax)  - list of nodal amplitude factors
#                                 for the standard list of constituents
#
    rad = np.pi/180.
    deg = 180.0/np.pi

    pw = p*rad
    nw = en*rad

    w1 = np.cos(nw)
    w2 = np.cos(2.0*nw)
    w3 = np.cos(3.0*nw)
    w4 = np.sin(nw)
    w5 = np.sin(2.0*nw)
    w6 = np.sin(3.0*nw)

    a1 = pw-nw
    a2 = 2.0*pw
    a3 = a2-nw
    a4 = a2-2.0*nw

# *** u's are computed in radians
# jelt: note matlab indexing starts from 1. Python starts from 0
# I.e. indices offset: matlab(i) --> python(i-1)
# To keep the code tracable fill [f,u] with indices 1:nconstituents (as in matlab).
# Then adjust as the end of the function.

    u = np.zeros((len(p),120+1)) # Add extra (121st) slot so indexing goes to 120. I.e. 0,1,..,120
    f = np.zeros((len(p),120+1))

    u[:,3]   =  0.0
    f[:,3]   =  1.0    -0.1300*w1 +0.0013*w2

    u[:,5]   =           -0.4143*w4 +0.0468*w5 -0.0066*w6
    f[:,5]   =  1.0429 +0.4135*w1 -0.004*w2

    u[:,10]  =            0.1885*w4 -0.0234*w5 +0.0033*w6
    f[:,10]  =  1.0089 +0.1871*w1 -0.0147*w2 +0.0014*w3

    x = 2.0*np.cos(pw)+0.4*np.cos(a1)
    y = np.sin(pw)+0.2*np.sin(a1)

    u[:,12]  = [math.atan2(y[i],x[i]) for i in range(len(p))]
# atan2 replaces need to locate quadrant as below
#      u(:,12]  = atan (y./x);
#      I = x < 0.0;
#      u(I,12] = u(I,12) + pi;

    f[:,12]  = np.sqrt(x**2 + y**2)

    u[:,17]  =           -0.1546*w4 +0.0119*w5 -0.0012*w6
    f[:,17]  =  1.0060 +0.1150*w1 -0.0088*w2 +0.0006*w3

    u[:,21]  =           -0.2258*w4 +0.0234*w5 -0.0033*w6
    f[:,21]  =  1.0129 +0.1676*w1 -0.0170*w2 +0.0016*w3

    f[:,23]  =  1.1027 +0.6504*w1 +0.0317*w2 -0.0014*w3
    u[:,23]  =           -0.6402*w4 +0.0702*w5 -0.0099*w6

    u[:,31]  =           -0.0374*w4
    f[:,31]  =  1.0004 -0.0373*w1 +0.0002*w2

    x  = 1.0-0.2505*np.cos(a2)-0.1102*np.cos(a3)-0.0156*np.cos(a4)-0.037*w1
    y  =    -0.2505*np.sin(a2)-0.1102*np.sin(a3)-0.0156*np.sin(a4)-0.037*w4

    u[:,34]  = [math.atan2(y[i],x[i]) for i in range(len(p))]
# atan2 replaces need to locate quadrant as below
#            u(:,34]  = atan (y./x);
#            I = x < 0.0;
#            u(I,34] = u(I,34) + pi;
    f[:,34]  = np.sqrt(x**2 + y**2)

    u[:,38]  =           -0.3096*w4 +0.0119*w5 -0.0007*w6
    f[:,38]  =  1.0241   +0.2863*w1 +0.0083*w2 -0.0015*w3

    u[:,1]   =  0.0
    u[:,2]   =  0.0
    u[:,4]   = -u[:,31]
    u[:,6]   =  u[:,10]
    u[:,7]   =  u[:,10]
    u[:,8]   =  u[:,10]
    u[:,9]   =  u[:,10]
    u[:,11]  =  u[:,31]
    u[:,13]  =  u[:,21]
    u[:,14]  =  0.0
    u[:,15]  =  0.0
    u[:,16]  =  0.0
    u[:,18]  =  0.0
    u[:,19]  =  0.0
    u[:,20]  =  u[:,21]
    u[:,22]  = -u[:,10]
    u[:,24]  =  2.0*u[:,10]
    u[:,25]  =  2.0*u[:,31]
    u[:,26]  =  u[:,31]
    u[:,27]  =  u[:,31]
    u[:,28]  =  u[:,31]
    u[:,29]  =  u[:,31]
    u[:,30]  =  u[:,10]
    u[:,32]  =  u[:,31]+u[:,38]
    u[:,33]  =  u[:,31]
    u[:,35]  =  0.0
    u[:,36]  =  0.0
    u[:,37]  =  0.0
    u[:,39]  =  0.0
    u[:,40]  =  u[:,17]+u[:,21]
    u[:,41]  =  u[:,4]
    u[:,42]  =  u[:,31]+u[:,10]
    u[:,43]  =  u[:,31]*1.5
    u[:,44]  =  u[:,10]
    u[:,45]  =  u[:,31]+u[:,17]
    u[:,46]  =  u[:,17]
    u[:,47]  =  u[:,25]
    u[:,48]  =  u[:,25]
    u[:,49]  =  u[:,31]
    u[:,50]  =  u[:,31]
    u[:,51]  =  u[:,32]
    u[:,52]  =  0.0
    u[:,53]  =  u[:,38]
    u[:,54]  =  u[:,25]+u[:,31]
    u[:,55]  =  u[:,54]
    u[:,56]  =  u[:,25]
    u[:,57]  =  u[:,25]
    u[:,58]  =  u[:,25]+u[:,38]
    u[:,59]  =  u[:,31]
    u[:,60]  =  u[:,32]
    u[:,61]  =  0.0
    u[:,62]  =  u[:,54]-u[:,38]
    u[:,63]  =  u[:,54]
    u[:,64]  =  u[:,58]
    u[:,65]  =  u[:,31]-u[:,38]
    u[:,66]  = -u[:,38]
    u[:,67]  =  0.0
    u[:,68]  =  u[:,42]
    u[:,69]  =  u[:,25]
    u[:,70]  =  u[:,25]-u[:,10]
    u[:,71]  =  u[:,54]-u[:,38]
    u[:,72]  =  u[:,54]
    u[:,73]  =  u[:,25]-u[:,38]
    u[:,74]  =  u[:,54]-u[:,17]
    u[:,75]  =  2.5*u[:,31]
    u[:,76]  =  u[:,54]-u[:,10]
    u[:,77]  =  2.0*u[:,25]
    u[:,78]  =  u[:,77]
    u[:,79]  =  u[:,77]-u[:,38]
    u[:,80]  =  u[:,77]
    u[:,81]  =  u[:,71]
    u[:,82]  =  u[:,54]
    u[:,83]  =  u[:,71]
    u[:,84]  =  u[:,54]
    u[:,85]  =  u[:,25]
    u[:,86]  =  u[:,51]+u[:,34]
    u[:,87]  =  u[:,77]
    u[:,88]  =  u[:,77]
    u[:,89]  =  u[:,77]
    u[:,90]  =  u[:,54]
    u[:,91]  =  u[:,54]
    u[:,92]  =  u[:,54]+u[:,38]
    u[:,93]  =  u[:,58]
    u[:,94]  =  u[:,25]
    u[:,95]  =  u[:,58]
    u[:,96]  =  u[:,77]
    u[:,97]  =  u[:,54]
    u[:,98]  =  5.0*u[:,31]
    u[:,99]  =  u[:,98]
    u[:,100] =  u[:,77]
    u[:,101] =  u[:,25]
    u[:,102] =  u[:,73]
    """
    % *** u(103),u(104) are changed according to
    % *** dr cartwright's notes of nov 15,1977
    % to reflect the fact that annual modulation of tides is due to radiation
    % rather than gravitational forcing
    """
    u[:,103] =  0.0
    u[:,104] =  0.0
    u[:,105] =  0.0
    u[:,106] = -u[:,65]
    u[:,107] =  u[:,54]
    u[:,108] =  u[:,25]
    u[:,109] =  0.0
    u[:,110] =  u[:,31]
    u[:,111] =  0.0
    u[:,112] =  0.0
    u[:,113] =  u[:,42]
    u[:,114] =  u[:,45]
    u[:,115] =  u[:,31]
    u[:,116] =  0.0
    u[:,117] =  0.0
    u[:,118] =  0.0
    u[:,119] =  0.0
    u[:,120] =  0.0

# *** convert into degrees

    u = np.mod(u*deg,360)
    u[u < 0.0] += 360

    f[:,1]   =  1.0
    f[:,2]   =  1.0
    f[:,4]   =  f[:,31]
    f[:,6]   =  f[:,10]
    f[:,7]   =  f[:,10]
    f[:,8]   =  f[:,10]
    f[:,9]   =  f[:,10]
    f[:,11]  =  f[:,31]
    f[:,13]  =  f[:,21]
    f[:,14]  =  1.0
    f[:,15]  =  1.0
    f[:,16]  =  1.0
    f[:,18]  =  1.0
    f[:,19]  =  1.0
    f[:,20]  =  f[:,21]
    f[:,22]  =  f[:,10]
    f[:,24]  =  f[:,10]**2
    f[:,25]  =  f[:,31]**2
    f[:,26]  =  f[:,31]
    f[:,27]  =  f[:,31]
    f[:,28]  =  f[:,31]
    f[:,29]  =  f[:,31]
    f[:,30]  =  f[:,10]
    f[:,32]  =  f[:,31]*f[:,38]
    f[:,33]  =  f[:,31]
    f[:,35]  =  1.0
    f[:,36]  =  1.0
    f[:,37]  =  1.0
    f[:,39]  =  f[:,25]
    f[:,40]  =  f[:,17]*f[:,21]
    f[:,41]  =  f[:,31]
    f[:,42]  =  f[:,31]*f[:,10]
    f[:,43]  =  f[:,31]**1.5
    f[:,44]  =  f[:,10]
    f[:,45]  =  f[:,31]*f[:,17]
    f[:,46]  =  f[:,17]
    f[:,47]  =  f[:,25]
    f[:,48]  =  f[:,25]
    f[:,49]  =  f[:,31]
    f[:,50]  =  f[:,31]
    f[:,51]  =  f[:,32]
    f[:,52]  =  1.0
    f[:,53]  =  f[:,38]
    f[:,54]  =  f[:,25]*f[:,31]
    f[:,55]  =  f[:,54]
    f[:,56]  =  f[:,25]
    f[:,57]  =  f[:,25]
    f[:,58]  =  f[:,25]*f[:,38]
    f[:,59]  =  f[:,31]
    f[:,60]  =  f[:,32]
    f[:,61]  =  1.0
    f[:,62]  =  f[:,54]*f[:,38]
    f[:,63]  =  f[:,54]
    f[:,64]  =  f[:,58]
    f[:,65]  =  f[:,32]
    f[:,66]  =  f[:,38]
    f[:,67]  =  f[:,25]*f[:,25]
    f[:,68]  =  f[:,42]
    f[:,69]  =  f[:,25]
    f[:,70]  =  f[:,25]*f[:,10]
    f[:,71]  =  f[:,54]*f[:,38]
    f[:,72]  =  f[:,54]
    f[:,73]  =  f[:,58]
    f[:,74]  =  f[:,54]*f[:,17]
    f[:,75]  =  0.5*(f[:,25]+f[:,54])
    f[:,76]  =  f[:,54]*f[:,10]
    f[:,77]  =  f[:,67]
    f[:,78]  =  f[:,67]
    f[:,79]  =  f[:,67]*f[:,38]
    f[:,80]  =  f[:,67]
    f[:,81]  =  f[:,71]
    f[:,82]  =  f[:,54]
    f[:,83]  =  f[:,71]
    f[:,84]  =  f[:,54]*f[:,25]
    f[:,85]  =  f[:,67]
    f[:,86]  =  f[:,51]*f[:,34]
    f[:,87]  =  f[:,67]
    f[:,88]  =  f[:,67]
    f[:,89]  =  f[:,67]
    f[:,90]  =  f[:,54]
    f[:,91]  =  f[:,54]
    f[:,92]  =  f[:,71]
    f[:,93]  =  f[:,58]
    f[:,94]  =  f[:,25]
    f[:,95]  =  f[:,58]
    f[:,96]  =  f[:,67]
    f[:,97]  =  f[:,54]
    f[:,98]  =  f[:,84]
    f[:,99]  =  f[:,84]
    f[:,100] =  f[:,67]
    f[:,101] =  f[:,25]
    f[:,102] =  f[:,58]
    """
    % *** f[:,103],f[104] are changed according to
    % *** dr cartwright's notes of nov 15,1977,
    % to reflect the fact that annual modulation of tides is due to radiation
    % rather than gravitational forcing
    """
    f[:,103] =  1.0
    f[:,104] =  1.0
    f[:,105] =  1.0
    f[:,106] =  f[:,32]
    f[:,107] =  f[:,54]
    f[:,108] =  f[:,25]
    f[:,109] =  1.0
    f[:,110] =  f[:,54]
    f[:,111] =  1.0
    f[:,112] =  1.0
    f[:,113] =  f[:,42]
    f[:,114] =  f[:,45]
    f[:,115] =  f[:,54]
    f[:,116] =  0.0
    f[:,117] =  0.0
    f[:,118] =  0.0
    f[:,119] =  0.0
    f[:,120] =  0.0

# Now shift the indexing from 1,..,nconstituents to 0,..,nconstituents-1
# and drop the last index
    u = np.roll(u, -1, axis=1)[:,:-1]
    f = np.roll(f, -1, axis=1)[:,:-1]
    return u,f
##########################################################################
def vsetfast(s,h,p,p1):
    """
    %
    % Computes mean phase (v, degrees), ignoring nodal adjustment factors,
    % for the standard list of tidal constituents at 00:00, given ecliptic
    % longitudes at that time of the moon (s), sun (h), lunar perigee (p),
    % and perihelion (p1)
    %
    % Input:
    %
    %  s (double) - ecliptic longitude of moon (degrees)
    %  h (double) - ecliptic longitude of sun (degrees)
    %  p (double) - ecliptic longitude of lunar perigee (degrees)
    % p1 (double) - ecliptic longitude of perihelion (degrees)
    %
    % Output:
    %
    % v (double array length ncmax) - list of phases at 00:00 (degrees)
    %                                 for the standard list of constituents
    %

    %     integer ncmax,k
    %     parameter (ncmax=120)
    %     real*8 s,p,p1,h,v(ncmax),h2,h3,h4,p2,s2,s3,s4
    """
    h2 = h+h
    h3 = h2+h
    h4 = h3+h

    s2 = s+s
    s3 = s2+s
    s4 = s3+s

    p2 = p+p

# *** v's computed in degrees.

    v = np.zeros((len(p),120+1)) # Extra column to allow for indexing: 0,..,120

    v[:,1]   =  h
    v[:,2]   =  h2
    v[:,3]   =  s-p
    v[:,4]   =  s2-h2
    v[:,5]   =  s2
    v[:,6]   =  h-s4+p2+270.0
    v[:,7]   =  h3-s4+270.0
    v[:,8]   =  h-s3+p+270.0
    v[:,9]   =  h3-s3-p+270.0
    v[:,10]  =  h-s2+270.0
    v[:,11]  =  h3-s2+90.0
    v[:,12]  =  h-s+90.0
    v[:,13]  =  h3-s-p+90.0
    v[:,14]  =  p1-h2+270.0
    v[:,15]  =  270.0-h
    v[:,16]  =  180.0
    v[:,17]  =  h+90.0
    v[:,18]  =  h2-p1+90.0
    v[:,19]  =  h3+90.0
    v[:,20]  =  s-h+p+90.0
    v[:,21]  =  s+h-p+90.0
    v[:,23]  =  s2+h+90.0
    v[:,26]  =  h2-s4+p2
    v[:,27]  =  h4-s4
    v[:,28]  =  h2-s3+p
    v[:,29]  =  h4-s3-p
    v[:,31]  =  h2-s2
    v[:,32]  =  h4-s2
    v[:,33]  =  p-s+180.0
    v[:,34]  =  h2-s-p+180.0
    v[:,35]  =  p1-h
    v[:,36]  =  0.0
    v[:,37]  =  h-p1+180.0
    v[:,38]  =  h2

    v[:,22]  = -v[:,10]
    v[:,24]  =  v[:,10]+v[:,8]
    v[:,25]  =  v[:,31]+v[:,28]
    v[:,30]  =  v[:,10]+v[:,15]
    v[:,39]  =  v[:,31]-v[:,28]
    v[:,40]  =  v[:,17]+v[:,21]
    v[:,41]  = -v[:,31]
    v[:,42]  =  v[:,31]+v[:,10]
    v[:,43]  =  h3-s3+180.0
    v[:,44]  =  v[:,10]
    v[:,45]  =  v[:,31]+v[:,17]
    v[:,46]  =  v[:,17]
    v[:,47]  =  v[:,25]
    v[:,48]  =  v[:,31]+v[:,31]
    v[:,49]  =  v[:,28]
    v[:,50]  =  v[:,31]
    v[:,51]  =  v[:,31]+v[:,38]
    v[:,52]  =  0.0
    v[:,53]  =  v[:,38]
    v[:,54]  =  v[:,48]+v[:,28]
    v[:,55]  =  v[:,48]+v[:,31]
    v[:,56]  =  v[:,47]
    v[:,57]  =  v[:,48]
    v[:,58]  =  v[:,48]+v[:,38]
    v[:,59]  =  v[:,31]
    v[:,60]  =  v[:,51]
    v[:,61]  =  v[:,54]
    v[:,62]  =  v[:,55]-v[:,38]
    v[:,63]  =  v[:,55]
    v[:,64]  =  v[:,25]+v[:,38]
    v[:,65]  =  v[:,28]-v[:,38]
    v[:,66]  =  -v[:,38]
    v[:,67]  =  v[:,48]-v[:,28]-v[:,28]
    v[:,68]  =  v[:,31]+v[:,8]
    v[:,69]  =  v[:,48]-v[:,15]
    v[:,70]  =  v[:,48]-v[:,8]
    v[:,71]  =  v[:,62]
    v[:,72]  =  v[:,55]
    v[:,73]  =  v[:,48]-v[:,38]
    v[:,74]  =  v[:,55]-v[:,17]
    v[:,75]  =  v[:,31]+v[:,43]
    v[:,76]  =  v[:,55]-v[:,10]
    v[:,77]  =  v[:,54]+v[:,28]
    v[:,78]  =  v[:,55]+v[:,28]

    v[:,89]  =  v[:,48]+v[:,48]

    v[:,79]  =  v[:,89]-v[:,38]
    v[:,80]  =  v[:,89]
    v[:,81]  =  v[:,54]-v[:,38]
    v[:,82]  =  v[:,48]+v[:,29]
    v[:,83]  =  v[:,62]
    v[:,84]  =  v[:,89]-v[:,28]
    v[:,85]  =  v[:,55]-v[:,28]
    v[:,86]  =  v[:,51]+v[:,34]
    v[:,87]  =  v[:,77]
    v[:,88]  =  v[:,78]

    v[:,90]  =  v[:,54]
    v[:,91]  =  v[:,55]
    v[:,92]  =  v[:,55]+v[:,38]
    v[:,93]  =  v[:,64]
    v[:,94]  =  v[:,48]
    v[:,95]  =  v[:,58]
    v[:,96]  =  v[:,89]
    v[:,97]  =  v[:,55]
    v[:,98]  =  v[:,89]+v[:,28]
    v[:,99]  =  v[:,89]+v[:,31]
    v[:,100] =  v[:,89]
    v[:,101] =  v[:,31]+v[:,29]
    v[:,102] =  v[:,73]

# *** v[:,103],v[:,104] are changed according to dr cartwrights
# *** notes of 15 nov,1977

    v[:,103] =  v[:,31]-h
    v[:,104] =  v[:,31]+h
    v[:,105] =  v[:,31]-v[:,29]
    v[:,106] =  v[:,38]-v[:,31]
    v[:,107] =  v[:,54]
    v[:,108] =  v[:,101]
    v[:,109] =  v[:,85]
    v[:,110] =  v[:,34]
    v[:,111] =  v[:,28]+v[:,35]
    v[:,112] =  v[:,28]-v[:,35]
    v[:,113] =  v[:,42]-270.0
    v[:,114] =  v[:,45]-90.0
#%     v[:,115] =  v[:,48]-v[:,28]
    v[:,115] =  v[:,34]
    v[:,116] =  0.0
    v[:,117] =  0.0
    v[:,118] =  0.0
    v[:,119] =  0.0
    v[:,120] =  0.0

    v = np.mod(v,360)
    v[v < 0.0] += 360.0

# Now shift the indexing from 1,..,nconstituents to 0,..,nconstituents-1
# and drop the last index
    v = np.roll(v, -1, axis=1)[:,:-1]
    return v

##########################################################################
def get_port(site_file_name="glad.txt"):
    """
    Obtain a reduced number of Gladstone harmonics from a file
    Returns: lat, lon [floats]
        z0 (ref height) [float]
        data (amp, pha, doo, labels) [pandas dataframe]

    $ head glad.txt

    Port Name: ENGLAND, WEST COAST $ LIVERPOOL (GLADSTONE DOCK)
    53 27.0 N 03 01.1 W
    z0=  5.249   OD= -4.930
    3.03800 320.72000   31 M2
    0.97800   4.70000   36 S2
    ...
    """

    fname = path.abspath(path.join(path.dirname(__file__), "anyTide_Cwrapper", site_file_name))
    if site_file_name == "glad.txt":
        data = pd.read_csv(fname, header=2, names=['amp', 'pha', 'doo', 'lab'], delimiter=r"\s+")
        lat = +(53 + 27.0 / 60)
        lon = -(3 + 1.1 / 60)
        z0 = +5.249

    elif site_file_name == "gladstone.txt":
        data = pd.read_csv(fname, header=2, names=['amp', 'pha', 'doo', 'lab'], delimiter=r"\s+")
        lat = +(53+27.0/60)
        lon = -( 3+ 1.1/60)
        z0  = +5.249

    elif site_file_name == "lowerlargo.txt":
        data = pd.read_csv(fname, header=2, names=['amp', 'pha', 'doo', 'lab'], delimiter=r"\s+")
        lat = 56.2
        lon = -2.925
        z0  = +3.1475328

    elif site_file_name == "dartmouth.txt":
        data = pd.read_csv(fname, header=3, delimiter=r"\s+")
        data.columns = ['amp','pha', 'doo']
        data["lab"] = "empty"
        lat = +(50+21./60)
        lon = -(3+34./60)
        z0  = +2.930

    elif site_file_name == "exmouth.txt":
        data = pd.read_csv(fname, header=3, delimiter=r"\s+")
        data.columns = ['amp','pha', 'doo']
        data["lab"] = "empty"
        lat = 50.6167
        lon = -1.830
        z0  = +2.106

    elif site_file_name == "poole.txt":
        data = pd.read_csv(fname, header=3, delimiter=r"\s+")
        data.columns = ['amp','pha', 'doo']
        data["lab"] = "empty"
        lat =+(50+43./60)
        lon = -(1+59./60)
        z0  = +1.50

    elif site_file_name == "weymouth.txt":
        data = pd.read_csv(fname, header=3, delimiter=r"\s+")
        data.columns = ['amp','pha', 'doo', 'label']
        data["lab"] = "empty"
        lat = +(50+36.5/60)
        lon = -( 2+ 26.9/60)
        z0  = +1.166


    else:
        print(f'Do not recognise file: {site_file_name}')

    return lat,lon, z0, data
##########################################################################
def get_coord_indices(ycoords,xcoords,lats,lons):
	# Find indices for specified coordinates
	[J_ll,I_ll] = findJI(min(ycoords), min(xcoords), lats, lons)  # Simple routine to find the nearest J,I coordinates for given lat lon
	[J_ur,I_ur] = findJI(max(ycoords), max(xcoords), lats, lons)  # Simple routine to find the nearest J,I coordinates for given lat lon
	J1 = min(J_ll,J_ur)
	J2 = max(J_ll,J_ur)+1
	I1 = min(I_ll,I_ur)
	I2 = max(I_ll,I_ur)+1

	return J1,J2,I1,I2
##########################################################################
def get_harmonic_arr(varstr='SSH',xcoords=[-3.1, -3.1],ycoords=[53.5, 53.5], coordsType='deg', dirname='/projectsa/pycnmix/jelt/AMM60/',filebase='AMM60_1d_20120801_20120831'):
	"""
	Get gridded harmonics and coordinate data.
	Get associated harmonic constituent labels and doodson numbers
	INPUT:
	varstr - name of harmonic variable to extract. STRING
	coordsType = 'deg' - input xcoords and ycoords are latitude and longitude
		e.g.:
		ycoords = [49.5, 51]; xcoords = [-3, 2] # Slice on Channel 49.5N : 51N, -3E : 2E
	        #ycoords = [53.5, 53.5]; xcoords = [-3.1, -3.1] # Nr Liverpool
		#ycoords = [43,63]; xcoords = [-13,13] # Whole domain
	coordsType = 'ind' - input xcoords and ycoords are indices

	RETURN:
	lat - array of latitudes [ny,nx]
	lon - array of longitudes [ny,nx]
	data - array of COMPLEX harmonic constituents [nh,ny,nx]
	doodson_list - a list of Doodson number that are available [nh]
	constit_list - the corresponding list of harmonic constituent labels [nh]

	In this example the data are stored in variables like M2x_SSH and M2y_SSH for the real and imaginary parts of the M2 SSH harmonic
	"""
	#dirname = '/projectsa/pycnmix/jelt/AMM60/'
	[ constit_list, period_list, doodson_list ] = harmonictable(dirname+'../harmonics_list.txt', doodson=True)
	#[ constit_list, period_list ] = harmonictable(dirname+'../harmonics_list.txt')
	print(doodson_list)
	print(constit_list)
	nh = len(constit_list)
	try:
		fD1 = Dataset(dirname + filebase + '_D1_Tides.nc')
		fD2 = Dataset(dirname + filebase + '_D2_Tides.nc')
		fD4 = Dataset(dirname + filebase + '_D4_Tides.nc')
	except:# AMM7 does not have separate files for harmonic groups
		fD1 = Dataset(dirname + filebase + '_Tides.nc')
		fD2 = fD1
		fD4 = fD1

	# Find indices for specified coordinates. Need full domain to find indices
	lats_full = fD2.variables['nav_lat_grid_T'][:]
	lons_full = fD2.variables['nav_lon_grid_T'][:]

	if coordsType == 'deg':
		[J1,J2,I1,I2] = get_coord_indices(ycoords,xcoords,lats_full,lons_full)
	else:
		J1,J2 = ycoords
		I1,I2 = xcoords

	# Load in subdomain
	lat_arr = fD2.variables['nav_lat_grid_T'][J1:J2,I1:I2]
	lon_arr = fD2.variables['nav_lon_grid_T'][J1:J2,I1:I2]

	[ny,nx] = np.shape(lat_arr)

	# Test for the dimensionality of the requested data (could pass as a variable). Initialise target array
	var_full_shape =  np.shape(fD2.variables['M2x_' + varstr][:])
	if len(var_full_shape) == 3:
		nz = var_full_shape[0]
		data_arr =  np.zeros((nh,nz,ny,nx) ,dtype=complex)
	elif len(var_full_shape) == 2:
		data_arr =  np.zeros((nh,ny,nx) ,dtype=complex)
	else:
		print('Panic!!')

	#for iconst in range(6,7): # M2 only
	for iconst in range(nh):

		# Get the harmonic file handle - verbose method for transparancy.
		if constit_list[iconst][-1] == '1':
			fileh = fD1
		elif constit_list[iconst][-1] == '2':
			fileh = fD2
		elif constit_list[iconst][-1] == '4':
			fileh = fD4
		else:
			print('{}: Not ready for that harmonic species band'.format(constit_list[iconst]))


		print('available: ', constit_list[iconst])
		constit = constit_list[iconst]
		#tmp_arr  = fileh.variables[constit+'x_' + varstr][...,ny,nx] + 1.j*fileh.variables[constit+'y_' + varstr][...,ny,nx]

		if len(var_full_shape) == 3: # 3D data
			data_arr[iconst,:,:,:]  =  fileh.variables[constit+'x_' + varstr][:,J1:J2,I1:I2] + 1.j*fileh.variables[constit+'y_' + varstr][:,J1:J2,I1:I2]
		if len(var_full_shape) == 2: # 2D data
			data_arr[iconst,:,:]  =  fileh.variables[constit+'x_' + varstr][J1:J2,I1:I2] + 1.j*fileh.variables[constit+'y_' + varstr][J1:J2,I1:I2]
	print('size of data: {}'.format( len(np.shape(data_arr)) ))
	try:
		fD1.close()
		fD2.close()
		fD4.close()
	except:
		print('Close harmonic files. Assuming that there was only one file')
		pass
	return lat_arr, lon_arr, data_arr, doodson_list, constit_list

##########################################################################

def UtcNow():
    now = datetime.datetime.utcnow()
    return now

def date2mjd(dates):
    """
    Convert datetime into Modified Julian Date (float)

    This is a float data reference to 2000-01-01T00:00:00 which has a mjd of 51544
    INPUT: datetime array
    OUTPUT: mjd array
    """

    mjd = [  (d - datetime.datetime(2000,1,1)).days     \
           + (d - datetime.datetime(2000,1,1)).seconds/86400. + 51544 for d in dates]

    return mjd
##########################################################################
def test_port(mjd, site_file_name="glad.txt"):
	"""
	Demonstration to load, reconstruct and plot port data.
	"""
	rad = np.pi / 180.
	# Obtain data
	lat, lon, z0, data = get_port(site_file_name)
	ha = data['amp'][:]
	ga = data['pha'][:]
	doo = data['doo'][:]
	lab = data['lab'][:]
	nh = len(ha)

	# Initialise output variable
	npred = np.shape(mjd)[0]            # vector: [npred]
	pred = np.zeros(npred)

	# Get the nodal corrections
	mjdns = np.floor(mjd)       # vector: [npred]
	hrs = 24 * (mjd - mjdns)   # hrs are the time variable required for the reconstruction. vector: [npred]
	[f,v] = phamp0fast(mjdns)  # [npred,120]
	print('shape of f,v {}'.format(np.shape(f)))
	print('shape of hrs time series {}'.format(np.shape(hrs)))

	# load the phase speeds and harmonic names for 120 harmonics
	names, sig0 = set_names_phases()

	#
	# Sum constant + constituents to give prediction
	#

	for j in range(nh):
		"""
		Need to match up input harmonics with tabulated harmonic. Can pair by names,
		or by doodson number (i.e. index). The former is preferable but there are
		issues with spelling conventions.
		"""
		#k = names.index(lab[j].upper()) # This is the index of the harmonic from vector sig0, of speeds
		k = doo[j]-1 # offset since indexing starts from zero
		try:
			print(names[k], lab[j])
		except:
			print("Issue with displaying name and label")
		## Port
		pred = pred + ha[j] * f[:,k] * np.cos(rad * ( sig0[k]*hrs + v[:,k] - ga[j] ))
		## Map
		#    pred = pred + ha(j+1) * f[:,k] * np.cos(rad * ( sig0[k]*hrs + v[:,k] - ga[j] ))

	print('shape of prediction {}'.format(np.shape(pred)))
	ssh = pred + z0

	print('prediction values: {}'.format(pred))

	return ssh


##########################################################################
def test_dataarray(xcoords, ycoords, mjd):
	"""
	INPUTS: xcoords : [49.5, 51]
		ycoords : [-3, 2]
		can be single value arrays for a point location
	"""
	# Obtain data
	lat_sub, lon_sub, data_sub, doodson_list, constit_list = get_harmonic_arr('SSH',xcoords,ycoords)
	[nh, ny, nx] = np.shape(data_sub)

	Ha = np.abs(data_sub)
	Ga = np.angle(data_sub, deg=True)

	# Reconstruct sea level
	# Initialise output variable
	npred = np.shape(mjd)[0]            # vector: [npred]
	pred = np.zeros((npred, ny, nx)) # npred, ny, nx

	# Get the nodal corrections
	mjdns = np.floor(mjd)       # vector: [npred]
	hrs = 24 * (mjd - mjdns)   # hrs are the time variable required for the reconstruction. vector: [npred]
	[f,v] = phamp0fast(mjdns)  # [npred,120]
	print('shape of f,v {}'.format(np.shape(f)))
	print('shape of hrs time series {}'.format(np.shape(hrs)))

	# load the phase speeds and harmonic names for 120 harmonics
	names, sig0 = set_names_phases()

	#
	# Sum constant + constituents to give prediction
	#

	for j in range(nh):
	    """
	    Need to match up input harmonics with tabulated harmonic. Can pair by names,
	    or by doodson number (i.e. index). The former is preferable but there are
	    issues with spelling conventions.
	    """
	    #k = names.index(lab[j].upper()) # This is the index of the harmonic from vector sig0, of speeds
	    k = doodson_list[j]-1 # offset since indexing starts from zero
	    print(names[k]) #, constit_list[j]
	    ## Port
	    #pred = pred + ha[j] * f[:,k] * np.cos(rad * ( sig0[k]*hrs + v[:,k] - ga[j] ))
	    ## Map
	    Ha_arr = np.tile(Ha[j,:,:],(npred,1,1))
	    f_arr = np.tile(f[:,k],(ny,nx,1)).transpose(2,0,1)
	    sig_hrs_arr = np.tile(sig0[k]*hrs,(ny,nx,1)).transpose(2,0,1)
	    v_arr = np.tile(v[:,k],(ny,nx,1)).transpose(2,0,1)
	    Ga_arr = np.tile(Ga[j,:,:],(npred,1,1))
	    pred = pred + Ha_arr * f_arr * np.cos( rad * (sig_hrs_arr + v_arr - Ga_arr) )

	ssh = np.squeeze(pred)

	print(np.shape(ssh))



	return lat_sub, lon_sub, ssh

##########################################################################
def plot_map(dates, lat_sub, lon_sub, ssh):

	ssh = np.ma.masked_where( ssh > 1E6, ssh)

	# fix for spurious lat value at zero
	lat_sub[ lat_sub == 0 ] = np.nan

	plt.figure()
	plt.pcolormesh(lon_sub,lat_sub,ssh)
	plt.clim([0, np.floor(np.nanmax(ssh)) ])
	plt.xlim(xcoords)
	plt.ylim(ycoords)
	plt.colorbar(extend='max')
	plt.title('max SSH over 24 hours [m]')
	plt.xlabel('longitude'), plt.ylabel('latitude')

	print(datetime.datetime.now() - startTime)

	plt.show()
	return

##########################################################################
def plot_port(dates, ssh):
    ssh = np.ma.masked_where( ssh > 1E6, ssh)

    # Plot sea level time series
    plt.figure()
    plt.plot(dates,[ssh[i] for i in range(len(dates))],'+-')
    plt.ylabel('Height (m)')
    plt.xlabel('Time since '+dates[0].strftime("%Y-%m-%d"))
    plt.show()

    return

##################################################################################
##################################################################################
## Now do the main routine stuff
##################################################################################
##################################################################################
if __name__ == '__main__':

    # Settings
    rad    = np.pi/180
    deg    = 1.0 / rad

    startTime = datetime.datetime.now() # For a timing test


    # Set the dates
    startdate = UtcNow()
    npred = 24
    dates = [startdate + datetime.timedelta(hours=hh) for hh in range(0, npred)]
    mjd = date2mjd( dates ) # convert to modified julian dates


    ## Compute reconstuction on port data.
    #####################################
    ssh = test_port(mjd, site_file_name="glad.txt")

    print('plot time series reconstruction of port data')
    plot_port(dates, ssh)


    try:
        ## Compute reconstruction on model data (1D or 2D)
        ##################################################
        ycoords = [49.5, 51]; xcoords = [-3, 2] # Slice on Channel 49.5N : 51N, -3E : 2E
        #ycoords = [53.5, 53.5]; xcoords = [-3.1, -3.1] # Nr Liverpool
        #ycoords = [43,63]; xcoords = [-13,13] # Whole domain

        [lat_sub, lon_sub, ssh] = test_dataarray(xcoords, ycoords, mjd)
        if len(np.shape(ssh)) > 1: # mapa
            print('plot max over 24hours')
            ssh = np.squeeze(np.nanmax(ssh, axis=0))
            plot_map(dates, lat_sub, lon_sub, ssh)
        else:
            print('plot time series reconstruction of point data')
            plot_port(dates, ssh)
    except:
        pass
