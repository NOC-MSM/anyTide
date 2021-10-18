#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Oban_demo.py

Harmonic reconstruction.
Reconstructs for a vector of times.

Usage: python Oban_demo.py

Author: jelt@noc.ac.uk
Date: 1 Oct 2019
"""
# For point reconstruction
import numpy as np
import math # atan2
import pandas as pd
import datetime
from os import path # fix the path
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

from NOCtidepred import UtcNow
from NOCtidepred import date2mjd
from NOCtidepred import phamp0fast
from NOCtidepred import set_names_phases

##########################################################################
def get_port():

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

        fname =  "Oban.txt"

        data = pd.read_csv(fname, header=3, delimiter=r"\s+")
        data.columns = ['amp','pha', 'doo']
        lat = +(56+25/60)
        lon = -( 5+29/60)
        z0  = +2.402

        return lat,lon, z0, data

##########################################################################
def test_port(mjd):
    """
    Demonstration to load, reconstruct and plot port data.
    """
    rad = np.pi / 180.
    # Obtain data
    lat, lon, z0, data = get_port()
    ha = data['amp'][:]
    ga = data['pha'][:]
    doo = data['doo'][:]
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
    ## Port
        pred = pred + ha[j] * f[:,k] * np.cos(rad * ( sig0[k]*hrs + v[:,k] - ga[j] ))
    ## Map
    #    pred = pred + ha(j+1) * f[:,k] * np.cos(rad * ( sig0[k]*hrs + v[:,k] - ga[j] ))

    print('shape of prediction {}'.format(np.shape(pred)))
    ssh = pred + z0

    print('prediction values: {}'.format(pred))

    return ssh

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
    # Create a vector of predictions times. E.g. 24 hourly instants. 
    startdate = UtcNow()
    npred = 24
    dates = [startdate + datetime.timedelta(hours=hh) for hh in range(0, npred)]
    mjd = date2mjd( dates ) # convert to modified julian dates


    ## Compute reconstuction on port data.
    #####################################
    ssh = test_port(mjd) # reconstuct ssh for the time vector
    print('plot time series reconstruction of port data')

    ssh = np.ma.masked_where( ssh > 1E6, ssh) # get rid of nasties

    # Plot sea level time series
    fig, ax = plt.subplots()
    ax.plot(np.array(dates),[ssh[i] for i in range(len(dates))],'+-')
    ax.set_ylabel('Height (m)')
    ax.set_xlabel('Hours since '+dates[0].strftime("%Y-%m-%d"))
    ax.set_title('Oban tide prediction')
    
    # Pain plotting time on the x-axis
    myFmt = mdates.DateFormatter('%H')
    ax.xaxis.set_major_formatter(myFmt)

    plt.show()
