#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Use anyTide compiled C-code to reconstruct tidal sea level from harmonics, for point locations.
Iterating over areas to produce maps is slow and the C needs to be modified.

Starting point: anyTide/tileserver/generate_velocity2/velocity4.py
The revelation was how to compile C and call it as a library in python.

Author: jelt
Date: 12 July 2018
"""

# For point reconstruction
import datetime
import polpredict
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Additional for map reconstruction
from netCDF4 import Dataset
import sys
sys.path.insert(0, "/login/jelt/python/ipynb/NEMO/")
import internaltideharmonics_NEMO as ITh


def UtcNow():
    now = datetime.datetime.utcnow()
    return now


def get_port():

	"""
	Obtain a reduced number of Gladstone harmonics from a file

	$ head glad.txt

	Port Name: ENGLAND, WEST COAST $ LIVERPOOL (GLADSTONE DOCK)
	53 27.0 N 03 01.1 W
	z0=  5.249   OD= -4.930
	3.03800 320.72000   31 M2
	0.97800   4.70000   36 S2
	...
	"""

	data = pd.read_csv('glad.txt', header=2, delimiter=r"\s+")
	data.columns = ['amp','pha', 'doo', 'lab']
	lat = +(53+27.0/60)
	lon = -(03+01.1/60)
	z0  = +5.249

	return lat,lon, z0, data


def get_map():
	dirname = '/projectsa/pycnmix/jelt/AMM60/'
	[ constit_list, period_list, doodson_list ] = ITh.harmonictable(dirname+'../harmonics_list2.txt', doodson=True)
	#[ constit_list, period_list ] = ITh.harmonictable(dirname+'../harmonics_list.txt')
	print doodson_list
	print constit_list
	nh = len(constit_list)
	f = Dataset(dirname + 'AMM60_1d_20120801_20120831_D2_Tides.nc')

	lat_arr = f.variables['nav_lat_grid_T'][:]
	lon_arr = f.variables['nav_lon_grid_T'][:]
	[ny,nx] = np.shape(lat_arr)

	data_arr = np.zeros((nh,ny,nx),dtype=complex)

	for iconst in range(6,7):
	#	print 'available: ', constit_list[iconst]
		constit = constit_list[iconst]
		data_arr[iconst,:,:]  = f.variables[constit+'x_SSH'][:] + 1.j*f.variables[constit+'y_SSH'][:]
	return lat_arr, lon_arr, data_arr, doodson_list


def reconstructPointSealevel(date, lat, lon, df):
    """
    Reconstruct a 24h (hourly) timeseries for whole input DAY, not time.

    Need to convert data to row of amplitude, phases, doodson number strings, with 40 sets
    AMP1, PHA1, DOOD1, AMP2, PHA2, DOOD2, ...,AMP40, PHA40, DOOD40
    """
    date_str = date.strftime("%Y-%m-%d 00:00:00")

    str1 = ""

    for x in range(len(df['amp'])): #range(0, 40):
        vv3H  = df["amp"][x]
        vv3G  = df["pha"][x]
        vv3K  = float(df["doo"][x])
        str1 += str(vv3H) + "," + str(vv3G) + "," + str(vv3K)
    if x < 39: #len(df['amp'])-1: # 39
        str1 += ",";

    # Pad out data to 40 constituents, because the C-code expects
    if len(df['amp']) < 40:
    	for x in range(len(df['amp']),40):
    		str1 += "0,0,0" # Add zero amplitudes onto an insgignificant Doodson number
    		if x< 39:
    			str1 += ",";
        ssh = polpredict.predictor(date_str, str(lat), str(lon), str1)
    return ssh

def reconstructMapSealevel(date, lat_arr, lon_arr, data_arr, doodson_list):

        date_str = date.strftime("%Y-%m-%d 00:00:00")

	[nh,ny,nx] = np.shape(data_arr)
	Npts = len(lat_arr.flatten())
	print Npts
	print nx,ny,nx*ny

	res = np.zeros(Npts)

	for count in range(Npts-1):
		#print count
		lats = lat_arr.flatten()
		lons = lon_arr.flatten()
		lat = lats[count]
		lon = lons[count]
		str1 = ""

		for iconst in range(nh):
			data = data_arr[iconst,:,:].flatten()
			vv3H = np.abs(data[count])
			vv3G = np.angle(data[count], deg=True)
			vv3K = doodson_list[iconst]
			str1 += str(vv3H) + "," + str(vv3G) + "," + str(vv3K)
			if iconst < 39:
				str1 += ",";

		# Pad out data to 40 constituents, because the C-code expects
		if nh < 40:
			for x in range(nh,40):
				str1 += "0,0,0" # Add zero amplitudes onto an insgignificant Doodson number
				if x< 39:
					str1 += ",";
		#print date_str
		#print lat
		#print lon
		#print str1
		res[count] = np.max( polpredict.predictor(date_str, str(lat), \
				str(lon), str1) )

	return res.reshape((ny, nx))



def test_port():

	# Set the date
	date_str = UtcNow()

	# Obtain data
	lat, lon, z0, data = get_port()

	# Reconstruct sea level
	ssh = reconstructPointSealevel(date_str, lat, lon, data) # build 24h (hourly) timeseries for whole input DAY, not time.

	# Plot sea level time series
	plt.figure()
	plt.plot([z0 + ssh[i] for i in range(len(ssh))],'+-')
	plt.ylabel('Height (m)')
	plt.xlabel('Hour since midnight')
	plt.show()

def test_map():

	# Set the date
	date_str = UtcNow()

	# Obtain the data
	lats, lons, data, doodson_list = get_map()
	nh,ny,nx = np.shape(data)

	# Slice on Channel 49.5N : 51N, -3E : 2E
	ycoords = [49.5, 51]
	xcoords = [-3, 2]
	[J_ll,I_ll] = ITh.findJI(min(ycoords), min(xcoords), lats, lons)  # Simple routine to find the nearest J,I coordinates for given lat lon
	[J_ur,I_ur] = ITh.findJI(max(ycoords), max(xcoords), lats, lons)  # Simple routine to find the nearest J,I coordinates for given lat lon
	J1 = min(J_ll,J_ur)
	J2 = max(J_ll,J_ur)+1
	I1 = min(I_ll,I_ur)
	I2 = max(I_ll,I_ur)+1

	lat_sub = lats[J1:J2,I1:I2]
	lon_sub = lons[J1:J2,I1:I2]
	data_sub = data[:,J1:J2,I1:I2]

	ssh = reconstructMapSealevel(date_str, lat_sub, lon_sub, data_sub, doodson_list)
	# unsophisticated managment of unphysical land points
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

	print datetime.datetime.now() - startTime
	plt.show()

##################################################################################
##################################################################################
## Now do the main routine stuff
##################################################################################
##################################################################################
if __name__ == '__main__':

	startTime = datetime.datetime.now()

	## Demonstration reconstruction of sea level for a port.
	#test_port()

	## Demonstration reconstruction of sea level for a port.
	########################################################################
	test_map()
