"""
extract_harmonics_from_database.py

Extract harmonics for single location from private coastal harmonic data product.
Then write harmonics to text file (for use with e.g. NOCtidepred.py)

Currently, harmonics need to be copied to anyTide_Cwrapper for use in NOCtidepred.py

jp
6 Oct 2023
"""
import numpy as np
import xarray as xr
import pandas as pd
from pathlib import Path

from typing import Union
from pathlib import Path
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from NOCtidepred import set_names_phases



class dataset:
    """
    Class for reading anyTide dataset from CVS file into an xarray object

            columns : ['lat', 'lon', 'depth',
                     'amp1', 'pha1', 'dood1',
                     'amp2', 'pha2', 'dood2',
                     ...
                     'amp40', 'pha40', 'dood40']
            rows : stations

    """

    def __init__(self, file_path: str = "/Users/jelt/DATA/anyTide/all_constit.txt"):
        """Init anyTide data object. Load data file

        Args:
            file_path (str) : path of data file
        """
        print(f"Creating a new {print(self)}")
        self.read_data(file_path)
        print(f"Loaded {file_path}.\nStations: {len(self.lat)}")


    def read_data(self, file_path: str) -> None:
        """Read the data file

        Expected format and variable names are

            lat lon depth amp(1) pha(1) dood(1)  ... amp(15) pha(15) dood(15) \n%% ...\n%%\n%% Note that the lat, lon, depth and dood(1:15) for the ith row is the same in U, V and Z tables

        xarray dataset to have dimension as time and coordinates as time, latitude and longitude

        Args:
            file_path (str) : path of data file
        """
        txt = np.loadtxt(file_path)

        # Save transposed np array
        self.lat = txt[:, 0].astype(float)
        self.lon = txt[:, 1].astype(float)
        self.dep = txt[:, 2].astype(float)
        self.amp = txt[:, 3::3].astype(float)
        self.pha = txt[:, 4::3].astype(float)
        self.dood = txt[:, 5::3].astype(int)

    def doodson_to_name_str(self, index) -> str:

        # load the phase speeds and harmonic names for 120 harmonics

        # The following phase speeds and names are indexed by what I'm called the
        # Doodson number, which start at one: Doodson(SA)=1.
        # Note that for python indexing these are called with indexing starting from
        # zero so that names[ doodson(SA)-1 ] = SA

        names, sig0 = set_names_phases()
        if index != 0: return names[index-1]
        else: return 0  # Don't return a name for zero index


    def plot_database(self, ds = None, pt = None, index = None):
        """ plot map of database locations. Highlight pt if specified """
        if ds == None:
            ds = self

        fig, ax = plt.subplots(1, 1, figsize=[5.5, 5.5])

        # Create inset of width 1.3 inches and height 0.9 inches
        # at the default upper right location
        axins = inset_axes(ax, width=1.3, height=0.9)

        ax.plot( ds.lon, ds.lat, 'k.')
        if index is not None:
            ax.plot(ds.lon[index], ds.lat[index], 'ro')

        # if point is given, plot zoom and show whole region as inset
        if pt is not None:
            ax.plot(pt[0], pt[1], 'r+')
            ax.set_xlim([pt[0]-0.1, pt[0]+0.1])
            ax.set_ylim([pt[1]-0.1, pt[1]+0.1])
            # Add inset for whole region
            axins.plot( ds.lon, ds.lat, 'k.')
            axins.plot(pt[0], pt[1], 'r+')
            # Turn ticklabels of inset off
            axins.tick_params(labelleft=False, labelbottom=False)

        plt.show()


    def find_nearest(self, lon: float = -3.0183, lat: float = 53.45):
        """
        find nearest location, in coordinate space
        args
            lon: float (deg) default Gladstone Dock
            lat: float (deg) default Gladstone Dock
        returns
         distance: float -- The distances to the nearest neighbors
         index: int -- The locations of the neighbors

        # pt = [-3.0183, 53.45]  # <-- Gladstone Dock
        # pt = [-2.88, 53.18]  # <-- Chester
         """
        A = np.array([self.lon, self.lat]).T
        pt = [lon, lat]

        distance, index = spatial.KDTree(A).query(pt, k=1)
        return distance, index

    def save_station_harmonics_txt(self, index, label: str = None):
        """
        Format of csv file:
            Port Name: ENGLAND, WEST COAST $ LIVERPOOL (GLADSTONE DOCK)
            53 27.0 N 03 01.1 W
            z0=  5.249   OD= -4.930
            3.03800 320.72000   31 M2
            0.97800   4.70000   36 S2
            ...
        Can be read by anyTide/NOCtidepred.py


        args:
            index: int - row for station in question
            label: str - name to be inserted into txt file, also file name
        """

        if label == None or type(label) is not str:
            label = input("Enter station name: ")

        # Creating pandas dataframe from numpy array
        df = pd.DataFrame({'A': self.amp[index, :],
                            'G': self.pha[index, :],
                            'K': self.dood[index, :],
                            #'lat': tt.lat[index],
                            #'lon': tt.lon[index],
                            #'dep': tt.dep[index]
                           })

        # construct a list of constituent names from Doodson numbers. To be added as a column.
        names = []
        for i in range(np.shape(self.dood[index, :])[0]):
            names.append(self.doodson_to_name_str( self.dood[index, i] ))

        df['names'] = names

        # Sort pd.dataframe by amplitude
        df = df.sort_values('A', ascending=False)
        # Remove rows with all zeros
        df = df.loc[(df != 0).any(axis=1)]

        # Write data to file. Then add header info
        #filepath = Path('folder/subfolder/out.txt')
        filepath = Path(f'{label}.txt')
        filepath.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(filepath, sep=' ', index=False, header=False)

        with open(filepath, 'r') as original: data = original.read()
        with open(filepath, 'w') as modified:
            modified.write(f'{label} src:extract_harmonics_from_database.py\n{self.lat[index]}N {self.lon[index]}E\nz0= {self.dep[index]}\n' + data)

        return df


##################################################################################
##################################################################################
## Now do the main routine stuff
##################################################################################
##################################################################################
from scipy import spatial
from collections import OrderedDict

if __name__ == '__main__':

    # Settings
    pt = [-3.0183, 53.45] # gladstone
    pt = [-(2+(56 +6.93/60)/60), 56+(12+44.20/60)/60 ] # lower largo

    # initialise dataset with data loaded
    tt = dataset()

    # find nearest location
    distance, index = tt.find_nearest(lon=pt[0], lat=pt[1])
    print(f'Input location: {pt[1]}N, {pt[0]}E'.format('%.1f'))
    print(f'Nearest distance: {distance}'.format('0.1f'))  # <-- The distances to the nearest neighbors
    print(f'has index: {index}')  # <-- The locations of the neighbors

    # Plot data base + station
    tt.plot_database(pt=pt, index=index)

    # Write txt file
    stn = tt.save_station_harmonics_txt(index)
    print(stn)






