import numpy as np
import scipy
import scipy.interpolate as interp
import scipy.stats as st
import scipy.signal as sig
import pandas as pd

import pickle
import os
import datetime
import pandas as pd

import astropy as ap
from astropy.io import fits
from astropy import table as t
from astropy.table import Table
from astropy import wcs

from aperture_spec import spec_measurements, aperture_measurements

print(datetime.datetime.now())

drpall = t.Table.read('/Volumes/Nitya/Data/mpl8/drpall-v2_5_3.fits')
drpall.add_index('plateifu')
index1 = np.where(drpall['srvymode']=='MaNGA dither')[0]
drpall1 = drpall[index1]
index2 = np.where(drpall1['nsa_z']>0)[0]
drpall = drpall1[index2]

afile = open(r'filename_dict','rb')
w = pickle.load(afile, encoding = 'latin1')
afile.close()

z = drpall['z']
plate_ifu = drpall['plateifu']

fits_files = [w[plateifu] for plateifu in plate_ifu]

"""
our "z_array" for each z is going to be
a bunch of z's below it...
"""
z_array = [0.5]

print(len(z))
m = 0
n = len(z)
sample = []
c = 1
for z_new in z_array:
    print('shifting galaxies to  '+ str(z_new))
    start = datetime.datetime.now()
    sample = []
    for i in range(m,n):
        redshift = z[i]
        galaxy = aperture_measurements(fits_files[i],redshift)
        if redshift < 0.02:
            a = galaxy.get_spec_measure(3,0.5,'all')
            measures = [a[0],a[1],z_new]
        else:
            measures = [np.nan,np.nan,z_new]
        sample.append(measures)
    end = datetime.datetime.now()
    elapsed = end - start
    print('time elapsed in seconds', elapsed.seconds)
    print('.....................time now', end)
    c += 1
    array = list(zip(drpall['plateifu'],drpall['z'],sample))
    filename = 'codetesting'+str(z_new)
    afile = open(filename, 'wb')
    pickle.dump(array,afile)
