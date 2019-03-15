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
z_array = np.linspace(0,0.15,16)
new_array = []
for redshift in z:
    new = [x for x in z_array if x < redshift]
    new.append(redshift)
    new_array.append(new)
z_dict = dict(zip(z,new_array))
print(len(z))
m = 0
n = len(z)
sample = []
c = 1

for i in range(m,n):
    start = datetime.datetime.now()
    redshift = z[i]
    galaxy = aperture_measurements(fits_files[i],redshift)
    measures = [galaxy.get_spec_measure(3,z_new,'all')
            for z_new in z_dict[redshift]]
    print(measures)
    sample.append(measures)
    end = datetime.datetime.now()
    elapsed = end - start
    print('time elapsed in seconds', elapsed.seconds)
    print('Executed '+str(c) + 'th galaxy.....' )
    print(fits_files[i])
    print('.....................time now', end)
    c += 1

print(datetime.datetime.now())
array = list(zip(drpall['plateifu'],drpall['z'],new_array,sample))
afile = open(r'redshift_data.pkl', 'wb')
pickle.dump(array,afile)
