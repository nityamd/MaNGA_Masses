import numpy as np
import scipy
import scipy.interpolate as interp
import scipy.stats as st
import scipy.signal as sig
import pandas as pd

import pickle
import os
import datetime

import astropy as ap
from astropy.io import fits
from astropy import table as t
from astropy.table import Table
from astropy import wcs

from aperture_spec import spec_measurements, aperture_measurements

print(datetime.datetime.now())

os.chdir("/Volumes/Nitya/Data/mpl8")
drpall = t.Table.read('drpall-v2_5_3.fits')
drpall.add_index('plateifu')
index1 = np.where(drpall['srvymode']=='MaNGA dither')[0]
drpall1 = drpall[index1]
index2 = np.where(drpall1['z']>0)[0]
drpall = drpall1[index2]

z = drpall['z']
plate_ifu = drpall['plateifu']
fits_files = ['manga-'+ str(x) + '-LOGCUBE.fits' for x in plate_ifu]

print('... before the calculation (5 arcsec) ...', datetime.datetime.now())

m = 4000
n = 6491
sample = []
c = 1
for i in range(m,n):
    start = datetime.datetime.now()
    galaxy = aperture_measurements(fits_files[i],z[i])
    measures = galaxy.get_spec_measure(5,'all')
    sample.append(measures)
    end = datetime.datetime.now()
    elapsed = end - start
    print('time elapsed in seconds', elapsed.seconds)
    print('Executed '+str(c) + 'th galaxy.....' )
    print(fits_files[i])
    print('.....................time now', end)
    c += 1

print(datetime.datetime.now())

#numpy.array([a, b, c],dtype=[('a','f8'),('b','f8'),('c','f8')])

sample_ifu = np.array(list(zip(plate_ifu[m:n],z[m:n],np.array(sample)[:,0], np.array(sample)[:,1])),
                        dtype = [('plate_ifu', 'U25'),('z', 'f8'), ('hdelta', 'f8'), ('dn4000', 'f8')])
afile = open(r'/Volumes/500GB/nmd299/Data/mpl8/5arcsec_c.pkl', 'wb')
pickle.dump(sample_ifu,afile)

print('... before the calculation (5 arcsec) ...', datetime.datetime.now())


sample = []
c = 1
for i in range(m,n):
    start = datetime.datetime.now()
    galaxy = aperture_measurements(fits_files[i],z[i])
    measures = galaxy.get_spec_measure(7,'all')
    sample.append(measures)
    end = datetime.datetime.now()
    elapsed = end - start
    if i%10==0:
        print('time elapsed in seconds', elapsed.seconds)
        print('Executed '+str(c) + 'th galaxy.....' )
        print(fits_files[i])
        print('.....................time now', end)
    c += 1

print(datetime.datetime.now())
sample_ifu = np.array(list(zip(plate_ifu[m:n],z[m:n],np.array(sample)[:,0], np.array(sample)[:,1])),
                        dtype = [('plate_ifu', 'U25'),('z', 'f8'), ('hdelta', 'f8'), ('dn4000', 'f8')])
afile = open(r'/Volumes/500GB/nmd299/Data/mpl8/7arcsec_c.pkl', 'wb')
pickle.dump(sample_ifu,afile)
