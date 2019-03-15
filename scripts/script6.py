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

dim_x = []
dim_y = []
dim_l = []
for file in fits_files:
        drp_logcube = fits.open(file)
        NL, NY, NX = (drp_logcube['FLUX']).data.shape
        dim_x.append(NX)
        dim_y.append(NY)
        dim_l.append(NL)


sample = np.array(list(zip(plate_ifu,z,dim_x,dim_y,dim_l)),
                        dtype = [('plate_ifu', 'U25'),('z', 'f8'),
                                 ('xdim', 'f8'), ('ydim', 'f8'),
                                 ('NL', 'f8')])
afile = open(r'dataset.pkl', 'wb')
pickle.dump(sample,afile)
