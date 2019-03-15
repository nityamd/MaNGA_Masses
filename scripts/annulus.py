import numpy as np
import scipy
from scipy.stats import binned_statistic
from scipy.stats import binned_statistic_2d
import scipy.interpolate as interp
import scipy.stats as st
import scipy.signal as sig
from scipy.optimize import curve_fit


import h5py
import gc
import sys
import subprocess
import pickle
import os


import astropy as ap
from astropy.io import fits
from astropy import table as t
from astropy.table import Table
from astropy import wcs
from astropy.cosmology import WMAP7


def get_dn4000(wave,spec):
    interp_spec = interp.interp1d(wave,spec)
    blue_wav = np.linspace(3850,3950,100)
    red_wav = np.linspace(4000,4100,100)
    d4000 = np.sum(interp_spec(red_wav)) / np.sum(interp_spec(blue_wav))
    return d4000

def get_HdA(wave,spec):
    bandw_HdA = np.logical_and(wave > 4083.500, wave < 4122.250)   # analogous to MPA-JHU
    bandw_HdA_blueside = np.logical_and(wave > 4041.600, wave < 4079.750)
    bandw_HdA_redside = np.logical_and(wave > 4128.500, wave < 4161.000)
    spec_av_blueside = np.sum(spec[bandw_HdA_blueside])/len(spec[bandw_HdA_blueside])
    spec_av_redside = np.sum(spec[bandw_HdA_redside])/len(spec[bandw_HdA_redside])
    a_spec = (spec_av_redside - spec_av_blueside)/((4161.000+4128.500)/2 - (4079.750+4041.600)/2)
    b_spec = spec_av_blueside - a_spec * (4079.750+4041.600)/2
    spec_cont_HdA = wave[bandw_HdA] * a_spec + b_spec
    HdA = tsum(wave[bandw_HdA],np.divide((spec_cont_HdA - spec[bandw_HdA]), spec_cont_HdA))
    return HdA

def tsum(xin, yin):
    tsum = np.sum(np.abs((xin[1:]-xin[:-1]))*(yin[1:]+yin[:-1])/2. )
    return tsum


def read_file(filename,)
dat = Table.read('manga-8455-3701-LOGCUBE.fits', format='fits')
drpall = t.Table.read('drpall-v2_0_1.fits')
drpall.add_index('plateifu')
drpall.loc['8455-3701']

drp_logcube = fits.open('manga-8455-3701-LOGCUBE.fits')

