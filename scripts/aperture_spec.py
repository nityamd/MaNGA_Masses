import numpy as np
import scipy
import scipy.interpolate as interp
import scipy.stats as st
import scipy.signal as sig

import pickle
import os

import astropy as ap
from astropy.io import fits
from astropy import table as t
from astropy.table import Table
from astropy import wcs
from astropy.cosmology import WMAP9 as cosmo


class spec_measurements():
    """
    Spectral measurements class initialized
    on a datacube instance.
    Currently measures the HdA and Dn4000
    (as estimated in the MPA-JHU)
    """
    def __init__(self,drp_logcube,z):
        """
        wave: set of wavelengths the spectrograph spans
        """

        self.wave = (drp_logcube['WAVE']).data/(1+z)
        self.bandw_HdA = np.logical_and(self.wave > 4083.500,
                                        self.wave < 4122.250)
        self.bandw_HdA_blueside = np.logical_and(self.wave > 4041.600,
                                                 self.wave < 4079.750)
        self.bandw_HdA_redside = np.logical_and(self.wave > 4128.500,
                                                self.wave < 4161.000)
        self.blues = len(self.wave[self.bandw_HdA_blueside])
        self.reds = len(self.wave[self.bandw_HdA_redside])

        self.blue_wav = np.linspace(3850,3950,100)
        self.red_wav = np.linspace(4000,4100,100)

    def tsum(self, xin, yin):
        """
        Trapezoidal Sum to estimate
        area under curve
        """
        tsum = np.sum(np.abs((xin[1:]-xin[:-1]))*(yin[1:]+yin[:-1])/2. )
        return tsum

    def dn4000_red(self,spec):
        interp_spec = interp.interp1d(self.wave,spec)
        d4000_r = np.sum(interp_spec(self.red_wav))
        return d4000_r

    def dn4000_blue(self,spec):
        interp_spec = interp.interp1d(self.wave,spec)
        d4000_b =  np.sum(interp_spec(self.blue_wav))
        return d4000_b

    def dn_4000(self,specs):
        dn_4000 = np.sum([self.dn4000_red(spec) for spec in specs])/np.sum(
            [self.dn4000_blue(spec) for spec in specs])
        return dn_4000

    def HdA(self,specs):

        spec_av_blueside = np.sum([np.sum(spec[self.bandw_HdA_blueside]) for
                                   spec in specs])/(self.blues*len(specs))
        spec_av_redside = np.sum([np.sum(spec[self.bandw_HdA_redside]) for
                                  spec in specs])/(self.reds*len(specs))

        a_spec = (spec_av_redside - spec_av_blueside)/(
                (4161.000+4128.500)/2.0 - (4079.750+4041.600)/2.0)
        b_spec = spec_av_blueside - a_spec * (4079.750+4041.600)/2

        spec_cont_HdA = self.wave[self.bandw_HdA] * a_spec + b_spec
        mean_dip = np.mean([spec[self.bandw_HdA] for spec in specs])

        HdA = self.tsum(self.wave[self.bandw_HdA],
                        np.divide(spec_cont_HdA - mean_dip, spec_cont_HdA))
        return HdA


class aperture_measurements():
    """
    Aperture/annuli measurements class
    initialized on a fits file.
    """

    def __init__(self, file, z):
        """
        For each datacube, we have:
        NX,NY: the dimensions of the IFU
        flux: in the form of [:,NX,NY]
        x,y: recreating the grid
        radii: gives the radius for any spaxel
        """
        self.drp_logcube = fits.open(file)
        self.z = z
        self.NL, self.NY, self.NX = (self.drp_logcube['FLUX']).data.shape
        self.flux = (self.drp_logcube['FLUX']).data

        """
        IMPORTANT SEMANTICS: The arrays here are in
        spaxel space, i.e.
        1 unit = 0.5" or
        a radius of 3 units is really 1.5" (or 3" aperture)
        """

        self.y = np.outer((np.arange(0, self.NX) + 0.5) - (self.NX/2.0),
                          np.ones(self.NX))
        self.x = np.transpose(self.y)
        self.grid = np.sqrt((self.y*self.y + self.x*self.x))
        self.radii = np.ravel(self.grid)

    def get_radii(self,z_new):
        #gets you the new radii for any redshift
        #radii_new = (1+z_new)*comdis(z_obs)*radii/((1+z_obs)*comdis(z_new))
        radii_new = ((1+z_new)*cosmo.comoving_distance(self.z)*
                     self.radii)/((1+self.z)*cosmo.comoving_distance(z_new))
        return np.array(radii_new)

    def get_spaxels(self, aperture, z_new):
        """
        Getting spaxels within an aperture
        """
        spectra = []
        for i in range(self.NX):
            for j in range(self.NY):
                thing = self.flux[:,i,j]
                spectra.append(thing)
        spectra = np.array(spectra)

        #In case we wanted all the spaxels
        if aperture == -1:
            return spectra
        else:
            index = np.where(self.get_radii(z_new)<= aperture)[0]
            spaxels_within = spectra[index]
            return spaxels_within

    def get_spaxels_annulus(self,aperture,z_new,ring_width):
        """
        Getting spaxels within an annulus
        defined by aperture + ring_width
        """
        spectra = []
        for i in range(self.NX):
            for j in range(self.NY):
                thing = self.flux[:,i,j]
                spectra.append(thing)
        spectra = np.array(spectra)
        index = np.where(self.get_radii(z_new)<= aperture)[0]
        spaxels_within = spectra[index]
        new_radii =self.get_radii(z_new)[index]
        index2 = np.where(new_radii<=aperture+ring_width)[0]
        return spaxels_within[index2]

    def get_spec_measure(self,aperture,z_new,spectral_measure):
        """
        Hda/Dn4000/Halpha/(?) within aperture
        """
        specs = self.get_spaxels(aperture,z_new)
        if spectral_measure == 'hdelta':
            return spec_measurements(self.drp_logcube,self.z).HdA(specs)
        elif spectral_measure == 'dn4000':
            return spec_measurements(self.drp_logcube,self.z).dn_4000(specs)
        elif spectral_measure == 'all':
            return [spec_measurements(self.drp_logcube,self.z).HdA(specs),
                    spec_measurements(self.drp_logcube,self.z).dn_4000(specs)]
        else:
            return False

    def get_spec_ann_measure(self,aperture,z_new,ring_width,spectral_measure):
        """
        Hda/Dn4000/Halpha/(?) within annulus
        """
        specs = self.get_spaxels_annulus(aperture,z_new,ring_width)
        if spectral_measure == 'hdelta':
            return spec_measurements(self.drp_logcube,self.z).HdA(specs)
        elif spectral_measure == 'dn4000':
            return spec_measurements(self.drp_logcube,self.z).dn_4000(specs)
        elif spectral_measure == 'all':
            return [spec_measurements(self.drp_logcube,self.z).HdA(specs),
                    spec_measurements(self.drp_logcube,self.z).dn_4000(specs)]
        else:
            return False
