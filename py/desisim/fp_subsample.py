"""
desisim.fp_subsample
=============

Functions and methods for mimicking a DESI footprint subsample survey with adequate exposure times.

"""
from __future__ import division, print_function
import numpy as np
import os
import healpy
from pkg_resources import resource_filename
import healpy as hp
from astropy.io import fits
from astropy.table import Table

def create_subsample_file(cat,outfile=None,nside=16,nest=True):
    pixels = hp.ang2pix(nside, np.pi/2.-cat['DEC']*np.pi/180.,cat['RA']*np.pi/180.,nest=nest)
    unique_pixels,qsocount = np.unique(pixels,return_counts=True)
    maxnumobs=max(cat['NUMOBS'])
    data = np.zeros((unique_pixels.size,2+maxnumobs))
    data[:,0] = unique_pixels
    pixarea = hp.pixelfunc.nside2pixarea(nside, degrees=True)
    density = qsocount/pixarea
    data[:,1] = density
    bins = np.arange(1,2+maxnumobs)
    for i,pix in zip(range(unique_pixels.size),unique_pixels):
        data[i,2:] = np.histogram(cat['NUMOBS'][pixels==pix],density=True,bins=bins)[0]
    header = 'DESI YEAR 5 HIGH Z QSOS FOOTPRINT \n\n'
    header+= 'NSIDE = {} \n'.format(nside)
    header+= 'NEST = {} \n'.format(nest)
    header+= f'pixel_index density [nb/deg^{2}] numobs_probability {list(bins[:-1])}\n'
    if outfile:
        np.savetxt(outfile,data,fmt='%u '+'%.5e '+'%.8e '*maxnumobs, header=header)
        print(f'Saved {outfile} file')
        return
    else:
        return data
    
class footprint_subsample:
    def __init__(self,fname,nside):
        self.nside = nside
        print('DESI subsample, nside',nside)
        self.subsample_pix,self.subsample_dens,self.numobs_prob = self.read_footprint_subsample(fname,self.nside)
        print('got data from file')
    
    def read_footprint_subsample(self,fname):
        print('Reading subsample footprint from file',fname)
        if not os.path.isfile(fname):
            print('DESI subsample file',fname)
            raise ValueError('file with DESI subsample footprint does not exist')
        data = np.loadtxt(fname)
        pix = data[:,0].astype(int)
        dens = data[:,1]
        numobs_prob = data[:,2:]
        return pix,dens,numobs_prob
    
    def lyaqsos_density(self,ra,dec):
        pixs = hp.ang2pix(self.nside, np.pi/2.-dec*np.pi/180.,ra*np.pi/180.,nest=True)
        thispix = np.unique(pixs)
        index = np.where(self.subsample_pix==thispix)[0]
        dens = self.subsample_dens[index]
        return dens
    
    def lyaqsos_obsprob(self,ra,dec):
        pixs = hp.ang2pix(self.nside, np.pi/2.-dec*np.pi/180.,ra*np.pi/180.,nest=True)
        thispix = np.unique(pixs)
        index = np.searchsorted(self.subsample_pix,thispix)[0]
        probs = self.numobs_prob[index]
        return probs

def dataset_subsample(ra,dec,input_density,subsample_footprint,nside):
    """ Downsample input list of angular positions based on DESI's year 5 footprint
        and input density of all quasars .
    Args:
        ra (ndarray): Right ascension (degrees)
        dec (ndarray): Declination (degrees)
        input_density (float): Input density of all quasars per sq.deg.
    Returns:
        selection (ndarray): mask to apply to downsample input list
    """
    # figure out expected DESI subsample footprint density, 
    # in quasars / sq.deg., at high-z
    N=len(ra)
    density = subsample_footprint.lyaqsos_density(ra,dec)
    pixarea = hp.pixelfunc.nside2pixarea(nside, degrees=True)
    N_desired = int(density*pixarea)
    index=np.arange(N)
    selection = np.random.choice(index,size=N_desired,replace=False)
    print(len(selection),'selected out of',N)
    return selection

def dataset_exptime(ra,dec,subsample_footprint):
    """ Array of random observation times based on probabilities on file.
    Args:
        ra (ndarray): Right ascension (degrees)
        dec (ndarray): Declination (degrees)
    Returns:
        exptime (ndarray): array of exposure time for each quasar.
    """
    N=len(ra)
    obs_prob = subsample_footprint.lyaqsos_obsprob(ra,dec)
    numobs = np.arange(1,len(obs_prob)+1)
    exptime = 1000*np.random.choice(numobs,size=N,p=obs_prob)
    print(f'Generated {N} qsos with {numobs} exposures with probabilities {obs_prob}')
    return exptime