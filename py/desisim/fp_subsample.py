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
    maxnumobs=max(cat['NUMOBS'])
    pixels = hp.ang2pix(nside, np.pi/2.-cat['DEC']*np.pi/180.,cat['RA']*np.pi/180.,nest=nest)
    pixarea = hp.pixelfunc.nside2pixarea(nside, degrees=True)
    unique_pixels = np.unique(pixels)
    
    data = np.zeros((unique_pixels.size,3+maxnumobs))
    data[:,0] = unique_pixels
    lowz = (cat['TRUEZ']>=1.8)&(cat['TRUEZ']<2.1)
    highz= cat['TRUEZ']>2.1
    
    bins = np.arange(1,2+maxnumobs)
    for i,pix in zip(range(unique_pixels.size),unique_pixels):
        lowz_qsos=np.count_nonzero(pixels[lowz]==pix)
        highz_qsos=np.count_nonzero(pixels[highz]==pix)
        if(lowz_qsos==0 or highz_qsos==0):
            print(lowz_qsos,highz_qsos,pix)
        data[i,1]=lowz_qsos/pixarea
        data[i,2]=highz_qsos/pixarea
        data[i,3:] = np.histogram(cat['NUMOBS'][highz][pixels[highz]==pix],density=True,bins=bins)[0]
    header = 'DESI YEAR 5 HIGH Z QSOS FOOTPRINT \n\n'
    header+= 'NSIDE = {} \n'.format(nside)
    header+= 'NEST = {} \n'.format(nest)
    header+= f'pixel_index lowz_density [nb/deg^{2}] highz_density [nb/deg^{2}] numobs_probability {list(bins[:-1])}\n'
    if outfile:
        np.savetxt(outfile,data,fmt='%u '+'%.5e '+'%.5e '+'%.8e '*maxnumobs, header=header)
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
    
    def read_footprint_subsample(self,fname,nside):
        print('Reading subsample footprint from file',fname)
        if not os.path.isfile(fname):
            print('DESI subsample file',fname)
            raise ValueError('file with DESI subsample footprint does not exist')
        data = np.loadtxt(fname)
        pix = data[:,0].astype(int)
        dens = {'LOWZ':data[:,1],'HIGHZ':data[:,2]}
        numobs_prob = data[:,3:]
        return pix,dens,numobs_prob
    
    def lyaqsos_density(self,ra,dec):
        pixs = hp.ang2pix(self.nside, np.pi/2.-dec*np.pi/180.,ra*np.pi/180.,nest=True)
        thispix = np.unique(pixs)
        index = np.where(self.subsample_pix==thispix) 
        dens=dict()
        dens['LOWZ']=self.subsample_dens['LOWZ'][index]
        dens['HIGHZ']=self.subsample_dens['HIGHZ'][index]
        return dens
    
    def lyaqsos_obsprob(self,ra,dec):
        pixs = hp.ang2pix(self.nside, np.pi/2.-dec*np.pi/180.,ra*np.pi/180.,nest=True)
        thispix = np.unique(pixs)
        index = np.searchsorted(self.subsample_pix,thispix)[0]
        probs = self.numobs_prob[index]
        return probs

def dataset_subsample(ra,dec,Z,subsample_footprint,nside):
    """ Downsample input list of angular positions based on DESI's year 5 footprint
        and input density of all quasars .
    Args:
        ra (ndarray): Right ascension (degrees)
        dec (ndarray): Declination (degrees)
        Z (ndarray): Redshift
    Returns:
        selection (ndarray): mask to apply to downsample input list
    """
    # figure out expected DESI subsample footprint density, 
    N=len(ra)
    density = subsample_footprint.lyaqsos_density(ra,dec)
    pixarea = hp.pixelfunc.nside2pixarea(nside, degrees=True)
    selection=dict()
    index={'LOWZ':np.where(Z<2.1)[0],'HIGHZ':np.where(Z>=2.1)[0]}
    for whichz in ['LOWZ','HIGHZ']:
        if len(density[whichz])==0: # Return empty array because there is no intersection with subsample footprint
            return np.array([])
        thisN=int(density[whichz]*pixarea)
        if len(index[whichz])<thisN:
            thisN=len(index[whichz])
            print(f'Not enough {whichz} quasars in metadata, taking {thisN} availables')
            
        whichIndeces=np.random.choice(index[whichz],size=thisN,replace=False)
        selection[whichz]=whichIndeces
        print(len(selection[whichz]),f'{whichz} quasars selected out of',len(index[whichz]))
    return np.concatenate((selection['LOWZ'],selection['HIGHZ']))

def dataset_exptime(ra,dec,Z,subsample_footprint):
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
    exptime[Z<2.1]=1000.
    print(f'Generated {len(exptime[Z<2.1])} Low-z qsos with 1000 exposure time')
    print(f'Generated {len(exptime[Z>=2.1])} High-z qsos with {numobs} exposures with probabilities {obs_prob}')
    return exptime