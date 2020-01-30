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
from astropy.table import Table, Column

def create_subsample_file(cat,outfile=None,nside=16,nest=True):
    maxnumobs=max(cat['NUMOBS'])
    pixels = hp.ang2pix(nside, np.pi/2.-cat['DEC']*np.pi/180.,cat['RA']*np.pi/180.,nest=nest)
    pixarea = hp.pixelfunc.nside2pixarea(nside, degrees=True)
    unique_pixels = np.unique(pixels)
    
    data = np.zeros((unique_pixels.size,3+maxnumobs))
    
    # REDSHIFT MASKS
    lowz = cat['TRUEZ']<1.8
    midz = (cat['TRUEZ']>=1.8)&(cat['TRUEZ']<2.1)
    highz= cat['TRUEZ']>=2.1
    
    bins = np.arange(1,2+maxnumobs)
    for i,pix in enumerate(unique_pixels):
        lowz_qsos = np.count_nonzero(pixels[lowz]==pix)
        midz_qsos=np.count_nonzero(pixels[midz]==pix)
        highz_qsos=np.count_nonzero(pixels[highz]==pix)
        data[i,0]=lowz_qsos/pixarea
        data[i,1]=midz_qsos/pixarea
        data[i,2]=highz_qsos/pixarea
        data[i,3:] = np.histogram(cat['NUMOBS'][highz][pixels[highz]==pix],density=True,bins=bins)[0]
        
    # CREATE A TABLE FROM DATA
    colnames = ['LOWZ_DENS','MIDZ_DENS','HIGHZ_DENS']+[f'PROB_NUMOBS_{i}' for i in range(1,maxnumobs+1)]
    tbl = Table(data,names=colnames)
    tbl.add_column(unique_pixels.astype(int),index=0,name='HPXPIXEL')
    tbl.meta['name'] = 'DESI YEAR 5 HIGH Z QSOS FOOTPRINT'
    tbl.meta['NSIDE'] = nside
    tbl.meta['NEST'] = nest
    if outfile:
        tbl.write(outfile,overwrite=True)
        print(f'Saved {outfile} file')
        return
    else:
        return tbl
    
class footprint_subsample:
    def __init__(self,fname,nside):
        self.nside = nside
        print('DESI subsample, nside',nside)
        self.subsample_pix,self.subsample_dens,self.numobs_prob = self.read_footprint_subsample(fname)
        print('got data from file')
    
    def read_footprint_subsample(self,fname):
        print('Reading subsample footprint from file',fname)
        if not os.path.isfile(fname):
            print('DESI subsample file',fname)
            raise ValueError('file with DESI subsample footprint does not exist')
        data = Table.read(fname)
        pix = data['HPXPIXEL']
        dens = data['LOWZ_DENS','MIDZ_DENS','HIGHZ_DENS']
        prob_numobs = data['PROB_NUMOBS']
        return pix,dens,prob_numobs
    
    def lyaqsos_density(self,ra,dec):
        pixs = hp.ang2pix(self.nside, np.pi/2.-dec*np.pi/180.,ra*np.pi/180.,nest=True)
        thispix = np.unique(pixs)
        index = np.where(self.subsample_pix==thispix)[0] 
        dens=self.subsample_dens[index]
        return dens
    
    def lyaqsos_obsprob(self,ra,dec):
        pixs = hp.ang2pix(self.nside, np.pi/2.-dec*np.pi/180.,ra*np.pi/180.,nest=True)
        thispix = np.unique(pixs)
        index = np.searchsorted(self.subsample_pix,thispix)[0]
        probs = self.numobs_prob[index]
        return probs

def dataset_subsample(ra,dec,Z,subsample_footprint):
    """ Downsample input list of angular positions based on DESI's year 5 footprint
        and input density of all quasars.
    Args:
        ra (ndarray): Right ascension (degrees)
        dec (ndarray): Declination (degrees)
        Z (ndarray): Redshift
    Returns:
        selection (ndarray): mask to apply to downsample input list
    """
    # figure out expected DESI subsample footprint density, 
    N=len(ra)
    nside = subsample_footprint.nside
    density = subsample_footprint.lyaqsos_density(ra,dec)
    
    if len(density)==0: # Return empty array because there is no intersection with subsample footprint
        return np.array([])
        
    pixarea = hp.pixelfunc.nside2pixarea(nside, degrees=True)
    selection=dict()
    index={'LOWZ':np.where(Z<1.8)[0],'MIDZ':np.where((Z>=1.8)&(Z<2.1))[0],'HIGHZ':np.where(Z>=2.1)[0]}
    for whichz in ['LOWZ','MIDZ','HIGHZ']:
        if len(index[whichz])==0:
            selection[whichz]= np.array([],dtype=int)
            continue
        thisN=int(density[f'{whichz}_DENS']*pixarea)
        if len(index[whichz])<thisN:
            thisN=len(index[whichz])
            print(f'Not enough {whichz} quasars in metadata, taking {thisN} availables')
            
        whichIndeces=np.random.choice(index[whichz],size=thisN,replace=False)
        selection[whichz]=whichIndeces
        print(len(selection[whichz]),f'{whichz} quasars selected out of',len(index[whichz]))
    return np.concatenate((selection['LOWZ'],selection['MIDZ'],selection['HIGHZ']))

def dataset_exptime(ra,dec,Z,subsample_footprint):
    """ Array of random observation times based on probabilities on file.
    Args:
        ra (ndarray): Right ascension (degrees)
        dec (ndarray): Declination (degrees)
    Returns:
        exptime (ndarray): array of exposure time for each quasar.
    """
    N=len(ra)
    prob_numobs = subsample_footprint.lyaqsos_obsprob(ra,dec)
    numobs = np.arange(1,len(prob_numobs)+1)
    exptime = 1000*np.random.choice(numobs,size=N,p=prob_numobs)
    exptime[Z<2.1]=1000.
    print(f'Generated {len(exptime[Z<1.8])} Low-z qsos with 1000 exposure time')
    print(f'Generated {len(exptime[(Z>=1.8)&(Z<2.1)])} Mid-z qsos with 1000 exposure time')
    print(f'Generated {len(exptime[Z>=2.1])} High-z qsos with {numobs} exposures with probabilities {prob_numobs}')
    return exptime