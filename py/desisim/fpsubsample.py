"""
desisim.fpsubsample
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

def create_density_file(cat,outfile=None,nside=16,nest=True):
    """Create a density file containing HPXPIXEL, quasar density,
    and number of observation probability.
    

    Args:
        cat (path): Input quasar catalog
        
    Optional:
        outfile (path): Output density file path, (default==None).
        nside (int): NSIDE for the sky (default==16)
        nest (bool): NEST for the sky (default==True)
        
    Returns:
        If outfile option is set None: Table object with the density 
        and number of observations probability information.
        If outfile option is given: None
    """
    maxnumobs=max(cat['NUMOBS'])
    pixels = hp.ang2pix(nside, np.pi/2.-cat['DEC']*np.pi/180.,cat['RA']*np.pi/180.,nest=nest)
    pixarea = hp.pixelfunc.nside2pixarea(nside, degrees=True)
    unique_pixels = np.unique(pixels)
    
    # GENERATE TABLE AND TABLE META
    colnames = ['HPXPIXEL','LOWZ_DENS','MIDZ_DENS','HIGHZ_DENS']
    tbl = Table(names=colnames,dtype=('i8','f8','f8','f8'))
    tbl.meta['name'] = 'DESI YEAR 5 QSOS FOOTPRINT'
    tbl.meta['comments'] =['This table contains information about the QSO density for each pixel']
    tbl.meta['NSIDE'] = nside
    tbl.meta['NEST'] = nest
    
    # REDSHIFT MASKS
    lowz = cat['TRUEZ']<1.8
    midz = (cat['TRUEZ']>=1.8)&(cat['TRUEZ']<2.1)
    highz= cat['TRUEZ']>=2.1

    
    bins = np.arange(1,2+maxnumobs)
    prob_numobs=np.zeros((unique_pixels.size,maxnumobs))
    for i,pix in enumerate(unique_pixels):
        lowz_dens = np.count_nonzero(pixels[lowz]==pix)/pixarea
        midz_dens = np.count_nonzero(pixels[midz]==pix)/pixarea
        highz_dens= np.count_nonzero(pixels[highz]==pix)/pixarea
        
        prob_numobs[i]=np.histogram(cat['NUMOBS'][highz][pixels[highz]==pix],density=True,bins=bins)[0]
        tbl.add_row((pix,lowz_dens,midz_dens,highz_dens))
    tbl.add_column(Column(prob_numobs,name='PROB_NUMOBS'))  

    if outfile:
        tbl.write(outfile,overwrite=True)
        print(f'Saved {outfile} file')
        return
    else:
        return tbl
    
class footprint_subsample:
    """Class to store useful functions to create a subsample footprint with\
    multiple exposure times.
    
    Args:
        fname (path): Input density file path.
        pixel (int): Input HPXPIXEL to get density and number of observations from.
        nside (int): NSIDE for the sky.
        hpxnest (int): NEST for the sky.
    """
    def __init__(self,fname,pixel,nside,hpxnest):
        self.nside = nside
        self.pixel = pixel
        self.hpxnest = hpxnest
        print('DESI subsample, nside',nside)
        self.subsample_pix,self.subsample_dens,self.numobs_prob = self.read_footprint_subsample(fname)
        print('got data from file')
    
    def read_footprint_subsample(self,fname):
        print('Reading subsample footprint from file',fname)
        if not os.path.isfile(fname):
            print('DESI subsample file',fname)
            raise ValueError('density file with DESI subsample footprint does not exist')
        data = Table.read(fname)
        
        if self.nside != data.meta['NSIDE']:
            raise ValueError('nside from density file does not match nside from transmission file')
        if self.hpxnest != data.meta['NEST']:
            raise ValueError('NEST option from density file does not match NEST from transmission file')
            
        pix = data['HPXPIXEL']
        dens = data['LOWZ_DENS','MIDZ_DENS','HIGHZ_DENS']
        prob_numobs = data['PROB_NUMOBS']
        return pix,dens,prob_numobs
    
    def lyaqsos_density(self):
        thispix = self.pixel
        dens = self.subsample_dens[self.subsample_pix==thispix]
        return dens
    
    def lyaqsos_obsprob(self):
        thispix = self.pixel
        probs = self.numobs_prob[self.subsample_pix==thispix][0]
        return probs

def dataset_subsample(Z,subsample_footprint):
    """ Sample list of quasars based on input density file.
    Args:
        Z (ndarray): Redshift array of objects.
        subsample_footprint: Class object created by footprint_subsample. 
    Returns:
        selection (ndarray): mask to apply to downsample input list
    """
    # figure out expected DESI subsample footprint density, 
    N=len(Z)
    nside = subsample_footprint.nside
    density = subsample_footprint.lyaqsos_density()
    
    if len(density)==0: # Return empty array because there is no intersection with subsample footprint
        return np.array([])
        
    pixarea = hp.pixelfunc.nside2pixarea(nside, degrees=True)
    selection=dict()
    index={'LOWZ':np.where(Z<1.8)[0],'MIDZ':np.where((Z>=1.8)&(Z<2.1))[0],'HIGHZ':np.where(Z>=2.1)[0]}
    for whichz in ['LOWZ','MIDZ','HIGHZ']:
        if len(index[whichz])==0:
            selection[whichz]= np.array([],dtype=int)
            continue
        thisN=np.round(density[f'{whichz}_DENS']*pixarea).astype(int)
        if len(index[whichz])<thisN:
            thisN=len(index[whichz])
            print(f'Not enough {whichz} quasars in metadata, taking {thisN} availables')
            
        whichIndeces=np.random.choice(index[whichz],size=thisN,replace=False)
        selection[whichz]=whichIndeces
        print(len(selection[whichz]),f'{whichz} quasars selected out of',len(index[whichz]))
    return np.concatenate((selection['LOWZ'],selection['MIDZ'],selection['HIGHZ']))

def dataset_exptime(Z,subsample_footprint):
    """ Array of random observation times based on probabilities on density file.
    Args:
        Z (ndarray): Redshift
        subsample_footprint: Class object created by footprint_subsample.
    Returns:
        exptime (ndarray): array of exposure time for each quasar.
    """
    N=len(Z)
    prob_numobs = subsample_footprint.lyaqsos_obsprob()
    numobs = np.arange(1,len(prob_numobs)+1)
    exptime = 1000*np.random.choice(numobs,size=N,p=prob_numobs)
    exptime[Z<2.1]=1000.
    print(f'Generated {len(exptime[Z<1.8])} Low-z qsos with 1000 exposure time')
    print(f'Generated {len(exptime[(Z>=1.8)&(Z<2.1)])} Mid-z qsos with 1000 exposure time')
    print(f'Generated {len(exptime[Z>=2.1])} High-z qsos with {numobs} exposures with probabilities {prob_numobs}')
    return exptime