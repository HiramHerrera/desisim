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
        self.fname = fname
        self.nside = nside
        self.pixel = pixel
        self.hpxnest = hpxnest
        print('DESI subsample, nside',nside)
        self.subsample_pix,self.subsample_dens,self.numobs_prob = self.read_footprint_subsample()
        print('got data from file')
    
    def read_footprint_subsample(self):
        print('Reading subsample footprint from file',self.fname)
        if not os.path.isfile(self.fname):
            print('DESI subsample file',self.fname)
            raise ValueError('density file with DESI subsample footprint does not exist')
        data = Table.read(self.fname)
        
        if self.nside != data.meta['NSIDE']:
            raise ValueError('nside from density file does not match nside from transmission file')
        if self.hpxnest != data.meta['NEST']:
            raise ValueError('NEST option from density file does not match NEST from transmission file')
            
        pix = data['HPXPIXEL']
        dens = data['MIDZ_DENS','HIGHZ_DENS']
        prob_numobs = data['PROB_NUMOBS']
        return pix,dens,prob_numobs
    
    def density(self):
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
    density = subsample_footprint.density()
    
    pixarea = hp.pixelfunc.nside2pixarea(nside, degrees=True)
    selection=np.array([],dtype=int)
    index={'MIDZ':np.where((Z>=1.8)&(Z<2.1))[0],'HIGHZ':np.where(Z>=2.1)[0]}
    if len(density)!=0:
        for whichz in ['MIDZ','HIGHZ']:
            thisN=np.round(density[f'{whichz}_DENS']*pixarea).astype(int)
            if len(index[whichz])==0 or thisN==0:
                continue   
            if len(index[whichz])<thisN:
                thisN=len(index[whichz])
                print(f'Not enough {whichz} quasars in metadata, taking {thisN} availables')

            whichIndices=np.random.choice(index[whichz],size=thisN,replace=False)
            print(f'{len(whichIndices)} {whichz} quasars selected out of {len(index[whichz])}')
            selection = np.concatenate((selection,whichIndices))
    return selection

def dataset_exptime(Z,subsample_footprint):
    """ Array of random observation times based on probabilities on density file.
    Args:
        Z (ndarray): Redshift
        subsample_footprint: Class object created by footprint_subsample.
    Returns:
        exptime (ndarray): array of exposure time for each quasar.
    """
    N=len(Z)
    exptime = np.full(len(Z),1000)
    if np.count_nonzero(Z>2.1)!=0:
        prob_numobs = subsample_footprint.lyaqsos_obsprob()
        numobs = np.arange(1,len(prob_numobs)+1)
        N_highz= np.count_nonzero(Z>2.1)
        print(f'Assigning {numobs} exposures with probabilities {prob_numobs} to {N_highz} High-z qsos')
        exp_highz = 1000*np.random.choice(numobs,size=N_highz,p=prob_numobs)
        
        
        # HARDCODED FOR MAIN SURVEY
        #w=exp_highz==1000
        #exptime_shift=np.random.choice([600.,  800., 1000., 1200., 1400., 1600., 1800.],size=sum(w),p=[0.04, 0.27, 0.29, 0.22, 0.11, 0.06, 0.01])
        #exp_highz[w]=exptime_shift
        
        #w=exp_highz==2000
        #exptime_shift=np.random.choice([1400., 1600., 1800., 2000., 2200., 2400., 2600., 2800.],size=sum(w),p=[0.12, 0.07, 0.08, 0.30, 0.15, 0.14, 0.09, 0.05])
        #exp_highz[w]=exptime_shift
        
        exptime[Z>2.1]=exp_highz
    return exptime