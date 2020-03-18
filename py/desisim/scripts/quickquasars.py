from __future__ import absolute_import, division, print_function

import sys, os
import argparse
import time

import numpy as np
from scipy.constants import speed_of_light
from scipy.stats import cauchy
from astropy.table import Table,Column
import astropy.io.fits as pyfits
import multiprocessing
import healpy
import yaml

from desiutil.log import get_logger
from desispec.io.util import write_bintable
from desispec.io.fibermap import read_fibermap
from desisim.simexp import reference_conditions
from desisim.templates import SIMQSO, QSO
from desisim.scripts.quickspectra import sim_spectra
from desisim.lya_spectra import read_lya_skewers , apply_lya_transmission, apply_metals_transmission, lambda_RF_LYA
from desisim.dla import dla_spec,insert_dlas
from desisim.bal import BAL
from desisim.io import empty_metatable
from desisim.fpsubsample import footprint_subsample, dataset_subsample, dataset_exptime
from desispec.interpolation import resample_flux

from desimodel.io import load_pixweight
from desimodel import footprint
from speclite import filters
from desitarget.cuts import isQSO_colors
from desiutil.dust import SFDMap, ext_odonnell

try:
    c = speed_of_light/1000. #- km/s
except TypeError:
    #
    # This can happen in documentation builds.
    #
    c = 299792458.0/1000.0

def parse(options=None):
    parser=argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Fast simulation of QSO Lya spectra into the final DESI format\
            (Spectra class) that can be directly used as an input to the redshift fitter\
            (redrock) or correlation function code (picca). The input file is a Lya\
            transmission skewer fits file which format is described in\
            https://desi.lbl.gov/trac/wiki/LymanAlphaWG/LyaSpecSim.")

    #- Required
    parser.add_argument('-i','--infile', type=str, nargs= "*", required=True, help="Input skewer healpix fits file(s)")

    parser.add_argument('-o','--outfile', type=str, required=False, help="Output spectra (only used if single input file)")
    
    parser.add_argument('-c', '--config', default='input.yaml')

    parser.add_argument('--outdir', type=str, default=".", required=False, help="Output directory")

    parser.add_argument('--nproc', type=int, default=1,help="number of processors to run faster")

    parser.add_argument('--overwrite', action = "store_true" ,help="rerun if spectra exists (default is skip)")

<<<<<<< HEAD
<<<<<<< HEAD
=======
    parser.add_argument('--nmax', type=int, default=None, help="Max number of QSO per input file, for debugging")


>>>>>>> 7b9207638d289ad316642d77f9eeeb75a7fb8227
=======
    parser.add_argument('--nmax', type=int, default=None, help="Max number of QSO per input file, for debugging")

>>>>>>> 907d7673a9ab14cbc84164294d045203281f4b87
    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)

    return args


def mod_cauchy(loc,scale,size,cut):
    samples=cauchy.rvs(loc=loc,scale=scale,size=3*size)
    samples=samples[abs(samples)<cut]
    if len(samples)>=size:
       samples=samples[:size]
    else:
        samples=mod_cauchy(loc,scale,size,cut)   ##Only added for the very unlikely case that there are not enough samples after the cut.
    return samples

def get_spectra_filename(args,nside,pixel):
    if args.outfile :
        return args.outfile
    filename="{}/{}/spectra-{}-{}.fits".format(pixel//100,pixel,nside,pixel)
    return os.path.join(args.outdir,filename)


def get_zbest_filename(params,pixdir,nside,pixel):
    if params.get('zbest') :
        return os.path.join(pixdir,"zbest-{}-{}.fits".format(nside,pixel))
    return None


def get_truth_filename(pixdir,nside,pixel):
    return os.path.join(pixdir,"truth-{}-{}.fits".format(nside,pixel))


def is_south(dec):
    """Identify which QSOs are in the south vs the north, since these are on
    different photometric systems. See
    https://github.com/desihub/desitarget/issues/353 for details.

    """
    return dec <= 32.125 # constant-declination cut!


def get_healpix_info(ifilename):
    """Read the header of the tranmission file to find the healpix pixel, nside
    and if we are lucky the scheme. If it fails, try to guess it from the
    filename (for backward compatibility).

    Args:
        ifilename: full path to input transmission file
    Returns:
        healpix: HEALPix pixel corresponding to the file
        nside: HEALPix nside value
        hpxnest: Whether HEALPix scheme in the file was nested
    """

    log = get_logger()

    print('ifilename',ifilename)

    healpix=-1
    nside=-1
    hpxnest=True

    hdulist=pyfits.open(ifilename)
    if "METADATA" in hdulist :
        head=hdulist["METADATA"].header
        for k in ["HPXPIXEL","PIXNUM"] :
            if k in head :
                healpix=int(head[k])
                log.info("healpix={}={}".format(k,healpix))
                break
        for k in ["HPXNSIDE","NSIDE"] :
            if k in head :
                nside=int(head[k])
                log.info("nside={}={}".format(k,nside))
                break
        for k in ["HPXNEST","NESTED","SCHEME"] :
            if k in head :
                if k == "SCHEME" :
                    hpxnest=(head[k]=="NEST")
                else :
                    hpxnest=bool(head[k])
                log.info("hpxnest from {} = {}".format(k,hpxnest))
                break

    if healpix >= 0 and nside < 0 :
        log.error("Read healpix in header but not nside.")
        raise ValueError("Read healpix in header but not nside.")

    if healpix < 0 :
        vals = os.path.basename(ifilename).split(".")[0].split("-")
        if len(vals)<3 :
            error_msg="Could not guess healpix info from {}".format(ifilename)
            log.error(error_msg)
            raise ValueError(error_msg)
        try :
            healpix=int(vals[-1])
            nside=int(vals[-2])
        except ValueError:
            error_msg="Could not guess healpix info from {}".format(ifilename)
            log.error(error_msg)
            raise ValueError(error_msg)
        log.warning("Guessed healpix and nside from filename, assuming the healpix scheme is 'NESTED'")

    return healpix, nside, hpxnest


def get_pixel_seed(pixel, nside, global_seed):
    if global_seed is None:
        # return a random seed
        return np.random.randint(2**32, size=1)[0]
    npix=healpy.nside2npix(nside)
    np.random.seed(global_seed)
    seeds = np.unique(np.random.randint(2**32, size=10*npix))[:npix]
    pixel_seed = seeds[pixel]
    return pixel_seed


def simulate_one_healpix(ifilename,args,params,model,obsconditions,decam_and_wise_filters,
                         bassmzls_and_wise_filters,bal=None,overwrite=False):
    log = get_logger()
    
    
    if 'zfit' in params.keys():
        zbest = params['zfit'].get('zbest')
        sigma_kms_fog = params['zfit'].get('sigma_kms_fog')
        if not sigma_kms_fog:
            sigma_kms_fog=150.
        shift_kms_los = params['zfit'].get('shift-kms-los')
        if not shift_kms_los: 
            shift_kms_los=0. 
        gamma_kms_zfit= params['zfit'].get('gamma_kms_zfit')
        if gamma_kms_zfit and not isinstance(gamma_kms_zfit,(int,float)):
            gamma_kms_zfit=400.
    else:
        zbest=False
        sigma_kms_fog=150.
        shift_kms_los=0.
        gamma_kms_zfit=False
        
    if 'transmission' in params.keys():
        bbfluxparam = params['transmission'].get('bbflux')
        add_LYB = params['transmission'].get('add-LYB')
        metals = params['transmission'].get('metals')
        metals_from_file = params['transmission'].get('metals-from-file')
        dla =  params['transmission'].get('dla')
        balprob =  params['transmission'].get('balprob')
        no_transmission =  params['transmission'].get('no-transmission')
    else:
        bbflux = False
        add_LYB = False
        metals = False
        metals_from_file = False
        dla =  False
        balprob =  False
        no_transmission = False

    # open filename and extract basic HEALPix information
    pixel, nside, hpxnest = get_healpix_info(ifilename)

    # using global seed (could be None) get seed for this particular pixel
    global_seed = params.get('seed')
    seed = get_pixel_seed(pixel, nside, global_seed)
    # use this seed to generate future random numbers
    np.random.seed(seed)

    # get output file (we will write there spectra for this HEALPix pixel)
    ofilename = get_spectra_filename(args,nside,pixel)
    # get directory name (we will also write there zbest file)
    pixdir = os.path.dirname(ofilename)

    # get filename for truth file
    truth_filename = get_truth_filename(pixdir,nside,pixel)

    # get filename for zbest file
    zbest_filename = get_zbest_filename(params,pixdir,nside,pixel)

    if not overwrite :
        # check whether output exists or not
        if zbest :
            if os.path.isfile(ofilename) and os.path.isfile(zbest_filename) :
                log.info("skip existing {} and {}".format(ofilename,zbest_filename))
                return
        else : # only test spectra file
            if os.path.isfile(ofilename) :
                log.info("skip existing {}".format(ofilename))
                return

    # create sub-directories if required
    if len(pixdir)>0 :
        if not os.path.isdir(pixdir) :
            log.info("Creating dir {}".format(pixdir))
            os.makedirs(pixdir)

    log.info("Read skewers in {}, random seed = {}".format(ifilename,seed))

    # Read transmission from files. It might include DLA information, and it
    # might add metal transmission as well (from the HDU file).
    log.info("Read transmission file {}".format(ifilename))

    trans_wave, transmission, metadata, dla_info = read_lya_skewers(ifilename,read_dlas=(dla=='file'),add_metals=metals_from_file,add_lyb=add_LYB)

    ### Add Finger-of-God, before generate the continua
    log.info("Add FOG to redshift with sigma {} to quasar redshift".format(sigma_kms_fog))
    DZ_FOG = sigma_kms_fog/c*(1.+metadata['Z'])*np.random.normal(0,1,metadata['Z'].size)
    metadata['Z'] += DZ_FOG

    ### Select quasar within a given redshift range
    if 'redshift' in params.keys():
        zmin = params['redshift'].get('zmin')
        if not zmin:
            zmin=0.
        zmax = params['redshift'].get('zmax')
        if not zmax:
            zmax=10.
    else:
        zmin,zmax=0.,10.
    w = (metadata['Z']>=zmin) & (metadata['Z']<=zmax)
    transmission = transmission[w]
    metadata = metadata[:][w]
    DZ_FOG = DZ_FOG[w]
      
    fpsubsample = params.get('fpsubsample')
    exptimeparam = params.get('exptime')
    if fpsubsample:
        selection = dataset_subsample(metadata["RA"], metadata["DEC"],metadata["Z"],footprint_subsample(fpsubsample,nside))
        log.info("Select QSOs in DESI subsample footprint {} -> {}".format(transmission.shape[0],selection.size))
        if selection.size == 0 :
            log.warning("No intersection with DESI subsample footprint")
            return
        transmission = transmission[selection]
        metadata = metadata[:][selection]
        DZ_FOG = DZ_FOG[selection]
        
        if not exptimeparam: #ADDED FOR DEBUGGING, might leave it or not.
            exptime = dataset_exptime(metadata["RA"],metadata["DEC"],metadata["Z"],footprint_subsample(fpsubsample,nside))
            obsconditions['EXPTIME']=exptime

    nqso=transmission.shape[0]
    nmax = params.get('nmax')
    if nmax:
        if nmax < nqso :
            log.info("Limit number of QSOs from {} to nmax={} (random subsample)".format(nqso,nmax))
            # take a random subsample
            indices = np.random.choice(np.arange(nqso),nmax,replace=False)
            transmission = transmission[indices]
            metadata = metadata[:][indices]
            DZ_FOG = DZ_FOG[indices]
            nqso = nmax
            
            if fpsubsample:
                obsconditions['EXPTIME']=obsconditions['EXPTIME'][indices]                                      

    # In previous versions of the London mocks we needed to enforce F=1 for
    # z > z_qso here, but this is not needed anymore. Moreover, now we also
    # have metal absorption that implies F < 1 for z > z_qso
    #for ii in range(len(metadata)):
    #    transmission[ii][trans_wave>lambda_RF_LYA*(metadata[ii]['Z']+1)]=1.0

    # if requested, add DLA to the transmission skewers
    if dla:
        # if adding random DLAs, we will need a new random generator
        if dla=='random':
            log.info('Adding DLAs randomly')
            random_state_just_for_dlas = np.random.RandomState(seed)
        elif dla=='file':
            log.info('Adding DLAs from transmission file')
        else:
            log.error("Wrong option for dla: "+dla)
            sys.exit(1)

        # if adding DLAs, the information will be printed here
        dla_filename=os.path.join(pixdir,"dla-{}-{}.fits".format(nside,pixel))
        dla_NHI, dla_z, dla_qid,dla_id = [], [], [],[]

        # identify minimum Lya redshift in transmission files
        min_lya_z = np.min(trans_wave/lambda_RF_LYA - 1)

        # loop over quasars in pixel

        for ii in range(len(metadata)):

            # quasars with z < min_z will not have any DLA in spectrum
            if min_lya_z>metadata['Z'][ii]: continue

            # quasar ID
            idd=metadata['MOCKID'][ii]
            dlas=[]

            if dla=='file':
                for dla in dla_info[dla_info['MOCKID']==idd]:

                    # Adding only DLAs with z < zqso
                    if dla['Z_DLA_RSD']>=metadata['Z'][ii]: continue
                    dlas.append(dict(z=dla['Z_DLA_RSD'],N=dla['N_HI_DLA'],dlaid=dla['DLAID']))
                transmission_dla = dla_spec(trans_wave,dlas)

            elif dla=='random':
                dlas, transmission_dla = insert_dlas(trans_wave, metadata['Z'][ii], rstate=random_state_just_for_dlas)
                for idla in dlas:
                   idla['dlaid']+=idd*1000      #Added to have unique DLA ids. Same format as DLAs from file.

            # multiply transmissions and store information for the DLA file
            if len(dlas)>0:
                transmission[ii] = transmission_dla * transmission[ii]
                dla_z += [idla['z'] for idla in dlas]
                dla_NHI += [idla['N'] for idla in dlas]
                dla_id += [idla['dlaid'] for idla in dlas]
                dla_qid += [idd]*len(dlas)
        log.info('Added {} DLAs'.format(len(dla_id)))
        # write file with DLA information
        if len(dla_id)>0:
            dla_meta=Table()
            dla_meta['NHI'] = dla_NHI
            dla_meta['Z_DLA'] = dla_z  #This is Z_DLA_RSD in transmision.
            dla_meta['TARGETID']=dla_qid
            dla_meta['DLAID'] = dla_id
            hdu_dla = pyfits.convenience.table_to_hdu(dla_meta)
            hdu_dla.name="DLA_META"
            del(dla_meta)
            log.info("DLA metadata to be saved in {}".format(truth_filename))
        else:
            hdu_dla=pyfits.PrimaryHDU()
            hdu_dla.name="DLA_META"

    # if requested, extend transmission skewers to cover full spectrum
    if bbfluxparam:
        wanted_min_wave = 3329. # needed to compute magnitudes for decam2014-r (one could have trimmed the transmission file ...)
        wanted_max_wave = 55501. # needed to compute magnitudes for wise2010-W2

        if trans_wave[0]>wanted_min_wave :
            log.info("Increase wavelength range from {}:{} to {}:{} to compute magnitudes".format(int(trans_wave[0]),int(trans_wave[-1]),int(wanted_min_wave),int(trans_wave[-1])))
            # pad with ones at short wavelength, we assume F = 1 for z <~ 1.7
            # we don't need any wavelength resolution here
            new_trans_wave = np.append([wanted_min_wave,trans_wave[0]-0.01],trans_wave)
            new_transmission = np.ones((transmission.shape[0],new_trans_wave.size))
            new_transmission[:,2:] = transmission
            trans_wave   = new_trans_wave
            transmission = new_transmission

        if trans_wave[-1]<wanted_max_wave :
            log.info("Increase wavelength range from {}:{} to {}:{} to compute magnitudes".format(int(trans_wave[0]),int(trans_wave[-1]),int(trans_wave[0]),int(wanted_max_wave)))
            # pad with ones at long wavelength because we assume F = 1
            coarse_dwave = 2. # we don't care about resolution, we just need a decent QSO spectrum, there is no IGM transmission in this range
            n = int((wanted_max_wave-trans_wave[-1])/coarse_dwave)+1
            new_trans_wave = np.append(trans_wave,np.linspace(trans_wave[-1]+coarse_dwave,trans_wave[-1]+coarse_dwave*(n+1),n))
            new_transmission = np.ones((transmission.shape[0],new_trans_wave.size))
            new_transmission[:,:trans_wave.size] = transmission
            trans_wave   = new_trans_wave
            transmission = new_transmission

    # whether to use QSO or SIMQSO to generate quasar continua.  Simulate
    # spectra in the north vs south separately because they're on different
    # photometric systems.
    no_simqso = params.get('no-simqso')
    south = np.where( is_south(metadata['DEC']) )[0]
    north = np.where( ~is_south(metadata['DEC']) )[0]
    meta, qsometa = empty_metatable(nqso, objtype='QSO', simqso=not no_simqso)
    if no_simqso:
        log.info("Simulate {} QSOs with QSO templates".format(nqso))
        tmp_qso_flux = np.zeros([nqso, len(model.eigenwave)], dtype='f4')
        tmp_qso_wave = np.zeros_like(tmp_qso_flux)
    else:
        log.info("Simulate {} QSOs with SIMQSO templates".format(nqso))
        tmp_qso_flux = np.zeros([nqso, len(model.basewave)], dtype='f4')
        tmp_qso_wave = model.basewave

    for these, issouth in zip( (north, south), (False, True) ):

        # number of quasars in these
        nt = len(these)
        if nt<=0: continue

        
        _tmp_qso_flux, _tmp_qso_wave, _meta, _qsometa \
            = model.make_templates(nmodel=nt,
                redshift=metadata['Z'][these],
                lyaforest=False, nocolorcuts=True,
                noresample=True, seed=seed, south=issouth)

        _meta['TARGETID'] = metadata['MOCKID'][these]
        _qsometa['TARGETID'] = metadata['MOCKID'][these]
        meta[these] = _meta
        qsometa[these] = _qsometa
        tmp_qso_flux[these, :] = _tmp_qso_flux

        if no_simqso:
            tmp_qso_wave[these, :] = _tmp_qso_wave

    log.info("Resample to transmission wavelength grid")
    qso_flux=np.zeros((tmp_qso_flux.shape[0],trans_wave.size))
    if no_simqso:
        for q in range(tmp_qso_flux.shape[0]) :
            qso_flux[q]=np.interp(trans_wave,tmp_qso_wave[q],tmp_qso_flux[q])
    else:
        for q in range(tmp_qso_flux.shape[0]) :
            qso_flux[q]=np.interp(trans_wave,tmp_qso_wave,tmp_qso_flux[q])

    tmp_qso_flux = qso_flux
    tmp_qso_wave = trans_wave

    # if requested, add BAL features to the quasar continua
    if balprob:
        if balprob<=1. and balprob >0:
            log.info("Adding BALs with probability {}".format(balprob))
            # save current random state
            rnd_state = np.random.get_state()
            tmp_qso_flux,meta_bal=bal.insert_bals(tmp_qso_wave,tmp_qso_flux, metadata['Z'],
                                                  balprob=balprob,seed=seed)
            # restore random state to get the same random numbers later
            # as when we don't insert BALs
            np.random.set_state(rnd_state)
            meta_bal['TARGETID'] = metadata['MOCKID']
            w = meta_bal['TEMPLATEID']!=-1
            meta_bal = meta_bal[:][w]
            hdu_bal=pyfits.convenience.table_to_hdu(meta_bal); hdu_bal.name="BAL_META"
            del meta_bal
        else:
            balstr=str(balprob)
            log.error("BAL probability is not between 0 and 1 : "+balstr)
            sys.exit(1)

    # Multiply quasar continua by transmitted flux fraction
    # (at this point transmission file might include Ly-beta, metals and DLAs)
    log.info("Apply transmitted flux fraction")
    if not no_transmission:
        tmp_qso_flux = apply_lya_transmission(tmp_qso_wave,tmp_qso_flux,
                            trans_wave,transmission)

    # if requested, compute metal transmission on the fly
    # (if not included already from the transmission file)
    if metals:
        if metals_from_file :
            log.error('you cannot add metals twice')
            raise ValueError('you cannot add metals twice')
        if no_transmission:
            log.error('you cannot add metals if asking for no-transmission')
            raise ValueError('can not add metals if using no-transmission')
        lstMetals = ''
        for m in metals: lstMetals += m+', '
        log.info("Apply metals: {}".format(lstMetals[:-2]))

        tmp_qso_flux = apply_metals_transmission(tmp_qso_wave,tmp_qso_flux,
                            trans_wave,transmission,metals)

    # if requested, compute magnitudes and apply target selection.  Need to do
    # this calculation separately for QSOs in the north vs south.
    if bbfluxparam:
        bands=['FLUX_G','FLUX_R','FLUX_Z', 'FLUX_W1', 'FLUX_W2']
        bbflux=dict()
        bbflux['SOUTH'] = is_south(metadata['DEC'])
        for band in bands:
            bbflux[band] = np.zeros(nqso)
        # need to recompute the magnitudes to account for lya transmission
        log.info("Compute QSO magnitudes")

        for these, filters in zip( (~bbflux['SOUTH'], bbflux['SOUTH']),
                                   (bassmzls_and_wise_filters, decam_and_wise_filters) ):
            if np.count_nonzero(these) > 0:
                maggies = filters.get_ab_maggies(1e-17 * tmp_qso_flux[these, :], tmp_qso_wave)
                for band, filt in zip( bands, maggies.colnames ):
                    bbflux[band][these] = np.ma.getdata(1e9 * maggies[filt]) # nanomaggies
                    
    log.info("Resample to a linear wavelength grid (needed by DESI sim.)")
    # careful integration of bins, not just a simple interpolation
    if 'wavelength' in params.keys():
        wmin = params['wavelength'].get('wmin',3500.)
        wmax = params['wavelength'].get('wmax',10000.)
        dwave = params['wavelength'].get('wmin',0.2)
    else:
        wmin,wmax,dwave = 3500.,10000.,0.2
    qso_wave=np.linspace(wmin,wmax,int((wmax-wmin)/dwave)+1)
    qso_flux=np.zeros((tmp_qso_flux.shape[0],qso_wave.size))
    for q in range(tmp_qso_flux.shape[0]) :
        qso_flux[q]=resample_flux(qso_wave,tmp_qso_wave,tmp_qso_flux[q])

    log.info("Simulate DESI observation and write output file")
    if "MOCKID" in metadata.dtype.names :
        #log.warning("Using MOCKID as TARGETID")
        targetid=np.array(metadata["MOCKID"]).astype(int)
    elif "ID" in metadata.dtype.names :
        log.warning("Using ID as TARGETID")
        targetid=np.array(metadata["ID"]).astype(int)
    else :
        log.warning("No TARGETID")
        targetid=None

    specmeta={"HPXNSIDE":nside,"HPXPIXEL":pixel, "HPXNEST":hpxnest}

    if bbfluxparam:
        fibermap_columns = dict(
            FLUX_G = bbflux['FLUX_G'],
            FLUX_R = bbflux['FLUX_R'],
            FLUX_Z = bbflux['FLUX_Z'],
            FLUX_W1 = bbflux['FLUX_W1'],
            FLUX_W2 = bbflux['FLUX_W2'],
            )
        photsys = np.full(len(bbflux['FLUX_G']), 'N', dtype='S1')
        photsys[bbflux['SOUTH']] = b'S'
        fibermap_columns['PHOTSYS'] = photsys
    else :
        fibermap_columns=None
    skyerr = params.get('skyerr')
    if not skyerr:
        skyerr=0.
    program = params.get('program')
    sim_spectra(qso_wave,qso_flux,program, obsconditions=obsconditions,spectra_filename=ofilename,
                sourcetype="qso", skyerr=skyerr,ra=metadata["RA"],dec=metadata["DEC"],targetid=targetid,
                meta=specmeta,seed=seed,fibermap_columns=fibermap_columns,use_poisson=False) # use Poisson = False to get reproducible results.

    ### Keep input redshift
    Z_spec = metadata['Z'].copy()
    Z_input = metadata['Z'].copy()-DZ_FOG

    ### Add a shift to the redshift, simulating the systematic imprecision of redrock
    DZ_sys_shift = shift_kms_los/c*(1.+Z_input)
    log.info('Added a shift of {} km/s to the redshift'.format(shift_kms_los))
    meta['REDSHIFT'] += DZ_sys_shift
    metadata['Z'] += DZ_sys_shift

    ### Add a shift to the redshift, simulating the statistic imprecision of redrock
    if gamma_kms_zfit:
        log.info("Added zfit error with gamma {} to zbest".format(gamma_kms_zfit))
        DZ_stat_shift = mod_cauchy(loc=0,scale=gamma_kms_zfit,size=nqso,cut=3000)/c*(1.+Z_input)
        meta['REDSHIFT'] += DZ_stat_shift
        metadata['Z'] += DZ_stat_shift

    ## Write the truth file, including metadata for DLAs and BALs
    log.info('Writing a truth file  {}'.format(truth_filename))
    meta.rename_column('REDSHIFT','Z')
    meta.add_column(Column(Z_spec,name='TRUEZ'))
    meta.add_column(Column(Z_input,name='Z_INPUT'))
    meta.add_column(Column(DZ_FOG,name='DZ_FOG'))
    meta.add_column(Column(DZ_sys_shift,name='DZ_SYS'))
    if fpsubsample:
        if exptimeparam:#Added for debugging and comparision
            exptime=exptime*np.ones(len(metadata['Z']))
        meta.add_column(Column(exptime,name='EXPTIME'))
    if gamma_kms_zfit:
        meta.add_column(Column(DZ_stat_shift,name='DZ_STAT'))
    if 'Z_noRSD' in metadata.dtype.names:
        meta.add_column(Column(metadata['Z_noRSD'],name='Z_NORSD'))
    else:
        log.info('Z_noRSD field not present in transmission file. Z_NORSD not saved to truth file')

    #Save global seed and pixel seed to primary header
    hdr=pyfits.Header()
    hdr['GSEED']=global_seed
    hdr['PIXSEED']=seed
    hdu = pyfits.convenience.table_to_hdu(meta)
    hdu.header['EXTNAME'] = 'TRUTH'
    hduqso=pyfits.convenience.table_to_hdu(qsometa)
    hduqso.header['EXTNAME'] = 'QSO_META'
    hdulist=pyfits.HDUList([pyfits.PrimaryHDU(header=hdr),hdu,hduqso])
    if dla:
        hdulist.append(hdu_dla)
    if balprob:
        hdulist.append(hdu_bal)
    hdulist.writeto(truth_filename, overwrite=True)
    hdulist.close()




    if zbest :
        log.info("Read fibermap")
        fibermap = read_fibermap(ofilename)
        log.info("Writing a zbest file {}".format(zbest_filename))
        columns = [
            ('CHI2', 'f8'),
            ('COEFF', 'f8' , (4,)),
            ('Z', 'f8'),
            ('ZERR', 'f8'),
            ('ZWARN', 'i8'),
            ('SPECTYPE', (str,96)),
            ('SUBTYPE', (str,16)),
            ('TARGETID', 'i8'),
            ('DELTACHI2', 'f8'),
            ('BRICKNAME', (str,8))]
        zbest = Table(np.zeros(nqso, dtype=columns))
        zbest['CHI2'][:] = 0.
        zbest['Z'][:] = metadata['Z']
        zbest['ZERR'][:] = 0.
        zbest['ZWARN'][:] = 0
        zbest['SPECTYPE'][:] = 'QSO'
        zbest['SUBTYPE'][:] = ''
        zbest['TARGETID'][:] = metadata['MOCKID']
        zbest['DELTACHI2'][:] = 25.
        hzbest = pyfits.convenience.table_to_hdu(zbest); hzbest.name='ZBEST'
        hfmap  = pyfits.convenience.table_to_hdu(fibermap);  hfmap.name='FIBERMAP'
        hdulist =pyfits.HDUList([pyfits.PrimaryHDU(),hzbest,hfmap])
        hdulist.writeto(zbest_filename, overwrite=True)
        hdulist.close() # see if this helps with memory issue

def _func(arg) :
    """ Used for multiprocessing.Pool """
    return simulate_one_healpix(**arg)

def main(args=None):

    log = get_logger()
    if isinstance(args, (list, tuple, type(None))):
        args = parse(args)

    if args.outfile is not None and len(args.infile)>1 :
        log.error("Cannot specify single output file with multiple inputs, use --outdir option instead")
        return 1

    if not os.path.isdir(args.outdir) :
        log.info("Creating dir {}".format(args.outdir))
        os.makedirs(args.outdir)
        
    log.info("Reading config file {}".format(args.config))
    with open(args.config, 'r') as cfile:
        params = yaml.safe_load(cfile)
    
    #- Generate obsconditions with args.program, then override as needed
    if 'program' not in params.keys():
        params['program'] = 'DARK'
    obsconditions = reference_conditions[params['program'].upper()]

    if params.get('no-simqso'):
        log.info("Load QSO model")
        model=QSO()
    else:
        log.info("Load SIMQSO model")
        #lya_simqso_model.py is located in $DESISIM/py/desisim/scripts/.
        #Uses a different emmision lines model than the default SIMQSO. 
        #We will update this soon to match with the one used in select_mock_targets. 
        model=SIMQSO(nproc=1,sqmodel='lya_simqso_model')
    decam_and_wise_filters = None
    bassmzls_and_wise_filters = None
    if params['transmission'].get('bbflux') :
        log.info("Load DeCAM and WISE filters for target selection sim.")
        # ToDo @moustakas -- load north/south filters
        decam_and_wise_filters = filters.load_filters('decam2014-g', 'decam2014-r', 'decam2014-z',
                                                      'wise2010-W1', 'wise2010-W2')
        bassmzls_and_wise_filters = filters.load_filters('BASS-g', 'BASS-r', 'MzLS-z',
                                                     'wise2010-W1', 'wise2010-W2')
        
    if 'zfit' in params.keys():
        if params['zfit'].get('gamma_kms_zfit') and not params['zfit'].get('zbest'):
            log.info("Setting --zbest to true as required by --gamma_kms_zfit")
            params['zbest'] = True

    if params['transmission'].get('balprob'):
        bal=BAL()
    else:
        bal=None

    if args.nproc > 1:
        func_args = [ {"ifilename":filename , \
                       "args":args,"args":args,"params":params, "model":model , \
                       "obsconditions":obsconditions , \
                       "decam_and_wise_filters": decam_and_wise_filters , \
                       "bassmzls_and_wise_filters": bassmzls_and_wise_filters , \
                       "bal":bal,"overwrite":args.overwrite\
                   } for i,filename in enumerate(args.infile) ]
        pool = multiprocessing.Pool(args.nproc)
        pool.map(_func, func_args)
    else:
        for i,ifilename in enumerate(args.infile) :
            simulate_one_healpix(ifilename,args,params,model,obsconditions,
                    decam_and_wise_filters,bassmzls_and_wise_filters,
                    bal=bal,overwrite=args.overwrite)
