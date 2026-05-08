import sys, os
import argparse
import time
import warnings

import numpy as np
from scipy.constants import speed_of_light
from astropy.table import Table,Column
import astropy.io.fits as pyfits
import multiprocessing
import healpy

import desisim
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
from desisim.survey_release import SurveyRelease
from desispec.interpolation import resample_flux

from desimodel.io import load_pixweight
from speclite import filters

try:
    c_kms = speed_of_light/1000. #- km/s
except TypeError:
    #
    # This can happen in documentation builds.
    #
    c_kms = 299792458.0/1000.0

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

    parser.add_argument('--outdir', type=str, default=".", required=False, help="Output directory")

    parser.add_argument('--seed', type=int, default=None, required = False, help="Global random seed (will be used to generate a seed per each file")

    parser.add_argument('--sigma_kms_fog',type=float,default=150, help="Adds a gaussian error to the quasar \
        redshift that simulate the fingers of god effect")

    parser.add_argument('--add-LYB', action='store_true', help = "Add LYB absorption from transmision file")

    parser.add_argument('--metals', type=str, default=None, required=False, help = "list of metals to get the\
        transmission from, if 'all' runs on all metals", nargs='*')

    parser.add_argument('--metals-from-file',type=str,const='all',help = "list of metals,'SI1260,SI1207' etc, to get from HDUs in file. \
Use 'all' or no argument for mock version < 7.3 or final metal runs. ",nargs='?')

    parser.add_argument('--dla',type=str,required=False, help="Add DLA to simulated spectra either randonmly\
        (--dla random) or from transmision file (--dla file)")

    parser.add_argument('--dlaplus',action='store_true', help="Add absorption in higher order Lyman lines for DLAs")

    parser.add_argument('--balprob',type=float,required=False, help="To add BAL features with the specified probability\
        (e.g --balprob 0.5). Expect a number between 0 and 1 ")

    parser.add_argument('--no-simqso',action = "store_true", help="Does not use desisim.templates.SIMQSO\
        to generate templates, and uses desisim.templates.QSO instead.")

    parser.add_argument('--save-continuum',action = "store_true", help="Save true continum to file")

    parser.add_argument('--no-transmission',action = 'store_true', help='Do not multiply continuum\
        by transmission, use F=1 everywhere')

    parser.add_argument('--nproc', type=int, default=1,help="number of processors to run faster")

    parser.add_argument('--overwrite', action = "store_true" ,help="rerun if spectra exists (default is skip)")

    parser.add_argument('--nmax', type=int, default=None, help="Max number of QSO per input file, for debugging")

    parser.add_argument('--source-contr-smoothing', type=float, default=10., \
        help="When this argument > 0 A, source electrons' contribution to the noise is smoothed " \
        "by a Gaussian kernel using FFT. Pipeline does this by 10 A. " \
        "Larger smoothing might be needed for better decoupling. Does not apply to eBOSS mocks.")
    
    parser.add_argument('--from-catalog', type=str, default=None, help="Input catalog of mock objects to simulate")
    
    parser.add_argument('--metal-strengths-label', type=str, default=None, choices=["lyacolore", "saclay", "lyacolore-2lpt"],
                        help="Use specific metal tuning (ignored if --metal-strengths is provided)")

    parser.add_argument('--metal-strengths', type=float, default=None, required=False, help = "list of strengths to appply\
        to metals. Should correspond to the --metals flag", nargs='*')


    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)

    return args


def get_spectra_filename(args,nside,pixel):
    if args.outfile :
        return args.outfile
    filename="{}/{}/spectra-{}-{}.fits".format(pixel//100,pixel,nside,pixel)
    return os.path.join(args.outdir,filename)


def get_zbest_filename(pixdir,nside,pixel):
    return os.path.join(pixdir,"zbest-{}-{}.fits".format(nside,pixel))


def get_truth_filename(pixdir,nside,pixel):
    return os.path.join(pixdir,"truth-{}-{}.fits".format(nside,pixel))


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
    hdulist.close()

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


def simulate_one_healpix(ifilename, args, model, bal=None):
    log = get_logger()

    # open filename and extract basic HEALPix information
    pixel, nside, hpxnest = get_healpix_info(ifilename)

    # using global seed (could be None) get seed for this particular pixel
    global_seed = args.seed
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
    zbest_filename = get_zbest_filename(pixdir,nside,pixel)

    if not args.overwrite :
        # check whether output exists or not
        if os.path.isfile(ofilename) and os.path.isfile(zbest_filename) :
            log.info("skip existing {} and {}".format(ofilename,zbest_filename))
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

    trans_wave, transmission, metadata, dla_info = read_lya_skewers(
            ifilename,
            read_dlas=(args.dla=='file'),
            add_metals=args.metals_from_file,
            add_lyb=args.add_LYB
    )

    ### Add Finger-of-God, before generate the continua
    log.info("Add FOG to redshift with sigma {} to quasar redshift".format(args.sigma_kms_fog))
    DZ_FOG = args.sigma_kms_fog/c_kms*(1.+metadata['Z'])*np.random.normal(0,1,metadata['Z'].size)
    metadata['Z'] += DZ_FOG

    #- Reference observing conditions for DARK program
    obsconditions = reference_conditions["DARK"]

    if args.from_catalog is None:
        raise ValueError("Quickquasars needs an input catalog now")
    else:
        log.info(f"Getting objects from catalog {args.from_catalog}")
        catalog = Table.read(args.from_catalog)
        # Get mockobjs in pixel that are in the catalog.
        selection = np.isin(metadata['MOCKID'],catalog['MOCKID'])
        if selection.sum()==0:
            log.warning(f'No intersection with catalog')
            return
        log.info(f'Catalog has {selection.sum()} QSOs in pixel {pixel}')
        transmission = transmission[selection]
        metadata = metadata[:][selection]
        DZ_FOG = DZ_FOG[selection] 
        
        this_pixel_targets = np.isin(catalog['MOCKID'],metadata['MOCKID'])
        catalog = catalog[this_pixel_targets]
        ids_index = [np.where(catalog['MOCKID']==mockid)[0][0] for mockid in metadata['MOCKID']]
        # Ensure catalog and metadata MOCKIDS are in the same order
        catalog = catalog[ids_index]
        if 'EXPTIME' in catalog.colnames:
            exptime=np.array(catalog['EXPTIME'])
            obsconditions['EXPTIME']=exptime
        else:
            raise ValueError("Input catalog in --from-catalog must have EXPTIME column")

        # Get magnitude from catalog 
        if 'FLUX_R' in catalog.colnames:
            mags = 22.5-2.5*np.log10(catalog['FLUX_R'])
        elif 'MAG_R' in catalog.colnames:
            mags=catalog['MAG_R']
        else:
            raise ValueError("Input catalog in --from-catalog must have FLUX_R or MAG_R column")

    nqso=transmission.shape[0]
    
    if args.nmax is not None :
        if args.nmax < nqso :
            log.info("Limit number of QSOs from {} to nmax={} (random subsample)".format(nqso,args.nmax))
            # take a random subsample
            indices = np.random.choice(np.arange(nqso),args.nmax,replace=False)
            transmission = transmission[indices]
            metadata = metadata[:][indices]
            DZ_FOG = DZ_FOG[indices]
            nqso = args.nmax

    # Generate quasar continua with QSO or SIMQSO
    log.info("Simulate {} quasar continua".format(nqso))
    tmp_qso_flux, tmp_qso_wave, meta, qsometa = model.make_templates(
            nmodel=nqso,
            redshift=metadata['Z'],
            mag=mags,
            lyaforest=False,
            nocolorcuts=True,
            noresample=True,
            seed=seed,
            south=True
    )
    meta['TARGETID'] = metadata['MOCKID']
    qsometa['TARGETID'] = metadata['MOCKID']

    # resample continua into Lya transmission skewers
    log.info("Resample to transmission wavelength grid")
    qso_flux=np.zeros((tmp_qso_flux.shape[0],trans_wave.size))
    if args.no_simqso:
        for q in range(tmp_qso_flux.shape[0]) :
            qso_flux[q]=np.interp(trans_wave,tmp_qso_wave[q],tmp_qso_flux[q])
    else:
        for q in range(tmp_qso_flux.shape[0]) :
            qso_flux[q]=np.interp(trans_wave,tmp_qso_wave,tmp_qso_flux[q])

    tmp_qso_flux = qso_flux
    tmp_qso_wave = trans_wave

    if args.save_continuum :
        true_wmin=3500.0
        true_wmax=10000.0
        true_dwave=2.0
        true_wave=np.linspace(true_wmin, true_wmax, int((true_wmax - true_wmin)/true_dwave)+1)
        true_flux=np.zeros((tmp_qso_flux.shape[0],true_wave.size))
        for q in range(tmp_qso_flux.shape[0]) :
            true_flux[q]=resample_flux(true_wave,tmp_qso_wave,tmp_qso_flux[q])
        continum_meta=Table()
        continum_meta['TARGETID'] = qsometa['TARGETID']
        continum_meta['TRUE_CONT'] = true_flux
        hdu_trueCont = pyfits.convenience.table_to_hdu(continum_meta)
        hdu_trueCont.name = "TRUE_CONT"
        hdu_trueCont.header['wmin'] = true_wmin
        hdu_trueCont.header['wmax'] = true_wmax
        hdu_trueCont.header['dwave'] = true_dwave

        del(continum_meta,true_wave,true_flux)
        log.info("True continum to be saved in {}".format(truth_filename))

    # if requested, add BAL features to the quasar continua
    if args.balprob:
        if args.balprob <= 1. and args.balprob > 0:
            from desisim.io import find_basis_template
            log.info("Adding BALs with probability {}".format(args.balprob))
            # save current random state
            rnd_state = np.random.get_state()
            tmp_qso_flux,meta_bal = bal.insert_bals(tmp_qso_wave, tmp_qso_flux, metadata['Z'],
                                                  balprob= args.balprob, seed=seed, qsoid=metadata['MOCKID'])
            # restore random state to get the same random numbers later
            # as when we don't insert BALs
            np.random.set_state(rnd_state)
            w = np.isin(qsometa['TARGETID'], meta_bal['TARGETID'])
            qsometa['BAL_TEMPLATEID'][w] = meta_bal['BAL_TEMPLATEID']
            hdu_bal=pyfits.convenience.table_to_hdu(meta_bal); hdu_bal.name="BAL_META"
            #Trim to only show the version, assuming it is located in os.environ['DESI_BASIS_TEMPLATES']
            hdu_bal.header["BALTEMPL"]=find_basis_template(objtype='BAL').split('basis_templates/')[1]
            del meta_bal
        else:
            balstr=str(args.balprob)
            log.error("BAL probability is not between 0 and 1 : "+balstr)
            sys.exit(1)

    if not args.no_transmission:
        # Do not add Lya absorption 
        log.info("Applying Lya transmission")
        tmp_qso_flux = apply_lya_transmission(
            tmp_qso_wave, tmp_qso_flux, trans_wave, transmission)
    else:
        log.info("Skipping Lya transmission (--no-transmission mode)")

    # Apply metals if requested
    if args.metals is not None:
        if args.metals_from_file is not None:
            raise ValueError("Cannot use --metals and --metals-from-file at the same time.")
        
        if args.no_transmission:
            log.warning("--metals with --no-transmission is not supported. Skipping metals.")
        else:
            lstMetals = ', '.join(args.metals)
            log.info("Applying metals: {}".format(lstMetals))
            tmp_qso_flux = apply_metals_transmission(
                tmp_qso_wave, tmp_qso_flux,
                trans_wave, transmission,
                args.metals,
                strengths_label=args.metal_strengths_label,
                strengths=args.metal_strengths
            )
            
    # if requested, add DLA to the transmission skewers
    if args.dla is not None :
        # Initialize HCD-only transmission array (ones = no absorption by default)
        transmission_hcds = np.ones_like(tmp_qso_flux)
        # if adding random DLAs, we will need a new random generator
        if args.dla=='random':
            log.info('Adding DLAs randomly')
            random_state_just_for_dlas = np.random.RandomState(seed)
        elif args.dla=='file':
            log.info('Adding DLAs from transmission file')
        else:
            log.error("Wrong option for args.dla: "+args.dla)
            sys.exit(1)

        # if adding DLAs, the information will be printed here
        dla_filename=os.path.join(pixdir,"dla-{}-{}.fits".format(nside,pixel))
        dla_NHI, dla_z, dla_qid,dla_id = [], [], [],[]

        # Identify the minimum redshift at which a DLA can be detected in the spectrum, given the wavelength range of the spectrum.
        wmin = np.min(tmp_qso_wave)
        min_lya_z = np.min(wmin/lambda_RF_LYA - 1)

        # loop over quasars in pixel
        for ii in range(len(metadata)):
            # quasars with z < min_z will not have any DLA in spectrum
            if min_lya_z>metadata['Z'][ii]: continue
            # quasar ID
            idd=metadata['MOCKID'][ii]
            dlas=[]
            if args.dla=='file':
                for dla in dla_info[dla_info['MOCKID']==idd]:
                    # Adding only DLAs with z < zqso
                    if dla['Z_DLA_RSD']>=metadata['Z'][ii]: continue
                    dlas.append(dict(z=dla['Z_DLA_RSD'],N=dla['N_HI_DLA'],dlaid=dla['DLAID']))
                transmission_dla = dla_spec(tmp_qso_wave,dlas,dlaplus=args.dlaplus)

            elif args.dla=='random':
                dlas, transmission_dla = insert_dlas(tmp_qso_wave, metadata['Z'][ii], rstate=random_state_just_for_dlas, dlaplus=args.dlaplus)
                for idla in dlas:
                   idla['dlaid']+=idd*1000      #Added to have unique DLA ids. Same format as DLAs from file.

            # multiply transmissions and store information for the DLA file
            if len(dlas)>0:
                transmission_hcds[ii] = transmission_dla  # store DLA-only transmission
                dla_z += [idla['z'] for idla in dlas]
                dla_NHI += [idla['N'] for idla in dlas]
                dla_id += [idla['dlaid'] for idla in dlas]
                dla_qid += [idd]*len(dlas)
        tmp_qso_flux *= transmission_hcds  # apply DLA absorption to the flux
        log.info('Added {} DLAs'.format(len(dla_id)))
        if args.dlaplus:
            log.info('Added higher order lines to DLAs')
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
        
        # Add comment to DLA_META extension if dlaplus is enabled
        if args.dlaplus:
            hdu_dla.header['COMMENT'] = 'Added higher order lines to DLAs'
                         
    log.info("Simulate DESI observation and write output file")
    if "MOCKID" in metadata.dtype.names :
        targetid=np.array(metadata["MOCKID"]).astype(int)
    elif "ID" in metadata.dtype.names :
        log.warning("Using ID as TARGETID")
        targetid=np.array(metadata["ID"]).astype(int)
    else :
        log.warning("No TARGETID")
        targetid=None

    log.info("Resample to a linear wavelength grid (needed by DESI sim.)")
    # careful integration of bins, not just a simple interpolation (hardcoded values)
    dwave=0.2
    wmin=3500.0
    wmax=10000.0
    qso_wave=np.linspace(wmin, wmax, int((wmax-wmin)/dwave)+1)
    qso_flux=np.zeros((tmp_qso_flux.shape[0],qso_wave.size))
    for q in range(tmp_qso_flux.shape[0]) :
        qso_flux[q]=resample_flux(qso_wave,tmp_qso_wave,tmp_qso_flux[q])

    ### use Poisson = False to get reproducible results.
    specmeta={"HPXNSIDE":nside,"HPXPIXEL":pixel, "HPXNEST":hpxnest}
    resolution=sim_spectra(qso_wave, qso_flux, program="DARK", obsconditions=obsconditions, spectra_filename=ofilename,
                           sourcetype="qso", ra=metadata["RA"], dec=metadata["DEC"], targetid=targetid,
                           meta=specmeta, seed=seed, fibermap_columns=None, use_poisson=False,
                           specsim_config_file='desi', dwave_out=0.8,
                           save_resolution=False, source_contribution_smoothing=args.source_contr_smoothing)

    ### Keep input redshift
    Z_spec = metadata['Z'].copy()
    Z_input = metadata['Z'].copy()-DZ_FOG

    ## Write the truth file, including metadata for DLAs and BALs
    log.info('Writing a truth file  {}'.format(truth_filename))
    meta.rename_column('REDSHIFT','Z')
    meta.add_column(Column(Z_spec,name='TRUEZ'))
    meta.add_column(Column(Z_input,name='Z_INPUT'))
    meta.add_column(Column(DZ_FOG,name='DZ_FOG'))
    if 'Z_noRSD' in metadata.dtype.names:
        meta.add_column(Column(metadata['Z_noRSD'],name='Z_NORSD'))
    else:
        log.info('Z_noRSD field not present in transmission file. Z_NORSD not saved to truth file')

    #Save global seed and pixel seed to primary header
    hdr=pyfits.Header()
    hdr['GSEED']=global_seed
    hdr['PIXSEED']=seed
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=".*nanomaggies.*")
        hdu = pyfits.convenience.table_to_hdu(meta)

    hdu.header['EXTNAME'] = 'TRUTH'
    hduqso=pyfits.convenience.table_to_hdu(qsometa)
    hduqso.header['EXTNAME'] = 'TRUTH_QSO'
    hdulist=pyfits.HDUList([pyfits.PrimaryHDU(header=hdr),hdu,hduqso])

    # save DLA, BAL and true continuum metadata
    if args.dla :
        hdulist.append(hdu_dla)
    if  args.balprob :
        hdulist.append(hdu_bal)
    if args.save_continuum :
        hdulist.append(hdu_trueCont)

    # Save one resolution matrix per camera to the truth file instead of one per quasar to the spectra files.
    for band in resolution.keys():
        hdu = pyfits.ImageHDU(name="{}_RESOLUTION".format(band.upper()))
        hdu.data = resolution[band].astype("f4")
        hdulist.append(hdu)

    hdulist.writeto(truth_filename, overwrite=True)
    hdulist.close()

    # write zbest file
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

    if args.no_simqso:
        log.info("Load QSO model")
        model=QSO()
    else:
        log.info("Load SIMQSO model")
        #lya_simqso_model.py is located in $DESISIM/py/desisim/scripts/.
        #Uses a different emmision lines model than the default SIMQSO. 
        #We will update this soon to match with the one used in select_mock_targets. 
        model=SIMQSO(nproc=1,sqmodel='lya_simqso_model')

    if args.balprob:
        bal=BAL()
    else:
        bal=None

    if args.nproc > 1:
        func_args = [ {"ifilename":filename, "args":args, "model":model, \
                       "bal":bal} for i,filename in enumerate(args.infile) ]
        with multiprocessing.Pool(args.nproc) as pool:
            pool.map(_func, func_args)
    else:
        for i,ifilename in enumerate(args.infile) :
            simulate_one_healpix(ifilename, args, model, bal=bal)
