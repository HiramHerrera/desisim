from desisim.io import empty_metatable
from speclite import filters
from copy import copy
from speclite.filters import load_filter
from redrock.utils import native_endian
from desiutil.log import get_logger
from desispec.interpolation import resample_flux
import numpy as np
from astropy.io import fits
from astropy import units

fluxunits = 1e-17 * units.erg / units.s / units.cm**2 / units.Angstrom


# AR read redrock template
def read_rrtemplate(fn, row):
    """
    Read a redrock-formatted template.

    Args:
        fn: template file name (str)
        row: zero-indexed row of the template to extract (int)

    Returns:
        ws: wavelengths (1d np array of floats)
        fs: fluxes (1d np array of floats)
    """
    fx = fits.open(fn)
    hdr = fx["BASIS_VECTORS"].header
    ws = np.asarray(
        hdr["CRVAL1"] + hdr["CDELT1"] * np.arange(hdr["NAXIS1"]), dtype=np.float64
    )
    if "LOGLAM" in hdr and hdr["LOGLAM"] != 0:
        ws = 10**ws
    fs = np.asarray(native_endian(fx["BASIS_VECTORS"].data), dtype=np.float64)
    fs = fs[row, :]

    # DISCARD invalid flux
    mask =~ np.isnan(fs)
    mask &= fs>0
    ws = ws[mask]
    fs = fs[mask]
    return ws, fs


def rescale_template2mag(ws, fs, mag, mag_band):
    """
    Return spectra normalized to input magnitudes.

    Args:
        ws: wavelengths (1d array of floats)
        fs: fluxes (1d array of floats)
        mag: desired magnitude (float)
        band: a band loaded with speclite.filters.load_filter() (str)
        verbose (optional, defaults to False): print infos on the prompt (bool)

    Returns:
        rescale_fs: rescaled fluxes

    Notes:
        Similar to desispec.fluxcalibration.normalize_templates,
            but handles padding and possibility to use other filter..
    """
    # AR filter
    filter_response = load_filter(mag_band)
    # AR zero-padding spectrum so that it covers the filter response
    pad_ws, pad_fs = ws.copy(), fs.copy()
    if (pad_ws.min() > filter_response.wavelength.min()) | (
        pad_ws.max() < filter_response.wavelength.max()
    ):
        pad_fs, pad_ws = filter_response.pad_spectrum(pad_fs, pad_ws, method="zero")
    # AR now scale
    orig_mag = filter_response.get_ab_magnitude(pad_fs, pad_ws)
    scalefac = 10 ** ((orig_mag - mag) / 2.5)
    rescale_fs = fs * scalefac
    return rescale_fs

def template_rf2z(rf_ws, rf_fs, z, mag, mag_band):
    """
    Redshift + rescale a template spectrum.

    Args:
        rf_ws: rest-frame wavelengths (1d array of floats)
        rf_fs: rest-frame fluxes (1d array of floats)
        cameras_ws: list of cameras (list of strs)
        z: redshift to redshift to
        mag: requested magnitude (float)
        mag_band: band for the requested magnitude (str)

    Returns:
        cameras_fs: dictionary, with the redshifted+rescaled flux
            (1d array of floats) for each camera (dictionary)
    """
    # AR first rescale to the requested magnitude
    ws = (1 + z) * rf_ws
    fs = rescale_template2mag(ws, rf_fs, mag, mag_band)
    return ws,fs

class LBG(object):
    def __init__(self, objtype='LBG', minwave=3600.0, maxwave=10000.0, cdelt=0.2, wave=None,
                 template_fn='/global/cfs/cdirs/desi/users/cmv/lbglae/RRTemplates/lbgV33sp-templates/rrtemplate-lbgV33sp.fits', template_id=2,
                normfilter_north='BASS-r', normfilter_south='decam2014-r'):
        self.objtype=objtype.upper()
        self.normfilter_north = normfilter_north
        self.normfilter_south = normfilter_south


        self.template_fn = template_fn
        self.template_id = template_id
        basewave, baseflux = read_rrtemplate(self.template_fn, row=template_id)
        self.baseflux = baseflux
        self.basewave = basewave
        self.basemeta = None

        # already provided.
        if wave is None:
            npix = int(round((maxwave-minwave) / cdelt))+1
            wave = np.linspace(minwave, maxwave, npix)
        self.wave = wave

        # Initialize the filter profiles.
        self.normfilt_north = filters.load_filters(self.normfilter_north)
        self.normfilt_south = filters.load_filters(self.normfilter_south)
        self.sdssrfilt = filters.load_filters('sdss2010-r')
        self.decamwise = filters.load_filters('decam2014-g', 'decam2014-r', 'decam2014-z')
        self.bassmzlswise = filters.load_filters('BASS-g', 'BASS-r', 'MzLS-z')


    def make_templates(self, nmodel, redshift, mag,  south=True, restframe=False, noresample=False):
    # Not random templates. i.e all objects get the same template, just re-scaled and shifted
    # TODO: More realistic approach using a random template?
        log = get_logger()
        # Read the rest-frame continuum basis, if not specified.
        meta, objmeta = empty_metatable(nmodel=nmodel,objtype=self.objtype,simqso=False)
                
        if len(redshift) != nmodel:
            log.fatal('Redshift must be an nmodel-length array')
            raise ValueError

        if len(mag) != nmodel:
            log.fatal('Mag must be an nmodel-length array')
            raise ValueError

        if noresample or restframe:
            outflux = np.zeros([nmodel, len(self.basewave)])
            outwave = self.basewave.copy()
        else:
            outflux = np.zeros([nmodel, len(self.wave)])
            outwave = self.wave.copy()

        for ii,(z,r) in enumerate(zip(redshift,mag)):
            if restframe:
                _, flux = template_rf2z(self.basewave, self.baseflux, z, r, 'lsst2023-r')
                zwave=self.basewave.copy()
            else:
                zwave, flux = template_rf2z(self.basewave, self.baseflux, z, r, 'lsst2023-r')
                # Shift template to correct redshift position in basewave 
                flux = resample_flux(self.basewave, zwave, flux, extrapolate=True)


            if noresample or restframe:
                outflux[ii, :] = flux
            else:
                outflux[ii, :] = resample_flux(self.wave, self.basewave, flux, extrapolate=True)
        meta['FLUX_R'] = 10**((22.5-mag)/2.5)
        meta['TEMPLATEID'] = self.template_id
        meta['REDSHIFT'] = redshift
        return 1e17 * outflux, outwave, meta, objmeta
        # Synthesize photometry to determine which models will pass the
        # color cuts.
        print(self.basewave)
        if south:
            magfilt = self.normfilt_south
            magfilter = self.normfilter_south
            maggies = self.decamwise.get_ab_maggies(outflux, self.basewave.copy(), mask_invalid=True)
        else:
            magfilt = self.normfilt_north
            magfilter = self.normfilter_north
            maggies = self.bassmzlswise.get_ab_maggies(outflux, self.basewave.copy(), mask_invalid=True)
        print(maggies)
        normmaggies = np.array(magfilt.get_ab_maggies(outflux, self.basewave.copy(), mask_invalid=True)[magfilter])
        print(normmaggies)
        assert(np.all(normmaggies > 0))

        synthnano = dict()
        for key in maggies.columns:
            synthnano[key] = 1E9 * maggies[key]

        if south:
            gflux, rflux, zflux, w1flux, w2flux = np.ma.getdata(synthnano['decam2014-g']), \
              np.ma.getdata(synthnano['decam2014-r']), np.ma.getdata(synthnano['decam2014-z']), \
              np.ma.getdata(synthnano['wise2010-W1']), np.ma.getdata(synthnano['wise2010-W2'])
        else:
            gflux, rflux, zflux, w1flux, w2flux = np.ma.getdata(synthnano['BASS-g']), \
              np.ma.getdata(synthnano['BASS-r']), np.ma.getdata(synthnano['MzLS-z']), \
              np.ma.getdata(synthnano['wise2010-W1']), np.ma.getdata(synthnano['wise2010-W2'])
        meta['TEMPLATEID'][ii] = template_id
        meta['REDSHIFT'][ii] = z
        meta['MAG'][ii] = -2.5 * np.log10(normmaggies)
        meta['FLUX_G'][ii] = gflux
        meta['FLUX_R'][ii] = rflux
        meta['FLUX_Z'][ii] = zflux
        meta['FLUX_W1'][ii] = w1flux
        meta['FLUX_W2'][ii] = w2flux

        return 1e17 * outflux, outwave, meta, objmeta
    