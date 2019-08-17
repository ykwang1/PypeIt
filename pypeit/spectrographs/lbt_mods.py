""" Module for LBT/MODS specific codes
"""
import numpy as np
from astropy.io import fits

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph
from pypeit.core import parse

# ToDo: test MODS1B and MODS2B
# ToDo: write a reader to read different quadrants. currently assuming one giant detector with one amplifier

class LBTMODSSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Shane/Kast specific code
    """
    def __init__(self):
        # Get it started
        super(LBTMODSSpectrograph, self).__init__()
        self.spectrograph = 'lbt_mods_base'
        self.telescope = telescopes.LBTTelescopePar()
        self.timeunit = 'isot'

    @staticmethod
    def default_pypeit_par():
        """
        Set default parameters for Shane Kast reductions.
        """
        par = pypeitpar.PypeItPar()
        # Frame numbers
        par['calibrations']['standardframe']['number'] = 1
        par['calibrations']['biasframe']['number'] = 5
        par['calibrations']['pixelflatframe']['number'] = 5
        par['calibrations']['traceframe']['number'] = 5
        par['calibrations']['arcframe']['number'] = 1

        # Scienceimage default parameters
        par['scienceimage'] = pypeitpar.ScienceImagePar()
        # Always flux calibrate, starting with default parameters
        par['fluxcalib'] = pypeitpar.FluxCalibrationPar()
        # Always correct for flexure, starting with default parameters
        par['flexure'] = pypeitpar.FlexurePar()
        # Set the default exposure time ranges for the frame typing
        par['calibrations']['biasframe']['exprng'] = [None, 1]
        par['calibrations']['darkframe']['exprng'] = [999999, None]     # No dark frames
        par['calibrations']['pinholeframe']['exprng'] = [999999, None]  # No pinhole frames
        par['calibrations']['pixelflatframe']['exprng'] = [0, None]
        par['calibrations']['traceframe']['exprng'] = [0, None]
        par['calibrations']['arcframe']['exprng'] = [None, None]
        par['calibrations']['standardframe']['exprng'] = [1, 200]
        par['scienceframe']['exprng'] = [200, None]
        return par

    def init_meta(self):
        """
        Generate the meta data dict
        Note that the children can add to this

        Returns:
            self.meta: dict (generated in place)

        """
        self.meta = {}
        # Required (core)
        self.meta['ra'] = dict(ext=0, card='OBJRA')
        self.meta['dec'] = dict(ext=0, card='OBJDEC')
        self.meta['target'] = dict(ext=0, card='OBJECT')
        self.meta['decker'] = dict(ext=0, card='MASKNAME')
        self.meta['binning'] = dict(card=None, compound=True)
        self.meta['mjd'] = dict(ext=0, card='MJD-OBS')
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')
        self.meta['dispname'] = dict(ext=0, card='GRATNAME')
        self.meta['dichroic'] = dict(ext=0, card='FILTNAME')
        self.meta['idname'] = dict(ext=0, card='IMAGETYP')


    def compound_meta(self, headarr, meta_key):
        """

        Args:
            headarr: list
            meta_key: str

        Returns:
            value

        """
        if meta_key == 'binning':
            binspatial, binspec = parse.parse_binning(np.array([headarr[0]['CCDXBIN'], headarr[0]['CCDYBIN']]))
            binning = parse.binning2string(binspatial, binspec)
            return binning
        msgs.error("Not ready for this compound meta")

    def configuration_keys(self):
        """
        Set the configuration keys

        Returns:
            cfg_keys: list

        """
        # decker is not included because standards are usually taken with a 5" slit and arc using 0.8" slit
        return ['dispname', 'binning' ]

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.
        """
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)
        if ftype in ['science', 'standard']:
            return good_exp & (fitstbl['idname'] == 'OBJECT') & (fitstbl['ra'] != 'none') \
                   & (fitstbl['dispname'] != 'Flat')
        if ftype == 'bias':
            return good_exp  & (fitstbl['idname'] == 'BIAS')
        if ftype in ['pixelflat', 'trace']:
            # Flats and trace frames are typed together
            return good_exp  & (fitstbl['idname'] == 'FLAT') & (fitstbl['decker'] != 'Imaging')
        if ftype in ['pinhole', 'dark']:
            # Don't type pinhole or dark frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['arc', 'tilt']:
            return good_exp & (fitstbl['idname'] == 'COMP') & (fitstbl['dispname'] != 'Flat')

        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)


class LBTMODS1RSpectrograph(LBTMODSSpectrograph):
    """
    Child to handle LBT/MODS1R specific code
    """
    def __init__(self):
        # Get it started
        super(LBTMODS1RSpectrograph, self).__init__()
        self.spectrograph = 'lbt_mods1r'
        self.camera = 'MODS1R'
        self.detector = [
                # Detector 1
                pypeitpar.DetectorPar(
                            dataext         = 0,
                            specaxis        = 1,
                            specflip        = True,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.123,
                            darkcurr        = 0.4,
                            saturation      = 65535.,
                            nonlinear       = 0.99,
                            numamplifiers   = 1,
                            gain            = 2.5,
                            ronoise         = 4.2,
                            datasec='[:, 49:8240]',
                            oscansec='[:, 8240:]',
                            suffix          = '_mods1r'
                            )]
        self.numhead = 1


    def default_pypeit_par(self):
        """
        Set default parameters for LBT MODS1R reductions.
        """
        par = LBTMODSSpectrograph.default_pypeit_par()
        par['rdx']['spectrograph'] = 'lbt_mods1r'

        # 1D wavelength solution
        par['calibrations']['wavelengths']['sigdetect'] = 5.
        par['calibrations']['wavelengths']['rms_threshold'] = 1.0
        par['calibrations']['wavelengths']['fwhm'] = 10.
        #par['calibrations']['wavelengths']['lamps'] = ['XeI','ArII','ArI','NeI','KrI']]
        par['calibrations']['wavelengths']['lamps'] = ['ArI','NeI','KrI','XeI']
        #par['calibrations']['wavelengths']['lamps'] = ['OH_MODS']
        par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        par['calibrations']['wavelengths']['n_first'] = 3
        #par['calibrations']['wavelengths']['n_final'] = 4
        par['calibrations']['wavelengths']['match_toler'] = 2.5

        # slit
        par['calibrations']['slits']['sigdetect'] = 1000

        # Set wave tilts order
        par['calibrations']['tilts']['spat_order'] = 5
        par['calibrations']['tilts']['spec_order'] = 5
        par['calibrations']['tilts']['maxdev_tracefit'] = 0.02
        par['calibrations']['tilts']['maxdev2d'] = 0.02

        # reidentification stuff
        #par['calibrations']['wavelengths']['method'] = 'reidentify'
        #par['calibrations']['wavelengths']['reid_arxiv'] = 'shane_kast_blue_600_4310_d55.json'

        #par['calibrations']['biasframe']['useframe'] = 'bias'


        return par


    def bpm(self, filename, det, shape=None):
        """ Generate a BPM

        Parameters
        ----------
        det : int, REQUIRED
        **null_kwargs:
           Captured and never used

        Returns
        -------
        badpix : ndarray

        """

        # Get the empty bpm: force is always True
        bpm_img = self.empty_bpm(shape=shape, filename=filename, det=det)

        msgs.info("Using hard-coded BPM for  MODS1R")

        # TODO: Fix this
        # Get the binning
        hdu = fits.open(filename)
        header = hdu[0].header
        xbin, ybin = header['CCDXBIN'], header['CCDYBIN']
        hdu.close()

        # Apply the mask
        bpm_img[6278//xbin:6289//xbin, 1544//ybin:1634//ybin] = 1
        bpm_img[4202//xbin:4204//xbin, 1474//ybin:1544//ybin] = 1
        bpm_img[3551//xbin:3558//xbin, 2391//ybin:2903//ybin] = 1
        bpm_img[3553//xbin:3558//xbin, 1454//ybin:1544//ybin] = 1

        bpm_img[5650//xbin, 1280//ybin:1544//ybin] = 1
        bpm_img[4780//xbin, 1406//ybin:1536//ybin] = 1
        bpm_img[3554//xbin, 1544//ybin:2392//ybin] = 1
        bpm_img[163//xbin, 1544//ybin:1963//ybin] = 1

        return bpm_img


#    def check_headers(self, headers):
#        """
#        Check headers match expectations for a LBT MODS1R exposure.
#
#        See also
#        :func:`pypeit.spectrographs.spectrograph.Spectrograph.check_headers`.
#
#        Args:
#            headers (list):
#                A list of headers read from a fits file
#        """
#        expected_values = {   '0.NAXIS': 2,
#                            '0.INSTRUME': 'MODS1R' }
#        super(LBTMODS1RSpectrograph, self).check_headers(headers,
#                                                             expected_values=expected_values)

class LBTMODS1BSpectrograph(LBTMODSSpectrograph):
    """
    Child to handle LBT/MODS1R specific code
    """
    def __init__(self):
        # Get it started
        super(LBTMODS1BSpectrograph, self).__init__()
        self.spectrograph = 'lbt_mods1b'
        self.camera = 'MODS1B'
        self.detector = [
                # Detector 1
                pypeitpar.DetectorPar(
                            dataext         = 0,
                            specaxis        = 1,
                            specflip        = True,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.120,
                            darkcurr        = 0.5,
                            saturation      = 65535.,
                            nonlinear       = 0.99,
                            numamplifiers   = 1,
                            gain            = 2.1,
                            ronoise         = 3.0,
                            datasec='[:, 49:8240]',
                            oscansec='[:, 8240:]',
                            suffix          = '_mods1b'
                            )]
        self.numhead = 1


    def default_pypeit_par(self):
        """
        Set default parameters for LBT MODS1B reductions.
        """
        par = LBTMODSSpectrograph.default_pypeit_par()
        par['rdx']['spectrograph'] = 'lbt_mods1b'

        # 1D wavelength solution
        par['calibrations']['wavelengths']['sigdetect'] = 5.
        par['calibrations']['wavelengths']['rms_threshold'] = 0.20
        par['calibrations']['wavelengths']['lamps'] = ['XeI','ArII','ArI','NeI','KrI']
        par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        par['calibrations']['wavelengths']['n_first'] = 1

        # slit
        par['calibrations']['slits']['sigdetect'] = 1000

        # Set wave tilts order
        par['calibrations']['tilts']['spat_order'] = 3
        par['calibrations']['tilts']['spec_order'] = 5
        par['calibrations']['tilts']['maxdev_tracefit'] = 0.02
        par['calibrations']['tilts']['maxdev2d'] = 0.02

        # reidentification stuff
        #par['calibrations']['wavelengths']['method'] = 'reidentify'
        #par['calibrations']['wavelengths']['reid_arxiv'] = 'shane_kast_blue_600_4310_d55.json'

        par['calibrations']['biasframe']['useframe'] = 'bias'


        return par

    def bpm(self, filename, det, shape=None):
        """ Generate a BPM

        Parameters
        ----------
        det : int, REQUIRED
        **null_kwargs:
           Captured and never used

        Returns
        -------
        badpix : ndarray

        """
        # ToDo: check whether MODS1B has the same flip direction with MODS1R
        #       and modify the BPM accordingly
        # Get the empty bpm: force is always True
        bpm_img = self.empty_bpm(shape=shape, filename=filename, det=det)

        msgs.info("Using hard-coded BPM for  MODS1B")

        # Get the binning
        hdu = fits.open(filename)
        header = hdu[0].header
        xbin, ybin = header['CCDXBIN'], header['CCDYBIN']
        hdu.close()

        # Apply the mask
        bpm_img[6390//xbin:6392//xbin, 1262//ybin:1545//ybin] = 1

        bpm_img[3064//xbin, 1437//ybin:1937//ybin] = 1
        bpm_img[6490//xbin, 1161//ybin:1545//ybin] = 1
        bpm_img[7306//xbin, 783//ybin:1531//ybin] = 1

        return bpm_img



#    def check_headers(self, headers):
#        """
#        Check headers match expectations for a LBT MODS1B exposure.
#
#        See also
#        :func:`pypeit.spectrographs.spectrograph.Spectrograph.check_headers`.
#
#        Args:
#            headers (list):
#                A list of headers read from a fits file
#        """
#        expected_values = {   '0.NAXIS': 2,
#                            '0.INSTRUME': 'MODS1B' }
#        super(LBTMODS1BSpectrograph, self).check_headers(headers,
#                                                             expected_values=expected_values)


class LBTMODS2RSpectrograph(LBTMODSSpectrograph):
    """
    Child to handle LBT/MODS1R specific code
    """
    def __init__(self):
        # Get it started
        super(LBTMODS2RSpectrograph, self).__init__()
        self.spectrograph = 'lbt_mods2r'
        self.camera = 'MODS2R'
        self.detector = [
                # Detector 1
                pypeitpar.DetectorPar(
                            dataext         = 0,
                            specaxis        = 1,
                            specflip        = True,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.123,
                            darkcurr        = 0.4,
                            saturation      = 65535.,
                            nonlinear       = 0.99,
                            numamplifiers   = 1,
                            gain            = 1.7,
                            ronoise         = 2.8,
                            datasec='[:, 49:8240]',
                            oscansec='[:, 8240:]',
                            suffix          = '_mods2r'
                            )]
        self.numhead = 1


    def default_pypeit_par(self):
        """
        Set default parameters for LBT MODS2R reductions.
        """
        par = LBTMODSSpectrograph.default_pypeit_par()
        par['rdx']['spectrograph'] = 'lbt_mods2r'

        # 1D wavelength solution
        par['calibrations']['wavelengths']['sigdetect'] = 5.
        par['calibrations']['wavelengths']['rms_threshold'] = 1.0
        par['calibrations']['wavelengths']['fwhm'] = 10.
        #par['calibrations']['wavelengths']['lamps'] = ['XeI','ArII','ArI','NeI','KrI']]
        par['calibrations']['wavelengths']['lamps'] = ['ArI','NeI','KrI','XeI']
        #par['calibrations']['wavelengths']['lamps'] = ['OH_MODS']
        par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        par['calibrations']['wavelengths']['n_first'] = 3
        #par['calibrations']['wavelengths']['n_final'] = 4
        par['calibrations']['wavelengths']['match_toler'] = 2.5

        # slit
        par['calibrations']['slits']['sigdetect'] = 1000

        # Set wave tilts order
        par['calibrations']['tilts']['spat_order'] = 3
        par['calibrations']['tilts']['spec_order'] = 5
        par['calibrations']['tilts']['maxdev_tracefit'] = 0.02
        par['calibrations']['tilts']['maxdev2d'] = 0.02

        # reidentification stuff
        #par['calibrations']['wavelengths']['method'] = 'reidentify'
        #par['calibrations']['wavelengths']['reid_arxiv'] = 'shane_kast_blue_600_4310_d55.json'

        par['calibrations']['biasframe']['useframe'] = 'bias'


        return par

    def bpm(self, filename, det, shape=None):
        """ Generate a BPM

        Parameters
        ----------
        det : int, REQUIRED
        **null_kwargs:
           Captured and never used

        Returns
        -------
        badpix : ndarray

        """
        # Get the empty bpm: force is always True
        bpm_img = self.empty_bpm(shape=shape, filename=filename, det=det)

        msgs.info("Using hard-coded BPM for  MODS2R")

        # Get the binning
        hdu = fits.open(filename)
        header = hdu[0].header
        xbin, ybin = header['CCDXBIN'], header['CCDYBIN']
        hdu.close()

        # Apply the mask
        bpm_img[6148//xbin:6150//xbin, 1333//ybin:1544//ybin] = 1
        bpm_img[6207//xbin:6209//xbin, 1396//ybin:1544//ybin] = 1

        bpm_img[6101//xbin, 1342//ybin:1544//ybin] = 1
        bpm_img[6159//xbin, 1399//ybin:1544//ybin] = 1
        bpm_img[6189//xbin, 1316//ybin:1544//ybin] = 1
        bpm_img[7552//xbin, 1544//ybin:2771//ybin] = 1
        bpm_img[7504//xbin, 1544//ybin:2774//ybin] = 1
        bpm_img[4203//xbin, 0:1544//ybin] = 1
        bpm_img[4155//xbin, 0:1544//ybin] = 1

        return bpm_img



#    def check_headers(self, headers):
#        """
#        Check headers match expectations for a LBT MODS2R exposure.
#
#        See also
#        :func:`pypeit.spectrographs.spectrograph.Spectrograph.check_headers`.
#
#        Args:
#            headers (list):
#                A list of headers read from a fits file
#        """
#        expected_values = {   '0.NAXIS': 2,
#                            '0.INSTRUME': 'MODS2R' }
#        super(LBTMODS2RSpectrograph, self).check_headers(headers,
#                                                             expected_values=expected_values)

class LBTMODS2BSpectrograph(LBTMODSSpectrograph):
    """
    Child to handle LBT/MODS1R specific code
    """
    def __init__(self):
        # Get it started
        super(LBTMODS2BSpectrograph, self).__init__()
        self.spectrograph = 'lbt_mods2b'
        self.camera = 'MODS2B'
        self.detector = [
                # Detector 1
                pypeitpar.DetectorPar(
                            dataext         = 0,
                            specaxis        = 1,
                            specflip        = True,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.120,
                            darkcurr        = 0.5,
                            saturation      = 65535.,
                            nonlinear       = 0.99,
                            numamplifiers   = 1,
                            gain            = 2.0,
                            ronoise         = 3.7,
                            datasec='[:, 49:8240]',
                            oscansec='[:, 8240:]',
                            suffix          = '_mods2b'
                            )]
        self.numhead = 1


    def default_pypeit_par(self):
        """
        Set default parameters for LBT MODS2B reductions.
        """
        par = LBTMODSSpectrograph.default_pypeit_par()
        par['rdx']['spectrograph'] = 'lbt_mods2b'

        # 1D wavelength solution
        par['calibrations']['wavelengths']['sigdetect'] = 5.
        par['calibrations']['wavelengths']['rms_threshold'] = 0.20
        par['calibrations']['wavelengths']['lamps'] = ['XeI','ArII','ArI','NeI','KrI']
        par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        par['calibrations']['wavelengths']['n_first'] = 1

        # slit
        par['calibrations']['slits']['sigdetect'] = 1000

        # Set wave tilts order
        par['calibrations']['tilts']['spat_order'] = 3
        par['calibrations']['tilts']['spec_order'] = 5
        par['calibrations']['tilts']['maxdev_tracefit'] = 0.02
        par['calibrations']['tilts']['maxdev2d'] = 0.02

        # reidentification stuff
        #par['calibrations']['wavelengths']['method'] = 'reidentify'
        #par['calibrations']['wavelengths']['reid_arxiv'] = 'shane_kast_blue_600_4310_d55.json'

        par['calibrations']['biasframe']['useframe'] = 'bias'


        return par

    def bpm(self, filename, det, shape=None):
        """ Generate a BPM

        Parameters
        ----------
        det : int, REQUIRED
        **null_kwargs:
           Captured and never used

        Returns
        -------
        badpix : ndarray

        """
        # ToDo: check whether MODS2B has the same flip direction with MODS1R
        #       and modify the BPM accordingly

        # Get the empty bpm: force is always True
        bpm_img = self.empty_bpm(shape=shape, filename=filename, det=det)

        msgs.info("Using hard-coded BPM for  MODS2B")

        # Get the binning
        hdu = fits.open(filename)
        header = hdu[0].header
        xbin, ybin = header['CCDXBIN'], header['CCDYBIN']
        hdu.close()

        # Apply the mask
        bpm_img[5176//xbin:5179//xbin, 549//ybin:1544//ybin] = 1
        bpm_img[5176//xbin:5179//xbin, 1544//ybin:1628//ybin] = 1
        bpm_img[4408//xbin:4410//xbin, 1544//ybin:2661//ybin] = 1
        bpm_img[4408//xbin:4411//xbin, 2660//ybin:2663//ybin] = 1
        bpm_img[2495//xbin:2499//xbin, 1326//ybin:1544//ybin] = 1
        bpm_img[2391//xbin:2394//xbin, 1048//ybin:1051//ybin] = 1
        bpm_img[1974//xbin:1980//xbin, 806//ybin:1544//ybin] = 1
        bpm_img[1975//xbin:1980//xbin, 1544//ybin:1607//ybin] = 1
        bpm_img[1972//xbin:1974//xbin, 1587//ybin:1589//ybin] = 1
        bpm_img[274//xbin:278//xbin, 1341//ybin:1544//ybin] = 1
        bpm_img[275//xbin:278//xbin, 1251//ybin:1341//ybin] = 1
        bpm_img[276//xbin:278//xbin, 1242//ybin:1251//ybin] = 1
        bpm_img[274//xbin:277//xbin, 1544//ybin:3066//ybin] = 1

        bpm_img[2392//xbin, 1051//ybin:1544//ybin] = 1
        bpm_img[275//xbin, 1220//ybin:1242//ybin] = 1

        return bpm_img


#    def check_headers(self, headers):
#        """
#        Check headers match expectations for a LBT MODS2B exposure.
#
#        See also
#        :func:`pypeit.spectrographs.spectrograph.Spectrograph.check_headers`.
#
#        Args:
#            headers (list):
#                A list of headers read from a fits file
#        """
#        expected_values = {   '0.NAXIS': 2,
#                            '0.INSTRUME': 'MODS2B' }
#        super(LBTMODS2BSpectrograph, self).check_headers(headers,
#                                                             expected_values=expected_values)

