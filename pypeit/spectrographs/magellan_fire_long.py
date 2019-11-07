""" Module for Magellan/FIRE specific codes
"""
import numpy as np

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph

class MagellanFIRELONGSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Magellan/FIRE specific code

    .. note::
        For FIRE longslit, standard star and
        calibrations are usually use Fowler 1 read mode in which case
        the read noise is ~20 electron.

    """
    def __init__(self):
        # Get it started
        super(MagellanFIRELONGSpectrograph, self).__init__()
        self.spectrograph = 'magellan_fire_long'
        self.telescope = telescopes.MagellanTelescopePar()
        self.camera = 'FIRE'
        self.numhead = 1
        self.detector = [
                # Detector 1
                pypeitpar.DetectorPar(
                            dataext         = 0,
                            specaxis        = 0,
                            specflip        = False,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.15,
                            darkcurr        = 0.01,
                            saturation      = 320000., #32000 for low gain, I set to a higher value to keep data in K-band
                            nonlinear       = 0.875,
                            numamplifiers   = 1,
                            gain            = 3.8,
                            ronoise         = 20.0,
                            datasec         = '[5:2044, 900:1250]',
                            oscansec        = '[:5, 900:1250]'
                            )]

    def default_pypeit_par(self):
        """
        Set default parameters.
        """
        par = pypeitpar.PypeItPar()
        par['rdx']['spectrograph'] = 'magellan_fire_long'
        # Frame numbers
        par['calibrations']['standardframe']['number'] = 1
        par['calibrations']['biasframe']['number'] = 0
        par['calibrations']['pixelflatframe']['number'] = 5
        par['calibrations']['traceframe']['number'] = 5
        par['calibrations']['arcframe']['number'] = 1
        # No overscan
        for key in par['calibrations'].keys():
            if 'frame' in key:
                par['calibrations'][key]['process']['overscan'] = 'none'
        # Wavelengths
        # 1D wavelength solution with OH lines
        par['calibrations']['wavelengths']['rms_threshold'] = 1.0
        par['calibrations']['wavelengths']['sigdetect']=3
        par['calibrations']['wavelengths']['fwhm'] = 20
        par['calibrations']['wavelengths']['n_first']=2
        par['calibrations']['wavelengths']['n_final']=4
        par['calibrations']['wavelengths']['lamps'] = ['ArI', 'ArII', 'ThAr', 'NeI']
        par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        par['calibrations']['wavelengths']['method'] = 'full_template'
        par['calibrations']['wavelengths']['reid_arxiv'] = 'magellan_fire_long.fits'
        par['calibrations']['wavelengths']['match_toler']=5.0

        # Set slits and tilts parameters
#        par['calibrations']['tilts']['order'] = 2
        par['calibrations']['tilts']['tracethresh'] = 5
        par['calibrations']['slits']['trace_npoly'] = 5
        par['calibrations']['slits']['sigdetect'] = 10
        par['calibrations']['slits']['maxshift'] = 0.5
#        par['calibrations']['slits']['pcatype'] = 'pixel'
        # Scienceimage parameters
        par['scienceimage']['sig_thresh'] = 2
        par['scienceimage']['maxnumber'] = 2
        par['scienceimage']['find_trim_edge'] = [5,5]
        # Always flux calibrate, starting with default parameters
        par['fluxcalib'] = pypeitpar.FluxCalibrationPar()
        # Do not correct for flexure
        par['flexure'] = None
        # Set the default exposure time ranges for the frame typing
        par['calibrations']['standardframe']['exprng'] = [None, 60]
        par['calibrations']['arcframe']['exprng'] = [1, 50]
        par['calibrations']['darkframe']['exprng'] = [20, None]
        par['scienceframe']['exprng'] = [20, None]
        return par

    def init_meta(self):
        """
        Generate the meta data dict
        Note that the children can add to this

        Returns:
            self.meta: dict (generated in place)

        """
        meta = {}
        # Required (core)
        meta['ra'] = dict(ext=0, card='RA')
        meta['dec'] = dict(ext=0, card='DEC')
        meta['target'] = dict(ext=0, card='OBJECT')
        meta['decker'] = dict(ext=0, card=None, default='default')
        meta['dichroic'] = dict(ext=0, card=None, default='default')
        meta['binning'] = dict(ext=0, card=None, default='1,1')

        meta['mjd'] = dict(ext=0, card='ACQTIME')
        meta['exptime'] = dict(ext=0, card='EXPTIME')
        meta['airmass'] = dict(ext=0, card='AIRMASS')
        # Extras for config and frametyping
        meta['dispname'] = dict(ext=0, card='GRISM')
        meta['idname'] = dict(ext=0, card='OBSTYPE')

        # Ingest
        self.meta = meta

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.
        """
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)
        if ftype in ['pinhole', 'bias']:
            # No pinhole or bias frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['pixelflat', 'trace']:
            return good_exp & (fitstbl['idname'] == 'PixFlat')
        if ftype == 'standard':
            return good_exp & (fitstbl['idname'] == 'Telluric')
        if ftype == 'science':
            return good_exp & (fitstbl['idname'] == 'Science')
        if ftype in ['arc', 'tilt']:
            return good_exp & (fitstbl['idname'] == 'Arc')
        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)

    def bpm(self, filename, det, shape=None):
        """
        Override parent bpm function with BPM specific to X-Shooter VIS.

        .. todo::
            Allow for binning changes.

        Parameters
        ----------
        det : int, REQUIRED
        **null_kwargs:
            Captured and never used

        Returns
        -------
        bpix : ndarray
          0 = ok; 1 = Mask

        """
        msgs.info("Custom bad pixel mask for NIRES")
        bpm_img = self.empty_bpm(filename, det, shape=shape)

        return bpm_img
