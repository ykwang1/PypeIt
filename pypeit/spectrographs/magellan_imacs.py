""" Implements DEIMOS-specific functions, including reading in slitmask design files.
"""

import glob
import re
import os
import numpy as np
import warnings

from scipy import interpolate

from astropy.io import fits
from astropy.time import Time

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import parse
from pypeit.core import framematch
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph
from pypeit.images import detector_container

from pypeit.utils import index_of_x_eq_y

from pypeit.spectrographs.slitmask import SlitMask
from pypeit.spectrographs.opticalmodel import ReflectionGrating, OpticalModel, DetectorMap
from IPython import embed

class MagellanIMACSSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Keck/DEIMOS specific code
    """
    ndet = 8

    def __init__(self):
        # Get it started
        super(MagellanIMACSSpectrograph, self).__init__()
        self.spectrograph = 'magellan_imacs'
        self.telescope = telescopes.MagellanTelescopePar()

        # Don't instantiate these until they're needed
        self.grating = None
        self.optical_model = None
        self.detector_map = None

    def get_detector_par(self, hdu, det):
        """
        Return a DectectorContainer for the current image

        Args:
            hdu (`astropy.io.fits.HDUList`):
                HDUList of the image of interest.
                Ought to be the raw file, or else..
            det (int):

        Returns:
            :class:`pypeit.images.detector_container.DetectorContainer`:

        """
        # Binning
        binning = self.get_meta_value(self.get_headarr(hdu), 'binning')  # Could this be detector dependent??

        # Mosaic3 + Slow -- http://www.lco.cl/telescopes-information/magellan/instruments/imacs/user-manual/the-imacs-user-manual
        # Detector 1
        detector_dict1 = dict(
            binning         = binning,
            det             = 1,
            dataext         = 1,
            specaxis        = 0,
            specflip        = False,
            spatflip        = False,
            platescale      = 0.1185,
            saturation      = 65535., # ADU
            nonlinear       = 0.95,
            mincounts       = -1e10,
            numamplifiers   = 1,
            gain            = np.atleast_1d(0.55),
            ronoise         = np.atleast_1d(2.8),
            )
        # Detector 2
        detector_dict2 = detector_dict1.copy()
        detector_dict2.update(dict(
            det=2,
            dataext=2,
            gain=np.atleast_1d(0.55),
        ))
        # Detector 3
        detector_dict3 = detector_dict1.copy()
        detector_dict3.update(dict(
            det=3,
            dataext=3,
            gain=np.atleast_1d(0.60),
        ))
        # Detector 4
        detector_dict4 = detector_dict1.copy()
        detector_dict4.update(dict(
            det=4,
            dataext=4,
            gain=np.atleast_1d(0.57),
            ronoise=np.atleast_1d(3.3),
        ))
        # Detector 5
        detector_dict5 = detector_dict1.copy()
        detector_dict5.update(dict(
            det=5,
            dataext=5,
            gain=np.atleast_1d(0.57),
            ronoise=np.atleast_1d(2.9),
        ))
        # Detector 6
        detector_dict6 = detector_dict1.copy()
        detector_dict6.update(dict(
            det=6,
            dataext=6,
            gain=np.atleast_1d(0.58),
            ronoise=np.atleast_1d(3.0),
        ))
        # Detector 7
        detector_dict7 = detector_dict1.copy()
        detector_dict7.update(dict(
            det=7,
            dataext=7,
            gain=np.atleast_1d(0.56),
            ronoise=np.atleast_1d(3.1),
        ))
        # Detector 8
        detector_dict8 = detector_dict1.copy()
        detector_dict8.update(dict(
            det=8,
            dataext=8,
            gain=np.atleast_1d(0.57),
            ronoise=np.atleast_1d(3.2),
        ))
        detectors = [detector_dict1, detector_dict2, detector_dict3, detector_dict4,
                     detector_dict5, detector_dict6, detector_dict7, detector_dict8]
        # Return
        return detector_container.DetectorContainer(**detectors[det-1])


    def default_pypeit_par(self):
        """
        Set default parameters for Keck DEIMOS reductions.
        """
        par = pypeitpar.PypeItPar()
        par['rdx']['spectrograph'] = 'keck_deimos'
        par['flexure']['spec_method'] = 'boxcar'
        # Set wave tilts order
        par['calibrations']['slitedges']['edge_thresh'] = 50.
        par['calibrations']['slitedges']['fit_order'] = 3
        par['calibrations']['slitedges']['minimum_slit_gap'] = 0.25
        par['calibrations']['slitedges']['minimum_slit_length'] = 4.
        par['calibrations']['slitedges']['sync_clip'] = False

        # 1D wavelength solution
        par['calibrations']['wavelengths']['lamps'] = ['ArI','HeI']
        #par['calibrations']['wavelengths']['nonlinear_counts'] \
        #        = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        par['calibrations']['wavelengths']['n_first'] = 3
        par['calibrations']['wavelengths']['match_toler'] = 2.5

        # Alter the method used to combine pixel flats
        par['calibrations']['pixelflatframe']['process']['combine'] = 'median'
        par['calibrations']['pixelflatframe']['process']['sig_lohi'] = [10.,10.]

        # Set the default exposure time ranges for the frame typing
        par['calibrations']['biasframe']['exprng'] = [None, 2]
        par['calibrations']['darkframe']['exprng'] = [999999, None]     # No dark frames
        par['calibrations']['pinholeframe']['exprng'] = [999999, None]  # No pinhole frames
        par['calibrations']['pixelflatframe']['exprng'] = [None, 29]
        par['calibrations']['traceframe']['exprng'] = [None, 29]
        par['scienceframe']['exprng'] = [29, None]
        
        # LACosmics parameters
        #par['scienceframe']['process']['sigclip'] = 4.0
        #par['scienceframe']['process']['objlim'] = 1.5

        return par

    def config_specific_par(self, scifile, inp_par=None):
        """
        Modify the PypeIt parameters to hard-wired values used for
        specific instrument configurations.

        .. todo::
            Document the changes made!
        
        Args:
            scifile (str):
                File to use when determining the configuration and how
                to adjust the input parameters.
            inp_par (:class:`pypeit.par.parset.ParSet`, optional):
                Parameter set used for the full run of PypeIt.  If None,
                use :func:`default_pypeit_par`.

        Returns:
            :class:`pypeit.par.parset.ParSet`: The PypeIt paramter set
            adjusted for configuration specific parameter values.
        """
        par = self.default_pypeit_par() if inp_par is None else inp_par

        headarr = self.get_headarr(scifile)

        '''
        # Longslit?
        if ('Long' in self.get_meta_value(headarr, 'decker')) or (
                'LVMslit' in self.get_meta_value(headarr, 'decker')):
            par['calibrations']['slitedges']['sync_predict'] = 'nearest'
        '''

        # Templates
        #if self.get_meta_value(headarr, 'dispname') == '600ZD':
        #    par['calibrations']['wavelengths']['method'] = 'full_template'
        #    par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_deimos_600.fits'
        #    par['calibrations']['wavelengths']['lamps'] += ['CdI', 'ZnI', 'HgI']
        #elif self.get_meta_value(headarr, 'dispname') == '830G':
        #    par['calibrations']['wavelengths']['method'] = 'full_template'
        #    par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_deimos_830G.fits'
        #elif self.get_meta_value(headarr, 'dispname') == '1200G':
        #    par['calibrations']['wavelengths']['method'] = 'full_template'
        #    par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_deimos_1200G.fits'

        # Return
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
        meta['ra'] = dict(ext=0, card='RA-D')
        meta['dec'] = dict(ext=0, card='DEC-D')
        meta['target'] = dict(ext=0, card='OBJECT')
        meta['decker'] = dict(ext=0, card='SLITMASK')
        meta['binning'] = dict(card=None, compound=True)

        meta['mjd'] = dict(ext=0, card=None, compound=True)
        meta['exptime'] = dict(ext=0, card='EXPTIME')
        meta['airmass'] = dict(ext=0, card='AIRMASS')
        meta['dispname'] = dict(ext=0, card='DISPERSR')
        # Extras for config and frametyping
        meta['dispangle'] = dict(ext=0, card='G-ANGLE', rtol=1e-3)
        # Image type
        meta['idname'] = dict(ext=0, card='EXPTYPE')
        # Lamps
        meta['lampstat01'] = dict(ext=0, card='LAMPS')

        # Ingest
        self.meta = meta

    def compound_meta(self, headarr, meta_key):
        """

        Args:
            headarr: list
            meta_key: str

        Returns:
            value

        """
        if meta_key == 'binning':
            binspatial, binspec = parse.parse_binning(headarr[0]['BINNING'])
            binning = parse.binning2string(binspec, binspatial)
            return binning
        elif meta_key == 'mjd':
            time = headarr[0]['UT-DATE']+'T'+headarr[0]['UT-TIME']
            ttime = Time(time, format='isot')
            return ttime.mjd
        else:
            msgs.error("Not ready for this compound meta")

    def configuration_keys(self):
        """
        Return the metadata keys that defines a unique instrument
        configuration.

        This list is used by :class:`pypeit.metadata.PypeItMetaData` to
        identify the unique configurations among the list of frames read
        for a given reduction.

        Returns:
            list: List of keywords of data pulled from meta
        """
        return ['dispname', 'decker', 'binning', 'dispangle']


    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.
        """
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)
        if ftype == 'science':
            #return good_exp & (fitstbl['lampstat01'] == 'Off') & (fitstbl['hatch'] == 'open')
            return good_exp & (fitstbl['idname'] == 'Object')
        if ftype == 'bias':
            return good_exp & (fitstbl['idname'] == 'Bias')
        if ftype in ['pixelflat', 'trace', 'illumflat']:
            # Flats and trace frames are typed together
            is_flat = [True if 'flat' in tname.lower() else False for tname in fitstbl['target']]
            return good_exp & np.array(is_flat)
        if ftype in ['pinhole', 'dark']:
            # Don't type pinhole or dark frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['arc', 'tilt']:
            is_arc = [True if 'lamp' in tname.lower() else False for tname in fitstbl['target']]
            return good_exp & (fitstbl['idname'] == 'Object') & is_arc

        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)


    def bpm(self, filename, det, shape=None, msbias=None):
        """
        Override parent bpm function with BPM specific to DEIMOS.

        .. todo::
            Allow for binning changes.

        Parameters
        ----------
        det : int, REQUIRED
        msbias : numpy.ndarray, required if the user wishes to generate a BPM based on a master bias
        **null_kwargs:
            Captured and never used

        Returns
        -------
        bpix : ndarray
          0 = ok; 1 = Mask

        """
        bpm_img = self.empty_bpm(filename, det, shape=shape)

        # Fill in bad pixels if a master bias frame is provided
        if msbias is not None:
            return self.bpm_frombias(msbias, det, bpm_img)

        '''
        if det == 1:
            bpm_img[:,1052:1054] = 1
        elif det == 2:
            bpm_img[:,0:4] = 1
            bpm_img[:,376:381] = 1
            bpm_img[:,489] = 1
            bpm_img[:,1333:1335] = 1
            bpm_img[:,2047] = 1
        elif det == 3:
            bpm_img[:,0:4] = 1
            bpm_img[:,221] = 1
            bpm_img[:,260] = 1
            bpm_img[:,366] = 1
            bpm_img[:,816:819] = 1
            bpm_img[:,851] = 1
            bpm_img[:,940] = 1
            bpm_img[:,1167] = 1
            bpm_img[:,1280] = 1
            bpm_img[:,1301:1303] = 1
            bpm_img[:,1744:1747] = 1
            bpm_img[:,-4:] = 1
        elif det == 4:
            bpm_img[:,0:4] = 1
            bpm_img[:,47] = 1
            bpm_img[:,744] = 1
            bpm_img[:,790:792] = 1
            bpm_img[:,997:999] = 1
        elif det == 5:
            bpm_img[:,25:27] = 1
            bpm_img[:,128:130] = 1
            bpm_img[:,1535:1539] = 1
        elif det == 7:
            bpm_img[:,426:428] = 1
            bpm_img[:,676] = 1
            bpm_img[:,1176:1178] = 1
        elif det == 8:
            bpm_img[:,440] = 1
            bpm_img[:,509:513] = 1
            bpm_img[:,806] = 1
            bpm_img[:,931:934] = 1
        '''

        return bpm_img

