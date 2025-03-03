"""
Main driver class for object finding, global skysubtraction and skymask construction

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""

import inspect
import numpy as np
import os

from astropy import stats
from abc import ABCMeta

from scipy import interpolate
from scipy.optimize import least_squares

from pypeit import specobjs
from pypeit import msgs, utils
from pypeit import masterframe, flatfield
from pypeit.display import display
from pypeit.core import skysub, pixels, qa, parse, flat
from pypeit.core import procimg
from pypeit.images import buildimage
from pypeit.core import findobj_skymask

from IPython import embed


class FindObjects:
    """
    This class will organize and run actions related to
    finding objects, global sky subtraction, and skymask construction for
    a Science or Standard star exposure

    Args:
        sciImg (:class:`~pypeit.images.pypeitimage.PypeItImage`):
            Image to reduce.
        spectrograph (:class:`~pypeit.spectrographs.spectrograph.Spectrograph`):
        par (:class:`~pypeit.par.pypeitpar.PypeItPar`):
        caliBrate (:class:`~pypeit.calibrations.Calibrations`):
        objtype (:obj:`str`):
           Specifies object being reduced 'science' 'standard' 'science_coadd2d'
        bkg_redux (:obj:`bool`, optional):
            If True, the sciImg has been subtracted by
            a background image (e.g. standard treatment in the IR)
        show (:obj:`bool`, optional):
           Show plots along the way?
        manual (:class:`~pypeit.manual_extract.ManualExtractObj`, optional):

    Attributes:
        ivarmodel (`numpy.ndarray`_):
            Model of inverse variance
        objimage (`numpy.ndarray`_):
            Model of object
        skyimage (`numpy.ndarray`_):
            Final model of sky
        initial_sky (`numpy.ndarray`_):
            Initial sky model after first pass with global_skysub()
        global_sky (`numpy.ndarray`_):
            Fit to global sky
        skymask (`numpy.ndarray`_):
            Mask of the sky fit
        outmask (`numpy.ndarray`_):
            Final output mask
        extractmask (`numpy.ndarray`_):
            Extraction mask
        slits (:class:`pypeit.slittrace.SlitTraceSet`):
        sobjs_obj (:class:`pypeit.specobjs.SpecObjs`):
            Objects found
        spat_flexure_shift (:obj:`float`):
        tilts (`numpy.ndarray`_):
            WaveTilts images generated on-the-spot
        waveimg (`numpy.ndarray`_):
            WaveImage image generated on-the-spot
        slitshift (`numpy.ndarray`_):
            Global spectral flexure correction for each slit (in pixels)
        vel_corr (:obj:`float`):
            Relativistic reference frame velocity correction (e.g. heliocentyric/barycentric/topocentric)

    """

    __metaclass__ = ABCMeta

    # Superclass factory method generates the subclass instance
    @classmethod
    def get_instance(cls, sciImg, spectrograph, par, caliBrate, objtype, bkg_redux=False,
                     find_negative=False, std_redux=False, show=False, basename=None, manual=None):
        """
        Instantiate the Reduce subclass appropriate for the provided
        spectrograph.

        The class must be subclassed from Reduce.  See :class:`Reduce` for
        the description of the valid keyword arguments.

        Args:
            sciImg (:class:`~pypeit.images.scienceimage.ScienceImage`):
                Image to reduce.
            spectrograph (:class:`~pypeit.spectrographs.spectrograph.Spectrograph`):
            par (:class:`~pypeit.par.pyepeitpar.PypeItPar`):
            caliBrate (:class:`~pypeit.calibrations.Calibrations`):
            objtype (:obj:`str`):
                Specifies object being reduced 'science' 'standard'
                'science_coadd2d'.  This is used only to determine the
                spat_flexure_shift and ech_order for coadd2d.
            bkg_redux (:obj:`bool`, optional):
                If True, the sciImg has been subtracted by
                a background image (e.g. standard treatment in the IR)
            find_negative (:obj:`bool`, optional):
                If True, the negative objects are found
            std_redux (:obj:`bool`, optional):
                If True the object being extracted is a standards star
                so that the reduction parameters can be adjusted accordingly.
            manual (:class:`~pypeit.manual_extract.ManualExtractObj`, optional):
                Class with info guiding the manual extraction
            basename (str, optional):
                Output filename used for spectral flexure QA
            show (:obj:`bool`, optional):
                Show plots along the way?
            **kwargs
                Passed to Parent init

        Returns:
            :class:`~pypeit.find_objects.FindObjects`:
        """
        return next(c for c in utils.all_subclasses(FindObjects)
                    if c.__name__ == (spectrograph.pypeline + 'FindObjects'))(
                            sciImg, spectrograph, par, caliBrate, objtype, bkg_redux=bkg_redux,
                            find_negative=find_negative, std_redux=std_redux, show=show,
                            basename=basename, manual=manual)

    def __init__(self, sciImg, spectrograph, par, caliBrate,
                 objtype, bkg_redux=False, find_negative=False, std_redux=False, show=False,
                 basename=None, manual=None):

        # Setup the parameters sets for this object. NOTE: This uses objtype, not frametype!

        # Instantiation attributes for this object
        self.sciImg = sciImg
        self.spectrograph = spectrograph
        self.objtype = objtype
        self.par = par
        self.caliBrate = caliBrate
        self.scaleimg = np.array([1.0], dtype=np.float)  # np.array([1]) applies no scale
        self.basename = basename
        self.manual = manual
        # Parse
        # Slit pieces
        #   WARNING -- It is best to unpack here then pass around self.slits
        #      Otherwise you have to keep in mind flexure, tweaking, etc.

        # TODO: The spatial flexure is not copied to the PypeItImage object if
        # the image (science or otherwise) is from a combination of multiple
        # frames.  Is that okay for this usage?
        # Flexure
        self.spat_flexure_shift = None
        if (objtype == 'science' and self.par['scienceframe']['process']['spat_flexure_correct']) or \
           (objtype == 'standard' and self.par['calibrations']['standardframe']['process']['spat_flexure_correct']):
            self.spat_flexure_shift = self.sciImg.spat_flexure
        elif objtype == 'science_coadd2d':
            self.spat_flexure_shift = None

        # Initialise the slits
        msgs.info("Initialising slits")
        self.initialise_slits()

        # Internal bpm mask
        # We want to keep the 'BOXSLIT', which has bpm=2. But we don't want to keep 'BOXSLIT'
        # with other bad flag (for which bpm>2)
        self.reduce_bpm = (self.slits.mask > 2) & (np.invert(self.slits.bitmask.flagged(
                        self.slits.mask, flag=self.slits.bitmask.exclude_for_reducing)))
        self.reduce_bpm_init = self.reduce_bpm.copy()

        # These may be None (i.e. COADD2D)
        self.waveTilts = caliBrate.wavetilts
        self.wv_calib = caliBrate.wv_calib

        # Load up other input items
        self.bkg_redux = bkg_redux
        self.find_negative = find_negative

        self.std_redux = std_redux
        self.det = caliBrate.det
        self.binning = caliBrate.binning
        self.pypeline = spectrograph.pypeline
        self.findobj_show = show

        self.steps = []

        # Key outputs images for extraction
        self.ivarmodel = None
        self.objimage = None
        self.skyimage = None
        self.initial_sky = None
        self.skymask = None
        self.outmask = None
        self.extractmask = None
        # SpecObjs object
        self.sobjs_obj = None
        self.slitshift = np.zeros(self.slits.nslits)  # Global spectral flexure slit shifts (in pixels) that are applied to all slits.
        self.vel_corr = None

        # Show?
        if self.findobj_show:
            self.show('image', image=sciImg.image, chname='processed', slits=True, clear=True)


    def create_skymask(self, sobjs_obj):
        r"""
        Creates a skymask from a SpecObjs object

        Args:
            sobjs_obj (:class:`pypeit.specobjs.SpecObjs`):
                Objects for which you would like to create the mask

        Returns:
            `numpy.ndarray`_: Boolean image with shape :math:`(N_{\rm spec},
            N_{\rm spat})` indicating which pixels are usable for global sky
            subtraction.  True = usable for sky subtraction, False = should be
            masked when sky subtracting.
        """
        # Masking options
        boxcar_rad_pix = None

        skymask = np.ones_like(self.sciImg.image, dtype=bool)
        gdslits = np.where(np.invert(self.reduce_bpm))[0]
        if sobjs_obj.nobj > 0:
            for slit_idx in gdslits:
                slit_spat = self.slits.spat_id[slit_idx]
                qa_title ="Generating skymask for slit # {:d}".format(slit_spat)
                msgs.info(qa_title)
                thismask = self.slitmask == slit_spat
                this_sobjs = sobjs_obj.SLITID == slit_spat
                # Boxcar mask?
                if self.par['reduce']['skysub']['mask_by_boxcar']:
                    boxcar_rad_pix = self.par['reduce']['extraction']['boxcar_radius'] / \
                                     self.get_platescale(slitord_id=self.slits.slitord_id[slit_idx])
                # Do it
                skymask[thismask] = findobj_skymask.create_skymask(sobjs_obj[this_sobjs], thismask,
                                                                   self.slits_left[:,slit_idx],
                                                                   self.slits_right[:,slit_idx],
                                                                   box_rad_pix=boxcar_rad_pix,
                                                                   trim_edg=self.par['reduce']['findobj']['find_trim_edge'])
        # Return
        return skymask

    def initialise_slits(self, initial=False):
        """
        Gather all the :class:`SlitTraceSet` attributes
        that we'll use here in :class:`FindObjects`

        Args:
            initial (:obj:`bool`, optional):
                Use the initial definition of the slits. If False,
                tweaked slits are used.
        """
        # Slits
        self.slits = self.caliBrate.slits
        # Select the edges to use
        self.slits_left, self.slits_right, _ \
            = self.slits.select_edges(initial=initial, flexure=self.spat_flexure_shift)

        # Slitmask
        self.slitmask = self.slits.slit_img(initial=initial, flexure=self.spat_flexure_shift,
                                            exclude_flag=self.slits.bitmask.exclude_for_reducing+['BOXSLIT'])
        # Now add the slitmask to the mask (i.e. post CR rejection in proc)
        # NOTE: this uses the par defined by EdgeTraceSet; this will
        # use the tweaked traces if they exist
        self.sciImg.update_mask_slitmask(self.slitmask)
#        # For echelle
#        self.spatial_coo = self.slits.spatial_coordinates(initial=initial, flexure=self.spat_flexure_shift)

    def run(self, std_trace=None, show_peaks=False):
        """
        Primary code flow for object finding in PypeIt reductions

        *NOT* used by COADD2D

        Parameters
        ----------
        std_trace : `numpy.ndarray`_, optional
            Trace of the standard star
        show_peaks : :obj:`bool`, optional
            Show peaks in find_objects methods

        Returns
        -------
        initial_sky : `numpy.ndarray`_
            Initial global sky model
        sobjs_obj : :class:`~pypeit.specobjs.SpecObjs`
            List of objects found
        """

        # Deal with dynamic calibrations
        # Tilts
        self.waveTilts.is_synced(self.slits)
        #   Deal with Flexure
        if self.par['calibrations']['tiltframe']['process']['spat_flexure_correct']:
            _spat_flexure = 0. if self.spat_flexure_shift is None else self.spat_flexure_shift
            # If they both shifted the same, there will be no reason to shift the tilts
            tilt_flexure_shift = _spat_flexure - self.waveTilts.spat_flexure
        else:
            tilt_flexure_shift = self.spat_flexure_shift
        msgs.info("Generating tilts image")
        self.tilts = self.waveTilts.fit2tiltimg(self.slitmask, flexure=tilt_flexure_shift)
        #

        # First pass object finding
        sobjs_obj, self.nobj = \
            self.find_objects(self.sciImg.image, std_trace=std_trace,
                              show_peaks=show_peaks,
                              show=self.findobj_show and not self.std_redux,
                              save_objfindQA=self.par['reduce']['findobj']['skip_second_find'] | self.std_redux)

        # create skymask using first pass sobjs_obj
        skymask_init = self.create_skymask(sobjs_obj)

        # Check if the user wants to overwrite the skymask with a pre-defined sky regions file.
        skymask_init, usersky = self.load_skyregions(skymask_init)

        # Global sky subtract (self.global_sky is also generated here)
        initial_sky = self.global_skysub(skymask=skymask_init).copy()

        # Second pass object finding on sky-subtracted image
        if (not self.std_redux) and (not self.par['reduce']['findobj']['skip_second_find']):
            sobjs_obj, self.nobj = self.find_objects(self.sciImg.image - initial_sky,
                                                     std_trace=std_trace, show=self.findobj_show,
                                                     show_peaks=show_peaks)
        else:
            msgs.info("Skipping 2nd run of finding objects")

        return initial_sky, sobjs_obj

    def find_objects(self, image, std_trace=None,
                     show_peaks=False, show_fits=False,
                     show_trace=False, show=False, save_objfindQA=True,
                     manual_extract_dict=None, debug=False):
        """
        Single pass at finding objects in the input image

        If self.find_negative is True, do a search for negative objects too

        Parameters
        ----------
        image : `numpy.ndarray`_
            Image to search for objects from. This floating-point image has shape
            (nspec, nspat) where the first dimension (nspec) is
            spectral, and second dimension (nspat) is spatial.
        std_trace : `numpy.ndarray`_, optional
            This is a one dimensional float array with shape = (nspec,) containing the standard star
            trace which is used as a crutch for tracing. If the no
            standard star is provided the code uses the the slit
            boundaries as the crutch.
        show_peaks : :obj:`bool`, optional
            Generate QA showing peaks identified by object finding
        show_fits : :obj:`bool`, optional
            Generate QA  showing fits to traces
        show_trace : :obj:`bool`, optional
            Generate QA  showing traces identified. Requires an open ginga RC
            modules window
        show : :obj:`bool`, optional
            Show all the QA
        save_objfindQA : :obj:`bool`, optional
            Save to disk (png file) QA showing the object profile
        manual_extract_dict : :obj:`dict`, optional
            This is only used by 2D coadd
        debug : :obj:`bool`, optional
            Show debugging plots?

        Returns
        -------
        sobjs_obj_single : :class:`~pypeit.specobjs.SpecObjs`
            Objects found
        nobj_single : :obj:`int`
            Number of objects found
        """
        # Positive image
        if manual_extract_dict is None:
            manual_extract_dict= self.manual.dict_for_objfind(self.det, neg=False) if self.manual is not None else None

        sobjs_obj_single, nobj_single = \
            self.find_objects_pypeline(image,
                                       std_trace=std_trace,
                                       show_peaks=show_peaks, show_fits=show_fits,
                                       show_trace=show_trace, save_objfindQA=save_objfindQA,
                                       manual_extract_dict=manual_extract_dict, 
                                       neg=False, debug=debug)

        # Find negative objects
        if self.find_negative:
            msgs.info("Finding objects in the negative image")
            # Parses
            manual_extract_dict = self.manual.dict_for_objfind(self.det, neg=True) if self.manual is not None else None
            sobjs_obj_single_neg, nobj_single_neg = \
                self.find_objects_pypeline(-image, std_trace=std_trace,
                                           show_peaks=show_peaks, show_fits=show_fits,
                                           show_trace=show_trace, save_objfindQA=save_objfindQA,
                                           manual_extract_dict=manual_extract_dict, neg=True,
                                           debug=debug)
            # Add (if there are any)
            sobjs_obj_single.append_neg(sobjs_obj_single_neg)

        if show:
            gpm = self.sciImg.select_flag(invert=True)
            self.show('image', image=image*gpm.astype(float), chname='objfind',
                      sobjs=sobjs_obj_single, slits=True)

        # For nobj we take only the positive objects
        return sobjs_obj_single, nobj_single

    # TODO maybe we don't need parent and children for this method. But IFU has a bunch of extra methods.
    def find_objects_pypeline(self, image, std_trace=None,
                              show_peaks=False, show_fits=False, show_trace=False,
                              show=False, save_objfindQA=False, neg=False, debug=False,
                              manual_extract_dict=None):

        """
         Dummy method for object finding. Overloaded by class specific object finding.

         Returns:

         """
        return None, None

    def get_platescale(self, slitord_id=None):
        """
        Return the platescale in binned pixels for the current detector/echelle order

        Over-loaded by the children

        Args:
            slitord_id (:obj:`int`, optional):
                slit spat_id (MultiSlit, IFU) or ech_order (Echelle) value

        Returns:
            :obj:`float`: plate scale in binned pixels

        """
        pass

    def global_skysub(self, skymask=None, update_crmask=True, trim_edg=(3,3),
                      previous_sky=None, show_fit=False, show=False, show_objs=False):
        """
        Perform global sky subtraction, slit by slit

        Wrapper to skysub.global_skysub

        Args:
            skymask (`numpy.ndarray`_, None):
                A 2D image indicating sky regions (1=sky)
            update_crmask (bool, optional):
            show_fit (bool, optional):
            show (bool, optional):
            show_objs (bool, optional):
            previous_sky (`numpy.ndarray`_, optional):
                Sky model estimate from a previous run of global_sky
                Used to generate an improved estimated of the variance

        Returns:
            `numpy.ndarray`_: image of the the global sky model

        """
        # Prep
        global_sky = np.zeros_like(self.sciImg.image)
        # Parameters for a standard star
        if self.std_redux:
            sigrej = 7.0
            update_crmask = False
            if not self.par['reduce']['skysub']['global_sky_std']:
                msgs.info('Skipping global sky-subtraction for standard star.')
                return global_sky
        else:
            sigrej = 3.0

        # We use this tmp bpm so that we exclude the BOXSLITS during the global_skysub
        tmp_bpm = (self.slits.mask > 0) & \
                          (np.invert(self.slits.bitmask.flagged(self.slits.mask,
                                                                flag=self.slits.bitmask.exclude_for_reducing)))
        gdslits = np.where(np.invert(tmp_bpm))[0]

        # Mask objects using the skymask? If skymask has been set by objfinding, and masking is requested, then do so
        skymask_now = skymask if (skymask is not None) else np.ones_like(self.sciImg.image, dtype=bool)

        # Allow for previous sky to better estimate ivar
        #  Unless we used a background image (i.e. bkg_redux=True)
        if (previous_sky is not None) and (not self.bkg_redux):
            # Estimate the variance using the input sky model
            var = procimg.variance_model(self.sciImg.base_var,
                                          counts=previous_sky,
                                          count_scale=self.sciImg.img_scale,
                                          noise_floor=self.sciImg.noise_floor)
            skysub_ivar = utils.inverse(var)
        else:
            skysub_ivar = self.sciImg.ivar


        # Loop on slits
        for slit_idx in gdslits:
            slit_spat = self.slits.spat_id[slit_idx]
            msgs.info("Global sky subtraction for slit: {:d}".format(slit_spat))
            thismask = self.slitmask == slit_spat
            inmask = self.sciImg.select_flag(invert=True) & thismask & skymask_now
            # All masked?
            if not np.any(inmask):
                msgs.warn("No pixels for fitting sky.  If you are using mask_by_boxcar=True, your radius may be too large.")
                self.reduce_bpm[slit_idx] = True
                continue

            # Find sky
            global_sky[thismask] = skysub.global_skysub(
                self.sciImg.image, skysub_ivar, self.tilts, thismask,
                self.slits_left[:,slit_idx], self.slits_right[:,slit_idx],
                inmask=inmask, sigrej=sigrej,
                bsp=self.par['reduce']['skysub']['bspline_spacing'],
                no_poly=self.par['reduce']['skysub']['no_poly'],
                pos_mask=(not self.bkg_redux), show_fit=show_fit)

            # Mask if something went wrong
            if np.sum(global_sky[thismask]) == 0.:
                msgs.warn("Bad fit to sky.  Rejecting slit: {:d}".format(slit_idx))
                self.reduce_bpm[slit_idx] = True

        if update_crmask and self.par['scienceframe']['process']['mask_cr']:
            # Find CRs with sky subtraction
            # TODO: Shouldn't the saturation flagging account for the
            # subtraction of the sky?
            self.sciImg.build_crmask(self.par['scienceframe']['process'],
                                     subtract_img=global_sky)
            # Update the fullmask
            self.sciImg.update_mask_cr(self.sciImg.crmask)

        # Step
        self.steps.append(inspect.stack()[0][3])

        if show:
            sobjs_show = None if show_objs else self.sobjs_obj
            # Global skysub is the first step in a new extraction so clear the channels here
            self.show('global', global_sky=global_sky, slits=True, sobjs=sobjs_show, clear=False)


        # Return
        return global_sky

    def load_skyregions(self, skymask_init):
        """
        Load or generate the sky regions

        Parameters
        ----------
        skymask_init :  `numpy.ndarray`_
            A boolean array of sky pixels (True is pixel is a sky region)

        Returns
        -------
        skymask_init :  `numpy.ndarray`_
            A boolean array of sky pixels (True is pixel is a sky region)
        usersky : bool
            If the user has defined the sky, set this variable to True (otherwise False).
        """
        usersky = False
        if self.par['reduce']['skysub']['load_mask']:
            # Check if a master Sky Regions file exists for this science frame
            file_base = os.path.basename(self.sciImg.files[0])
            prefix = os.path.splitext(file_base)
            if prefix[1] == ".gz":
                sciName = os.path.splitext(prefix[0])[0]
            else:
                sciName = prefix[0]

            # Setup the master frame name
            master_dir = self.caliBrate.master_dir
            master_key = self.caliBrate.fitstbl.master_key(0, det=self.det) + "_" + sciName

            regfile = masterframe.construct_file_name(buildimage.SkyRegions,
                                                      master_key=master_key,
                                                      master_dir=master_dir)
            # Check if a file exists
            if os.path.exists(regfile):
                msgs.info("Loading SkyRegions file for: {0:s} --".format(sciName) + msgs.newline() + regfile)
                skyreg = buildimage.SkyRegions.from_file(regfile)
                skymask_init = skyreg.image.astype(np.bool)
                usersky = True
            else:
                msgs.warn("SkyRegions file not found:" + msgs.newline() + regfile)
        elif self.par['reduce']['skysub']['user_regions'] is not None and \
                len(self.par['reduce']['skysub']['user_regions']) != 0:
            skyregtxt = self.par['reduce']['skysub']['user_regions']
            if type(skyregtxt) is list:
                skyregtxt = ",".join(skyregtxt)
            msgs.info("Generating skysub mask based on the user defined regions   {0:s}".format(skyregtxt))
            # The resolution probably doesn't need to be a user parameter
            maxslitlength = np.max(self.slits_right-self.slits_left)
            # Get the regions
            status, regions = skysub.read_userregions(skyregtxt, self.slits.nslits, maxslitlength)
            # Generate image
            skymask_init = skysub.generate_mask(self.pypeline, regions, self.slits, self.slits_left,
                                                self.slits_right, spat_flexure=self.spat_flexure_shift)
            usersky = True
        return skymask_init, usersky

    def show(self, attr, image=None, global_sky=None, showmask=False, sobjs=None,
             chname=None, slits=False,clear=False):
        """
        Show one of the internal images

        .. todo::
            Should probably put some of these in ProcessImages

        Parameters
        ----------
        attr : str
          global -- Sky model (global)
          sci -- Processed science image
          rawvar -- Raw variance image
          modelvar -- Model variance image
          crmasked -- Science image with CRs set to 0
          skysub -- Science image with global sky subtracted
          image -- Input image
        display : str, optional
        image : ndarray, optional
          User supplied image to display

        Returns
        -------

        """

        if showmask:
            mask_in = self.sciImg.fullmask
            bitmask_in = self.sciImg.bitmask
        else:
            mask_in = None
            bitmask_in = None

        img_gpm = self.sciImg.select_flag(invert=True)
        detname = self.spectrograph.get_det_name(self.det)

        if attr == 'global' and all([a is not None for a in [self.sciImg.image, global_sky, self.sciImg.fullmask]]):
            # global sky subtraction
            # sky subtracted image
            image = (self.sciImg.image - global_sky) * img_gpm.astype(float)
            mean, med, sigma = stats.sigma_clipped_stats(image[img_gpm], sigma_lower=5.0,
                                                         sigma_upper=5.0)
            cut_min = mean - 1.0 * sigma
            cut_max = mean + 4.0 * sigma
            ch_name = chname if chname is not None else f'global_sky_{detname}'
            viewer, ch = display.show_image(image, chname=ch_name, bitmask=bitmask_in,
                                            mask=mask_in, clear=clear, wcs_match=True)
                                          #, cuts=(cut_min, cut_max))
        elif attr == 'image':
            ch_name = chname if chname is not None else 'image'
            viewer, ch = display.show_image(image, chname=ch_name, clear=clear, wcs_match=True)
        else:
            msgs.warn("Not an option for show")

        if sobjs is not None:
            for spec in sobjs:
                color = 'magenta' if spec.hand_extract_flag else 'orange'
                display.show_trace(viewer, ch, spec.TRACE_SPAT, spec.NAME, color=color)

        if slits and self.slits_left is not None:
            display.show_slits(viewer, ch, self.slits_left, self.slits_right)

    def __repr__(self):
        txt = '<{:s}: nimg={:d}'.format(self.__class__.__name__,
                                        self.nsci)
        if len(self.steps) > 0:
            txt+= ' steps: ['
            for step in self.steps:
                txt += '{:s}, '.format(step)
            txt = txt[:-2]+']'  # Trim the trailing comma
        txt += '>'
        return txt


class MultiSlitFindObjects(FindObjects):
    """
    Child of Reduce for Multislit and Longslit reductions

    See parent doc string for Args and Attributes

    """
    def __init__(self, sciImg, spectrograph, par, caliBrate, objtype, **kwargs):
        super().__init__(sciImg, spectrograph, par, caliBrate, objtype, **kwargs)

    def get_platescale(self, slitord_id=None):
        """
        Return the platescale in binned pixels for the current detector/echelle order

        Args:
            slitord_id (:obj:`int`, optional):
                slit spat_id (MultiSlit, IFU) or ech_order (Echelle) value

        Returns:
            :obj:`float`: plate scale in binned pixels

        """
        bin_spec, bin_spat = parse.parse_binning(self.binning)
        return self.sciImg.detector.platescale * bin_spec

    def find_objects_pypeline(self, image, std_trace=None,
                              manual_extract_dict=None,
                              show_peaks=False, show_fits=False, show_trace=False,
                              show=False, save_objfindQA=False, neg=False, debug=False):
        """
        Pipeline specific find objects routine

        Parameters
        ----------

        image : `numpy.ndarray`_
            Image to search for objects from. This floating-point image has shape
            (nspec, nspat) where the first dimension (nspec) is
            spectral, and second dimension (nspat) is spatial.
        std_trace : `numpy.ndarray`_, optional
            This is a one dimensional float array with shape = (nspec,) containing the standard star
            trace which is used as a crutch for tracing. If the no
            standard star is provided the code uses the the slit
            boundaries as the crutch.
        manual_extract_dict : :obj:`dict`, optional
            Dict guiding the manual extraction
        show_peaks : :obj:`bool`, optional
            Generate QA showing peaks identified by object finding
        show_fits : :obj:`bool`, optional
            Generate QA  showing fits to traces
        show_trace : :obj:`bool`, optional
            Generate QA  showing traces identified. Requires an open ginga RC
            modules window
        show : :obj:`bool`, optional
            Show all the QA
        save_objfindQA : :obj:`bool`, optional
            Save to disk (png file) QA showing the object profile
        neg : :obj:`bool`, optional
            Is this a negative image?
        debug : :obj:`bool`, optional
            Show debugging plots?

        Returns
        -------
        specobjs : :class:`~pypeot.specobjs.Specobjs`
            Container holding Specobj objects
        nobj : :obj:`int`
            Number of objects identified
        """
        gdslits = np.where(np.invert(self.reduce_bpm))[0]

        # Instantiate the specobjs container
        sobjs = specobjs.SpecObjs()

        # Loop on slits
        for slit_idx in gdslits:
            slit_spat = self.slits.spat_id[slit_idx]
            qa_title ="Finding objects on slit # {:d}".format(slit_spat)
            msgs.info(qa_title)
            thismask = self.slitmask == slit_spat
            inmask = self.sciImg.select_flag(invert=True) & thismask
            specobj_dict = {'SLITID': slit_spat,
                            'DET': self.sciImg.detector.name,
                            'OBJTYPE': self.objtype,
                            'PYPELINE': self.pypeline}

            # This condition allows to not use a threshold to find objects in alignment boxes
            # because these boxes are smaller than normal slits and the stars are very bright,
            # the detection threshold would be too high and the star not detected.
            if self.slits.bitmask.flagged(self.slits.mask[slit_idx], flag='BOXSLIT'):
                sig_thresh = 0.
            else:
                sig_thresh = self.par['reduce']['findobj']['sig_thresh']

            # Set objfind QA filename
            objfindQA_filename = None
            if save_objfindQA and (self.basename is not None):
                out_dir = os.path.join(self.par['rdx']['redux_path'], self.par['rdx']['qadir'])
                if self.find_negative:
                    basename = 'neg_' + self.basename if neg else 'pos_' + self.basename
                else:
                    basename = self.basename
                detname = self.spectrograph.get_det_name(self.det)
                objfindQA_filename = qa.set_qa_filename(basename, 'obj_profile_qa', slit=slit_spat,
                                                        det=detname, out_dir=out_dir)

            sobjs_slit = \
                    findobj_skymask.objs_in_slit(image, thismask,
                                self.slits_left[:,slit_idx],
                                self.slits_right[:,slit_idx],
                                inmask=inmask, has_negative=self.find_negative,
                                ncoeff=self.par['reduce']['findobj']['trace_npoly'],
                                std_trace=std_trace,
                                sig_thresh=sig_thresh,
                                cont_sig_thresh=self.par['reduce']['findobj']['cont_sig_thresh'],
                                hand_extract_dict=manual_extract_dict,
                                specobj_dict=specobj_dict, show_peaks=show_peaks,
                                show_fits=show_fits, show_trace=show_trace,
                                trim_edg=self.par['reduce']['findobj']['find_trim_edge'],
                                cont_fit=self.par['reduce']['findobj']['find_cont_fit'],
                                npoly_cont=self.par['reduce']['findobj']['find_npoly_cont'],
                                fwhm=self.par['reduce']['findobj']['find_fwhm'],
                                use_user_fwhm=self.par['reduce']['extraction']['use_user_fwhm'],
                                boxcar_rad=self.par['reduce']['extraction']['boxcar_radius'] / self.get_platescale(),  #pixels
                                maxdev=self.par['reduce']['findobj']['find_maxdev'],
                                find_min_max=self.par['reduce']['findobj']['find_min_max'],
                                qa_title=qa_title, nperslit=self.par['reduce']['findobj']['maxnumber'],
                                objfindQA_filename=objfindQA_filename,
                                debug_all=debug)
            # Record
            sobjs.add_sobj(sobjs_slit)

        # Steps
        self.steps.append(inspect.stack()[0][3])
        if show:
            gpm = self.sciImg.select_flag(invert=True)
            self.show('image', image=image*gpm.astype(float), chname='objfind', sobjs=sobjs,
                      slits=True)

        # Return
        return sobjs, len(sobjs)


class EchelleFindObjects(FindObjects):
    """
    Child of Reduce for Echelle reductions

    See parent doc string for Args and Attributes

    """
    def __init__(self, sciImg, spectrograph, par, caliBrate, objtype, **kwargs):
        super().__init__(sciImg, spectrograph, par, caliBrate, objtype, **kwargs)

        # JFH For 2d coadds the orders are no longer located at the standard locations
        self.order_vec = spectrograph.orders if 'coadd2d' in self.objtype \
                            else self.slits.ech_order
#                            else self.spectrograph.order_vec(self.spatial_coo)
        if self.order_vec is None:
            msgs.error('Unable to set Echelle orders, likely because they were incorrectly '
                       'assigned in the relevant SlitTraceSet.')

    def get_platescale(self, slitord_id=None):
        """
        Return the platescale in binned pixels for the current detector/echelle order

        Args:
            slitord_id (:obj:`int`, optional):
                slit spat_id (MultiSlit, IFU) or ech_order (Echelle) value

        Returns:
            :obj:`float`: plate scale in binned pixels

        """
        if slitord_id is None:
            msgs.error('slitord_id is missing. Plate scale for current echelle order cannot be determined.')
        return self.spectrograph.order_platescale(slitord_id, binning=self.binning)[0]


# TODO This does not appear to be used anywhere
#    def get_positive_sobj(self, specobjs, iord):
#        """
#        Return the current object from self.sobjs_obj
#
#        Args:
#            iord (int):
#                Echelle order index
#
#        Returns:
#            :class:`pypeit.specobj.SpecObj`:
#
#        """
#        # pos indices of objects for this slit
#        thisobj = (self.sobjs_obj.ech_orderindx == iord) & (self.sobjs_obj.ech_objid > 0)
#        return self.sobjs_obj[np.where(thisobj)[0][0]]

    def find_objects_pypeline(self, image, std_trace=None,
                              show=False, show_peaks=False, show_fits=False,
                              show_trace=False, save_objfindQA=False, neg=False, debug=False,
                              manual_extract_dict=None):
        """
        Pipeline specific find objects routine

        Parameters
        ----------
        image : `numpy.ndarray`_
            Image to search for objects from. This floating-point image has shape
            (nspec, nspat) where the first dimension (nspec) is
            spectral, and second dimension (nspat) is spatial.
        std_trace : `numpy.ndarray`_, optional
            This is a one dimensional float array with shape = (nspec,) containing the standard star
            trace which is used as a crutch for tracing. If the no
            standard star is provided the code uses the the slit
            boundaries as the crutch.
        manual_extract_dict : :obj:`dict`, optional
            Dict guiding the manual extraction
        show_peaks : :obj:`bool`, optional
            Generate QA showing peaks identified by object finding
        show_fits : :obj:`bool`, optional
            Generate QA  showing fits to traces
        show_trace : :obj:`bool`, optional
            Generate QA  showing traces identified. Requires an open ginga RC modules window
        save_objfindQA : :obj:`bool`, optional
            Save to disk (png file) QA showing the object profile
        neg : :obj:`bool`, optional
            Is this a negative image?
        show : :obj:`bool`, optional
        debug : :obj:`bool`, optional

        Returns
        -------
        specobjs : :class:`~pypeit.specobjs.Specobjs`
            Container holding Specobj objects
        nobj : :obj:`int`
            Number of objects identified
        """

        plate_scale = self.spectrograph.order_platescale(self.order_vec, binning=self.binning)
        inmask = self.sciImg.select_flag(invert=True)
        # Find objects
        # TODO -- Eliminate this specobj_dict thing
        # TODO: Not sure how this fairs if self.det is a tuple...
        specobj_dict = {'SLITID': 999, 'DET': self.sciImg.detector.name, 'OBJTYPE': self.objtype,
                        'PYPELINE': self.pypeline}

        # Set objfind QA filename
        objfindQA_filename = None
        if save_objfindQA and (self.basename is not None):
            out_dir = os.path.join(self.par['rdx']['redux_path'], self.par['rdx']['qadir'])
            if self.find_negative:
                basename = 'neg_' + self.basename if neg else 'pos_' + self.basename
            else:
                basename = self.basename
            detname = self.spectrograph.get_det_name(self.det)
            objfindQA_filename = qa.set_qa_filename(basename, 'obj_profile_qa', slit=999,
                                                    det=detname, out_dir=out_dir)

        sobjs_ech = findobj_skymask.ech_objfind(
            image, self.sciImg.ivar, self.slitmask, self.slits_left, self.slits_right,
            self.order_vec, self.reduce_bpm, det=self.det,
            spec_min_max=np.vstack((self.slits.specmin, self.slits.specmax)),
            inmask=inmask, has_negative=self.find_negative, ncoeff=self.par['reduce']['findobj']['trace_npoly'],
            hand_extract_dict=manual_extract_dict, plate_scale=plate_scale,
            std_trace=std_trace,
            specobj_dict=specobj_dict,
            sig_thresh=self.par['reduce']['findobj']['sig_thresh'],
            cont_sig_thresh=self.par['reduce']['findobj']['cont_sig_thresh'],
            show_peaks=show_peaks, show_fits=show_fits,
            trim_edg=self.par['reduce']['findobj']['find_trim_edge'],
            cont_fit=self.par['reduce']['findobj']['find_cont_fit'],
            npoly_cont=self.par['reduce']['findobj']['find_npoly_cont'],
            fwhm=self.par['reduce']['findobj']['find_fwhm'],
            use_user_fwhm=self.par['reduce']['extraction']['use_user_fwhm'],
            nperorder=self.par['reduce']['findobj']['maxnumber'],
            maxdev=self.par['reduce']['findobj']['find_maxdev'],
            max_snr=self.par['reduce']['findobj']['ech_find_max_snr'],
            min_snr=self.par['reduce']['findobj']['ech_find_min_snr'],
            nabove_min_snr=self.par['reduce']['findobj']['ech_find_nabove_min_snr'],
            box_radius=self.par['reduce']['extraction']['boxcar_radius'],  # arcsec
            show_trace=show_trace, objfindQA_filename=objfindQA_filename)

        # Steps
        self.steps.append(inspect.stack()[0][3])
        if show:
            gpm = self.sciImg.select_flag(invert=True)
            self.show('image', image=image*gpm.astype(float), chname='ech_objfind',
                      sobjs=sobjs_ech, slits=False)

        return sobjs_ech, len(sobjs_ech)


class IFUFindObjects(MultiSlitFindObjects):
    """
    Child of Reduce for IFU reductions

    See parent doc string for Args and Attributes

    """
    def __init__(self, sciImg, spectrograph, par, caliBrate, objtype, **kwargs):
        super().__init__(sciImg, spectrograph, par, caliBrate, objtype, **kwargs)
        self.initialise_slits(initial=True)

    def find_objects_pypeline(self, image, std_trace=None,
                              show_peaks=False, show_fits=False, show_trace=False,
                              show=False, save_objfindQA=False, neg=False, debug=False,
                              manual_extract_dict=None):
        """
        See MultiSlitReduce for slit-based IFU reductions
        """
        if self.par['reduce']['cube']['slit_spec']:
            return super().find_objects_pypeline(image, std_trace=std_trace,
                                                 show_peaks=show_peaks, show_fits=show_fits, show_trace=show_trace,
                                                 show=show, save_objfindQA=save_objfindQA, neg=neg,
                                                 debug=debug, manual_extract_dict=manual_extract_dict)
        return None, None, None

    def apply_relative_scale(self, scaleImg):
        """Apply a relative scale to the science frame (and correct the varframe, too)

         Args:
             scaleImg (`numpy.ndarray`_):
                scale image to divide the science frame by
        """
        # Check that scaleimg is set to the correct shape
        if self.scaleimg.size == 1:
            self.scaleimg = np.ones_like(self.sciImg.image)
        # Correct the relative illumination of the science frame
        msgs.info("Correcting science frame for relative illumination")
        self.scaleimg *= scaleImg.copy()
        self.sciImg.image, _bpm, varImg = flat.flatfield(self.sciImg.image, scaleImg,
                                                         varframe=utils.inverse(self.sciImg.ivar))
        if np.any(_bpm):
            self.sciImg.update_mask('BADSCALE', indx=_bpm)
        self.sciImg.ivar = utils.inverse(varImg)

    def illum_profile_spatial(self, skymask=None, trim_edg=(0, 0), debug=False):
        """
        Calculate the residual spatial illumination profile using the sky regions.

        The redisual is calculated using the differential:

        .. code-block:: python

            correction = amplitude * (1 + spatial_shift * (dy/dx)/y)

        where ``y`` is the spatial profile determined from illumflat, and
        spatial_shift is the residual spatial flexure shift in units of pixels.

         Args:
            skymask (`numpy.ndarray`_):
                Mask of sky regions where the spatial illumination will be determined
            trim_edg (:obj:`tuple`):
                A tuple of two ints indicated how much of the slit edges should be
                trimmed when fitting to the spatial profile.
            debug (:obj:`bool`):
                Show debugging plots?
        """

        msgs.info("Performing spatial sensitivity correction")
        # Setup some helpful parameters
        skymask_now = skymask if (skymask is not None) else np.ones_like(self.sciImg.image, dtype=bool)
        hist_trim = 0  # Trim the edges of the histogram to take into account edge effects
        gpm = self.sciImg.select_flag(invert=True)
        slitid_img_init = self.slits.slit_img(pad=0, initial=True, flexure=self.spat_flexure_shift)
        spatScaleImg = np.ones_like(self.sciImg.image)
        # For each slit, grab the spatial coordinates and a spline
        # representation of the spatial profile from the illumflat
        rawimg = self.sciImg.image.copy()
        numbins = int(np.max(self.slits.get_slitlengths(initial=True, median=True)))
        spatbins = np.linspace(0.0, 1.0, numbins + 1)
        spat_slit = 0.5 * (spatbins[1:] + spatbins[:-1])
        slitlength = np.median(self.slits.get_slitlengths(median=True))
        coeff_fit = np.zeros((self.slits.nslits, 2))
        for sl, slitnum in enumerate(self.slits.spat_id):
            msgs.info("Deriving spatial correction for slit {0:d}/{1:d}".format(sl + 1, self.slits.spat_id.size))
            # Get the initial slit locations
            onslit_b_init = (slitid_img_init == slitnum)

            # Synthesize ximg, and edgmask from slit boundaries. Doing this outside this
            # routine would save time. But this is pretty fast, so we just do it here to make the interface simpler.
            spatcoord, edgmask = pixels.ximg_and_edgemask(self.slits_left[:, sl], self.slits_right[:, sl],
                                                          onslit_b_init, trim_edg=trim_edg)

            # Make the model histogram
            xspl = np.linspace(0.0, 1.0, 10 * int(slitlength))  # Sub sample each pixel with 10 subpixels
            modspl = self.caliBrate.flatimages.illumflat_spat_bsplines[sl].value(xspl)[0]
            gradspl = interpolate.interp1d(xspl, np.gradient(modspl) / modspl, kind='linear', bounds_error=False,
                                           fill_value='extrapolate')

            # Ignore skymask
            coord_msk = onslit_b_init & gpm
            hist, _ = np.histogram(spatcoord[coord_msk], bins=spatbins, weights=rawimg[coord_msk])
            cntr, _ = np.histogram(spatcoord[coord_msk], bins=spatbins)
            hist_slit_all = hist / (cntr + (cntr == 0))
            histmod, _ = np.histogram(spatcoord[coord_msk], bins=spatbins, weights=gradspl(spatcoord[coord_msk]))
            hist_model = histmod / (cntr + (cntr == 0))

            # Repeat with skymask
            coord_msk = onslit_b_init & gpm & skymask_now
            hist, _ = np.histogram(spatcoord[coord_msk], bins=spatbins, weights=rawimg[coord_msk])
            cntr, _ = np.histogram(spatcoord[coord_msk], bins=spatbins)
            hist_slit = hist / (cntr + (cntr == 0))

            # Prepare for fit - take the non-zero elements and trim slit edges
            if hist_trim == 0:
                ww = (hist_slit != 0)
                xfit = spat_slit[ww]
                yfit = hist_slit_all[ww]
                mfit = hist_model[ww]
            else:
                ww = (hist_slit[hist_trim:-hist_trim] != 0)
                xfit = spat_slit[hist_trim:-hist_trim][ww]
                yfit = hist_slit_all[hist_trim:-hist_trim][ww]
                mfit = hist_model[hist_trim:-hist_trim][ww]

            # Fit the function
            spat_func = lambda par, ydata, model: par[0]*(1 + par[1] * model) - ydata
            res_lsq = least_squares(spat_func, [np.median(yfit), 0.0], args=(yfit, mfit))
            spatnorm = spat_func(res_lsq.x, 0.0, gradspl(spatcoord[onslit_b_init]))
            spatnorm /= spat_func(res_lsq.x, 0.0, gradspl(0.5))
            # Set the scaling factor
            spatScaleImg[onslit_b_init] = spatnorm
            coeff_fit[sl, :] = res_lsq.x

        if debug:
            from matplotlib import pyplot as plt
            xplt = np.arange(24)
            plt.subplot(121)
            plt.plot(xplt[0::2], coeff_fit[::2, 0], 'rx')
            plt.plot(xplt[1::2], coeff_fit[1::2, 0], 'bx')
            plt.subplot(122)
            plt.plot(xplt[0::2], coeff_fit[::2, 1]/10, 'rx')
            plt.plot(xplt[1::2], coeff_fit[1::2, 1]/10, 'bx')
            plt.show()
            plt.imshow(spatScaleImg, vmin=0.99, vmax=1.01)
            plt.show()
            plt.subplot(133)
            plt.plot(xplt[0::2], coeff_fit[::2, 2], 'rx')
            plt.plot(xplt[1::2], coeff_fit[1::2, 2], 'bx')
            plt.show()
        # Apply the relative scale correction
        self.apply_relative_scale(spatScaleImg)

    def illum_profile_spectral(self, global_sky, skymask=None):
        """Calculate the residual spectral illumination profile using the sky regions.
        This uses the same routine as the flatfield spectral illumination profile.

         Args:
             global_sky (`numpy.ndarray`_):
                Model of the sky
             skymask (`numpy.ndarray`_, optional):
                Mask of sky regions where the spatial illumination will be determined
        """
        trim = self.par['calibrations']['flatfield']['slit_trim']
        ref_idx = self.par['calibrations']['flatfield']['slit_illum_ref_idx']
        gpm = self.sciImg.select_flag(invert=True)
        scaleImg = flatfield.illum_profile_spectral(self.sciImg.image.copy(), self.waveimg, self.slits,
                                                    slit_illum_ref_idx=ref_idx, model=global_sky, gpmask=gpm,
                                                    skymask=skymask, trim=trim, flexure=self.spat_flexure_shift)
        # Now apply the correction to the science frame
        self.apply_relative_scale(scaleImg)

    def joint_skysub(self, skymask=None, update_crmask=True, trim_edg=(0,0),
                     show_fit=False, show=False, show_objs=False, adderr=0.01):
        """ Perform a joint sky model fit to the data. See Reduce.global_skysub()
        for parameter definitions.
        """
        msgs.info("Performing joint global sky subtraction")
        # Mask objects using the skymask? If skymask has been set by objfinding, and masking is requested, then do so
        skymask_now = skymask if (skymask is not None) else np.ones_like(self.sciImg.image, dtype=bool)
        global_sky = np.zeros_like(self.sciImg.image)
        thismask = (self.slitmask > 0)
        inmask = (self.sciImg.select_flag(invert=True) & thismask & skymask_now).astype(np.bool)
        # Convert the wavelength image to A/pixel, registered at pixel 0 (this gives something like
        # the tilts frame, but conserves wavelength position in each slit)
        wavemin = self.waveimg[self.waveimg != 0.0].min()
        tilt_wave = (self.waveimg - wavemin) / (self.waveimg.max() - wavemin)

        # Parameters for a standard star
        sigrej = 3.0
        if self.std_redux:
            sigrej = 7.0
            update_crmask = False
            if not self.par['reduce']['skysub']['global_sky_std']:
                msgs.info('Skipping global sky-subtraction for standard star.')
                return global_sky

        # Iterate to use a model variance image
        numiter = 4
        model_ivar = self.sciImg.ivar
        for nn in range(numiter):
            msgs.info("Performing iterative joint sky subtraction - ITERATION {0:d}/{1:d}".format(nn+1, numiter))
            global_sky[thismask] = skysub.global_skysub(self.sciImg.image, model_ivar, tilt_wave,
                                                             thismask, self.slits_left, self.slits_right, inmask=inmask,
                                                             sigrej=sigrej, trim_edg=trim_edg,
                                                             bsp=self.par['reduce']['skysub']['bspline_spacing'],
                                                             no_poly=self.par['reduce']['skysub']['no_poly'],
                                                             pos_mask=(not self.bkg_redux), show_fit=show_fit)
            # Update the ivar image used in the sky fit
            msgs.info("Updating sky noise model")
            # Choose the highest counts out of sky and object
            counts = global_sky# + np.clip(self.sciImg.image-self.global_sky, 0, None)
            _scale = None if self.sciImg.img_scale is None else self.sciImg.img_scale[thismask]
            # NOTE: darkcurr must be a float for the call below to work.
            var = procimg.variance_model(self.sciImg.base_var[thismask], counts=counts[thismask],
                                         count_scale=_scale, noise_floor=adderr)
            model_ivar[thismask] = utils.inverse(var)
            # var = np.abs(self.global_sky - np.sqrt(2.0) * np.sqrt(self.sciImg.rn2img)) + self.sciImg.rn2img
            # var = var + adderr ** 2 * (np.abs(self.global_sky)) ** 2
            # model_ivar = utils.inverse(var)
            # Redo the relative spectral illumination correction with the improved sky model
            if self.par['scienceframe']['process']['use_specillum']:
                self.illum_profile_spectral(global_sky, skymask=skymask)

        if update_crmask:
            # Find CRs with sky subtraction
            self.sciImg.build_crmask(self.par['scienceframe']['process'],
                                     subtract_img=global_sky)
            # Update the fullmask
            self.sciImg.update_mask_cr(self.sciImg.crmask)

        # Step
        self.steps.append(inspect.stack()[0][3])

        if show:
            sobjs_show = None if show_objs else self.sobjs_obj
            # Global skysub is the first step in a new extraction so clear the channels here
            self.show('global', global_sky=global_sky, slits=True, sobjs=sobjs_show, clear=False)
        return global_sky

    def global_skysub(self, skymask=None, update_crmask=True, trim_edg=(0,0),
                      previous_sky=None, show_fit=False, show=False, show_objs=False):
        """
        Perform global sky subtraction. This IFU-specific routine ensures that the
        edges of the slits are not trimmed, and performs a spatial and spectral
        correction using the sky spectrum, if requested. See Reduce.global_skysub()
        for parameter definitions.
        """
        # Generate a global sky sub for all slits separately
        global_sky_sep = super().global_skysub(skymask=skymask, update_crmask=update_crmask,
                                               trim_edg=trim_edg, show_fit=show_fit, show=show,
                                               show_objs=show_objs)
        # If the joint fit or spec/spat sensitivity corrections are not being performed, return the separate slits sky
        if not self.par['reduce']['skysub']['joint_fit'] and \
                not self.par['scienceframe']['process']['use_specillum'] and \
                not self.par['scienceframe']['process']['use_illumflat']:
            return global_sky_sep

        # Do the spatial scaling first
        # if self.par['scienceframe']['process']['use_illumflat']:
        #     # Perform the correction
        #     self.illum_profile_spatial(skymask=skymask)
        #     # Re-generate a global sky sub for all slits separately
        #     global_sky_sep = Reduce.global_skysub(self, skymask=skymask, update_crmask=update_crmask, trim_edg=trim_edg,
        #                                           show_fit=show_fit, show=show, show_objs=show_objs)

        # Wavelengths (on unmasked slits)
        msgs.info("Generating wavelength image")
        # It's needed in `illum_profile_spectral`
        # TODO maybe would be better to move it inside `illum_profile_spectral`
        self.waveimg = self.wv_calib.build_waveimg(self.tilts, self.slits, spat_flexure=self.spat_flexure_shift)

        if self.par['scienceframe']['process']['use_specillum']:
            self.illum_profile_spectral(global_sky_sep, skymask=skymask)

        # Fit to the sky
        if self.par['reduce']['skysub']['joint_fit']:
            # Use sky information in all slits to perform a joint sky fit
            global_sky = self.joint_skysub(skymask=skymask, update_crmask=update_crmask, trim_edg=trim_edg,
                                                show_fit=show_fit, show=show, show_objs=show_objs)
        else:
            # Re-run global skysub on individual slits, with the science frame now scaled
            global_sky = super().global_skysub(skymask=skymask, update_crmask=update_crmask,
                                                    trim_edg=trim_edg, show_fit=show_fit,
                                                    show=show, show_objs=show_objs)

        # TODO remove? This does not seem to be usable
        # debug = False
        # if debug:
        #     embed()
        #     wavefull = np.linspace(3950, 4450, 10000)
        #     import matplotlib.pylab as pl
        #     from matplotlib import pyplot as plt
        #     colors = pl.cm.jet(np.linspace(0, 1, gdslits.size))
        #     plt.subplot(121)
        #     for sl, slit_idx in enumerate(gdslits):
        #         slit_spat = self.slits.spat_id[slit_idx]
        #         thismask = self.slitmask == slit_spat
        #         wav = self.waveimg[thismask]
        #         flx = global_sky_sep[thismask]
        #         argsrt = np.argsort(wav)
        #         spl = interpolate.interp1d(wav[argsrt], flx[argsrt], bounds_error=False)
        #         if sl == 0:
        #             ref = spl(wavefull)
        #             plt.plot(wavefull, ref / np.nanmedian(ref), color=colors[sl], linestyle=':')
        #         plt.plot(wavefull, spl(wavefull) / ref, color=colors[sl])
        #     plt.subplot(122)
        #     for sl, slit_idx in enumerate(gdslits):
        #         slit_spat = self.slits.spat_id[slit_idx]
        #         thismask = self.slitmask == slit_spat
        #         wav = self.waveimg[thismask]
        #         flx = self.global_sky[thismask]
        #         argsrt = np.argsort(wav)
        #         spl = interpolate.interp1d(wav[argsrt], flx[argsrt], bounds_error=False)
        #         if sl == 0:
        #             ref = spl(wavefull)
        #             plt.plot(wavefull, ref / np.nanmedian(ref), color=colors[sl], linestyle=':')
        #         plt.plot(wavefull, spl(wavefull) / ref, color=colors[sl])
        #         print(sl, np.median(spl(wavefull) / ref))
        #         # plt.plot(wavefull, spl(wavefull), color=colors[sl])
        #
        #     plt.show()

        return global_sky



