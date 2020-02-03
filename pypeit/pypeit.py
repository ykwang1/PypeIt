"""
Main driver class for PypeIt run
"""
import time
import os
from collections import OrderedDict

from configobj import ConfigObj

import numpy as np

from astropy.io import fits

from pypeit import msgs
from pypeit import calibrations
from pypeit.images import scienceimage
from pypeit import ginga
from pypeit import reduce
from pypeit.core import qa
from pypeit.core import save
from pypeit import specobjs
from pypeit.spectrographs.util import load_spectrograph

from pypeit.par.util import parse_pypeit_file
from pypeit.par import PypeItPar
from pypeit.metadata import PypeItMetaData

from IPython import embed

class PypeIt:
    """
    This class runs the primary calibration and extraction in PypeIt

    .. todo::
        - Fill in list of attributes!
        - Explain the setup.
        - Give example usage.

    Args:
        pypeit_file (:obj:`str`):
            PypeIt filename.
        verbosity (:obj:`int`, optional):
            Verbosity level of system output.  Can be:

                - 0: No output
                - 1: Minimal output (default)
                - 2: All output

        overwrite (:obj:`bool`, optional):
            Flag to overwrite any existing files/directories.
        reuse_masters (:obj:`bool`, optional):
            Reuse any pre-existing calibration files
        logname (:obj:`str`, optional):
            The name of an ascii log file with the details of the
            reduction.
        show: (:obj:`bool`, optional):
            Show reduction steps via plots (which will block further
            execution until clicked on) and outputs to ginga. Requires
            remote control ginga session via ``ginga --modules=RC &``
        redux_path (:obj:`str`, optional):
            Over-ride reduction path in PypeIt file (e.g. Notebook usage)

    Attributes:
        pypeit_file (:obj:`str`):
            Name of the pypeit file to read.  PypeIt files have a
            specific set of valid formats. A description can be found
            :ref:`pypeit_file`.
        fitstbl (:obj:`pypit.metadata.PypeItMetaData`): holds the meta info

    """
    def __init__(self, pypeit_file, verbosity=2, overwrite=True, reuse_masters=False, logname=None,
                 show=False, redux_path=None):

        # Load
        cfg_lines, data_files, frametype, usrdata, setups \
                = parse_pypeit_file(pypeit_file, runtime=True)
        self.pypeit_file = pypeit_file

        # Spectrograph
        cfg = ConfigObj(cfg_lines)
        spectrograph_name = cfg['rdx']['spectrograph']
        self.spectrograph = load_spectrograph(spectrograph_name, ifile=data_files[0])
        msgs.info('Loaded spectrograph {0}'.format(self.spectrograph.spectrograph))

        # --------------------------------------------------------------
        # Get the full set of PypeIt parameters
        #   - Grab a file for configuration specific parameters. Allowed
        #     frame types for this are defined by the spectrograph
        #     ``config_par_frametype`` method. The order of the frame
        #     types in that list are treated as a prioritization.
        config_ex_file = None
        for f in self.spectrograph.config_par_frametypes():
            indx = [np.any(np.isin(d.split(','),[f])) for d in usrdata['frametype']]
            if not np.any(indx):
                continue
            config_ex_file = data_files[np.where(indx)[0][0]]
            break
        #   - Configuration specific parameters for the spectrograph
        if config_ex_file is not None:
            msgs.info('Setting configuration-specific parameters using {0}'.format(
                      os.path.split(config_ex_file)[1]))
        # NOTE: The call below will fault if the spectrograph requires
        # a file to set the configuration specific parameters
        spectrograph_cfg_lines = self.spectrograph.config_specific_par(config_ex_file).to_config()
        #   - Build the full set, merging with any user-provided
        #     parameters
        self.par = PypeItPar.from_cfg_lines(cfg_lines=spectrograph_cfg_lines, merge_with=cfg_lines)
        msgs.info('Built full PypeIt parameter set.')

        # Check the output paths are ready
        if redux_path is not None:
            self.par['rdx']['redux_path'] = redux_path

        # TODO: Write the full parameter set here?
        # --------------------------------------------------------------

        # --------------------------------------------------------------
        # Build the meta data
        #   - Re-initilize based on the file data
        msgs.info('Compiling metadata')
        self.fitstbl = PypeItMetaData(self.spectrograph, self.par, files=data_files,
                                      usrdata=usrdata, strict=True)
        #   - Interpret automated or user-provided data from the PypeIt
        #   file
        self.fitstbl.finalize_usr_build(frametype, setups[0])
        # --------------------------------------------------------------
        #   - Write .calib file (For QA naming amongst other things)
        calib_file = pypeit_file.replace('.pypeit', '.calib')
        self.fitstbl.write_calib(calib_file)

        # Other Internals
        self.logname = logname
        self.overwrite = overwrite

        # Currently the runtime argument determines the behavior for
        # reuse_masters.
        self.reuse_masters = reuse_masters
        self.show = show

        # TODO: I think this should go back to being an @property
        # method. I also think everything should always be relative to
        # redux_path.
        # Set paths
        if self.par['calibrations']['caldir'] == 'default':
            self.calibrations_path = os.path.join(self.par['rdx']['redux_path'], 'Masters')
        else:
            self.calibrations_path = self.par['calibrations']['caldir']

        # Report paths
        msgs.info('Setting reduction path to {0}'.format(self.par['rdx']['redux_path']))
        msgs.info('Master calibration data output to: {0}'.format(self.calibrations_path))
        msgs.info('Science data output to: {0}'.format(self.science_path))
        msgs.info('Quality assessment plots output to: {0}'.format(self.qa_path))

        # TODO: Should we have separate calibration and science QA
        # directories?  (Are there any "science" QA plots?)

        # Instantiate Calibrations class
        self.caliBrate \
            = calibrations.MultiSlitCalibrations(self.fitstbl, self.par['calibrations'],
                                                 self.spectrograph,
                                                 caldir=self.calibrations_path,
                                                 qadir=self.qa_path,
                                                 reuse_masters=self.reuse_masters,
                                                 show=self.show)
        # Init
        self.verbosity = verbosity
        # TODO: I don't think this ever used

        self.frame = None
        self.det = None

        self.tstart = None
        self.basename = None
        self.sciI = None
        self.obstime = None
        self.ir_redux = False

    @property
    def science_path(self):
        """Return the path to the science directory."""
        return os.path.join(self.par['rdx']['redux_path'], self.par['rdx']['scidir'])

    @property
    def qa_path(self):
        """Return the path to the top-level QA directory."""
        return None if self.par['rdx']['qadir'] is None else \
                    os.path.join(self.par['rdx']['redux_path'], self.par['rdx']['qadir'])

    def build_qa(self):
        """
        Generate QA HTML wrappers
        """
        if self.qa_path is not None:
            msgs.info('Generating QA HTML')
            qa.gen_mf_html(self.pypeit_file, self.qa_path)
            qa.gen_exp_html()

    # TODO: This should go in a more relevant place
    def spec_output_file(self, frame, twod=False):
        """
        Return the path to the spectral output data file.
        
        Args:
            frame (:obj:`int`):
                Frame index from :attr:`fitstbl`.
            twod (:obj:`bool`, optional):
                Name for the 2D output file; 1D file otherwise.
        
        Returns:
            :obj:`str`: The path for the output file
        """
        return os.path.join(self.science_path, 'spec{0}d_{1}.fits'.format('2' if twod else '1',
                                                    self.fitstbl.construct_basename(frame)))

    def outfile_exists(self, frame):
        """
        Check whether the 2D outfile of a given frame already exists.

        Args:
            frame (:obj:`int`):
                Frame index from fitstbl

        Returns:
            :obj:`bool`: True if the 2d file exists, False if it does
            not exist
        """
        return os.path.isfile(self.spec_output_file(frame, twod=True))

    def get_std_outfile(self, standard_frames):
        """
        Grab the output filename from an input list of standard_frame
        indices.

        If more than one index is provided, the first is taken.

        Args:
            standard_frames (:obj:`list`):
                List of indices corresponding to standard stars

        Returns:
            :obj:`str`: Full path to the standard spec1d output file
        """
        # TODO: Need to decide how to associate standards with
        # science frames in the case where there is more than one
        # standard associated with a given science frame.  Below, I
        # just use the first standard

        std_outfile = None
        std_frame = None if len(standard_frames) == 0 else standard_frames[0]
        # Prepare to load up standard?
        if std_frame is not None:
            std_outfile = self.spec_output_file(std_frame) \
                            if isinstance(std_frame, (int,np.integer)) else None
        if std_outfile is not None and not os.path.isfile(std_outfile):
            msgs.error('Could not find standard file: {0}'.format(std_outfile))
        return std_outfile

    def reduce_all(self):
        """
        Main driver to perform all reductions.

        The algorithm does the following:

            - Checks that the appropriate parameter sets are defined.
            - If no standard or science frames are present, all the
              calibration groups are reduced. See
              :func:`reduce_calibration_group`.
            - If standard frames are available, all of the
              combination groups with standard-star exposures are
              reduced. See :func:`reduce_combination_group`.
            - If science frames are available, all of the
              combination groups with science exposures are reduced.
              If standard frames were also reduced, the reduced
              standard output is used in the reduction. See
              :func:`reduce_combination_group`.

        """
        # Validate the parameter set
        # TODO: Is this still correct?
        # TODO: Can this be moved to the __init__?
        required = ['rdx', 'calibrations', 'scienceframe', 'reduce', 'flexure', 'fluxcalib']
        can_be_None = ['flexure', 'fluxcalib']
        self.par.validate_keys(required=required, can_be_None=can_be_None)

        self.tstart = time.time()

        # Find the standard and science frames
        is_standard = self.fitstbl.find_frames('standard')
        is_science = self.fitstbl.find_frames('science')

        # Frame indices
        frame_indx = np.arange(len(self.fitstbl))

        # Allow the code to proceed by only reducing the calibrations
        if np.sum(is_standard | is_science) == 0:
            msgs.warn('No science or standard-star frames found.  Only calibration frames will '
                      'be reduced.')
            for i in range(self.fitstbl.n_calib_groups):
                self.reduce_calibration_group(i)

        # Reduce any standards first
        for i in range(self.fitstbl.n_calib_groups * int(np.sum(is_standard) > 0)):
            # Reduce all the standard frames, loop on unique comb_id
            grp_standards = frame_indx[is_standard & self.fitstbl.find_calib_group(i)]
            for comb_id in enumerate(np.unique(self.fitstbl['comb_id'][grp_standards])):
                self.reduce_combination_group(comb_id)

        # Define the standard star to use during the science
        # reductions.
        # TODO: Allow for standards to be science-frame specific? If
        # there's more than one standard, this just uses the first one.
        stdfile = self.get_std_outfile(np.where(is_standard)[0]) if np.any(is_standard) else None

        # Reduce any science frames
        for i in range(self.fitstbl.n_calib_groups * int(np.sum(is_science) > 0)):
            # Reduce all the standard frames, loop on unique comb_id
            grp_science = frame_indx[is_science & self.fitstbl.find_calib_group(i)]
            for comb_id in enumerate(np.unique(self.fitstbl['comb_id'][grp_science])):
                self.reduce_combination_group(comb_id, stdfile=stdfile)

        # Done
        msgs.info('Data reduction complete')
        self.print_end_time()

    # TODO: Let the calibration group be something other than a
    # 0-indexed running number?
    def reduce_calibration_group(self, grp):
        """
        Only reduce the calibrations for the specified group.

        This method follows the same algorithm as
        :func:`reduce_exposure` without doing any science frame
        reductions.

        Args:
            grp (:obj:`int`):
                0-indexed calibration group to reduce.
        """
        # Check input
        if grp < 0 or grp >= self.fitstbl.n_calib_groups:
            msgs.error('Calibration group {0} is undefined!'.format(grp))

        # Find all the frames in this calibration group
        grp_frames = np.arange(len(self.fitstbl))[self.fitstbl.find_calib_group(grp)]
        if len(grp_frames) == 0:
            return

        # If show is set, clear the ginga channels
        if self.show:
            # TODO: Put this in a try/except block?
            ginga.clear_all()

        # Print status message
        msgs.info('Reducing calibration group {0}'.format(grp))

        # Find the detectors to reduce
        detectors = PypeIt.select_detectors(detnum=self.par['rdx']['detnum'],
                                            ndet=self.spectrograph.ndet)
        if len(detectors) != self.spectrograph.ndet:
            msgs.warn('Not reducing detectors: {0}'.format(' '.join([str(d) for d in 
                                set(np.arange(self.spectrograph.ndet))-set(detectors)])))

        # Loop on Detectors
        # TODO: Attempt to put in a multiprocessing call here?
        for self.det in detectors:
            msgs.info("Working on detector {0}".format(self.det))
            # Prep the calibrations for this group/detector
            self.caliBrate.set_config(grp_frames[0], self.det, self.par['calibrations'],
                                      calib_id=grp)
            # Calibrate
            self.caliBrate.run_the_steps()

    def reduce_combination_group(self, comb_id, stdfile=None):
        """
        Reduce a combination group.

        It is expected that science frames are reduced *after*
        standard frames. To reduce all combination frames for both
        standard and science frames in the correct order, use
        :func:`reduce_all`.

        Args:
            comb_id (:obj:`int`):
                The combination ID to reduce.
            stdfile (:obj:`str`, optional):
                If reducing a set of science frames, this is the name
                of a file with a standard star reduction to use in
                :func:`reduce_exposure`.
        """
        # Find the frames to combine and any associated background
        # frames
        frames, bg_frames = self.fitstbl.find_frame_combinations(comb_id)
        if not self.outfile_exists(frames[0]) or self.overwrite:
            std_dict = self.reduce_exposure(frames, bg_frames=bg_frames, std_outfile=stdfile)
            # TODO come up with sensible naming convention for
            # save_exposure for combined files
            self.save_exposure(frames[0], std_dict, self.basename)
        else:
            # TODO: Should this fault?
            msgs.warn('Output file: {:s} already exists'.format(
                        self.fitstbl.construct_basename(frames[0]))
                      + '. Set overwrite=True to recreate and overwrite.')

    # This is a static method to allow for use in coadding script 
    @staticmethod
    def select_detectors(detnum=None, ndet=1):
        """
        Return the 1-indexed list of detectors to reduce.

        Args:
            detnum (:obj:`int`, :obj:`list`, optional):
                One or more detectors to reduce.  If None, return the
                full list for the provided number of detectors (`ndet`).
            ndet (:obj:`int`, optional):
                The number of detectors for this instrument.  Only used
                if `detnum is None`.

        Returns:
            :obj:`list`: List of detectors to be reduced.
        """
        if detnum is None:
            return np.arange(1, ndet+1).tolist()
        return [detnum] if isinstance(detnum, int) else detnum

    # TODO: Need to rename this method.
    def reduce_exposure(self, frames, bg_frames=None, std_outfile=None):
        """
        Reduce a set of standard or science frames.

        The provided ``frames`` can be a single frame or a list of
        frames to be combined. (To reduce exposures based on their
        combination group, see :func:`reduce_combination_group`.) All
        of the frames should be standard or science frames (but this
        is not checked). First, the calibration groups associated
        with these frames are reduced (see
        :func:`pypeit.calibrations.Calibrations.run_the_steps`) and
        then the standard/science frames are reduced using
        :func:`extract_one`.

        Importantly, the *first* frame is passed to
        :func:`pypeit.calibrations.Calibrations.set_config`; see the
        description of that function for how this frame is used. For
        additional uses of only the *first* frame, see
        :func:`extract_one`.
        
        Args:
            frames (array-like):
                0-indexed row(s) in :attr:`fitstbl` with the frames
                to reduce. Note that the first frame in this list is
                passed to
                :func:`pypeit.calibrations.Calibrations.set_config`.
            bg_frames (array-like, optional):
                List of frame indices for the background. Can be None
                or an empty list.
            std_outfile (:obj:`str`, optional):
                File with a previously reduced standard spectrum from
                PypeIt.

        Returns:
            :obj:`dict`: The dictionary containing the primary
            outputs of extraction.
        """
        # Check the input
        target = self.fitstbl['target'][frames[0]].strip()
        if len(frames) > 0 \
                and np.any([self.fitstbl['target'][f].strip() != target for f in frames]):
            # TODO: Fault?
            msgs.warn('All of the frames in this group do not have the same target!')

        # If show is set, clear the ginga channels at the start of each new sci_ID
        if self.show:
            # TODO: Put this in a try/except block?
            ginga.clear_all()

        # Are background frames provided?
        # TODO: Do we need both this and self.ir_redux
        has_bg = bg_frames is not None and len(bg_frames) > 0

        # If so, assume this is an IR reduction?
        # TODO: Why specific to IR?
        self.ir_redux = has_bg

        # Print status message
        # TODO: Print these when the frames are actually combined,
        # backgrounds are used, etc?
        msgs_string = 'Reducing target {:s}'.format(target) + msgs.newline() \
                        + 'Combining frames:' + msgs.newline()
        for iframe in frames:
            msgs_string += '{0:s}'.format(self.fitstbl['filename'][iframe]) + msgs.newline()
        if has_bg:
            msgs_string += msgs.newline() + 'Using background from frames:' + msgs.newline()
            for iframe in bg_frames:
                msgs_string += '{0:s}'.format(self.fitstbl['filename'][iframe]) + msgs.newline()
        msgs.info(msgs_string)

        # Find the detectors to reduce
        detectors = PypeIt.select_detectors(detnum=self.par['rdx']['detnum'],
                                            ndet=self.spectrograph.ndet)
        if len(detectors) != self.spectrograph.ndet:
            msgs.warn('Not reducing detectors: {0}'.format(' '.join([str(d) for d in 
                                set(np.arange(self.spectrograph.ndet))-set(detectors)])))

        # Instantiate the dictionary with the science data
        # TODO: JFH Why does this need to be ordered?
        sci_dict = OrderedDict()  # This needs to be ordered
        sci_dict['meta'] = {}
        sci_dict['meta']['ir_redux'] = self.ir_redux

        # Loop on Detectors
        # TODO: Attempt to put in a multiprocessing call here?
        for self.det in detectors:
            msgs.info("Working on detector {0}".format(self.det))
            sci_dict[self.det] = {}

            # NOTE: From set_config, frame[0] is used to set the
            # calibration group, the master-frame key, and it
            # identifies the frame that is passed to
            # :func:`pypeit.spectrographs.spectrograph.Spectrograph.bpm`
            # to set the bad-pixel mask.

            # Calibrate
            self.caliBrate.set_config(frames[0], self.det, self.par['calibrations'])
            self.caliBrate.run_the_steps()

            # Extract
            # TODO: pass back the background frame, pass in background
            # files as an argument. extract one takes a file list as an
            # argument and instantiates science within
            sci_dict[self.det]['sciimg'], sci_dict[self.det]['sciivar'], \
                sci_dict[self.det]['skymodel'], sci_dict[self.det]['objmodel'], \
                sci_dict[self.det]['ivarmodel'], sci_dict[self.det]['outmask'], \
                sci_dict[self.det]['specobjs'], \
                        = self.extract_one(frames, self.det, bg_frames=bg_frames,
                                           std_outfile=std_outfile)
            # JFH TODO write out the background frame?

        # Return
        return sci_dict

    def get_sci_metadata(self, frame, det):
        """
        Grab the meta data for a given science frame and detector.

        Args:
            frame (:obj:`int`):
                0-indexed frame number.
            det (:obj:`int`):
                1-indexed detector number.

        Returns:
            :obj:`tuple`: Returns (1) the object type (either
            standard or science), (2) the setup name (see
            :func:`pypeit.metadata.PypeItMetaData.master_key`), (3)
            the time of the observation as an `astropy.time.Time`_
            object (see
            :func:`pypeit.metadata.PypeItMetaData.construct_obstime),
            (4) the "base name" for the output science files (see
            :func:`pypeit.metadata.PypeItMetaData.construct_basename),
            and (5) the detector binning (e.g., '1,1').
        """
        # Set binning, obstime, basename, and objtype
        binning = self.fitstbl['binning'][frame]
        obstime  = self.fitstbl.construct_obstime(frame)
        basename = self.fitstbl.construct_basename(frame, obstime=obstime)
        objtype  = self.fitstbl['frametype'][frame]
        if 'science' in objtype:
            objtype_out = 'science'
        elif 'standard' in objtype:
            objtype_out = 'standard'
        else:
            msgs.error('Unrecognized object type.  Must be science or standard.')
        setup = self.fitstbl.master_key(frame, det=det)
        return objtype_out, setup, obstime, basename, binning

    # TODO: Why can't std_redux be handled by the calling function?
    # What does the explanation of std_redux mean? Add the shape
    # information to the return docstring.
    def get_std_trace(self, std_redux, det, std_outfile):
        """
        Returns the trace of the standard if it is applicable to the
        current reduction.

        Args:
            std_redux (:obj:`bool`):
                If False, proceed
            det (:obj:`int`):
                1-indexed detector number.
            std_outfile (:obj:`str`):
                "spec1d" file  with the standard star spectrum.

        Returns:
            `numpy.ndarray`_: The trace of the standard star on the
            input detector.
        """
        # TODO: The nesting here is largely unnecessary. Use returns to
        # interupt the control flow.
        if std_redux is False and std_outfile is not None:
            sobjs = specobjs.SpecObjs.from_fitsfile(std_outfile)
            # Does the detector match?
            # TODO: Instrument specific logic here could be implemented
            # with the parset. For example LRIS-B or LRIS-R we we would
            # use the standard from another detector
            this_det = sobjs.DET == det
            if np.any(this_det):
                sobjs_det = sobjs[this_det]
                sobjs_std = sobjs_det.get_std()
                std_trace = sobjs_std.TRACE_SPAT
                # flatten the array if this multislit
                if 'MultiSlit' in self.spectrograph.pypeline:
                    std_trace = std_trace.flatten()
                elif 'Echelle' in self.spectrograph.pypeline:
                    std_trace = std_trace.T
                else:
                    msgs.error('Unrecognized pypeline')
            else:
                std_trace = None
        else:
            std_trace = None

        return std_trace

    # TODO: Need to rename this method.
    # TODO: Can "frames" be an int?
    def extract_one(self, frames, det, bg_frames=None, std_outfile=None):
        """
        Reduce a single detector for a set of frames.

        The calibrations must have already been completed before
        calling this method. To ensure the calibrations are reduced
        before this method is called, use :func:`reduce_exposure`. To
        reduce exposures based on their combination group, use
        :func:`reduce_combination_group`.

        The list of frames (or observing set) should contain a list
        of science (or standard) exposures in a single combination
        group. Importantly, the *first* frame is passed to
        :func:`get_sci_metadata` and used to set the relevant RA/DEC.

        Args:
            frames (:obj:`list`):
                List of frames to extract. The frames are stacked if
                more than one is provided.
            det (:obj:`int`):
                1-indexed detector.
            bg_frames (:obj:`list`, optional):
                List of frames to use as the background. Can be None
                or an empty list.
            std_outfile (:obj:`str`, optional):
                Filename for the standard star spec1d file. Passed
                directly to :func:`get_std_trace`.

        Returns:
            :obj:`tuple`: Returns six `numpy.ndarray`_ objects and a
            :class:`pypeit.specobjs.SpecObjs` object with the
            extracted spectra from this exposure/detector pair. The
            six `numpy.ndarray`_ objects are (1) the science image,
            (2) its inverse variance, (3) the sky model, (4) the
            object model, (5) the model inverse variance, and (6) the
            mask.

        """
        # Grab some meta-data needed for the reduction from the fitstbl
        self.objtype, self.setup, self.obstime, self.basename, self.binning \
                = self.get_sci_metadata(frames[0], det)
        # Is this a standard star?
        self.std_redux = 'standard' in self.objtype
        # Get the standard trace if need be
        std_trace = self.get_std_trace(self.std_redux, det, std_outfile)

        # Force the illumination flat to be None if the user doesn't
        # want to apply the correction.
        illum_flat = self.caliBrate.msillumflat \
                        if self.par['calibrations']['flatfield']['illumflatten'] else None

        # TODO: report if illum_flat is None? Done elsewhere, but maybe
        # want to do so here either only or also.

        # Build Science image
        sci_files = self.fitstbl.frame_paths(frames)
        self.sciImg = scienceimage.build_from_file_list(
            self.spectrograph, det, self.par['scienceframe']['process'],
            self.caliBrate.msbpm, sci_files, self.caliBrate.msbias,
            self.caliBrate.mspixelflat, illum_flat=illum_flat)

        # Background Image?
        if bg_frames is not None and len(bg_frames) > 0:
            bg_file_list = self.fitstbl.frame_paths(bg_frames)
            self.sciImg = self.sciImg - scienceimage.build_from_file_list(
                self.spectrograph, det, self.par['scienceframe']['process'],
                self.caliBrate.msbpm, bg_file_list, self.caliBrate.msbias,
                self.caliBrate.mspixelflat, illum_flat=illum_flat)

        # Update mask for slitmask; uses pad in EdgeTraceSetPar
        self.sciImg.update_mask_slitmask(self.caliBrate.slits.slit_img())

        # For QA on crash
        msgs.sciexp = self.sciImg

        # Instantiate Reduce object
        # TODO: Do we need these maskslits to be part of self?
        self.maskslits = self.caliBrate.slits.mask.copy()
        # Required for pypeline specific object
        # TODO -- caliBrate should be replaced by the ~3 primary Objects needed
        #   once we have the data models in place.
        self.redux = reduce.instantiate_me(self.sciImg, self.spectrograph,
                                           self.par, self.caliBrate,
                                           maskslits=self.maskslits,
                                           ir_redux=self.ir_redux,
                                           std_redux=self.std_redux,
                                           objtype=self.objtype,
                                           setup=self.setup,
                                           show=self.show,
                                           det=det, binning=self.binning)
        # Show?
        if self.show:
            self.redux.show('image', image=self.sciImg.image, chname='processed',
                            slits=True, clear=True)

        # Prep for manual extraction (if requested)
        manual_extract_dict = self.fitstbl.get_manual_extract(frames, det)

        self.skymodel, self.objmodel, self.ivarmodel, self.outmask, self.sobjs \
                = self.redux.run(std_trace=std_trace, manual_extract_dict=manual_extract_dict,
                                 show_peaks=self.show, basename=self.basename,
                                 ra=self.fitstbl["ra"][frames[0]],
                                 dec=self.fitstbl["dec"][frames[0]], obstime=self.obstime)

        # Return
        return self.sciImg.image, self.sciImg.ivar, self.skymodel, self.objmodel, \
                    self.ivarmodel, self.outmask, self.sobjs

    # TODO: Why not use self.frame?
    def save_exposure(self, frame, sci_dict, basename):
        """
        Save the extraction outputs.

        This method should be called after the data have been reduced
        using :func:`extract_one`. To reduce the relevant
        calibrations, run the science reduction, and then save the
        result, use :func:`reduce_exposure`.

        This is a basic wrapper for
        :func:`pypeit.core.save.save_all`.

        Args:
            frame (:obj:`int`):
                0-indexed row in the metadata table with the frame
                that has been reduced.
            sci_dict (:obj:`dict`):
                Dictionary containing the primary outputs of
                extraction
            basename (:obj:`str`):
                The root name for the output file.
        """
        # TODO: Need some checks here that the exposure has been reduced

        # Determine the headers
        head1d = self.fitstbl[frame]
        # Need raw file header information
        rawfile = self.fitstbl.frame_paths(frame)
        head2d = fits.getheader(rawfile, ext=self.spectrograph.primary_hdrext)
        refframe = 'pixel' if self.caliBrate.par['wavelengths']['reference'] == 'pixel' else \
            self.caliBrate.par['wavelengths']['frame']

        # Determine the paths/filenames
        save.save_all(sci_dict, self.caliBrate.master_key_dict, self.caliBrate.master_dir,
                      self.spectrograph, head1d, head2d, self.science_path, basename,
                      update_det=self.par['rdx']['detnum'], binning=self.fitstbl['binning'][frame])

    def msgs_reset(self):
        """
        Reset the pypeit logger.
        """
        # Reset the global logger
        msgs.reset(log=self.logname, verbosity=self.verbosity)
        msgs.pypeit_file = self.pypeit_file

    def print_end_time(self):
        """
        Print the elapsed time.

        This uses the current value of :attr:`tstart` (currently only
        defined in :func:`reduce_all`).
        """
        # Capture the end time and print it to user
        tend = time.time()
        codetime = tend-self.tstart
        if codetime < 60.0:
            msgs.info('Execution time: {0:.2f}s'.format(codetime))
        elif codetime/60.0 < 60.0:
            mns = int(codetime/60.0)
            scs = codetime - 60.0*mns
            msgs.info('Execution time: {0:d}m {1:.2f}s'.format(mns, scs))
        else:
            hrs = int(codetime/3600.0)
            mns = int(60.0*(codetime/3600.0 - hrs))
            scs = codetime - 60.0*mns - 3600.0*hrs
            msgs.info('Execution time: {0:d}h {1:d}m {2:.2f}s'.format(hrs, mns, scs))

    # TODO: Move this to fitstbl?
    def show_science(self):
        """
        Simple print of science frames.
        """
        indx = self.fitstbl.find_frames('science')
        print(self.fitstbl[['target','ra','dec','exptime','dispname']][indx])

    def __repr__(self):
        # Generate sets string
        return '<{:s}: pypeit_file={}>'.format(self.__class__.__name__, self.pypeit_file)


