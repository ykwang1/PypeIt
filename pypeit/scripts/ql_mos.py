#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-
"""
This script runs PypeIt on a set of MultiSlit images
"""
import argparse

from pypeit import msgs
from pypeit.scripts import run_pypeit

import warnings

def parser(options=None):

    args = argparse.ArgumentParser(description='Script to run PypeIt in QuickLook on a set of MOS files')
    args.add_argument('spectrograph', type=str, help='Name of spectograph, e.g. shane_kast_blue')
    args.add_argument('full_rawpath', type=str, help='Full path to the raw files')
    args.add_argument('arc', type=str, help='Arc frame filename')
    args.add_argument('flat', type=str, help='Flat frame filename')
    args.add_argument('science', type=str, help='Science frame filename')
    args.add_argument('-b', '--box_radius', type=float, help='Set the radius for the boxcar extraction (arcsec)')
    args.add_argument("--user_pixflat", type=str, help="Use a user-supplied pixel flat (e.g. keck_lris_blue)")
    args.add_argument('-d', '--detnum', metavar='D', type=int, nargs='+', default=None,
                      help='One or more detectors to reduce (1-indexed).  If not provided, all '
                           'detectors are reduced.  If the output files exist and -o is used, '
                           'the outputs for the input detector will be replaced.')
    args.add_argument('-i', '--ignore_bad_headers', default=False, action='store_true',
                      help='Ignore headers that do not have the required metadata.  WARNING: '
                           'Use at your own risk.')

    return args.parse_args() if options is None else args.parse_args(options), \
                run_pypeit.specified_args(args, options=options)


def main(args, specified):

    import os
    import numpy as np

    from IPython import embed

    from pypeit import pypeit
    from pypeit import pypeitsetup
    from pypeit.par.pypeitpar import ExecutionPar
    from pypeit.core import framematch

    spec = args.spectrograph

    # Config the run
    cfg_lines = ['[rdx]']
    cfg_lines += ['    spectrograph = {0}'.format(spec)]
    cfg_lines += ['    redux_path = {0}_A'.format(os.path.join(os.getcwd(),spec))]
    cfg_lines += ['[scienceframe]']
    cfg_lines += ['    [[process]]']
    cfg_lines += ['          cr_reject = False']
    if args.user_pixflat is not None:
        cfg_lines += ['[calibrations]']
        cfg_lines += ['    [[flatfield]]']
        cfg_lines += ['        frame = {0}'.format(args.user_pixflat)]
    cfg_lines += ['[reduce]']
    cfg_lines += ['    [[extraction]]']
    cfg_lines += ['         skip_optimal = True']
    if args.box_radius is not None: # Boxcar radius
        cfg_lines += ['    boxcar_radius = {0}'.format(args.box_radius)]
    cfg_lines += ['    [[findobj]]']
    cfg_lines += ['         skip_second_find = True']

    # Data files
    data_files = [os.path.join(args.full_rawpath, args.arc),
                  os.path.join(args.full_rawpath, args.flat),
                  os.path.join(args.full_rawpath,args.science)]

    # Setup
    ps = pypeitsetup.PypeItSetup(data_files, path='./', spectrograph_name=spec,
                                 cfg_lines=cfg_lines)
    ps.build_fitstbl()
    # TODO -- Get the type_bits from  'science'
    bm = framematch.FrameTypeBitMask()
    file_bits = np.zeros(3, dtype=bm.minimum_dtype())
    file_bits[0] = bm.turn_on(file_bits[0], ['arc', 'tilt'])
    file_bits[1] = bm.turn_on(file_bits[1],
                              ['pixelflat', 'trace'] if args.user_pixflat is None else 'trace')
    file_bits[2] = bm.turn_on(file_bits[2], 'science')

    # PypeItSetup sorts according to MJD
    #   Deal with this
    asrt = []
    for ifile in data_files:
        bfile = os.path.basename(ifile)
        idx = ps.fitstbl['filename'].data.tolist().index(bfile)
        asrt.append(idx)
    asrt = np.array(asrt)

    # Set bits
    ps.fitstbl.set_frame_types(file_bits[asrt])
    ps.fitstbl.set_combination_groups()
    # Extras
    ps.fitstbl['setup'] = 'A'

    # Write
    ofiles = ps.fitstbl.write_pypeit('', configs=['A'], write_bkg_pairs=True, cfg_lines=cfg_lines)
    if len(ofiles) > 1:
        msgs.error("Bad things happened..")

    # Instantiate the main pipeline reduction object
    execpar = ExecutionPar(**run_pypeit.exec_kwargs(args, specified))
    execpar['verbosity'] = 2
    execpar['reuse_masters'] = True
    execpar['overwrite'] = True
    execpar['logfile'] = 'mos.log'
    execpar['show'] = False             # False is the default...
    pypeIt = pypeit.PypeIt(ofiles[0], par=execpar)

    pypeIt.reduce_all()
    pypeIt.build_qa()

    return 0

