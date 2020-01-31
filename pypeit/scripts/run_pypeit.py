#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-
"""
This script runs PypeIt
"""
import argparse

from pypeit import msgs

import warnings

def parser(options=None):

    parser = argparse.ArgumentParser(description=msgs.usage('PypeIt'),
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('pypeit_file', type=str,
                        help='PypeIt reduction file (must have .pypeit extension)')
    parser.add_argument('-v', '--verbosity', type=int, default=2,
                        help='Verbosity level between 0 [none] and 2 [all]')
    # JFH TODO Are the -t and -r keyword still valid given that run_pypeit no longer runs setup?
#    parser.add_argument('-t', '--hdrframetype', default=False, action='store_true',
#                        help='Use file headers and the instument-specific keywords to determine'
#                             'the type of each frame')
#    parser.add_argument('-r', '--sort_dir', default=None,
#                        help='Directory used to store the sorted files.  Default is to omit '
#                             'writing these files.')
    parser.add_argument('-m', '--use_masters', default=False, action='store_true',
                        help='Load previously generated MasterFrames')

    # TODO: Does the ginga window still need to be opened before
    # calling `run_pypeit`?
    parser.add_argument('-s', '--show', default=False, action='store_true',
                        help='Show reduction steps via plots (which will block further execution '
                             'until clicked on) and outputs to ginga. Requires remote control '
                             'ginga session via "ginga --modules=RC &"')

    # JFH Should the default now be true with the new definition.
    # TODO: Is this propagated to all the relevant modules?
    parser.add_argument('-o', '--overwrite', default=False, action='store_true',
                        help='Overwrite any existing files/directories')

#    group = parser.add_mutually_exclusive_group()
#    group.add_argument('-p', '--prep_setup', default=False, action='store_true',
#                       help='Run pypeit to prepare the setup only')
#    group.add_argument('-c', '--calcheck', default=False, action='store_true',
#                       help='Run pypeit only as a check on the calibrations')

    parser.add_argument('-d', '--detector', default=None,
                        help='Detector to limit reductions on. If the output files exist and -o '
                             'is used, the outputs for the input detector will be replaced.')

#    parser.add_argument('-q', '--quick', default=False, help='Quick reduction',
#                        action='store_true')
#    parser.add_argument('-c', '--cpus', default=False, action='store_true',
#                         help='Number of CPUs for parallel processing')
#    parser.print_help()

    parser.add_argument('-q', '--noqa', dest='qa', default=True, action='store_false',
                        help='Do not produce any qa plots.  Can also set qadir=none (the little n '
                             'is important) in the pypeit file.')
    parser.add_argument('-l', '--nolog', dest='log', default=True, action='store_false',
                        help='By default, a log of the printed messages are written to a file '
                             'with the same name as your pypeit file but with a \'.log\' '
                             'extension.  This suppresses output to the log file.')

    return parser.parse_args() if options is None else parser.parse_args(options)


def main(args):

    import os
    import sys
    import traceback

    from pypeit import pypeit
    from pypeit import pypeitsetup
    from pypeit import debugger

    # Initiate logging for bugs and command line help
    # These messages will not be saved to a log file
    # Set the default variables
#    qck = False
#    cpu = 1
    #vrb = 2

    # Load options from command line
    splitnm = os.path.splitext(args.pypeit_file)
    if splitnm[1] != '.pypeit':
        msgs.error("Bad extension for PypeIt reduction file."+msgs.newline()+".pypeit is required")
    logname = splitnm[0] + ".log" if args.log else None

    # Instantiate the main pipeline reduction object
    pypeIt = pypeit.PypeIt(args.pypeit_file, verbosity=args.verbosity,
                           reuse_masters=args.use_masters, overwrite=args.overwrite,
                           logname=logname, show=args.show)

    # JFH I don't see why this is an optional argument here. We could
    # allow the user to modify an infinite number of parameters from
    # the command line? Why do we have the PypeIt file then? This
    # detector can be set in the pypeit file. Detector?

    # KBW: There are certain parameters that I think are useful to
    # modify from the command line so that you don't have to edit the
    # pypeit file. You're right that this could get unwieldy, but a
    # small subset of parameters make sense to me. Detector is one, QA
    # is another.
    if args.detector is not None:
        msgs.info("Restricting reductions to detector={}".format(args.detector))
        pypeIt.par['rdx']['detnum'] = int(args.detector)
    if not args.qa or pypeIt.par['rdx']['qadir'] is None:
        msgs.info('No QA plots will be produced.')
        pypeIt.par['rdx']['qadir'] = None

    # NOTE: Messages moved to these methods
    pypeIt.reduce_all()
    pypeIt.build_qa()

    return 0

