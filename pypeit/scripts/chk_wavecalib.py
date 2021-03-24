"""
This script displays the flat images in an RC Ginga window.
"""

def parse_args(options=None, return_parser=False):

    import argparse

    parser = argparse.ArgumentParser(description='Print QA on Wavelength Calib to the screen',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('master_file', type=str,
                        help='PypeIt MasterWaveCalib file [e.g. MasterWaveCalib_A_1_01.fits]')
    parser.add_argument('--plot', type=int,  
                        help='Plot input SpatID')

    if return_parser:
        return parser

    return parser.parse_args() if options is None else parser.parse_args(options)


def main(args):

    import numpy as np
    from pypeit import wavecalib
    from matplotlib import pyplot as plt
    
    from IPython import embed

    # Load
    waveCalib = wavecalib.WaveCalib.from_file(args.master_file)#, chk_version=(not args.try_old))

    # Do it
    waveCalib.print_diagnostics()

    # Plot?
    if args.plot is not None:
        plt.clf()
        ax = plt.gca()
        spec = waveCalib.wv_fits[args.plot].spec
        ax.plot(np.arange(spec.size), spec)
        plt.show()


def entry_point():
    main(parse_args())


if __name__ == '__main__':
    entry_point()
