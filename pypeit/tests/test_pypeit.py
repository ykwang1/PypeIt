"""
Module to run tests on arsave
"""
import os
import shutil

from IPython import embed

import numpy as np

import pytest

from pypeit.par.pypeitpar import ExecutionPar
from pypeit.par.util import make_pypeit_file
from pypeit import pypeitsetup
from pypeit.pypeit import PypeIt
from pypeit.scripts import run_pypeit

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_initialization():
    """ Load input PYPIT file
    """
    # Generate a PYPIT file
    pypeit_file = data_path('test.pypeit')
    make_pypeit_file(pypeit_file, 'shane_kast_blue', [data_path('b*fits.gz')], setup_mode=True)

    # Perform the setup
    setup = pypeitsetup.PypeItSetup.from_pypeit_file(pypeit_file)
    par, spectrograph, fitstbl = setup.run(sort_dir=data_path(''))

    # Test
    assert spectrograph.spectrograph == 'shane_kast_blue'
    assert len(fitstbl) == 2
    assert par['calibrations']['arcframe']['number'] == 1

    # Clean-up
    os.remove(data_path('test.calib'))
    os.remove(data_path('test.pypeit'))

def test_select_detectors():
    # Generate a PYPIT file
    pypeit_file = data_path('test.pypeit')
    make_pypeit_file(pypeit_file, 'shane_kast_blue', [data_path('b*fits.gz')], setup_mode=True)

    # Perform the setup
    setup = pypeitsetup.PypeItSetup.from_pypeit_file(pypeit_file)
    par, spectrograph, fitstbl = setup.run(sort_dir=data_path(''))

    assert PypeIt.select_detectors(detnum=par['rdx']['detnum'], ndet=spectrograph.ndet) == [1], \
            'Incorrect detectors selected.'

    # Clean-up
    os.remove(data_path('test.calib'))
    os.remove(data_path('test.pypeit'))

    assert np.array_equal(PypeIt.select_detectors(), [1]), 'Incorrect detectors selected.'
    assert np.array_equal(PypeIt.select_detectors(detnum=3, ndet=5), [3]), \
            'Incorrect detectors selected.'
    assert np.array_equal(PypeIt.select_detectors(ndet=5), [1,2,3,4,5]), \
            'Incorrect detectors selected.'
    assert np.array_equal(PypeIt.select_detectors(detnum=[1,3]), [1,3]), \
            'Incorrect detectors selected.'

def test_arg_input():

    cfg_lines = ['[rdx]']
    cfg_lines += ['    spectrograph = shane_kast_blue']
    cfg_lines += ['    detnum = 2, 4']
    cfg_lines += ['    qadir = none']
    cfg_lines += ['    overwrite = True']
    cfg_lines += ['    redux_path = files/shane_kast_blue_A']

    # Set file names and remove them if they exist
    pypeit_file = data_path('test.pypeit')
    if os.path.isfile(pypeit_file):
        os.remove(pypeit_file)
    calib_file = data_path('test.calib')
    if os.path.isfile(calib_file):
        os.remove(calib_file)
    setup_dir = data_path('shane_kast_blue_A')
    if os.path.isdir(setup_dir):
        shutil.rmtree(setup_dir)

    # Generate a PYPIT file
    make_pypeit_file(pypeit_file, 'shane_kast_blue', [data_path('b*fits.gz')], cfg_lines=cfg_lines,
                     setup_mode=True)

    # Perform the setup
    ps = pypeitsetup.PypeItSetup.from_pypeit_file(pypeit_file)
    ps.run(sort_dir=data_path(''))
    ofile = data_path('shane_kast_blue_A/shane_kast_blue_A.pypeit')
    if os.path.isfile(ofile):
        os.remove(ofile)
    ps.fitstbl.write_pypeit(pypeit_file, cfg_lines=ps.user_cfg, write_bkg_pairs=True)
    pypeit_file = ofile

    assert os.path.isfile(pypeit_file), 'pypeit file was not written'

    # Mimic the run_pypeit command-line arguments
    detector = [1,3]
    redux_path = os.getcwd()
    args, specified = run_pypeit.parser([pypeit_file, '-d'] + [str(d) for d in detector]
                                         + ['-p', redux_path, '-i'])
    execpar = ExecutionPar(**run_pypeit.exec_kwargs(args, specified))
    pypeIt = PypeIt(pypeit_file, par=execpar)

    assert pypeIt.par['rdx']['detnum'] == [1,3], 'Detectors not specified correctly'
    assert pypeIt.par['rdx']['overwrite'], 'Overwrite not set correctly'
    assert pypeIt.par['rdx']['redux_path'] == redux_path, 'Reduction path not set correctly'

    # Clean-up
    os.remove(data_path('test.calib'))
    os.remove(data_path('test.pypeit'))
    shutil.rmtree(setup_dir)
