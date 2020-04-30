.. highlight:: rest

*************
Magellan Mage
*************


Overview
========

This file summarizes several instrument specific
settings that are related to Magellan/Mage.

Setup
=====

You will need to run :ref:`pypeit_setup` with `-e c1.fits`
or `-e c1.fits.gz` to restrict the code to view only one of
the 8 detectors.  Here is an example call::

    pypeit_setup -s magellan_imacs -r /data/Magellan/IMACS/ut191017_18/ift0 -e c1.fits

Note, you can use any of the detectors with this, e.g. replace `c1` with `c3`.

The resulting :doc:`pypeit_file` will only show the files for the single
detector but the code will attempt to reduce them all.

Short slits
===========

There are several issues related to the very short
slits of Magellan/Mage  (34 pixels or 10" unbinned).

Find Objects
------------

To have enough slit to 'properly' find objects,
we restrict the find_trim_edge parameter, i.e.::

    par['scienceimage']['find_trim_edge'] = (4,4)    # Slit is too short to trim 5,5 especially with 2x binning

For spatial binning, we recommend you to further reduce
this by the binning factor.
