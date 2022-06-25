'''
2020-04-30, Dennis Alp, dalp@kth.se
python -W ignore xmatch.py
Cross-match observations with known SNe and GRBs.



Observation logs and catalogs are dictionaries that have the following keys:
obs: Observation ID (a string)
tt0: Time of observation (from astropy.time import Time)
exp: Exposure time in seconds
poi: Pointing (coordinates, to compute off-axis angle)
OR (if a catalog provides sources and not detections)
nam: Name of the source
coo: Coordinates

and the following optional keys:
tar: Target
coo: Coorindates (for detected sources, to compute offset)
pub: Date when the data is made public (else assuming 'yes'; a string)



SNe and GRB catalogs are dictionaries with the following keys:
nam: Name
coo: Coordinates (position of transient)
tt0: Time of transient

and the following optional keys:
hos: Host galaxy
typ: Subtype of SNe or GRB
'''

from __future__ import division, print_function
import os
from pdb import set_trace as db
import sys
from glob import glob
import time
from datetime import date

import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.io import fits
from astropy import units
from astropy.time import Time

from hlp import *




################################################################
xmm_log = get_xmm_log()
nus_log = get_nus_log()

xmm_cat = get_xmm_cat()
xrt_cat = get_xrt_cat()
cxo_cat = get_cxo_cat()

sne_cat = get_sne_cat()
grb_cat = get_grb_cat()
ext_cat = get_ext_cat()

################################################################
def xmm():
    lab = 'sne_observed_by_xmm'
    txt = 'This is a list of SNe within the FoV.'
    lim = 15*units.arcmin
    mk_xmatch(xmm_log, sne_cat, lab, txt, lim)
    
    lab = 'sne_detected_by_xmm'
    txt = 'This is a list of SNe matched with catalog detections (8 arcsec).'
    lim = 8*units.arcsec
    mk_xmatch(xmm_cat, sne_cat, lab, txt, lim)
    
    lab = 'grb_observed_by_xmm'
    txt = 'This is a list of GRBs within the FoV.'
    lim = 15*units.arcmin
    mk_xmatch(xmm_log, grb_cat, lab, txt, lim)
    
    lab = 'grb_detected_by_xmm'
    txt = 'This is a list of GRBs matched with catalog detections (20 arcsec).'
    lim = 20*units.arcsec
    mk_xmatch(xmm_cat, grb_cat, lab, txt, lim)

def nustar():
    lab = 'sne_observed_by_nustar'
    txt = 'This is a list of SNe within the FoV.'
    lim = 7*units.arcmin
    mk_xmatch(nus_log, sne_cat, lab, txt, lim)
    
    lab = 'grb_observed_by_nustar'
    txt = 'This is a list of GRBs within the FoV.'
    lim = 7*units.arcmin
    mk_xmatch(nus_log, grb_cat, lab, txt, lim)

    lab = 'ext_observed_by_nustar'
    txt = 'This is a list of GRBs within the extended FoV (0-reflection events, aperture photons).'
    lim = 8*units.deg
    mk_xmatch(nus_log, ext_cat, lab, txt, lim)
    
def swift():
    lab = 'sne_detected_by_xrt'
    txt = 'This is a list of SNe matched with catalog detections (8 arcsec).\nUsing 2MASS corrected astrometry when available.'
    lim = 8*units.arcsec
    mk_xmatch(xrt_cat, sne_cat, lab, txt, lim)
    
    lab = 'grb_detected_by_xrt'
    txt = 'This is a list of GRBs matched with catalog detections (20 arcsec).\nUsing 2MASS corrected astrometry when available.'
    lim = 20*units.arcsec
    mk_xmatch(xrt_cat, grb_cat, lab, txt, lim)

def chandra():
    lab = 'sne_detected_by_cxo'
    txt = 'This is a list of SNe matched with catalog detections (1 arcsec; CSC 2.0 limit).'
    lim = 1*units.arcsec
    mk_xmatch(cxo_cat, sne_cat, lab, txt, lim)
    
    lab = 'grb_detected_by_cxo'
    txt = 'This is a list of GRBs matched with catalog detections (1 arcsec; CSC 2.0 limit).'
    lim = 1*units.arcsec
    mk_xmatch(cxo_cat, grb_cat, lab, txt, lim)


xmm()
nustar()
swift()
chandra()
