#!/Users/silver/anaconda3/envs/astroconda/bin/python
'''
2019-09-17, Dennis Alp, dalp@kth.se

Get intra-observation XMM/EPIC X-ray light curves for a SIMBAD object
or provided coordinates.

Examples:
Either strict SIMBAD name:
util_xmm_lc.py 0854590801 1 10 100 0.3 10 "FRB 180916.J0158+65"
util_xmm_lc.py 0854590801 1 10 100   5 10 "FRB 180916.J0158+65"
util_xmm_lc.py 0804980201 1 10 100 0.3 10 "SN 1987A"
util_xmm_lc.py 0804980201 1 10 100   5 10 "SN 1987A"
util_xmm_lc.py 0650450301 1 10 100 0.3 10 "CXOU J232327.8+584842"
util_xmm_lc.py 0650450301 1 10 100   5 10 "CXOU J232327.8+584842"

Or any coordinate format understood by astropy:
util_xmm_lc.py 0107860301 1 10 100 0.3 10 255.303 64.1342

or:
05:36:10.123 -69:35:09.107
"23 23 53.316" "+58 40 12.06"

Pay attention to how Bash interprets spaces
'''

from __future__ import division, print_function
import os
from pdb import set_trace as db
import sys
import time
from glob import glob
import subprocess
import re

import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.wcs import WCS
from astropy.io import fits
from astropy import units
from astroquery.simbad import Simbad



################################################################
def get_src():
    try:
        print('Trying to interpret input as coordinates')
        ra = sys.argv[-2].replace(',', '')
        de = sys.argv[-1].replace(',', '')
        if re.search(r"[0-9].[0-9]", ra):
            print('Assuming degree angles')
            return SkyCoord(ra, de, unit=(u.deg, u.deg))
        else:
            print('Assuming hour angle')
            return SkyCoord(ra, de, unit=(u.hourangle, u.deg))
    except:
        print('Failed to interpret as coordinates, searching SIMBAD')
        try:
            res = Simbad.query_object(sys.argv[-1])
            ra = res['RA'][0]
            de = res['DEC'][0]
            return SkyCoord(ra, de, unit=(u.hourangle, u.deg))
        except:
            print('All interpretations failed.\nExiting\n')
            sys.exit()
        
class Xoi:
    def __init__(self):
        self.src_tt = np.empty(0)
        self.src_pi = np.empty(0)
        self.bkg_tt = np.empty(0)
        self.bkg_pi = np.empty(0)
        self.src_rad = 400
        self.bkg_rad = [800, 1200]
        self.backscal = self.src_rad**2/(self.bkg_rad[1]**2-self.bkg_rad[0]**2)

    def ext_evt(self, wcs):
        print('Extracting events')
        for ff in glob(cwd + obs + '/*FTZ'):
            self.cam = ff.split('/')[-1][11:13]
            xx, yy = np.array(wcs.all_world2pix(coo.ra, coo.dec, 1))
            self.coo = np.array([xx, yy])
            dd = fits.open(ff)[1].data

            if self.cam == 'PN':
                good = dd['FLAG'] <= XMMEA_EP
                good = good & (dd['PATTERN'] <= 4)
            elif self.cam == 'M1' or self.cam == 'M2':
                good = dd['FLAG'] <= XMMEA_EM
                good = good & (dd['PATTERN'] <= 12)

            good = good & (dd['PI'] > e_min) & (dd['PI'] < e_max)
            rad = (dd['X']-xx)**2+(dd['Y']-yy)**2

            src = rad < self.src_rad**2
            src = good & src
            bkg = (rad > self.bkg_rad[0]**2) & (rad < self.bkg_rad[1]**2)
            bkg = good & bkg

            self.src_tt = np.concatenate((self.src_tt, dd['TIME'][src]))
            self.src_pi = np.concatenate((self.src_pi, dd['PI'][src]/1e3))
            self.bkg_tt = np.concatenate((self.bkg_tt, dd['TIME'][bkg]))
            self.bkg_pi = np.concatenate((self.bkg_pi, dd['PI'][bkg]/1e3))

        if self.src_tt.size == 0:
            print('No source events\n')
            sys.exit(1)
        elif self.src_tt.size < 10:
            print('Only {0:d} events, skipping\n'.format(self.src_tt.size))
            sys.exit(1)

        order = np.argsort(self.src_tt)
        self.src_tt = self.src_tt[order]
        self.src_pi = self.src_pi[order]
        order = np.argsort(self.bkg_tt)
        self.bkg_tt = self.bkg_tt[order]
        self.bkg_pi = self.bkg_pi[order]

    def mk_lc(self, tb):
        self.tlc = np.arange(self.src_tt[0], self.src_tt[-1], tb)
        self.src_lc, _ = np.histogram(self.src_tt, self.tlc)
        self.src_lc = self.src_lc/np.diff(self.tlc)
        self.bkg_lc, _ = np.histogram(self.bkg_tt, self.tlc)
        self.bkg_lc = self.backscal*self.bkg_lc/np.diff(self.tlc)
        self.sub_lc = self.src_lc - self.bkg_lc

        tmp = np.empty(2*self.tlc.size-2)
        tmp[:-1:2] = self.tlc[:-1]
        tmp[1::2] = self.tlc[1:]
        shift = int(tmp[0]//1000*1000)
        tmp = tmp - shift

        fig = plt.figure(figsize=(5, 3.75))
        ax = plt.gca()
        plt.plot(tmp, self.src_lc.repeat(2), label='Source', color='greenyellow')
        plt.plot(tmp, self.bkg_lc.repeat(2), label='Background', color='gold')
        plt.plot(tmp, self.sub_lc.repeat(2), label='Subtracted', color='k')

        tmp = '{0:s}, {1:.0f} s, {2:.0f}-{3:.0f} keV'
        tmp = tmp.format(obs, tb, e_min/1.e3, e_max/1.e3)
        ax.set_title(tmp)

        ax.legend()
        ax.set_xlabel('Time-{0:d} (s)'.format(shift))
        ax.set_ylabel('Rate (cts/s)')
        ax.set_ylim(bottom=0)

        tmp = '{0:s}_{1:.0f}s_{2:.0f}-{3:.0f}kev.pdf'
        tmp = tmp.format(obs, tb, e_min/1.e3, e_max/1.e3)
        plt.savefig(tmp, bbox_inches='tight', pad_inches=0.1, dpi=300)
   
def get_wcs():
    ff = glob(cwd + obs + '/*FTZ')[0]
    dd = fits.open(ff)
    wcs = WCS()
    wcs.wcs.crpix = [dd[0].header['REFXCRPX'], dd[0].header['REFYCRPX']]
    wcs.wcs.cdelt = [dd[0].header['REFXCDLT'], dd[0].header['REFYCDLT']]
    wcs.wcs.crval = [dd[0].header['REFXCRVL'], dd[0].header['REFYCRVL']]
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    return wcs
    
def download_data():
    def dl_hlp(cmd):
        retcode = subprocess.call(cmd)
        gtar = glob('GUEST*.tar')
        gftz = glob('*.FTZ')
        file_missing = ' '.join(cmd).split('/')[-1]
        if len(gtar) == 1:
            retcode = subprocess.call(['tar', '-xvf', cwd + gtar[0], '-C', cwd], stderr=subprocess.DEVNULL)
            retcode = subprocess.call(['mv'] + glob(cwd + obs + '/pps/*') + [cwd + obs])
            retcode = subprocess.call(['rm', '-rf', cwd + obs + '/pps'])
            retcode = subprocess.call(['rm', cwd + gtar[0]])
            return True
        elif len(gftz) == 1:
            retcode = subprocess.call(['mkdir', '-p', cwd + obs])
            retcode = subprocess.call(['mv', cwd + gftz[0], cwd + obs])
            return True
        elif os.path.isfile(file_missing):
            retcode = subprocess.call(['rm', file_missing])
            return False

        print('ERROR Unknown data format delivered from XSA AIO:', obs)
        print(' '.join(cmd) + '\n')
        sys.exit(1)

    if os.path.isdir(cwd + obs): return True

    print('Downloading data')
    retcode = subprocess.call(['rm', '-rf', cwd + obs])
    retcode = subprocess.call(['rm'] + glob('GUEST*.tar'), stderr=subprocess.DEVNULL)
    retcode = subprocess.call(['rm'] + glob('*.FTZ'), stderr=subprocess.DEVNULL)

    tmp = dl_hlp(cmd + [purl.format(obs)])
    return dl_hlp(cmd + [murl.format(obs)]) or tmp




################################################################
cwd = '/Users/silver/Desktop/'
os.chdir(cwd)
cts_per_bin = 25 # For the light curves
cmd = ['curl', '-sS', '-O', '-J']
purl = 'http://nxsa.esac.esa.int/nxsa-sl/servlet/data-action-aio?level=PPS&extension=FTZ&name=PIEVLI&obsno={0:10s}'
murl = 'http://nxsa.esac.esa.int/nxsa-sl/servlet/data-action-aio?level=PPS&extension=FTZ&name=MIEVLI&obsno={0:10s}'
XMMEA_EP = 65584
XMMEA_EM = 65000


################################################################
obs = sys.argv[1]
tbin = [float(x) for x in sys.argv[2:5]]
e_min = float(sys.argv[5])*1.e3
e_max = float(sys.argv[6])*1.e3

tmp = '{0:s}, {1:.1f}-{2:.1f} keV'
tmp = tmp.format(obs, e_min/1.e3, e_max/1.e3)
print(tmp)
coo = get_src()
xoi = Xoi()

if not download_data():
    print('Download from XSA failed for {0:s}\n'.format(obs))
    sys.exit(1)

wcs = get_wcs()
xoi.ext_evt(wcs)
print('Making light curves')
for tb in tbin:
    xoi.mk_lc(tb)
# plt.show()
print('Done\n')
