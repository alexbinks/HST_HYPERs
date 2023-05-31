import os
import glob
import numpy as np
from warnings import filterwarnings

from hst_modules.source_detection import run_source_detection
from hst_modules.download_data import move_files
from hst_modules.download_data import get_hst_data
#from hst_modules.process_images import run_raw
#from hst_modules.process_images import run_cln
#from hst_modules.process_images import plot_bkg
from hst_modules.process_images import *

from hst_modules.artificial_phot import run_artificial_phot
from hst_modules.artificial_phot import artificial_recovery

flc_files = sorted(glob.glob('./data/iehc03/iehc03j8q_flc.fits'))
nthresh_arr = np.array([1, 2, 3, 5, 10])
mag1, mag2, n_mag = 27., 28., 2
mag_art_arr = np.round(np.linspace(mag1, mag2, n_mag))

im_cut = [[2000, 2200], [800, 1000]]
#im_cut = None

for flc_file in flc_files:
    d_sci = plot_bkg(flc_file) ### PLOT THE ORIGINAL FLC FILE, THE BACKGROUND, AND THE BACKGROUND-SUBTRACTED IMAGE
    d_sci = make_bkg_sub_fits(flc_file, d_sci) ### MAKE FITS FILES FOR TILE 1 AND TILE 2
    update_flc(flc_file, d_sci)
    d_sci = im_stats(d_sci)
    d_sci = wcs_trans(d_sci)
    d_sci = med_mad(d_sci)    
    print(d_sci.keys())
    fwhm_in = 10.
#    print(os.path.basename(flc_file), ' running DAOPhot')
    run_source_detection(nthresh_arr, flc_file, d_sci, [1, 4], 
                         is_art=0, scramble=0, fwhm=fwhm_in, im_cut=im_cut)

#    print(os.path.basename(flc_file), ' scrambling the arrays to check for false positives')
    run_source_detection(nthresh_arr, flc_file, d_sci, [1, 4], 
                         is_art=0, scramble=1, fwhm=fwhm_in, im_cut=im_cut)

#make some artificial star lists.
    if d_sci['head']['FILTER'] in ['F625W', 'F775W']:
#        for f, f_file in enumerate(d_sci["flc_n"]):
        run_artificial_phot(flc_file, d_sci, mag_art_arr, nthresh_arr, fwhm=fwhm_in)
        artificial_recovery(d_sci, mag_art_arr, nthresh_arr, flc_file)

#else:
#    print(os.path.basename(flc_file), ' processing already complete: no need to get coordinates')
