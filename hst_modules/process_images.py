import os
import glob
import subprocess

from astropy.stats import sigma_clipped_stats, sigma_clip, SigmaClip
from astropy.io import fits
from astropy.wcs import WCS

from astroscrappy import detect_cosmics

from photutils.background import Background2D
from photutils.background import MedianBackground

from scipy.stats import median_abs_deviation


from wfc3tools import calwf3
from acstools.calacs import calacs

import numpy as np
from collections import defaultdict

from .image_functions import make_image


def run_raw(raw_files, bkg_sigma=3., bkg_filt=3, bkg_kern=50, cr_clip=1., cr_objlim=30.):
    '''Perform the background subtraction and cosmic ray correction for each raw fits file.
    
    Save the file with suffix *_crz_raw.fits.
    '''
    bkg_clip = SigmaClip(sigma=bkg_sigma)
    rf_dir = os.path.dirname(raw_files[0])+'/'

    if len(glob.glob(f"{rf_dir}*_crz_raw.fits")) == 0:
        for raw_file in raw_files:
            print(f'processing {raw_file}')
            # Define the naming convention to store files
            rf_base = raw_file[:-9]
            rf_name = os.path.basename(raw_file)[:-9]

            # Run the background subtraction and cosmic ray corrections.
            try:
                with fits.open(raw_file, mode='readonly') as hdu:
                    head = hdu[0].header
                    
                    # The header in the fits file that contains the filter type (e.g., F625W, F775W, etc)
                    # is titled "FILTER1" in the ACS (prefix "j") files. Make a fits header titled "FILTER"
                    # to make processing easier later on.
                    if rf_name[0] == 'j':
                        head.set('FILTER', head['FILTER1'], 'filter name')

                    for i, ext in enumerate([1, 4]):
                        data = hdu[ext].data

                # Plot the raw fits files (before any processing steps)
                        make_image(head, data, f'{i+1}', rf_name, rf_dir, ext='raw')
                # Perform the background subtraction.
                        bkg = Background2D(data, bkg_kern, filter_size=bkg_filt, sigma_clip=bkg_clip, 
                                        bkg_estimator=MedianBackground())
                        bkg_corr = data - bkg.background

                # Do the cosmic ray zapping.
                        crm, cli = detect_cosmics(bkg_corr, sigclip=cr_clip, objlim=cr_objlim)

                # Save the background-subtracted, cosmic ray zapped images.
                        hdu[ext].data = cli

                    hdu.writeto(f'{rf_base}_crz_raw.fits')
                    hdu.close()
                print(f'{raw_file} done!')

            except:
                raise Exception(f'Unable to process {raw_file}')
    else:
        print(f'All initial processing has been done for raw files in the {rf_dir} directory.')
        
        
        
        
def run_cln(crz_files):
    '''Run the hst reduction pipeline to acquire flat-fielded, CTE-corrected images.
    '''
    for crz_file in crz_files:
        if len(glob.glob(f'{os.path.dirname(crz_file)}/*flc_crz.fits')) == 0:
            subprocess.check_output(f'crds bestrefs --update-bestrefs --sync-references=1 --files {crz_file}', \
                                    shell=True,\
                                    stderr=subprocess.DEVNULL)
            if os.path.basename(crz_file)[0] == 'j':
                calacs(crz_file)
            if os.path.basename(crz_file)[0] == 'i':
                calwf3(crz_file)
            base = f'{"_".join(crz_file.split("_")[:-2])}'
            subprocess.run(f'mv {base}_crz_flc.fits {base}_flc.fits', shell=True)
            subprocess.run(f'mv {base}_crz_flt.fits {base}_flt.fits', shell=True)
        else:
            print(f'{os.path.basename(crz_file)} processing already complete: no need to run calwfc')
            
            
            
            
def make_fits(f_in, images, im_names):
    '''Create a multi-extension fits file from a set of input images.
    '''
    fits_dict = {}
    fits_dict['hdu0'] = fits.PrimaryHDU()
    f = fits.open(f_in)

    for i, image in enumerate(images):
        fits_dict[f'hdu{i+1}'] = fits.ImageHDU(image)
    hlist = ['hdu'+str(i) for i in range(len(images)+1)]
    new_hdul = fits.HDUList([d_sci[h] for h in hlist])
    new_hdul[0].header = f[0].header
    for i, (image, im_name) in enumerate(zip(images, im_names)):
        new_hdul[i+1].header = f[i+1].header
        new_hdul[i+1].data = image
        new_hdul[i+1].name = im_name
    return new_hdul




def make_bkg_sub_fits(flc_file, d_sci):
    '''Save the background and background-subtracted images as fits files.
    
    This is important when dealing with the artificial photometry, which
    requires the fits files as single-tile images.
    '''

    # THE ARTIFICIAL PHOTOMETRY ONLY ALLOWS FILES THAT HAVE HST NOMENCLATURE
    # THEREFORE WE ARE RETURNING THE "SUB" FILE AS "FLC".
    d_sci['flc_dir'] = os.path.dirname(flc_file)+'/'
    d_sci['flc_base'] = flc_file[:-9]
    d_sci['flc_name'] = d_sci['flc_base'][-9:]
    d_sci['targ_and_filt'] = '_'.join([d_sci['flc_name'], d_sci['head']['FILTER']  ])
#flc_dir: ./data/iehc03/
#flc_base: ./data/iehc03/iehc03j8q
#flc_name: iehc03j8q
    for i in [0,1]:
        for ext, new in zip(['bkg', 'sub'], ['bkg', 'flc']):
            image_fits = make_fits(flc_file, [d_sci[ext][i], d_sci['err'][i], d_sci['dqa'][i]],
                                ['SCI', 'ERR', 'DQ'])
            i_n = f'{d_sci["flc_dir"]}{d_sci["targ_and_filt"]}_{i+1}_{new}.fits'
            d_sci[f'{new}_n'].append(i_n)
            image_fits.writeto(i_n)
    return d_sci



def update_flc(flc_file, d_sci):
    '''Do the background correction to the *flc.fits files and save.
    '''
    with fits.open(flc_file, mode='update') as hdul:
        head = hdul[0].header
        # Only do the background correction once!! The if statement ensures this.
        if "BKG_CORR" not in list(head.keys()):
            hdul[1].data = hdul[1].data - d_sci['bkg'][0]
            hdul[4].data = hdul[4].data - d_sci['bkg'][1]
            head.set('BKG_CORR', 1, 'is the image background-corrected?')
            hdul.flush()
            print('background has just been corrected - good job!')
        else:
            print('background has already been corrected - good job!')
    
    


def plot_flc_bkg_sub(flc_file, bkg_sigma=3., bkg_filt=3, bkg_kern=50):
    '''Make plots of the *flc.fits, *bkg.fits and *sub.fits images
    '''
    bkg_clip = SigmaClip(sigma=bkg_sigma)
    flc_base = flc_file[:-9]
    flc_name = os.path.basename(flc_file)[:-9]
    flc_dir = os.path.dirname(flc_file)+'/'

    with fits.open(flc_file) as hdu:
        d_sci = defaultdict(list)
        head = hdu[0].header
        if os.path.basename(flc_file)[0] == 'j':
            head.set('FILTER', head['FILTER1'], 'filter name')
        d_sci['head'] = head
        for i, ext in enumerate([1, 4]):
            d_sci['head_data'].append(hdu[ext].header)
            d_sci['flc'].append(hdu[ext].data)
            d_sci['err'].append(hdu[ext+1].data)
            d_sci['dqa'].append(hdu[ext+2].data)

            backg = Background2D(d_sci['flc'][i], bkg_kern, filter_size=bkg_filt, sigma_clip=bkg_clip,
                                  bkg_estimator=MedianBackground())
            d_sci['bkg'].append(backg.background)
            d_sci['sub'].append((d_sci['flc'][i]-backg.background).astype('float32'))
            for ext in ['flc', 'bkg', 'sub']:
                make_image(head, d_sci[ext][i], f'{i+1}', flc_name, flc_dir, ext=ext)
        return d_sci
        
        
        
        
        
def im_stats(d_sci, sig_val=3.0):
    '''Get image statistics for the background and background-subtracted image.
    
    Returns the mean, median and standard deviation. Probably best to use the median
    from the background image as a threshold for source detection.
    
    I PREFER TO MAKE MY OWN SIGMA CLIPPING AND USE THE MEDIAN AND MEDIAN-ABSOLUTE-DEVIATION
    '''
    for i in [0,1]:
        sc = sigma_clip(d_sci['sub'][i], sigma=sig_val, masked=False, axis=None)

        d_sci['stat_arr'].append(sigma_clipped_stats(d_sci['bkg'][i], sigma=sig_val))
        d_sci['subt_arr'].append(sigma_clipped_stats(d_sci['sub'][i], sigma=sig_val))
        d_sci['med_sc'].append(np.median(sc))
        d_sci['mad_sc'].append(median_abs_deviation(sc, scale='normal'))

    return d_sci
    
    
    
    
    
def wcs_trans(d_sci):
    '''Read in the WCS data for the pixel transformation (RA, DEC -> x, y)
    '''
    for i in [0,1]:
        d_sci['mywcs'].append(WCS(d_sci['head_data'][i]))
    return d_sci
    
