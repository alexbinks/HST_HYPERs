from astropy.io import fits, ascii
from astropy.table import Table
import pandas as pd
import matplotlib.pyplot as plt

import numpy as np
import os
import subprocess
import glob

from astropy.visualization import ImageNormalize, ZScaleInterval

from .source_detection import run_source_detection
from .source_detection import crossmatch
from .source_detection import make_histo_crossmatch
from .source_detection import make_mag_table_combined
from .image_functions import make_image

def make_artificial_stars(f_in, mag_art_arr=23.+0.5*np.arange(5), n_stars=200):
    '''Create a text file containing the X, Y and magnitudes for the artificial
       stars, which will be used in the Anderson+22 processing in Fortran.
       The function generates "Nflux" flux points, sampled evenly in logarithmic
       space, where each flux value is placed at "Nsamp" random X, Y position.
       Fluxes are converted to instrumental magnitude using the exposure time
       for each image, which are the necessary units for the Fortran code.

       For each of the two tiles in the images a text file is
       created called "art_stars#.XYM", where # refers to either
       tile "1" or "2".       
    '''
    f = fits.open(f_in)
    for i, ext in enumerate([1, 4]):
        head = f[ext].header
        STMAG = -2.5*np.log10(head['PHOTFLAM'])+head['PHOTZPT']
        for mag in mag_art_arr:
            mag_inst = mag-STMAG
            x = head["NAXIS1"]*np.random.uniform(low=0.0, high=1.0, size=n_stars)
            y = head["NAXIS2"]*np.random.uniform(low=0.0, high=1.0, size=n_stars)
            m_arr = np.repeat(mag_inst, n_stars)
            df = pd.DataFrame([x, y, m_arr])
            df.T.to_csv(f'art_{mag}_{i+1}.XYM', sep=' ', index=False, header=False)




    
def add_artificial_stars(d_sci):
    '''Here is where we actually run HST1PASS to generate the
       fits files containing the artificial stars. They are
       then moved to their respective subdirectories.
    '''
    dict_PSFs = {'F502N': 'F467M',
                 'F555W': 'F555W',
                 'F625W': 'F621M',
                 'F656N': 'F621M',
                 'F775W': 'F775W'}

    loc_files = os.getcwd()+"/andersonPSF/hst1pass/"

    filt_name = d_sci['head']['FILTER']
    filt_art = dict_PSFs[filt_name]
    image_dir = f'{d_sci["flc_dir"]}'
    image_base = f'{d_sci["flc_base"]}'
    image_name = f'{d_sci["flc_name"]}'

    for i in [0, 1]:
        image_path_full = f'{d_sci["flc_n"][i]}'

#image_path_full: ./data/iehc03/iehc03j8q_F625W_1_flc.fits
#image_dir: ./data/iehc03/
#image_base: ./data/iehc03/iehc03j8q
#image_name: iehc03j8q

        art_files = sorted(glob.glob(f'art*{i+1}.XYM'))
        fin_files = []
        for art_file in art_files:
            mag = art_file.split("_")[-2]
            fin_files.append(f'{image_dir}{d_sci["targ_and_filt"]}_{mag}_art_{i+1}.txt')

        for art_file, fin_file in zip(art_files, fin_files):
            f = open(fin_file, "w")
            x = subprocess.run([f'{loc_files}./hst1pass.e',
                            f'HMIN=5',
                            f'FMIN=100',
                            f'OUT=xymXYMUVrd',
                            f'ART_XYM={art_file}',
                            f'PSF={loc_files}STDPSF_WFC3UV_{filt_art}.fits',
                            f'REG=xy',
                            f'SHOW_ART+',
                            f'DOSATD-',
                            f'{image_path_full}'], stdout=f)
            mag_in = art_file.split("_")[-2]
            fx_s = glob.glob(f'{image_name}*.fits')
            for fx in fx_s:
                if filt_name in fx:
                    # I PROMISE TO TIDY THIS AT SOME STAGE TO MAKE IT MORE COHERENT!
                    # IT MOVES THE FILES STORED IN THE BASE DIRECTORY TO THE SUBDIRECTORY
                    # CONTAINING THE FITS FILES AND DATA FOR A SPECIFIC OBSERVATION.
                    ch_ext = f'mv {fx} {image_dir}{fx[:8]}q{fx[9:-5]}_{mag_in}{fx[-5:]}'
                else:
                    ch_ext = f'mv {fx} {image_dir}{fx[:8]}q_{filt_name}_{i+1}{fx[9:-5]}_{mag_in}{fx[-5:]}'
                print(ch_ext)
                subprocess.run(ch_ext, shell=True)    
            fx_s = glob.glob(f'{image_name}*')
            for fx in fx_s:
                ch_ext = 'rm ' + fx
                subprocess.run(ch_ext, shell=True)
        ch_ext = f'mv art_*_{i+1}.XYM {image_dir}'
        subprocess.run(ch_ext, shell=True)





def add_artifical_stars_to_image(im_file, d_sci, i, art_files):
    '''Add the artificial photometry to the data and save the fits file.
    '''
    f = fits.open(im_file)
    h0 = f[0].header
    df, hf = f[1].data, f[1].header
    for art_file in art_files:
        with fits.open(art_file) as a:
            da = a[0].data[0:hf['NAXIS2'],0:hf['NAXIS1']]
            da = da*h0['EXPTIME']
            d_tot = da+df
            x = a
            x[0].data = d_tot
            new_name = art_file.replace("add","com")
            x.writeto(new_name, output_verify='ignore', overwrite=True)
            make_image(h0, x[0].data, f'{i+1}', d_sci['flc_name'], d_sci['flc_dir'], ext='art')





def combine_art_fits(c1_file, c2_file):
    '''Merge the data+artificial images into one fits file.
    
    This makes the cross-matching function a lot easier later on.
    '''
    with fits.open(c1_file) as hdu1:
        c0_head = hdu1[0].header
        c1_data = hdu1[0].data
    with fits.open(c2_file) as hdu2:
        c2_data = hdu2[0].data

    hdu1 = fits.PrimaryHDU(header=c0_head)
    hdu2 = fits.ImageHDU(c1_data)
    hdu3 = fits.ImageHDU(c2_data)
    new_hdul = fits.HDUList([hdu1, hdu2, hdu3])
    
    i_n1 = "_".join(c1_file.split("_")[:2])
    i_n2 = "_".join(c1_file.split("_")[3:])
    i_n = f'{i_n1}_{i_n2}'
    new_hdul.writeto(i_n)
    subprocess.run(f'rm {c1_file} {c2_file}', shell=True)


        
        
        
        
def run_artificial_phot(flc_file, d_sci, mag_art_arr, nthresh_arr, fwhm=3.):
    '''Make the artificial photometry, add them to the data and run a crossmatch.
    '''
    make_artificial_stars(flc_file, mag_art_arr=mag_art_arr)    
    add_artificial_stars(d_sci)
    
    for i in [0, 1]:
        astars = sorted(glob.glob(f'{d_sci["flc_dir"]}*_{i+1}_add_*.fits'))
        add_artifical_stars_to_image(d_sci['flc_n'][i], d_sci, i, astars)
                    
    cstars1 = sorted(glob.glob(f'{d_sci["flc_dir"]}*_1_com_*.fits'))
    cstars2 = sorted(glob.glob(f'{d_sci["flc_dir"]}*_2_com_*.fits'))

    for cstar1, cstar2 in zip(cstars1, cstars2):
        combine_art_fits(cstar1, cstar2)
    
    cstars = sorted(glob.glob(f'{d_sci["flc_dir"]}*_com_*.fits'))
    for cstar in cstars:
        run_source_detection(nthresh_arr, cstar, d_sci, [1, 2], 
                             is_art=1, scramble=0, fwhm=fwhm)





def artificial_recovery(d_sci, mag_art_arr, nthresh_arr, flc_file):
    '''Determine how many artificial stars were recovered for a given magnitude/threshold
    '''
    table_name = f'{d_sci["flc_dir"]}{d_sci["targ_and_filt"]}_magnitude_recovery.csv'
    t = Table(names=['mag', 'extn', 't_lim', 'Nrecov', 'frac'])
    for mag in mag_art_arr:
        for i in [0, 1]:
            sf_t = ascii.read(f'{d_sci["flc_dir"]}{d_sci["targ_and_filt"]}_{mag}_s{i+1}.dat')
            af_t = ascii.read(f'{d_sci["flc_dir"]}art_{mag}_{i+1}.XYM', names=['X','Y','m_inst'])
            for nth in nthresh_arr:
                sn_t = sf_t
                s_min = (d_sci['med_sc'][i] + nth*d_sci['mad_sc'][i])/ \
                        (d_sci['med_sc'][i] + nthresh_arr[0]*d_sci['mad_sc'][i])
                sn_t = sf_t[sf_t["flux"]> s_min]

                sf_xy = np.array([sn_t['xcentroid'], sn_t['ycentroid']]).T
                af_xy = np.array([af_t['X']-1, af_t['Y']-1]).T
                x1, x2 = crossmatch(af_xy, sf_xy, 3.0)
                xn = np.array(x1)#.astype(int)
                make_histo_crossmatch(mag, f'{i+1}', sn_t, xn, d_sci['head']['FILTER'],
                                      d_sci['head_data'][0], d_sci['head_data'][1],
                                      nth, flc_file)
                t.add_row([mag, f'{i+1}', nth, len(x1), (1.0*len(x1))/(1.0*len(af_xy))])
    t.write(table_name, overwrite=True)
    x_t = make_mag_table_combined(table_name, d_sci)
