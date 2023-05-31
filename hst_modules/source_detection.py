import os
import numpy as np
from astropy.table import Table, unique
from astropy.io import fits, ascii
from datetime import datetime
from photutils.detection import DAOStarFinder
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
from .image_functions import make_image






def scramble_array(f, df, ext):
    rng = np.random.default_rng()
    sub_dim = f[ext].data.shape
    flat = f[ext].data.ravel()
    rng.shuffle(flat)
    scrambled_image = flat.reshape(sub_dim) 
    s = df(scrambled_image)
    return scrambled_image, s




def run_source_detection(nthresh_arr, f_in, d_sci, ext_names, 
                         is_art=0, scramble=0, fwhm=3., im_cut=None):
    '''A function to run the source detection algorithm.
    '''
# Flux values from DAOStarFinder are calculated as the peak density in the convolved image divided by the detection threshold.
# So the "lowest" accepted source in nthresh_arr[0] "just" passes the threshold and has flux ~1.000001. This means we can simply
# multiply the flux value by (med + nth*MAD)/(med + nthresh_arr[0]*MAD) to match the other thresholds.

    start = datetime.now()
    print(d_sci['flc_base'], d_sci['flc_name'], d_sci['head']['FILTER'])
    t_res = Table(names=['ObsID', 'Filter', 'ext', 'Nthresh', 'median', 'MAD', 'Nmatch'],
                  dtype=[str, str, int, int, np.float32, np.float32, int])

    for i, ext in enumerate(ext_names):
        f = fits.open(f_in)
        threshold = d_sci['med_sc'][i] + nthresh_arr[0]*d_sci['mad_sc'][i]
        f[ext].data = f[ext].data - d_sci['med_sc'][i]
        df = DAOStarFinder(threshold=threshold, fwhm=fwhm)
        if im_cut:
            f[ext].data = f[ext].data[im_cut[1][0]:im_cut[1][1]+1, im_cut[0][0]:im_cut[0][1]+1]

        if not is_art:
            if scramble:
                print("SCRAMBLING!")
                t_match = f'{d_sci["flc_dir"]}{d_sci["targ_and_filt"]}_s{i+1}_scr.dat'
                t_res_name = f'{d_sci["flc_dir"]}{d_sci["targ_and_filt"]}_scr.csv'
                scrambled_image, s = scramble_array(f, df, ext)
                make_image(f[0].header, scrambled_image, f'{i+1}', d_sci['flc_name'], d_sci['flc_dir'], sources=True, s_file=s, ext='scr')
            else:
                t_match = f'{d_sci["flc_dir"]}{d_sci["targ_and_filt"]}_s{i+1}_sou.dat'
                t_res_name = f'{d_sci["flc_dir"]}{d_sci["targ_and_filt"]}_sou.csv'
                s = df(f[ext].data)
        else:
            mag_bit = f_in.split("_")[3][:4]
            t_match = f'{d_sci["flc_dir"]}{d_sci["targ_and_filt"]}_{mag_bit}_s{i+1}.dat'
            t_res_name = f'{d_sci["flc_dir"]}{d_sci["targ_and_filt"]}_{mag_bit}_art.csv'
            s = df(f[ext].data)
        
        for col in s.colnames:
            s[col].info.format = '%.8g' # for consistent table output
            s.sort(['flux'], reverse=True)
        print(f'extn{i+1}: thresh0={threshold}, Nsources={len(s)}')

        if im_cut:
            s['xcentroid'] += im_cut[1][0]
            s['ycentroid'] += im_cut[0][0]
            
        radec = SkyCoord.from_pixel(s['xcentroid'], s['ycentroid'], d_sci['mywcs'][i], 0)
        s['RA'], s['DEC'] = radec.ra.deg, radec.dec.deg

        ascii.write(s, t_match, overwrite=True)

        for nth in nthresh_arr:            
            sn = s
            s_min = (d_sci['med_sc'][i] + nth*d_sci['mad_sc'][i])/ \
                    (d_sci['med_sc'][i] + nthresh_arr[0]*d_sci['mad_sc'][i])
            sn = s[s["flux"]> s_min]

            t_res.add_row([d_sci['flc_name'], d_sci['head']['FILTER'], i+1, nth, d_sci['med_sc'][i], d_sci['mad_sc'][i], len(sn)])

        finish = datetime.now()
        print(f'completed {d_sci["flc_name"]}, [ext{ext}] in {(finish-start).total_seconds()} secs')

    t_res.write(t_res_name, overwrite=True)
    print(f'total time = {(finish-start).total_seconds()} secs')
        
        
        
        
        
        


def crossmatch(pos1, pos2, max_d):
    '''Cross-matching function copied from Tara Murphy's data-driven astronomy course

    https://www.coursera.org/learn/data-driven-astronomy
    '''
    
    matches = []
    no_matches = []
    # sort the second catalog in order of "ycentroid"
    asc_y = np.argsort(pos2[:,1])
    pos2_sorted = pos2[asc_y]
    # make a column called y2_sorted with "ycentroid" in asc order.
    y2_sorted = pos2_sorted[:,1]
    for id1, (x1, y1) in enumerate(pos1):
        closest_dist = np.inf
        closest_id2 = None
        min_y = y1 - max_d
        max_y = y1 + max_d
        # Start and end indices of the search
        start = y2_sorted.searchsorted(min_y, side='left')
        end = y2_sorted.searchsorted(max_y, side='right')
        N_match = 0
        for s_id2, (x2, y2) in enumerate(pos2_sorted[start:end+1]):
            dist = np.sqrt((x1-x2)**2 + (y1-y2)**2)
            if dist < closest_dist:
                N_match += 1
                closest_dist = dist
                closest_id2 = start+s_id2
        
        # Ignore match if it's outside the maximum radius
        if closest_dist > max_d:
            no_matches.append(id1)
        else:
            closest_id2 = asc_y[closest_id2]
            matches.append([id1, closest_id2, closest_dist, N_match])
    return np.array(matches), np.array(no_matches)
    
        
    
    

def make_histo_crossmatch(mag, ext, sf_t, xn, filt_name,
                          head_data1, head_data2,
                          nth, flc_file):
    if (xn.shape[0]) != 0:                        
        x_list = xn[:,1].astype(int)
        fig, ax = plt.subplots(figsize=(10,5), nrows=1, ncols=2)
        ax[0].set_xlabel("$\Delta$mag", fontsize=14)
        ax[1].set_xlabel("$\Delta$pixel", fontsize=14)
        if ext == 1:
            ax[0].hist(-2.5*np.log10(head_data1['PHOTFLAM']) +
                       head_data1['PHOTZPT'] + 
                       sf_t["mag"][x_list] - mag)
            ax[1].hist(xn[:,2])
        else:
            ax[0].hist(-2.5*np.log10(head_data2['PHOTFLAM']) +
                       head_data2['PHOTZPT'] + 
                       sf_t["mag"][x_list] - mag)
            ax[1].hist(xn[:,2])
        fig.savefig(f'{os.path.dirname(flc_file)}/mag_histo_{filt_name}_{mag}_{ext}_nt_{nth:03d}.png', bbox_inches='tight', facecolor='white')
        fig.clf()





def make_mag_table_combined(file_in, d_sci):
    '''Make a new table which mean combines the number of stars and
    the fraction recovered from tiles 1 and 2. Construct a plot of
    recovery fraction vs magnitude for each detection threshold.
    '''
    filt_name = d_sci['head']['FILTER']
    f_in = Table.read(file_in)
    print(file_in)
    u_vals = unique(f_in, keys=['mag', 't_lim'])
    t_mean = Table(names=['mag', 't_lim', 'Nrecov', 'frac'],
                   dtype=[np.float64, np.float64, np.float64, np.float64])
    for r,rows in enumerate(u_vals):
        g = f_in[(f_in['mag'] == rows['mag']) & (f_in['t_lim'] == rows['t_lim'])]
        t_mean.add_row([g["mag"][0], g["t_lim"][0], np.mean(g["Nrecov"]), np.mean(g["frac"])])
    obs_by_name = t_mean.group_by("t_lim")
    x_marks = np.unique(t_mean["t_lim"])

    cmap = plt.cm.get_cmap('copper', len(x_marks))
    fig, ax = plt.subplots(figsize=(8,6))
    ax.grid()
    ax.set_xlabel('magnitude', fontsize=20)
    ax.set_ylabel('recovery fraction', fontsize=20)
    for i, (key, group) in enumerate(zip(obs_by_name.groups.keys, obs_by_name.groups)):
        ax = plt.plot(group["mag"], group["frac"], '--', lw=3, c=cmap(i+1), label=group['t_lim'][0])
    plt.legend(fontsize=15)
    plot_name = f'{d_sci["flc_dir"]}{d_sci["targ_and_filt"]}_magnitude_recovery.png'
    plt.savefig(plot_name, bbox_inches='tight', facecolor='white')
    plt.clf()
    return t_mean
