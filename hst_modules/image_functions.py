from PIL import Image
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy.visualization import ImageNormalize, ZScaleInterval
from mpl_toolkits.axes_grid1 import make_axes_locatable

def make_image(fh, fd, num, name, dty, sources=False, s_file=None, ext=''):
    '''Create 2D histograms using image data from fits files.
    '''
    fig, ax = plt.subplots(figsize=(20, 20), dpi=200)
#    file_type = fh['FILENAME'].split('_')[-1][:-5]
    n=ImageNormalize((fd), interval=ZScaleInterval())

    kwargs={'font.size': 40}
    mpl.rcParams.update(kwargs)
    
    ax.set_xlabel('X pixel')
    ax.set_ylabel('Y pixel')

    im = ax.imshow(fd, origin='lower', norm=n,
                   cmap='Blues_r', interpolation='nearest')
    
    ax.text(0, 0, fh['FILTER'], fontsize=40, color='yellow')
    if sources:
        ax.scatter(s_file['xcentroid'], s_file['ycentroid'], facecolors='none', edgecolors='r', s=30)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(im, cax=cax)
    cbar.set_label('counts ($e^{-}$/s)', fontsize=25)
    strlist = [name, num, fh['FILTER'], ext]
    plt.savefig(f'{dty}{"_".join(strlist)}.png', bbox_inches='tight', facecolor='white')
    plt.clf()
    mpl.rcParams.update(mpl.rcParamsDefault)
