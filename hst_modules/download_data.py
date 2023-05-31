import subprocess
import os
import glob
from astroquery.mast import Observations
from collections.abc import Iterable

def move_files(dl_files):
    '''Move the downloaded raw fits files to the correct directory.
    
    All raw fits files, when downloaded are stored in a directory titled "mastDownload"
    This function move the files from there to the appropriate directory for storing
    HST data.
    '''
    for dl_file in dl_files:
        dir_name = dl_file.split('/')[-2][:-3]
        if os.path.exists(f'./data/{dir_name}') == False:
            subprocess.run(f'mkdir data/{dir_name}', shell=True)
        subprocess.run(f'mv {dl_file} data/{dir_name}', shell=True)
    subprocess.run(f'rm -r mastDownload', shell=True)        
    
    
    
    
def get_hst_data(obs_ids):
    '''Download the raw fits files and store them in an appropriate directory.
    
    The MAST archive assigns a naming convention to files such that observations
    of the same field of view have the same 6 character prefix, and the following
    3 characters define the unique filter used in the observation.
    
    This function reads in either a single string, or an iterable of strings, each
    6 characters long, and downloads all raw fits files from the MAST archive whose
    prefixes match the string.
    '''

    # Make the data directory if it doesn't already exist.
    isExist = os.path.exists('./data')
    if not isExist:
        os.makedirs('./data')

    # Data can be either:
    # An iterable of strings.
    if not isinstance(obs_ids, str) and isinstance(obs_ids, Iterable):
        obs_ids = [obs_id+'*' for obs_id in obs_ids]
    # A single string
    elif isinstance(obs_ids, str):
        obs_ids = [obs_ids+'*']
    #...but nothing else...
    else:
        raise Exception('Input is neither a single string, or an iterable of strings')

    for obs_id in obs_ids:
        obs_table = Observations.query_criteria(dataproduct_type=["*"],
                                        proposal_id=["15888","16359"],
                                        obs_id=obs_id,
                                        calib_level=[1,2,3],
                                       )
        data = Observations.get_product_list(obs_table)
        raw_files = data[data['productSubGroupDescription']==['RAW']]
        raw_names = raw_files['obs_id']
        mv_file = False
        for raw_file, raw_name in zip(raw_files, raw_names):
            raw_name_local = f'./data/{raw_name[:-3]}/{raw_name}_raw.fits'
            if not os.path.exists(raw_name_local):
                manifest = Observations.download_products(raw_file, productSubGroupDescription=["RAW"])
                mv_file = True
            else:
                print(f'{raw_name_local} already exists!')
        if mv_file:
            dl_files = sorted(glob.glob(f'./mastDownload/**/*.fits', recursive=True))
            move_files(dl_files)
