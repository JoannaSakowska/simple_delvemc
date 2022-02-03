#!/usr/bin/env python
"""
Perform simple binning search
"""
__author__ = "Joanna Sakowska"

import os
import glob
import subprocess
import yaml
import healpy as hp

import simple_adl.search

#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~ DELVE-MC DWARF SEARCH 2021  ~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#

"""
Python parallelisation for general clusters
The csub package is retired in favour of internal python packages
namely multiprocessing

Joanna Sakowska
"""
import pandas as pd
import multiprocessing
from multiprocessing import Pool

# Enter number of threads to use 

n_threads = 3

#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Sub-routines ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
# JDS: To do: if local cat search, need to add 'cfg' argument I think
#

# Function to tell us job status
# Use if desired

#def info(title):
#    print(title)
#    print('Parent process:', os.getppid())
#    print('Process id:', os.getpid())


def submit_job(ra, dec):
    
    # Output
    outfile = 'results_{:0.2f}_{:0.2f}.csv'.format(ra, dec)

    # Commands
    command = 'python {}/search.py --ra {:0.2f} --dec {:0.2f} --outfile {}'.format(os.path.dirname(simple_adl.search.__file__), ra, dec, outfile)

    # Log
    # Change log directory here
    log_dir = '/home/js01093/dwarf/simple_adl/simple_adl/log_dir'
    log_name = '{}/results_{:0.2f}_{:0.2f}.log'.format(log_dir, ra, dec)

    # Main function

    print(command)

    log = open(log_name, "w")

    search_output = (
        subprocess.run(
            command.split(' '), stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True
        )
        .stdout.decode("utf-8")
    )

    print(search_output)
    
    log.write(str(search_output))
    
    log.close()

    return

#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Main search ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--config',type=str,required=False,default='config.yaml',
                        help='config file [.yaml]')
    args = parser.parse_args()

    # Reading config
    with open(args.config, 'r') as ymlfile:
        cfg = yaml.load(ymlfile, Loader=yaml.SafeLoader)

        # JDS: TO FIX
        # Implement below structure to remove path dependencies for future version

        #survey = cfg['survey']
        #simple_dir = cfg['setup']['simple_dir']
        #jobs = cfg['batch']['jobs']
        #nside = cfg[survey]['nside']
        #datadir = cfg[survey]['datadir']

    # JDS: To search locally
    #data_dir = '/home/js01093/dwarf/simple_adl/simple_adl/data_dir'
    #infiles = glob.glob('{}/*.fits'.format(cfg['catalog']['data_dir']))

    # Reading in coordinates
    search_coordinates = pd.read_csv('{}.csv'.format('search_coordinates'), delimiter=',', header=None)

    ra_search_all, dec_search_all = [], []

    results_dir = '/home/js01093/dwarf/simple_adl/simple_adl/results_dir/'

    ra_search_all.extend(search_coordinates.iloc[:,0])
    dec_search_all.extend(search_coordinates.iloc[:,1])

    ra_search, dec_search = [], []

    for ra, dec in zip(ra_search_all, dec_search_all):
        if os.path.exists(os.path.join(results_dir, 'results_{:0.2f}_{:0.2f}.csv'.format(ra, dec))):
            print('EXISTS results_{:0.2f}_{:0.2f}.csv'.format(ra, dec))
        else:
            ra_search.append(ra)
            dec_search.append(dec)

    # Clean memory
    ra_search_all, dec_search_all = [], []

    print('Ready to search the Magellanic Clouds!')

    # Zipping arguments for command
    search_arguments = [*zip(ra_search,dec_search)]

    with Pool(n_threads) as p:
        p.starmap(submit_job, search_arguments)

        
        
