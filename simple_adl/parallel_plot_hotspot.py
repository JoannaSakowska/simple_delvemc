#!/usr/bin/env python
"""
Python parallelisation to plot results from the candidate list 
in result_dir
Adapted for general clusters / personal machine
"""
__author__ = "Joanna Sakowska"

import yaml
import os
import glob
import subprocess

import simple_adl.search

import csv
import pandas as pd

import multiprocessing
from multiprocessing import Pool

#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~ DELVE-MC DWARF SEARCH 2021  ~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#

# Enter number of threads to use 

# Use $ lspcu to check how many threads per core
# e.g. Joanna's has 2 per core, so max 14 threads = 7 cores

n_threads = 14

#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Sub-routines ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#

def submit_plot_job(ra, dec, mod):
    
    # Output
    outfile = 'candidate_{:0.2f}_{:0.2f}.png'.format(ra, dec)
    
    # Commands
    command = 'python {}/plot_hotspot.py --ra {:0.2f} --dec {:0.2f} --mod {:0.2f} --outfile {}'.format(os.path.dirname(simple_adl.search.__file__), ra, dec, mod, outfile)
    print(command)
    print('Preparing plotting jobs...')
    subprocess.run(command.split(' '), shell=False)



    return

#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Main search ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#

if __name__ == '__main__':

    # Reading config
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--config',type=str,required=False,default='config.yaml',
                        help='config file [.yaml]')
    args = parser.parse_args()

    with open(args.config, 'r') as ymlfile:
        cfg = yaml.load(ymlfile, Loader=yaml.SafeLoader)


    # Reading in candidates
    # JDS: Temporarily setting result directory
    # JDS: Need to formalise this to use the config
    results_dir = '/home/js01093/dwarf/simple_adl/simple_adl/results_dir'
    plots_dir = '/home/js01093/dwarf/simple_adl/simple_adl/plot_dir/'


    # Reading in candidates from folder
    #candidates = pd.read_csv('{}/{}.csv'.format(results_dir, 'out'), delimiter=',', header=None)

    # To read in candidates after make_list.py
    #candidates = pd.read_csv('{}.csv'.format('candidate_list'), delimiter=',', header=None)

    # To read in candidates AFTER 5 sigma cut
    candidates = pd.read_csv('{}.csv'.format('candidate_list_5sig'), delimiter=',', header=None)

    ra_candidate_all, dec_candidate_all, mod_candidate_all = [], [], []

    ra_candidate_all.extend(candidates.iloc[:,1])
    dec_candidate_all.extend(candidates.iloc[:,2])
    mod_candidate_all.extend(candidates.iloc[:,3])

    ra_candidate, dec_candidate, mod_candidate = [], [], []

    # JDS: Method to skip overwriting plots in plot_dir
    for ra, dec, mod in zip(ra_candidate_all, dec_candidate_all, mod_candidate_all):
        if os.path.exists(os.path.join(plots_dir, 'candidate_{:0.2f}_{:0.2f}.png'.format(ra, dec))):
            print('EXISTS candidate_{:0.2f}_{:0.2f}.png'.format(ra, dec))
        else:
            ra_candidate.append(ra)
            dec_candidate.append(dec)
            mod_candidate.append(mod)

    # Clean memory
    ra_candidate_all, dec_candidate_all, mod_candidate_all = [], [], []

    print('Ready to plot candidates!')

    # Zipping arguments to feed into command
    plot_arguments = [*zip(ra_candidate,dec_candidate,mod_candidate)]

    with Pool(n_threads) as p:
        p.starmap(submit_plot_job, plot_arguments)
        

