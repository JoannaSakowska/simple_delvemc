# simple_adl
An adaptation of the simple binner dwarf galaxy search algorithm for use on generic clusters for the MC 2021 dwarf search. 

This version follows [`simple_adl`](https://github.com/sidneymau/simple_adl) which uses data pulled directly from the Astro Data Lab. Thanks Sidney and Will!

Simple_adl originates from [`simple`](https://github.com/DarkEnergySurvey/simple).

## Installation

Clone this repository to the machine you wish to run the code from.

Run the following command (consider adding to your `.bashrc`/`.zshrc`/etc.):
```
export PYTHONPATH=<path to simple_adl>:$PYTHONPATH
```
where `<path to simple_adl>` is the path to the inner simple_adl repository (ie, the full string should end with `*/simple_adl:$PYTHONPATH`).

Create the conda environment with
```
conda env create -f conda.yaml
```
This will create a conda environment called `simple_adl` that includes all of the necessary software.

## Use

Activate the conda environment through `conda activate simple_adl`. Note that this needs to be done from a directory that does not have a `simple_adl` folder.

Run `init_simple.py` to generate a default `config.yaml` configuration file.
This can be edited as needed; in particular, you may want to change the `profile` entry.

## New procedure for parallel searching

Update the search_coordinates.csv file with the coordinates to be searched. 

Open parallel_search.py. Enter the number of threads to run the search on in the `n_threads` parameter. On linux, use:
```
 $ lspcu 
```
to check how many threads per core your machine has (e.g., if there are 2 per core then a 7 core computer uses 14 threads). 

NOTE: as the algorithm is quite cpu intense, it is recommended to try a few test runs with less threads (monitoring your cpu status) and find
the best compromise.

Adapt the parallel_search.py filepaths to your machine.

Run on the command line:
```
python parallel_search.py 
```
Monitor the log_dir directory for output log information. 

NOTE: if any of the log files are empty, the search failed (most likely due to an interuption of the program). Run it again.

Once finished, run:
```
 python make_list.py
```
Inspect your candidate_list.csv.

Open parallel_plot_hotspot.py. You can either read in candidate_list.csv or, after
applying any cuts (such as 5 sigma), read in the new candidate_list file with a 
modified name (such as candidate_list_5sigma.csv) and keep your old one for reference. 

Set your number of threads and paths, run:
```
python parallel_plot_hotspot.py
```

And voila! Happy searching! 

## Example

To search for [DELVE 2](https://arxiv.org/abs/2009.08550), run the following:
```
python search.py --ra 28.77 --dec -68.25 
```
or
```
python search.py --nside 32 --ipix 11812
```

To plot DELVE 2, run
```
python plot_hotspot.py --ra 28.77 --dec -68.25 --mod 19.26
```
