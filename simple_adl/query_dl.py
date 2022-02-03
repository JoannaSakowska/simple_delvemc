#!/usr/bin/env python
"""
Generic python script.

Modified for DELVE-MC by 
Joanna Sakowska
"""
__author__ = "Sidney Mau"

import numpy as np
from dl import queryClient as qc

#-------------------------------------------------------------------------------

# DES Redenning coefficients
#R_g = 3.185
#R_r = 2.140
#R_i = 1.571

# Yumi:
# DES DR1 update

#R_u = 3.9631
R_g = 3.1863
R_r = 2.1401
R_i = 1.5690
R_z = 1.196
R_y = 1.048 

# For SMASH:
#R_u = 4.329
#R_g = 3.303
#R_r = 2.285
#R_i = 1.263

#--------------------------------------------------------------------------------

def query(profile, ra, dec, radius=1.0, gmax=24.2, stars=True, galaxies=True):
    """Return data queried from datalab
    Parameters
    ----------
    profile : Profile for data lab query [str]
    ra      : Right Ascension [deg]
    dec     : Declination [deg]
    radius  : radius around (ra, dec) [deg]

    Returns
    -------
    data : numpy recarray of data
    """


    qc.set_profile(profile)


    sql_stars = f'''
        SELECT ra,
               dec,
               gmag,
#               gmag-{R_g}*ebv AS gmag_dered, -- dereddend magnitude
               gerr,
#               rmag-{R_r}*ebv AS rmag_dered, -- dereddend magnitude
               rmag,
               rerr,
               ebv
        FROM delvemc_y2t2.object 
        WHERE q3c_radial_query(ra,dec,{ra},{dec},{radius})
              AND brickuniq = 1 -- NEW SG CUT
              AND gmag < 23.5      -- NEW SG CUT
              AND rmag < 23.5      -- NEW SG CUT
              AND abs(gmag - rmag) < 1  -- NEW SG CUT
              AND prob < 10 -- NEW SG CUT
              AND abs(sharp) < 0.8 -- for star-galaxy separation
              AND chi < 3 -- NEW SG CUT
    '''

 
    sql_galaxies = f'''
        SELECT ra,
               dec,
               gmag,
#               gmag-{R_g}*ebv AS gmag_dered, -- dereddend magnitude
               gerr,
#               rmag-{R_r}*ebv AS rmag_dered, -- dereddend magnitude
               rmag,
               rerr,
               ebv
        FROM delvemc_y2t2.object
        WHERE q3c_radial_query(ra,dec,{ra},{dec},{radius})
              AND brickuniq = 1  -- NEW SG CUT
              AND gmag < 23.5      -- NEW SG CUT
              AND rmag < 23.5      -- NEW SG CUT
              AND abs(gmag - rmag) < 1  -- NEW SG CUT
              AND prob < 10 -- NEW SG CUT
              AND abs(sharp) >= 0.8 -- NEW SG CUT
    '''


    if stars:
        query_stars = qc.query(sql=sql_stars,fmt='structarray')
    if galaxies:
        query_galaxies = qc.query(sql=sql_galaxies,fmt='structarray')


    if stars and not galaxies:
        return(query_stars)
    elif not stars and galaxies:
        return(query_galaxies)
    elif stars and galaxies:
        return(query_stars, query_galaxies)
    else:
        return(None)

#-------------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--profile',type=str,required=False,
                        help='Profile for data lab query [str]')
    parser.add_argument('--ra',type=float,required=True,
                        help='Right Ascension of target position [deg]')
    parser.add_argument('--dec',type=float,required=True,
                        help='Declination of target position [deg]')
    parser.add_argument('--radius',type=float,default=1.0,
                        help='Radius around target position [deg]')
    parser.add_argument('--gmax',type=float,default=23.5,
                        help='Maximum g-band magnitude [mag]')
    args = parser.parse_args()
    data = query(args.profile, args.ra, args.dec, args.radius, args.gmax)
