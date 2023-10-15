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

#R_u = 3.9631
R_g = 3.1863
R_r = 2.1401
R_i = 1.5690
R_z = 1.196
R_y = 1.048 


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
               gmag-{R_g}*ebv AS gmag_dered,
               gerr,
               rmag,
               rmag-{R_r}*ebv AS rmag_dered,
               rerr,
               ebv,
               ndetg,
               ndetr,
               chi,
               prob,
               sharp,
               brickuniq
        FROM delvemc_y4t2.object 
        WHERE q3c_radial_query(ra,dec,{ra},{dec},{radius})
              AND brickuniq = 1
              AND gmag-{R_g}*ebv < 24.2
              AND rmag-{R_r}*ebv < 24.2
              AND ndetg > 2
              AND ndetr > 2
              AND abs(gmag-{R_g}*ebv - rmag-{R_r}*ebv) < 1 
              AND prob >0.8
              AND abs(sharp) < 0.5
              AND chi < 3
    '''

 
    sql_galaxies = f'''
        SELECT ra,
               dec,
               gmag,
               gmag-{R_g}*ebv AS gmag_dered,
               gerr,
               rmag,
               rmag-{R_r}*ebv AS rmag_dered,
               rerr,
               ebv,
               ndetg,
               ndetr,
               chi,
               prob,
               sharp,
               brickuniq
        FROM delvemc_y4t2.object
        WHERE q3c_radial_query(ra,dec,{ra},{dec},{radius})
              AND brickuniq = 1
              AND gmag-{R_g}*ebv < 24.2 
              AND rmag-{R_r}*ebv < 24.2 
              AND ndetg > 2
              AND ndetr > 2
              AND abs(gmag-{R_g}*ebv - rmag-{R_r}*ebv) < 1
              AND prob < 0.5
              AND abs(sharp) >= 0.5
              AND chi > 3
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
