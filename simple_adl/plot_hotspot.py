#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "Sid Mau"

# Python libraries
import os
from sre_constants import AT_END_STRING
import yaml
import numpy as np
import healpy as hp
import scipy.ndimage

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import simple_adl.survey
import simple_adl.isochrone
import simple_adl.coordinate_tools
from simple_adl.search import cut_isochrone_path

#-------------------------------------------------------------------------------

props = dict(facecolor='white', edgecolor='black', linewidth=1) #dict(facecolor='white', edgecolor='none', alpha=0.7)

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--config',type=str,required=False,default='config.yaml',
                        help='Config file [.yaml]')
    parser.add_argument('--outfile',type=str,required=False,default='out.png',
                        help='Output file [.png/.pdf]')
    parser.add_argument('--ra',type=float,required=True,
                        help='Right Ascension of target position [deg]')
    parser.add_argument('--dec',type=float,required=True,
                        help='Declination of target position [deg]')
    parser.add_argument('--mod',type=float,required=False,
                        help='Distance modulus [mag]')
    args = parser.parse_args()

    with open(args.config, 'r') as ymlfile:
        cfg = yaml.load(ymlfile, Loader=yaml.SafeLoader)
        survey = simple_adl.survey.Survey(cfg)

    #---------------------------------------------------------------------------

    region = simple_adl.survey.Region(survey, args.ra, args.dec)
    print('Plot coordinates: (RA, Dec) = ({:0.2f}, {:0.2f})'.format(region.ra, region.dec))

    # Apply filters

    data = region.load_data(stars=True, galaxies=False) # Data is stars as default
    stars = data

    galaxies = region.load_data(stars=False, galaxies=True)

    print('Found {} objects'.format(len(region.data)))
    if (len(region.data) == 0):
        print('Ending search.')
        exit()



    #---------------------------------------------------------------------------

    iso = simple_adl.isochrone.Isochrone(survey=survey.isochrone['survey'],
                              band_1=survey.band_1.lower(),
                              band_2=survey.band_2.lower(),
                              age=12.0, #survey.isochrone['age'],
                              metallicity=0.00010) #survey.isochrone['metallicity']
    if args.mod:
        iso.distance_modulus = args.mod
    #iso_sep = iso.separation(data[mag_1], data[mag_2])

    
    # Isochrone for stars

    iso_filter = cut_isochrone_path(
        stars[survey.mag_dered_1], 
        stars[survey.mag_dered_2], 
        stars[survey.mag_err_1], 
        stars[survey.mag_err_2], 
        iso, 
        mag_max=survey.catalog['mag_max'], 
        radius=0.1, 
        return_all=False)

        
    
    # projection of image
    proj = region.proj
    x_stars, y_stars = proj.sphereToImage(stars[survey.catalog['basis_1']], stars[survey.catalog['basis_2']])
    x_galaxies, y_galaxies = proj.sphereToImage(galaxies[survey.catalog['basis_1']], galaxies[survey.catalog['basis_2']])

    # hess
    mag = stars[survey.mag_dered_1]
    color = stars[survey.mag_dered_1] - stars[survey.mag_dered_2]
    
    # Filters
    #extension = 0.05 # JS: Minimised the radius to show closer objects
    extension = 0.025
    r0 = 3.0 * extension # 3.0 # g-radius # can be minimised or made larger # IS R READ IN?
    r1 = 5.0 * extension # 5.0 # rnear
    r2 = np.sqrt(r0**2 + r1**2) # rfar
    angsep = simple_adl.coordinate_tools.angsep(args.ra, args.dec, stars[survey.catalog['basis_1']], stars[survey.catalog['basis_2']])
    inner = (angsep < r0)
    outer = ((angsep > r1) & (angsep < r2))
    background = (angsep > r2)

    #---------------------------------------------------------------------------

    fig, axs = plt.subplots(nrows=2,  ncols=3, figsize=(16, 8), dpi=600)
    fig.subplots_adjust(wspace=0.5)

    #---------------------------------------------------------------------------
    
    # Stellar histogram
    ax = axs[0,0]
    plt.sca(ax)

    #r0 = g_radius
    
    bound = 0.5
    steps = 100
    bins = np.linspace(-bound, bound, steps)
    signal = np.histogram2d(x_stars[iso_filter], y_stars[iso_filter], bins=[bins, bins])[0]
    sigma = 0.01 * (0.25 * np.arctan(0.25 * r0 * 60. - 1.5) + 1.3)
    convolution = scipy.ndimage.filters.gaussian_filter(signal, sigma/(bound/steps)).T
    pc = ax.pcolormesh(bins, bins, convolution, cmap='Greys', rasterized=True)
    
    # search kernel
    #x, y = ra_proj, dec_proj
    #delta_x = 0.01
    #area = delta_x**2
    #smoothing = 2. / 60. # Was 3 arcmin
    #bins = np.arange(-8., 8. + 1.e-10, delta_x)
    #centers = 0.5 * (bins[0: -1] + bins[1:])
    #yy, xx = np.meshgrid(centers, centers)
    #h = np.histogram2d(x, y, bins=[bins, bins])[0]
    #h_g = scipy.ndimage.filters.gaussian_filter(h, smoothing / delta_x)
    #pc = ax.pcolormesh(bins, bins, h_g.T, cmap='Greys', rasterized=True)
    
    ax.text(0.05, 0.95, 'Stars', transform=ax.transAxes, verticalalignment='top', bbox=props)
    
    ax.set_xlim(0.5, -0.5)
    ax.set_ylim(-0.5, 0.5)
    ax.set_xlabel(r'$\Delta {\rm RA}$ (deg)')
    ax.set_ylabel(r'$\Delta {\rm Dec}$ (deg)')
    
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0)
    cb = fig.colorbar(pc, cax=cax)
    #cb.set_label('Counts')

    #---------------------------------------------------------------------------
    
    # Galactic histogram
    ax = axs[0,1]
    plt.sca(ax)

    bound = 0.5
    steps = 100
    bins = np.linspace(-bound, bound, steps)
    signal = np.histogram2d(x_galaxies, y_galaxies, bins=[bins, bins])[0]
    sigma = 0.01 * (0.25 * np.arctan(0.25 * r0 * 60. - 1.5) + 1.3)
    convolution = scipy.ndimage.filters.gaussian_filter(signal, sigma/(bound/steps)).T
    pc = ax.pcolormesh(bins, bins, convolution, cmap='Greys', rasterized=True)
    
    ## search kernel
    ##x, y = ra_proj, dec_proj
    ##delta_x = 0.01
    ##area = delta_x**2
    ##smoothing = 2. / 60. # Was 3 arcmin
    ##bins = np.arange(-8., 8. + 1.e-10, delta_x)
    ##centers = 0.5 * (bins[0: -1] + bins[1:])
    ##yy, xx = np.meshgrid(centers, centers)
    ##h = np.histogram2d(x, y, bins=[bins, bins])[0]
    ##h_g = scipy.ndimage.filters.gaussian_filter(h, smoothing / delta_x)
    ##pc = ax.pcolormesh(bins, bins, h_g.T, cmap='Greys', rasterized=True)
    
    ax.text(0.05, 0.95, 'Galaxies', transform=ax.transAxes, verticalalignment='top', bbox=props)
    
    ax.set_xlim(0.5, -0.5)
    ax.set_ylim(-0.5, 0.5)
    ax.set_xlabel(r'$\Delta \alpha_{2000}$ (deg)')
    ax.set_ylabel(r'$\Delta \delta_{2000}$ (deg)')
    
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0)
    cb = fig.colorbar(pc, cax=cax)
    #cb.set_label('Counts')

    #---------------------------------------------------------------------------
    
    # Hess
    #ax = axs[1]
    ax = axs[0,2]
    plt.sca(ax)
    
    xbins = np.arange(-0.3, 1.1, 0.1)
    ybins = np.arange(16., 24.0 + 0.5, 0.5)
    foreground = np.histogram2d(color[inner], mag[inner], bins=[xbins, ybins])
    background = np.histogram2d(color[outer], mag[outer], bins=[xbins, ybins])
    fg = foreground[0].T
    bg = background[0].T
    fg_abs = np.absolute(fg)
    bg_abs = np.absolute(bg)
    mask_abs = fg_abs + bg_abs
    mask_abs[mask_abs == 0.] = np.nan # mask common zeroes
    signal = fg - bg
    signal = np.ma.array(signal, mask=np.isnan(mask_abs)) # mask nan
    pc = ax.pcolormesh(xbins, ybins, signal, cmap='viridis', rasterized=True)
    if args.mod:
        plt.plot(iso.color, iso.mag+iso.distance_modulus, lw=1, c='k')
    
    ax.set_xlim(-0.3, 1.0)
    ax.set_ylim(24.0, 16.0)
    ax.set_xlabel('${} - {}$'.format(survey.band_1.lower(), survey.band_2.lower()))
    ax.set_ylabel('${}$'.format(survey.band_1.lower()))
    
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0)
    cb = fig.colorbar(pc, cax=cax)
    #cb.set_label('Counts')

    #---------------------------------------------------------------------------
   
    #file_name = 'candidate_{:0.2f}_{:0.2f}'.format(args.ra, args.dec)
    #fig.savefig('./{}.pdf'.format(file_name), bbox_inches='tight')
    #fig.savefig(os.path.join(survey.output['save_dir'],args.outfile), bbox_inches='tight')
    #plt.close(fig)

    #---------------------------------------------------------------------------
    # TEST FIGS

    #fig2, axs = plt.subplots(1, 3, figsize=(16, 4))
    #fig2.subplots_adjust(wspace=0.5)

    #---------------------------------------------------------------------------

    """Color-magnitude plot STARS"""
    ax = axs[1,0]
    plt.sca(ax)


    #angsep = ugali.utils.projector.angsep(ra, dec, data[basis_1], data[basis_2])
    #annulus = (angsep > g_radius) & (angsep < 1.)

    #if type == 'stars':
    #    filter = simple.filters.star_filter(survey, data)
    #    #ax.set_title('Stellar Color-Magnitude')
    #    ax.text(0.05, 0.95, 'Stars', transform=ax.transAxes, verticalalignment='top')
    #elif type == 'galaxies':
    #    filter = simple.filters.galaxy_filter(survey, data)
    #    #ax.set_title('Galactic Color-Magnitude')
    #    ax.text(0.05, 0.95, 'Galaxies', transform=ax.transAxes, verticalalignment='top')

    #iso_filter = simple.simple_utils.cut_isochrone_path(data[mag_dered_1], data[mag_dered_2], data[mag_err_1], data[mag_err_2], iso, radius=0.1, return_all=False)
    
    angsep = angsep
    annulus = (angsep > r0) & (angsep < 1.)
    nbhd = (angsep < r0)

    ## Plot background objects
    ax.scatter(stars[survey.mag_dered_1][annulus & iso_filter] - stars[survey.mag_dered_2][annulus & iso_filter], stars[survey.mag_dered_1][annulus & iso_filter], c='k', alpha=0.1, edgecolor='none', s=1)

    # FOCUS ON STARS
    ax.scatter(stars[survey.mag_dered_1][annulus] - stars[survey.mag_dered_2][annulus], stars[survey.mag_dered_1][annulus], c='k', alpha=0.1, edgecolor='none', s=1)

    ## Plot isochrone
    ##ugali.utils.plotting.drawIsochrone(iso, lw=2, label='{} Gyr, z = {}'.format(iso.age, iso.metallicity))
    #ax.plot(iso.color, iso.mag_1 + iso.distance_modulus, c='k', lw=1)
    ax.plot(iso.color, iso.mag + iso.distance_modulus, c='k', lw=1)

    ## Plot objects in nbhd
    #ax.scatter(data[mag_dered_1][filter & nbhd] - data[mag_dered_2][filter & nbhd], data[mag_dered_1][filter & nbhd], c='g', s=5, label='r < {:.3f}$^\circ$'.format(g_radius))
    ax.scatter(stars[survey.mag_dered_1][nbhd] - stars[survey.mag_dered_2][nbhd], stars[survey.mag_dered_1][nbhd], c='g', s=5, label='r < {:.3f}$^\circ$'.format(r0))

    ## Plot objects in nbhd and near isochrone
    #ax.scatter(data[mag_dered_1][filter & nbhd & iso_filter] - data[mag_dered_2][filter & nbhd & iso_filter], data[mag_dered_1][filter & nbhd & iso_filter], c='r', s=5, label='$\Delta$CM < 0.1')

    ax.scatter(stars[survey.mag_dered_1][nbhd & iso_filter] - stars[survey.mag_dered_2][nbhd & iso_filter], stars[survey.mag_dered_1][nbhd & iso_filter], c='r', s=5, label='$\Delta$CM < 0.1')

    ax.legend(loc='upper right')

    ax.set_xlim(-0.5, 1)
    ax.set_ylim(survey.catalog['mag_max'], 16)
    ax.text(0.05, 0.95, 'Stars', transform=ax.transAxes, verticalalignment='top', bbox=props)

    #ax.set_xlabel('{} - {} (mag)'.format(band_1.lower(), band_2.lower()))
    #ax.set_ylabel('{} (mag)'.format(band_1.lower()))

    """Color-magnitude plot GALAXIES"""
    ax = axs[1,1]    
    plt.sca(ax)


    #angsep = ugali.utils.projector.angsep(ra, dec, data[basis_1], data[basis_2])
    #annulus = (angsep > g_radius) & (angsep < 1.)

    #if type == 'stars':
    #    filter = simple.filters.star_filter(survey, data)
    #    #ax.set_title('Stellar Color-Magnitude')
    #    ax.text(0.05, 0.95, 'Stars', transform=ax.transAxes, verticalalignment='top')
    #elif type == 'galaxies':
    #    filter = simple.filters.galaxy_filter(survey, data)
    #    #ax.set_title('Galactic Color-Magnitude')
    #    ax.text(0.05, 0.95, 'Galaxies', transform=ax.transAxes, verticalalignment='top')

    #iso_filter = simple.simple_utils.cut_isochrone_path(data[mag_dered_1], data[mag_dered_2], data[mag_err_1], data[mag_err_2], iso, radius=0.1, return_all=False)
    
    angsep2 = simple_adl.coordinate_tools.angsep(args.ra, args.dec, galaxies[survey.catalog['basis_1']], galaxies[survey.catalog['basis_2']])

    annulus = (angsep2 > r0) & (angsep2 < 1.)
    nbhd = (angsep2 < r0)

    iso_filter = cut_isochrone_path(
        galaxies[survey.mag_dered_1], 
        galaxies[survey.mag_dered_2], 
        galaxies[survey.mag_err_1], 
        galaxies[survey.mag_err_2], 
        iso, 
        mag_max=survey.catalog['mag_max'], 
        radius=0.1, 
        return_all=False)

    ## Plot background objects
    ax.scatter(galaxies[survey.mag_dered_1][annulus & iso_filter] - galaxies[survey.mag_dered_2][annulus & iso_filter], galaxies[survey.mag_dered_1][annulus & iso_filter], c='k', alpha=0.1, edgecolor='none', s=1)

    # FOCUS ON GALAXIES
    ax.scatter(galaxies[survey.mag_dered_1][annulus] - galaxies[survey.mag_dered_2][annulus], galaxies[survey.mag_dered_1][annulus], c='k', alpha=0.1, edgecolor='none', s=1)

    ## Plot isochrone
    ##ugali.utils.plotting.drawIsochrone(iso, lw=2, label='{} Gyr, z = {}'.format(iso.age, iso.metallicity))
    #ax.plot(iso.color, iso.mag_1 + iso.distance_modulus, c='k', lw=1)
    ax.plot(iso.color, iso.mag + iso.distance_modulus, c='k', lw=1)

    ## Plot objects in nbhd
    #ax.scatter(data[mag_dered_1][filter & nbhd] - data[mag_dered_2][filter & nbhd], data[mag_dered_1][filter & nbhd], c='g', s=5, label='r < {:.3f}$^\circ$'.format(g_radius))
    ax.scatter(galaxies[survey.mag_dered_1][nbhd] - galaxies[survey.mag_dered_2][nbhd], galaxies[survey.mag_dered_1][nbhd], c='g', s=5, label='r < {:.3f}$^\circ$'.format(r0))

    ## Plot objects in nbhd and near isochrone
    #ax.scatter(data[mag_dered_1][filter & nbhd & iso_filter] - data[mag_dered_2][filter & nbhd & iso_filter], data[mag_dered_1][filter & nbhd & iso_filter], c='r', s=5, label='$\Delta$CM < 0.1')

    ax.scatter(galaxies[survey.mag_dered_1][nbhd & iso_filter] - galaxies[survey.mag_dered_2][nbhd & iso_filter], galaxies[survey.mag_dered_1][nbhd & iso_filter], c='r', s=5, label='$\Delta$CM < 0.1')

    ax.legend(loc='upper right')

    ax.set_xlim(-0.5, 1)
    ax.set_ylim(survey.catalog['mag_max'], 16)
    #ax.set_xlabel('{} - {} (mag)'.format(band_1.lower(), band_2.lower()))
    #ax.set_ylabel('{} (mag)'.format(band_1.lower()))
    ax.text(0.05, 0.95, 'Galaxies', transform=ax.transAxes, verticalalignment='top', bbox=props)



    #---------------------------------------------------------------------------
   
    #file_name = 'candidate_{:0.2f}_{:0.2f}'.format(args.ra, args.dec)
    #fig.savefig('./{}.pdf'.format(file_name), bbox_inches='tight')
    fig.savefig(os.path.join(survey.output['save_dir'],args.outfile), bbox_inches='tight')
    plt.close(fig)





