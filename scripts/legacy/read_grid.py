from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import glob
from irspec.datacube import Datacube
from irspec.cubespec import CubeSpec
from photutils.aperture import CircularAperture
from regions import CircleSkyRegion
from astropy.visualization.wcsaxes import SphericalCircle
import astropy.units as u
from matplotlib.patches import Circle




def generate_grid(gridPointsX, gridPointsY, gridstep_dist, pixel_scale, wcs, cube, x_pix, y_pix, r_ap):
    NX = np.arange(0,gridPointsX)
    NY = np.arange(0,gridPointsY)
    
    gridstep_dist_pix = gridstep_dist / pixel_scale
    if r_ap == -1:
        print('Warning: The spaxel size is not defined. Using as default the distance between spaxels.') 
        r_ap = gridstep_dist / 2
        r_pix = gridstep_dist_pix / 2
    else:    
        r_pix = r_ap / pixel_scale
    
    print(r_pix, r_ap / pixel_scale, gridstep_dist_pix)
    sky = wcs.pixel_to_world(x_pix, y_pix)
    user_x, user_y = cube.wcs.world_to_pixel(sky)
    print((NX - (gridPointsX-1)/2), gridstep_dist_pix)
    grids_xs = user_x + (NX - (gridPointsX-1)/2) * gridstep_dist_pix
    grids_ys = user_y + (NY - (gridPointsY-1)/2) * gridstep_dist_pix
    
    sky_list = []
    pixels_list = []
    real_pixels_grid = []
    names = []
    sky_ra = []
    sky_dec = []
    for i in range(len(grids_xs)):
        for j in range(len(grids_ys)):
            sky = cube.wcs.pixel_to_world(grids_xs[i], grids_ys[j])
            real_pixels_grid.append([grids_xs[i], grids_ys[j]])   
            sky_list.append(sky)
            pixels_list.append([i,j])
            names.append(str(i)+"_"+str(j))
            sky_ra.append(sky.ra)
            sky_dec.append(sky.dec)
    return sky_list, pixels_list, real_pixels_grid, names, sky_ra, sky_dec


datacube = Datacube(glob.glob(f"./irspec/data/ch1-short_s3d.fits")[0], redshift=0.044601, verbose=True)
sky_list, pixels_list, real_pixels_grid, names, sky_ra, sky_dec = generate_grid(7, 7, 1.6, 0.196, datacube.wcs, datacube, 24, 27, 0.8)
def generate_param_files(sky_list, sky_ra, sky_dec):
    for idx, coord in enumerate(sky_list):
        ra, dec = coord.to_string('hmsdms').split(" ")
        with open(f"./irspec/param_files/manual_grid/{idx}_single_param.txt", 'w') as param_file:
            param_file.write("cubes = g395h, ch1-short, ch1-medium, ch1-long, ch2-short, ch2-medium, ch2-long, ch3-short, ch3-medium, ch3-long, ch4-short, ch4-medium, ch4-long\n")
            param_file.write("user_r_ap = 0.8                 # [arcsec] Aperture size\n")
            param_file.write(f"user_ra = {ra}         # RA of the extraction location\n")
            param_file.write(f"user_dec = {dec.replace('+', '')}        # Dec of the extraction location\n")
            param_file.write("point_source = False            # [True/False] True = Cone extraction. False = Cylinder extraction\n")
            param_file.write("lambda_ap = 5.4                 # [microns] The user-specified aperture size (user_r_ap) will be defined at this wavelength (ignored if point_source=False; cylinder extraction)\n")
            param_file.write("aperture_correction = False     # [True/False]\n")
            param_file.write("centering = True                # [True/False] At lambda_cent, this will update user_ra,user_dec to recenter the extraction location\n")
            param_file.write("lambda_cent = 5.4               # [microns] The wavelength at which to perform the centering (ignored if centering=False)\n")
            param_file.write("background_sub = False          # [True/False]\n")
            param_file.write("r_ann_in = 0.0                  # [arcsec] When background_sub=True, this defines the inner annulus size. This size will NOT scale with aperture size even when point_source=True\n")
            param_file.write("ann_width = 0.0                 # [arcsec] When background_sub=True, this defines the annulus width. This size will NOT scale with aperture size even when point_source=True")

#generate_param_files(sky_list, sky_ra, sky_dec)
fig, ax = plt.subplots()
for idx, pix in enumerate(real_pixels_grid):
    ap = CircularAperture(pix, 1/0.196 * 0.8)
    ap.plot(ax, label=idx)
ax.set_xlim(0, 50)
ax.set_ylim(0, 50)
plt.show()
"""fig = plt.figure()
ax = fig.add_subplot(projection=datacube.wcs)
for idx, coord in enumerate(sky_ra):
    #s = Circle((coord.value, sky_dec[idx].value), 0.196, edgecolor='black', facecolor='blue')
    #ax.add_patch(s)
    ax.scatter(coord.to(u.radian), sky_dec[idx].to(u.radian), marker="o", s=2, alpha=0.3)
plt.show()"""