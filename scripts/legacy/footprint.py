import glob 
from astropy.coordinates import SkyCoord
from irspec.datacube import Datacube
from regions import PolygonSkyRegion, PixCoord
import matplotlib.pyplot as plt


g395h_path = glob.glob(f"./irspec/data/g395h-f290lp_s3d.fits")[0]
g395hcube = Datacube(g395h_path, redshift=0.044601, verbose=True)
g395hregion_data = g395hcube.science_header["S_REGION"].split()
g395h_pixel_x = [6, -3, 82, 92]
g395h_pixel_y = [-3, 71, 83, 10]
g395h_pixcoord = PixCoord(g395h_pixel_x, g395h_pixel_y)
g395h_skycoord = g395h_pixcoord.to_sky(g395hcube.wcs)

ch1s_path = glob.glob(f"./irspec/data/ch1-short_s3d.fits")[0]
ch1scube = Datacube(ch1s_path, redshift=0.044601, verbose=True)
ch1sregion_data = ch1scube.science_header["S_REGION"].split()
ch1s_pixel_x = [26, 0, 22, 48]
ch1s_pixel_y = [0, 28, 48, 21]
ch1s_pixcoord = PixCoord(ch1s_pixel_x, ch1s_pixel_y)
ch1s_skycoord = ch1s_pixcoord.to_sky(ch1scube.wcs)

ch1m_path = glob.glob(f"./irspec/data/ch1-medium_s3d.fits")[0]
ch1mcube = Datacube(ch1m_path, redshift=0.044601, verbose=True)
ch1mregion_data = ch1mcube.science_header["S_REGION"].split()
ch1m_pixel_x = [24, 0, 22, 47]
ch1m_pixel_y = [0, 28, 49, 21]
ch1m_pixcoord = PixCoord(ch1m_pixel_x, ch1m_pixel_y)
ch1m_skycoord = ch1m_pixcoord.to_sky(ch1mcube.wcs)

ch1l_path = glob.glob(f"./irspec/data/ch1-long_s3d.fits")[0]
ch1lcube = Datacube(ch1l_path, redshift=0.044601, verbose=True)
ch1lregion_data = ch1lcube.science_header["S_REGION"].split()
ch1l_pixel_x = [24, -1, 22, 48]
ch1l_pixel_y = [0, 28, 49, 21]
ch1l_pixcoord = PixCoord(ch1l_pixel_x, ch1l_pixel_y)
ch1l_skycoord = ch1l_pixcoord.to_sky(ch1lcube.wcs)

ch2s_path = glob.glob(f"./irspec/data/ch2-short_s3d.fits")[0]
ch2scube = Datacube(ch2s_path, redshift=0.044601, verbose=True)
ch2sregion_data = ch2scube.science_header["S_REGION"].split()
ch2s_pixel_x = [25, 0, 20, 44]
ch2s_pixel_y = [0, 25, 45, 19]
ch2s_pixcoord = PixCoord(ch2s_pixel_x, ch2s_pixel_y)
ch2s_skycoord = ch2s_pixcoord.to_sky(ch2scube.wcs)

ch2m_path = glob.glob(f"./irspec/data/ch2-medium_s3d.fits")[0]
ch2mcube = Datacube(ch2m_path, redshift=0.044601, verbose=True)
ch2mregion_data = ch2mcube.science_header["S_REGION"].split()
ch2m_pixel_x = [25, 0, 20, 44]
ch2m_pixel_y = [0, 25, 45, 19]
ch2m_pixcoord = PixCoord(ch2m_pixel_x, ch2m_pixel_y)
ch2m_skycoord = ch2m_pixcoord.to_sky(ch2mcube.wcs)

ch2l_path = glob.glob(f"./irspec/data/ch2-long_s3d.fits")[0]
ch2lcube = Datacube(ch2l_path, redshift=0.044601, verbose=True)
ch2lregion_data = ch2lcube.science_header["S_REGION"].split()
ch2l_pixel_x = [25, 0, 20, 44]
ch2l_pixel_y = [0, 25, 45, 19]
ch2l_pixcoord = PixCoord(ch2l_pixel_x, ch2l_pixel_y)
ch2l_skycoord = ch2l_pixcoord.to_sky(ch2lcube.wcs)

ch3s_path = glob.glob(f"./irspec/data/ch3-short_s3d.fits")[0]
ch3scube = Datacube(ch3s_path, redshift=0.044601, verbose=True)
ch3sregion_data = ch3scube.science_header["S_REGION"].split()
ch3s_pixel_x = [24, -1, 24, 48]
ch3s_pixel_y = [-3, 26, 49, 21]
ch3s_pixcoord = PixCoord(ch3s_pixel_x, ch3s_pixel_y)
ch3s_skycoord = ch3s_pixcoord.to_sky(ch3scube.wcs)

ch3m_path = glob.glob(f"./irspec/data/ch3-medium_s3d.fits")[0]
ch3mcube = Datacube(ch3m_path, redshift=0.044601, verbose=True)
ch3mregion_data = ch3mcube.science_header["S_REGION"].split()
ch3m_pixel_x = [26, 0, 24, 48]
ch3m_pixel_y = [-1, 27, 49, 21]
ch3m_pixcoord = PixCoord(ch3m_pixel_x, ch3m_pixel_y)
ch3m_skycoord = ch3m_pixcoord.to_sky(ch3mcube.wcs)

ch3l_path = glob.glob(f"./irspec/data/ch3-long_s3d.fits")[0]
ch3lcube = Datacube(ch3l_path, redshift=0.044601, verbose=True)
ch3lregion_data = ch3lcube.science_header["S_REGION"].split()
ch3l_pixel_x = [26, 0, 23, 48]
ch3l_pixel_y = [-1, 27, 49, 21]
ch3l_pixcoord = PixCoord(ch3l_pixel_x, ch3l_pixel_y)
ch3l_skycoord = ch3l_pixcoord.to_sky(ch3lcube.wcs)

ch4s_path = glob.glob(f"./irspec/data/ch4-short_s3d.fits")[0]
ch4scube = Datacube(ch4s_path, redshift=0.044601, verbose=True)
ch4sregion_data = ch4scube.science_header["S_REGION"].split()
ch4s_pixel_x = [19, 0, 16, 35]
ch4s_pixel_y = [-1, 18, 34, 15]
ch4s_pixcoord = PixCoord(ch4s_pixel_x, ch4s_pixel_y)
ch4s_skycoord = ch4s_pixcoord.to_sky(ch4scube.wcs)

ch4m_path = glob.glob(f"./irspec/data/ch4-medium_s3d.fits")[0]
ch4mcube = Datacube(ch4m_path, redshift=0.044601, verbose=True)
ch4mregion_data = ch4mcube.science_header["S_REGION"].split()
ch4m_pixel_x = [19, -1, 15, 35]
ch4m_pixel_y = [-1, 19, 35, 15]
ch4m_pixcoord = PixCoord(ch4m_pixel_x, ch4m_pixel_y)
ch4m_skycoord = ch4m_pixcoord.to_sky(ch4mcube.wcs)

ch4l_path = glob.glob(f"./irspec/data/ch4-long_s3d.fits")[0]
ch4lcube = Datacube(ch4l_path, redshift=0.044601, verbose=True)
ch4region_data = ch4lcube.science_header["S_REGION"].split()
ch4l_pixel_x = [15, 35, 19, -1]
ch4l_pixel_y = [35, 15, -1, 19]
ch4l_pixcoord = PixCoord(ch4l_pixel_x, ch4l_pixel_y)
ch4l_skycoord = ch4l_pixcoord.to_sky(ch4lcube.wcs)

# NORTHERN CUTOUTS

nch1s_path = glob.glob(f"./irspec/north_data/ch1-short_s3d.fits")[0]
nch1scube = Datacube(nch1s_path, redshift=0.044601, verbose=True)
nch1sregion_data = nch1scube.science_header["S_REGION"].split()
nch1s_pixel_x = [26, 0, 22, 48]
nch1s_pixel_y = [0, 28, 48, 21]
nch1s_pixcoord = PixCoord(nch1s_pixel_x, nch1s_pixel_y)
nch1s_skycoord = nch1s_pixcoord.to_sky(nch1scube.wcs)

nch4l_path = glob.glob(f"./irspec/north_data/ch4-long_s3d.fits")[0]
nch4lcube = Datacube(nch4l_path, redshift=0.044601, verbose=True)
nch4region_data = nch4lcube.science_header["S_REGION"].split()
nch4l_pixel_x = [15, 35, 19, -1]
nch4l_pixel_y = [35, 15, -1, 19]
nch4l_pixcoord = PixCoord(nch4l_pixel_x, nch4l_pixel_y)
nch4l_skycoord = nch4l_pixcoord.to_sky(nch4lcube.wcs)

fig = plt.figure()
#ax = plt.subplot(projection=nch1scube.wcs)
ax = plt.subplot()
sky_region = PolygonSkyRegion(nch4l_skycoord)
ax.set_xlim(-10, 100)
ax.set_ylim(-10, 100)
ax.imshow(nch4lcube.science_data[0], origin="lower")
pixel_region = sky_region.to_pixel(nch4lcube.wcs)
pixel_region.plot(ax=ax, color="black", lw=2.0)
plt.show()