
'''Calcualtes mass from column density maps in a contoured radio region'''

import spectral_cube
from spectral_cube import SpectralCube
from astropy import units as u
import astropy.constants as con
import numpy as np
from astropy.io import fits
import pyfits
import pylab as pl
import math as ma
import matplotlib.pyplot as plt
from astropy.table import hstack
from mpl_toolkits.mplot3d.axes3d import Axes3D
from mpl_toolkits.axes_grid1.inset_locator import mark_inset,zoomed_inset_axes
#import aplpy, atpy
import matplotlib.cm as cm
from astropy import wcs
from astropy.coordinates import SkyCoord
import matplotlib.patches as patches
import matplotlib
import matplotlib.colors as colors
import mpmath

def mass_est(H2col_masked,cell,dist):
	#dist in kpc cell in arcsec
	col_data=H2col_masked[H2col_masked>0]
	numpix=1.#np.count_nonzero(col_data)
	numarea=numpix*(cell*(dist*1000*3e18)/206265)**2
	mass=np.nansum(col_data)*(2*1.00794*1.6e-24)*(numarea)/(1.989e33)
	print 9.11e20*(2*1.00794*1.6e-24)*(numarea)/(1.989e33)
	print np.count_nonzero(col_data)*(cell*(dist*1000*3e18)/206265)**3
	return(mass,numpix,numarea)

datadir = '/mnt/bigdata/tetarenk/ALMA_GRS1915_105/'

#Get rid of extra wcs axes
'''
fits_file1=datadir+'columns_plots/H2_column2.fits'
filename=fits_file1
fh = fits.open(filename)
data = fh[0].data.squeeze() # drops the size-1 axes
header = fh[0].header
mywcs = wcs.WCS(header)
new_header = mywcs.to_header()
new_fh = fits.PrimaryHDU(data=data, header=new_header)
new_fh.writeto(datadir+'columns_plots/GRSALMA_colH.fits')'''

fits_file1=datadir+'columns_plots/GRSALMA_colH.fits'
hdulist = fits.open('/mnt/bigdata/tetarenk/VLA_grs1915_images/GRSVLA.fits')[0]
hdulist.header.remove('CRPIX3')
hdulist.header.remove('CRVAL3')
hdulist.header.remove('CDELT3')
hdulist.header.remove('CUNIT3')
hdulist.header.remove('CTYPE3')
hdulist.header.remove('CRPIX4')
hdulist.header.remove('CRVAL4')
hdulist.header.remove('CDELT4')
#hdulist.header.remove('CUNIT4')
hdulist.header.remove('CTYPE4')
hdulist.header['WCSAXES']=2
data=hdulist.data
wmap=wcs.WCS(hdulist.header)
hdulist1 = fits.open(fits_file1)[0]
hdulist1.header.remove('CRPIX3')
hdulist1.header.remove('CRVAL3')
hdulist1.header.remove('CDELT3')
hdulist1.header.remove('CUNIT3')
hdulist1.header.remove('CTYPE3')
hdulist1.header.remove('CRPIX4')
hdulist1.header.remove('CRVAL4')
hdulist1.header.remove('CDELT4')
#hdulist.header.remove('CUNIT4')
hdulist1.header.remove('CTYPE4')
hdulist1.header['WCSAXES']=2
data1=hdulist1.data
wmap1=wcs.WCS(hdulist1.header)
proj=spectral_cube.lower_dimensional_structures.Projection.from_hdu(hdulist)
re=proj.reproject(hdulist1.header)
rad_alma=np.array(re)

mask=rad_alma>4*0.00005


alma_masked=np.nan_to_num(data1*mask)

mass,numpix,numarea=mass_est(alma_masked,0.2,8.6)

print 'The mass in the selected region is ', "%.2f" % mass,'solar masses'


