'''Plots regridded glimpse map'''

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

datadir = '/mnt/bigdata/tetarenk/ALMA_GRS1915_105/'

fits_file1='/mnt/bigdata/pafreema/GLM_04550-0075_mosaic_I4.fits'
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
proj=spectral_cube.lower_dimensional_structures.Projection.from_hdu(hdulist1)
re=proj.reproject(hdulist.header)


fits.writeto(filename='/mnt/bigdata/tetarenk/ALMA_GRS1915_105/other_data/glimpse8_new.fits',output_verify='ignore',\
	clobber=True,data=np.array(re),header=re.header)


fits_file1='/mnt/bigdata/tetarenk/ALMA_GRS1915_105/other_data/glimpse8_new.fits'
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
hdulist1 = fits.open(fits_file1)[0]#[1]
data1=hdulist1.data*(9.382805e-7)*1e3*(2.5**2)/(1.2**2)#*(6.6**2)/(3.2**2)
wmap1=wcs.WCS(hdulist1.header)
coord0=SkyCoord('19h15m51.0s','+10d38m24.2s',frame='icrs')
coord1=SkyCoord('19h15m26.9s','+10d44m12.68s',frame='icrs')
x1=float(wmap1.wcs_world2pix(coord0.ra.value,coord0.dec.value,1)[0])
y1=float(wmap1.wcs_world2pix(coord0.ra.value,coord0.dec.value,1)[1])
x2=float(wmap1.wcs_world2pix(coord1.ra.value,coord1.dec.value,1)[0])
y2=float(wmap1.wcs_world2pix(coord1.ra.value,coord1.dec.value,1)[1])
coord3=SkyCoord('19h15m37.3s','+10d41m01s',frame='icrs')
x3=float(wmap1.wcs_world2pix(coord3.ra.value,coord3.dec.value,1)[0])
y3=float(wmap1.wcs_world2pix(coord3.ra.value,coord3.dec.value,1)[1])


x=np.arange(0,len(data[0,:]))
y=np.arange(0,len(data[:,0]))
X, Y = np.meshgrid(x, y)
Z=data
levels=np.array([4,6,8,10,15,20,40,60])*0.00005
Z2=data1
levels2=np.array([0.35,0.5,0.6,0.8,1,2,4,6])*1.
fig=plt.figure()
plt.rcdefaults()
plt.rc('xtick.major', size=4)
#plt.rc('xtick', color='w', labelsize='large')
ax1 = fig.add_subplot(111, projection=wmap1.celestial)
im=plt.imshow(np.nan_to_num(data1[:,:]),origin="lower",cmap=cm.get_cmap('hot', 500),norm=colors.PowerNorm(gamma=0.45),vmin=0,vmax=10)#,vmax=0.8)#65,55,65,0.9,0.9,0.9,0.9/x,x,x,0.1,0.03,0.025,0.025,0.025,0.1
cbar=plt.colorbar(im, orientation='vertical',fraction=0.04,pad=0)
cbar.set_label('mJy/beam')
ax1.tick_params(axis='both', which='major', labelsize=15,width=3,length=7,color='k')
ax1.tick_params(axis='both', which='minor', labelsize=15,width=1,length=7,color='k')
ax1.coords['ra'].set_axislabel('Right Ascension')
ax1.coords['dec'].set_axislabel('Declination',minpad=-0.1)
ax1.coords['ra'].set_major_formatter('hh:mm:ss.s')
ax1.set_ylim(y1,y2)
ax1.set_xlim(x1,x2)
#ax1.text(470,550,'4-8 GHz',color='k')
#e1 = patches.Ellipse((x3,y3), 5.31, 4.50,angle=-50.1963, linewidth=2, fill=False,color='m')
#ax1.add_patch(e1)
#ae = AnchoredEllipse(ax1.transData, width=5.31, height=4.5, angle=-50.1963,loc=4, pad=0.5, borderpad=0.4, frameon=True)
#ax1.add_artist(ae)
#plt.contour(X2,Y2,Z2,levels2,colors='w')
plt.contour(X,Y,Z,levels,colors='k',transform=ax1.get_transform(wmap))
plt.savefig(datadir+'for_paper/glimpse_contour.pdf',bbox_inches='tight')
plt.show()
