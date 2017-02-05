from spectral_cube import SpectralCube
from astropy import units as u
import numpy as np
from astropy.io import fits
import pyfits
import pylab as pl
import math as ma
import matplotlib.pyplot as plt
from astropy.table import hstack
from mpl_toolkits.mplot3d.axes3d import Axes3D
#import aplpy, atpy
import matplotlib.cm as cm
from astropy import wcs
from astropy.coordinates import SkyCoord
import matplotlib.patches as patches
import matplotlib
import matplotlib.colors as colors


datadir = '/mnt/bigdata/tetarenk/ALMA_GRS1915_105/'
lines=['13CO','12CO','18CO','H2CO_303_202','H2CO_322_221','H2CO_321_220','SiO','N2Dp','H30a']

def make_moment_maps(fitsfile,mind,ddir,line,flag):
	header = fits.getdata(fitsfile, header=True)[1]
	X=SpectralCube.read(fitsfile)
	X.allow_huge_operations=True
	if flag=='K':
		X1=X.to(u.K, equivalencies=X.beam.jtok_equiv(X.header['RESTFRQ']*u.Hz))
		Xkms=X1.with_spectral_unit(u.km/u.s,velocity_convention='radio')
	else:
		Xkms=X.with_spectral_unit(u.km/u.s,velocity_convention='radio')
	mom=Xkms.moment(order=mind)
	mom_array=np.array(mom)
	fits.writeto(filename=ddir+'moment_maps/'+line+'_moment'+str(mind)+flag+'.fits',output_verify='ignore',\
	clobber=True,data=mom_array,header=header)


#for i in lines:
	#make_moment_maps(datadir+'alex_imaging_'+i+'/GRS1915_modelimg_'+i+'.image.pbcor.fits',0,datadir,i)
	#make_moment_maps(datadir+'alex_imaging_'+i+'/GRS1915_modelimg_'+i+'.image.pbcor.fits',1,datadir,i)
	#make_moment_maps(datadir+'alex_imaging_'+i+'/GRS1915_modelimg_'+i+'.image.pbcor.fits',2,datadir,i)

i='18CO'
flag='K'
make_moment_maps(datadir+'alex_imaging_'+i+'_fix_all/GRS1915_modelimg_'+i+'.image.pbcor.fits',0,datadir,i,flag)
#make_moment_maps(datadir+'alex_imaging_'+i+'/GRS1915_modelimg_'+i+'.image.pbcor.fits',1,datadir,i)
#make_moment_maps(datadir+'alex_imaging_'+i+'/GRS1915_modelimg_'+i+'.image.pbcor.fits',2,datadir,i)
#imview(datadir+'moment_maps/'+i+'_moment'+str(0)+'.fits')

fits_file1=datadir+'moment_maps/'+i+'_moment'+str(0)+flag+'.fits'
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
data1=hdulist1.data
wmap1=wcs.WCS(hdulist1.header)


x=np.arange(0,len(data[0,:]))
y=np.arange(0,len(data[:,0]))
X, Y = np.meshgrid(x, y)
Z=data#[0,0,490:550,470:550]
levels=np.array([1,2,3,4,5,6,7])*0.000345
#evels=np.array([1,2,3,4,5,6,7])*0.000345
fig=plt.figure()
#plt.rcdefaults()
plt.rc('xtick.major', size=4)
plt.rc('xtick', color='w', labelsize='large')
ax1 = fig.add_subplot(111, projection=wmap1.celestial)
im=plt.imshow(np.nan_to_num(data1),origin="lower",cmap=cm.get_cmap('jet', 500),norm=colors.PowerNorm(gamma=1),vmin=0.0)
cbar=plt.colorbar(im, orientation='vertical',fraction=0.04,pad=0)
#cbar.set_label('$N_{{\\rm H}}\\,\\times10^{22}\\, {\\rm cm}^{-2}$')
if flag=='K':
	cbar.set_label('K km/s')
else:
	cbar.set_label('Jy km/s')
ax1.tick_params(axis='both', which='major', labelsize=15,width=3,length=7,color='k')
ax1.tick_params(axis='both', which='minor', labelsize=15,width=1,length=7,color='k')
ax1.coords['ra'].set_axislabel('Right Ascension')
ax1.coords['dec'].set_axislabel('Declination',minpad=-0.1)
ax1.coords['ra'].set_major_formatter('hh:mm:ss.s')
ax1.set_ylim(150, 700)
ax1.set_xlim(100, 650)
plt.contour(X,Y,Z,levels,colors='k',transform=ax1.get_transform(wmap))
plt.savefig(datadir+'moment_maps/'+i+'_moment'+str(0)+'_contour.pdf',bbox_inches='tight')
plt.show()

import scipy.ndimage as nd
#rad=100#100,125,125,125,125,100
#data_mask=nd.morphology.binary_erosion(np.nan_to_num(data1), np.ones((rad, rad)))
#import os
#os.system('cp -r '+datadir+'moment_maps/'+i+'_moment'+str(0)+'_contour.pdf /home/ubuntu/Dropbox')



