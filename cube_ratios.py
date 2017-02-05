from spectral_cube import SpectralCube
from astropy.io import fits
import astropy.units as u
import pyfits
import numpy as np
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

def line_ratio(fitsfile1,fitsfile2,ddir,noisemult,noise1,noise2,line1,line2):
	header = fits.getdata(fitsfile1, header=True)[1]
	X=SpectralCube.read(fitsfile1)
	Xkms=X.with_spectral_unit(u.km/u.s,velocity_convention='radio')
	Y=SpectralCube.read(fitsfile2)
	Ykms=Y.with_spectral_unit(u.km/u.s,velocity_convention='radio')
	Y_interp=Ykms.spectral_interpolate(spectral_grid=Xkms.spectral_axis)
	val1=pyfits.open(fitsfile1)
	val2=pyfits.open(fitsfile2)
	mask1=val1[0].data>noisemult*noise1
	mask2=val2[0].data>noisemult*noise2
	badmask1=val1[0]<=0
	badmask2=val2[0]<=0
	keep1=mask1*(~badmask1)
	keep2=mask2*(~badmask2)
	X_masked=Xkms.with_mask(keep1[0,:,:,:])
	Y_masked=Y_interp.with_mask(keep2[0,0:60,:,:])
	ratio=X_masked/Y_masked
	fits.writeto(filename=ddir+'ratio_maps/'+line1+'_'+line2+'_ratiomap.fits',output_verify='ignore',\
	clobber=True,data=ratio,header=header)
def line_ratio_tmax(tmaxfile1,tmaxfile2,ddir,line1,line2):
	header = fits.getdata(tmaxfile1, header=True)[1]
	tmax2=pyfits.getdata(tmaxfile2)
	tmax1=pyfits.getdata(tmaxfile1)
	badmask2=tmax2<=0
	badmask1=tmax1<=0
	keep2=tmax2*~badmask2
	keep1=tmax1*~badmask1
	ratio=np.array(keep1)/np.array(keep2)
	fits.writeto(filename=ddir+'ratio_maps/'+line1+'_'+line2+'_ratiomap_tmax.fits',output_verify='ignore',\
	clobber=True,data=ratio,header=header)



#make same number of channels for both!!!!
line1='13CO'
line2='18CO'
line_ratio(datadir+'alex_imaging_'+line1+'_fix/GRS1915_modelimg_'+line1+'.image.pbcor.fits',\
datadir+'alex_imaging_'+line2+'_fix/GRS1915_modelimg_'+line2+'.image.pbcor.fits',datadir,3.,2.,2.,line1,line2)
line_ratio_tmax(datadir+'T_max_maps/'+line1+'_tmax.fits',\
datadir+'T_max_maps/'+line2+'_tmax.fits',datadir,line1,line2)
raw_input('go')
#view
fig=plt.figure()
line1='13CO'
line2='18CO'
filea=pyfits.getdata(datadir+'ratio_maps/'+line1+'_'+line2+'_ratiomap_tmax.fits')
a=plt.imshow(filea)
cbar=fig.colorbar(a)#ticks=[0,0.5,1,1.5,2.])
plt.clim(0,10)
cbar.set_label('13CO/18CO')
plt.gca().invert_yaxis()
plt.show()

#ra/dec version
fits_file1=datadir+'ratio_maps/'+line1+'_'+line2+'_ratiomap_tmax.fits'
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
im=plt.imshow(np.nan_to_num(data1),origin="lower",cmap=cm.get_cmap('jet', 500),norm=colors.PowerNorm(gamma=1.),vmin=0.0,vmax=10)
cbar=plt.colorbar(im, orientation='vertical',fraction=0.04,pad=0)
cbar.set_label(line1+'/'+line2)
ax1.tick_params(axis='both', which='major', labelsize=15,width=3,length=7,color='k')
ax1.tick_params(axis='both', which='minor', labelsize=15,width=1,length=7,color='k')
ax1.coords['ra'].set_axislabel('Right Ascension')
ax1.coords['dec'].set_axislabel('Declination',minpad=-0.1)
ax1.coords['ra'].set_major_formatter('hh:mm:ss.s')
ax1.set_ylim(150, 700)
ax1.set_xlim(100, 650)
plt.contour(X,Y,Z,levels,colors='k',transform=ax1.get_transform(wmap))
plt.savefig(datadir+'ratio_maps/'+line1+'_'+line2+'_tmax'+'_contour.pdf',bbox_inches='tight')
plt.show()

#import os
#os.system('cp -r '+datadir+'ratio_maps/'+line1+'_'+line2+'_tmax'+'_contour.pdf /home/ubuntu/Dropbox')
	
