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
	ratio=Y_masked/X_masked
	fits.writeto(filename=ddir+'ratio_maps/'+line2+'_'+line1+'_ratiomap.fits',output_verify='ignore',\
	clobber=True,data=ratio,header=header)
def line_ratio_tmax(file1,file2,ddir,line1,line2,chan1):
	header = fits.getdata(file1, header=True)[1]
	t2=pyfits.getdata(file2)
	t1=pyfits.getdata(file1)
	#badmask2=t2<=0
	#badmask1=t1<=0
	#keep2=t2*~badmask2
	#keep1=t1*~badmask1
	t1max=np.nanmax(t1)
	print t1max
	maxvchan=np.nanargmax(t1)
	maxind=np.unravel_index(maxvchan, t1.shape)[1]
	print maxind,t1max.shape
	plt.imshow(t1max)
	plt.show()
	t2max=t2[0,maxind,:,:]
	#t1max=t1[0,chan1,:,:]
	#t2max=t2[0,chan1,:,:]
	ratio=t2max/t1max
	fits.writeto(filename=ddir+'ratio_maps/'+line2+'_'+line1+'_ratiomap_tmax.fits',output_verify='ignore',\
	clobber=True,data=ratio,header=header)


def ratioss(fitsfile1,fitsfile2,ddir,line1,line2):
	header = fits.getdata(fitsfile1, header=True)[1]
	X=SpectralCube.read(fitsfile1)
	Y=SpectralCube.read(fitsfile2)
	Xkms=X.with_spectral_unit(u.km/u.s,velocity_convention='radio')
	Ykms=Y.with_spectral_unit(u.km/u.s,velocity_convention='radio')
	XK=Xkms.to(u.K, equivalencies=Xkms.beam.jtok_equiv(Xkms.header['RESTFRQ']*u.Hz))
	YK=Ykms.to(u.K, equivalencies=Ykms.beam.jtok_equiv(Ykms.header['RESTFRQ']*u.Hz))
	Xs=XK.spectral_slab(60*u.km/u.s,80*u.km/u.s)
	Ys=YK.spectral_slab(60*u.km/u.s,80*u.km/u.s)
	Yss=Ys.spectral_interpolate(spectral_grid=Xs.spectral_axis)
	tmax1=Xs.apply_numpy_function(np.nanmax,axis=0)
	plt.figure()	
	im0=plt.imshow(tmax1)
	plt.colorbar(im0)
	plt.show()
	masknum=raw_input('enter mask-->')
	X_masked=Xs.with_mask(Xs>float(masknum)*u.K)
	Xarr=np.array(X_masked)
	Yarr=np.array(Yss)
	ratio=np.empty((Xarr.shape[1],Xarr.shape[2]))
	for i in range(0,Xarr.shape[1]):
		for j in range(0,Xarr.shape[2]):
			data=Xarr[:,i,j]
			maxval1=np.nanmax(data)
			if maxval1==np.nan:
				maxind=0
			else:
				try:
					maxind=np.nanargmax(data)
				except ValueError:
					maxind=0
			maxval2=Yarr[maxind,i,j]
			ratio[i][j]=maxval2/maxval1
	fits.writeto(filename=ddir+'ratio_maps/'+line2+'_'+line1+'_ratiomap_tmax_matched_mask.fits',output_verify='ignore',\
	clobber=True,data=ratio,header=header)
#make same number of channels for both!!!!
line2='12CO'
line1='13CO'
#line_ratio(datadir+'alex_imaging_'+line1+'_fix/GRS1915_modelimg_'+line1+'.image.pbcor.fits',\
#datadir+'alex_imaging_'+line2+'_fix/GRS1915_modelimg_'+line2+'.image.pbcor.fits',datadir,3.,2.,2.,line1,line2)
#line_ratio_tmax(datadir+'alex_imaging_'+line1+'_fix/GRS1915_modelimg_'+line1+'.image.pbcor.fits',\
#datadir+'alex_imaging_'+line2+'_fix/GRS1915_modelimg_'+line2+'.image.pbcor.fits',datadir,line1,line2,24)
ratioss(datadir+'alex_imaging_'+line1+'_fix/GRS1915_modelimg_'+line1+'.image.pbcor.fits',\
datadir+'alex_imaging_'+line2+'_fix/GRS1915_modelimg_'+line2+'.image.pbcor.fits',datadir,line1,line2)
raw_input('go')
#view
fig=plt.figure()
filea=pyfits.getdata(datadir+'ratio_maps/'+line2+'_'+line1+'_ratiomap_tmax_matched_mask.fits')
a=plt.imshow(filea,vmin=0,vmax=10)
cbar=fig.colorbar(a)#ticks=[0,0.5,1,1.5,2.])
#plt.clim(0,2)
cbar.set_label(line2+'/'+line1)
plt.gca().invert_yaxis()
plt.show()
raw_input('go')
#ra/dec version
fits_file1=datadir+'ratio_maps/'+line2+'_'+line1+'_ratiomap_tmax_matched_mask.fits'
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
levels=np.array([4,6,8,10,15,20,40,60])*0.00005
#evels=np.array([1,2,3,4,5,6,7])*0.000345
coord0=SkyCoord('19h15m41.1s','+10d40m57s',frame='icrs')
coord1=SkyCoord('19h15m37.3s','+10d41m49s',frame='icrs')
x1=float(wmap1.wcs_world2pix(coord0.ra.value,coord0.dec.value,0,0,1)[0])
y1=float(wmap1.wcs_world2pix(coord0.ra.value,coord0.dec.value,0,0,1)[1])
x2=float(wmap1.wcs_world2pix(coord1.ra.value,coord1.dec.value,0,0,1)[0])
y2=float(wmap1.wcs_world2pix(coord1.ra.value,coord1.dec.value,0,0,1)[1])
fig=plt.figure()
plt.rcdefaults()
plt.rc('xtick.major', size=4)
#plt.rc('xtick', color='w', labelsize='large')
ax1 = fig.add_subplot(111, projection=wmap1.celestial)
im=plt.imshow(np.nan_to_num(data1),origin="lower",cmap=cm.get_cmap('hot_r', 500),\
norm=colors.PowerNorm(gamma=1.),vmin=1,vmax=3)
cbar=plt.colorbar(im, orientation='vertical',fraction=0.04,pad=0)
#cbar.set_label(line2+'/'+line1)
cbar.set_label('$^{12}$CO/$^{13}$CO')
ax1.tick_params(axis='both', which='major', labelsize=15,width=3,length=7,color='k')
ax1.tick_params(axis='both', which='minor', labelsize=15,width=1,length=7,color='k')
ax1.coords['ra'].set_axislabel('Right Ascension')
ax1.coords['dec'].set_axislabel('Declination',minpad=-0.1)
ax1.coords['ra'].set_major_formatter('hh:mm:ss.s')
#ax1.set_ylim(150, 700)
#ax1.set_xlim(100, 650)
ax1.set_ylim(y1,y2)
ax1.set_xlim(x1,x2)
plt.contour(X,Y,Z,levels,colors='k',transform=ax1.get_transform(wmap))
plt.savefig(datadir+'for_paper/'+line2+'_'+line1+'_tmax'+'_contour_zoomed_mask.pdf',bbox_inches='tight')
plt.show()

#import os
#os.system('cp -r '+datadir+'ratio_maps/'+line1+'_'+line2+'_tmax'+'_contour.pdf /home/ubuntu/Dropbox')
	
