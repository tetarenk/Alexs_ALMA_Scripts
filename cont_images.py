'''plot alma continuum images, primary beam response, and field map'''

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
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredEllipse


datadir = '/mnt/bigdata/tetarenk/ALMA_GRS1915_105/'

#12m
fits_file1=datadir+'alex_imaging_contspw_fix_redo/GRS1915_12m_cont.image.pbcor.fits'
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
coord0=SkyCoord('19h15m40.8s','+10d40m58s',frame='icrs')
coord1=SkyCoord('19h15m37.1s','+10d41m46s',frame='icrs')
x1=float(wmap1.wcs_world2pix(coord0.ra.value,coord0.dec.value,0,0,1)[0])
y1=float(wmap1.wcs_world2pix(coord0.ra.value,coord0.dec.value,0,0,1)[1])
x2=float(wmap1.wcs_world2pix(coord1.ra.value,coord1.dec.value,0,0,1)[0])
y2=float(wmap1.wcs_world2pix(coord1.ra.value,coord1.dec.value,0,0,1)[1])
coord3=SkyCoord('19h15m37.5s','+10d41m04s',frame='icrs')
x3=float(wmap1.wcs_world2pix(coord3.ra.value,coord3.dec.value,0,0,1)[0])
y3=float(wmap1.wcs_world2pix(coord3.ra.value,coord3.dec.value,0,0,1)[1])


x=np.arange(0,len(data[0,:]))
y=np.arange(0,len(data[:,0]))
X, Y = np.meshgrid(x, y)
Z=data
levels=np.array([4,6,8,10,15,20,40,60])*0.00005
fig=plt.figure()
#plt.rcdefaults()
plt.rc('xtick.major', size=4)
plt.rc('xtick', color='w', labelsize='large')
ax1 = fig.add_subplot(111, projection=wmap1.celestial)
im=plt.imshow(np.nan_to_num(data1[0,0,:,:])*1000.,origin="lower",cmap=cm.get_cmap('hot_r', 500),norm=colors.PowerNorm(gamma=0.7),vmin=0.0,vmax=2)#65,55,65,0.9,0.9,0.9,0.9/x,x,x,0.1,0.03,0.025,0.025,0.025,0.1
cbar=plt.colorbar(im, orientation='vertical',fraction=0.04,pad=0)
cbar.set_label('mJy/beam')
ax1.tick_params(axis='both', which='major', labelsize=15,width=3,length=7,color='k')
ax1.tick_params(axis='both', which='minor', labelsize=15,width=1,length=7,color='k')
ax1.coords['ra'].set_axislabel('Right Ascension')
ax1.coords['dec'].set_axislabel('Declination',minpad=-0.1)
ax1.coords['ra'].set_major_formatter('hh:mm:ss.s')
ax1.set_ylim(y1,y2)
ax1.set_xlim(x1,x2)
e1 = patches.Ellipse((x3,y3), 6.5, 6.0,angle=-145.2, linewidth=2, fill=False,color='m')
#ax1.add_patch(e1)
ae = AnchoredEllipse(ax1.transData, width=6.5, height=6.0, angle=-145.2,loc=4, pad=0.5, borderpad=0.4, frameon=True)
ax1.add_artist(ae)
plt.contour(X,Y,Z,levels,colors='k',transform=ax1.get_transform(wmap),lw=8)
plt.savefig(datadir+'other_data/cont12_contour.pdf',bbox_inches='tight')
plt.show()

#7m
fits_file1=datadir+'alex_imaging_contspw_fix_redo/GRS1915_7m_cont.image.pbcor.fits'
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
coord0=SkyCoord('19h15m40.8s','+10d40m58s',frame='icrs')
coord1=SkyCoord('19h15m37.1s','+10d41m46s',frame='icrs')
x1=float(wmap1.wcs_world2pix(coord0.ra.value,coord0.dec.value,0,0,1)[0])
y1=float(wmap1.wcs_world2pix(coord0.ra.value,coord0.dec.value,0,0,1)[1])
x2=float(wmap1.wcs_world2pix(coord1.ra.value,coord1.dec.value,0,0,1)[0])
y2=float(wmap1.wcs_world2pix(coord1.ra.value,coord1.dec.value,0,0,1)[1])
coord3=SkyCoord('19h15m37.5s','+10d41m04s',frame='icrs')
x3=float(wmap1.wcs_world2pix(coord3.ra.value,coord3.dec.value,0,0,1)[0])
y3=float(wmap1.wcs_world2pix(coord3.ra.value,coord3.dec.value,0,0,1)[1])

x=np.arange(0,len(data[0,:]))
y=np.arange(0,len(data[:,0]))
X, Y = np.meshgrid(x, y)
Z=data
levels=np.array([4,6,8,10,15,20,40,60])*0.00005
fig=plt.figure()
#plt.rcdefaults()
plt.rc('xtick.major', size=4)
plt.rc('xtick', color='w', labelsize='large')
ax1 = fig.add_subplot(111, projection=wmap1.celestial)
im=plt.imshow(np.nan_to_num(data1[0,0,:,:])*1000.,origin="lower",cmap=cm.get_cmap('hot_r', 500),norm=colors.PowerNorm(gamma=1),vmin=0.0)#65,55,65,0.9,0.9,0.9,0.9/x,x,x,0.1,0.03,0.025,0.025,0.025,0.1
cbar=plt.colorbar(im, orientation='vertical',fraction=0.04,pad=0)
cbar.set_label('mJy/beam')
ax1.tick_params(axis='both', which='major', labelsize=15,width=3,length=7,color='k')
ax1.tick_params(axis='both', which='minor', labelsize=15,width=1,length=7,color='k')
ax1.coords['ra'].set_axislabel('Right Ascension')
ax1.coords['dec'].set_axislabel('Declination',minpad=-0.1)
ax1.coords['ra'].set_major_formatter('hh:mm:ss.s')
ax1.set_ylim(y1,y2)
ax1.set_xlim(x1,x2)
e1 = patches.Ellipse((x3,y3), 8.61, 5.14,angle=-164.04, linewidth=2, fill=False,color='m')
#ax1.add_patch(e1)
ae = AnchoredEllipse(ax1.transData, width=8.61, height=5.14, angle=-164.04,loc=4, pad=0.5, borderpad=0.4, frameon=True)
ax1.add_artist(ae)
plt.contour(X,Y,Z,levels,colors='k',transform=ax1.get_transform(wmap))
plt.savefig(datadir+'other_data/cont7_contour.pdf',bbox_inches='tight')
plt.show()

#12m+7m
fits_file1=datadir+'alex_imaging_contspw_fix_redo/GRS1915_12m7m_cont.image.pbcor.fits'
fits_file1a=datadir+'alex_imaging_contspw_fix_redo/GRSALMA.fits'
#convert image so only 2 dimensional
'''filename=datadir+'alex_imaging_contspw_fix_redo/GRS1915_12m7m_cont.image.pbcor.fits'
fh = fits.open(filename)
data = fh[0].data.squeeze() # drops the size-1 axes
header = fh[0].header
mywcs = wcs.WCS(header)
new_header = mywcs.to_header()
new_fh = fits.PrimaryHDU(data=data, header=new_header)
new_fh.writeto(datadir+'alex_imaging_contspw_fix_redo/GRSALMA.fits')'''
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
hdulist1a = fits.open(fits_file1a)[0]
data1a=hdulist1a.data
hdulist1a.header.remove('CRPIX3')
hdulist1a.header.remove('CRVAL3')
hdulist1a.header.remove('CDELT3')
hdulist1a.header.remove('CUNIT3')
hdulist1a.header.remove('CTYPE3')
hdulist1a.header.remove('CRPIX4')
hdulist1a.header.remove('CRVAL4')
hdulist1a.header.remove('CDELT4')
#hdulist.header.remove('CUNIT4')
hdulist1a.header.remove('CTYPE4')
hdulist1a.header['WCSAXES']=2
wmap1a=wcs.WCS(hdulist1a.header)
coord0=SkyCoord('19h15m40.8s','+10d40m58s',frame='icrs')
coord1=SkyCoord('19h15m37.1s','+10d41m46s',frame='icrs')
x1=float(wmap1.wcs_world2pix(coord0.ra.value,coord0.dec.value,0,0,1)[0])
y1=float(wmap1.wcs_world2pix(coord0.ra.value,coord0.dec.value,0,0,1)[1])
x2=float(wmap1.wcs_world2pix(coord1.ra.value,coord1.dec.value,0,0,1)[0])
y2=float(wmap1.wcs_world2pix(coord1.ra.value,coord1.dec.value,0,0,1)[1])
coord3=SkyCoord('19h15m37.2s','+10d41m00s',frame='icrs')
x3=float(wmap1.wcs_world2pix(coord3.ra.value,coord3.dec.value,0,0,1)[0])
y3=float(wmap1.wcs_world2pix(coord3.ra.value,coord3.dec.value,0,0,1)[1])

x=np.arange(0,len(data[0,:]))
y=np.arange(0,len(data[:,0]))
X, Y = np.meshgrid(x, y)
Z=data
levels=np.array([4,6,8,10,15,20,40,60])*0.00005
x00=np.arange(0,len(data1a[0,:]))
y00=np.arange(0,len(data1a[:,0]))
X00, Y00 = np.meshgrid(x00, y00)
Z00=data1a
levels00=np.array([2,3,4])*0.00025
fig=plt.figure()
plt.rcdefaults()
plt.rc('xtick.major', size=4)
#plt.rc('xtick', color='w', labelsize='large')
ax1 = fig.add_subplot(111, projection=wmap1.celestial)
im=plt.imshow(np.nan_to_num(data1[0,0,:,:])*1000.,origin="lower",cmap=cm.get_cmap('hot_r', 500),norm=colors.PowerNorm(gamma=1),vmin=0.0,vmax=0.002*1000.)#65,55,65,0.9,0.9,0.9,0.9/x,x,x,0.1,0.03,0.025,0.025,0.025,0.1
cbar=plt.colorbar(im, orientation='vertical',fraction=0.04,pad=0)
cbar.set_label('mJy/beam')
ax1.tick_params(axis='both', which='major', labelsize=15,width=3,length=7,color='k')
ax1.tick_params(axis='both', which='minor', labelsize=15,width=1,length=7,color='k')
ax1.coords['ra'].set_axislabel('Right Ascension')
ax1.coords['dec'].set_axislabel('Declination',minpad=-0.1)
ax1.coords['ra'].set_major_formatter('hh:mm:ss.s')
ax1.set_ylim(y1,y2)
ax1.set_xlim(x1,x2)
e1 = patches.Ellipse((x3,y3), 6.78, 6.36,angle=-150.02356, linewidth=2, fill=False,color='m')
ax1.add_patch(e1)
ae = AnchoredEllipse(ax1.transData, width=6.78, height=6.36, angle=-150.02356,loc=4, pad=0.5, borderpad=0.4, frameon=True)
ax1.add_artist(ae)
plt.contour(X,Y,Z,levels,colors='k',transform=ax1.get_transform(wmap))
#plt.contour(X00,Y00,Z00,levels00,colors='c',transform=ax1.get_transform(wmap1a))
plt.savefig(datadir+'for_paper/cont127_contour.pdf',bbox_inches='tight')
plt.show()

#PB noise mp

fits_file1=datadir+'alex_imaging_12CO_fix/CLEAN12CO_combinewithtp_flux.fits'
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


coord0=SkyCoord('19h15m40.8s','+10d40m58s',frame='icrs')
coord1=SkyCoord('19h15m37.1s','+10d41m46s',frame='icrs')
x1=float(wmap1.wcs_world2pix(coord0.ra.value,coord0.dec.value,0,0,1)[0])
y1=float(wmap1.wcs_world2pix(coord0.ra.value,coord0.dec.value,0,0,1)[1])
x2=float(wmap1.wcs_world2pix(coord1.ra.value,coord1.dec.value,0,0,1)[0])
y2=float(wmap1.wcs_world2pix(coord1.ra.value,coord1.dec.value,0,0,1)[1])
coord3=SkyCoord('19h15m37.2s','+10d41m00s',frame='icrs')
x3=float(wmap1.wcs_world2pix(coord3.ra.value,coord3.dec.value,0,0,1)[0])
y3=float(wmap1.wcs_world2pix(coord3.ra.value,coord3.dec.value,0,0,1)[1])

x=np.arange(0,len(data[0,:]))
y=np.arange(0,len(data[:,0]))
X, Y = np.meshgrid(x, y)
Z=data#[0,0,490:550,470:550]
levels=np.array([4,6,8,10,15,20,40,60])*0.00005
fig=plt.figure()
plt.rcdefaults()
plt.rc('xtick.major', size=4)
ax1 = fig.add_subplot(111, projection=wmap1.celestial)
im=plt.imshow(np.nan_to_num(data1[0,0,:,:]),origin="lower",cmap=cm.get_cmap('hot_r', 500),norm=colors.PowerNorm(gamma=1),vmin=0.0,vmax=1)#65,55,65,0.9,0.9,0.9,0.9/x,x,x,0.1,0.03,0.025,0.025,0.025,0.1
cbar=plt.colorbar(im, orientation='vertical',fraction=0.04,pad=0)
cbar.set_label('Primary Beam Response')
ax1.tick_params(axis='both', which='major', labelsize=15,width=3,length=7,color='k')
ax1.tick_params(axis='both', which='minor', labelsize=15,width=1,length=7,color='k')
ax1.coords['ra'].set_axislabel('Right Ascension')
ax1.coords['dec'].set_axislabel('Declination',minpad=-0.1)
ax1.coords['ra'].set_major_formatter('hh:mm:ss.s')
ax1.set_ylim(150, 700)
ax1.set_xlim(100, 650)
plt.contour(X,Y,Z,levels,colors='w',transform=ax1.get_transform(wmap))
plt.savefig(datadir+'for_paper/PB_noise_map.pdf',bbox_inches='tight')
plt.show()


#ALMA field map

fits_file1=datadir+'alex_imaging_12CO_fix/CLEAN12CO_combinewithtp_flux.fits'
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


coord0=SkyCoord('19h15m40.8s','+10d40m58s',frame='icrs')
coord1=SkyCoord('19h15m37.1s','+10d41m46s',frame='icrs')
x1=float(wmap1.wcs_world2pix(coord0.ra.value,coord0.dec.value,0,0,1)[0])
y1=float(wmap1.wcs_world2pix(coord0.ra.value,coord0.dec.value,0,0,1)[1])
x2=float(wmap1.wcs_world2pix(coord1.ra.value,coord1.dec.value,0,0,1)[0])
y2=float(wmap1.wcs_world2pix(coord1.ra.value,coord1.dec.value,0,0,1)[1])
coord3=SkyCoord('19h15m37.2s','+10d41m00s',frame='icrs')
x3=float(wmap1.wcs_world2pix(coord3.ra.value,coord3.dec.value,0,0,1)[0])
y3=float(wmap1.wcs_world2pix(coord3.ra.value,coord3.dec.value,0,0,1)[1])

x=np.arange(0,len(data[0,:]))
y=np.arange(0,len(data[:,0]))
X, Y = np.meshgrid(x, y)
Z=data#[0,0,490:550,470:550]
levels=np.array([4,6,8,10,15,20,40,60])*0.00005
fig=plt.figure()
plt.rcdefaults()
plt.rc('xtick.major', size=4)
ax1 = fig.add_subplot(111, projection=wmap1.celestial)
im=plt.imshow(np.nan_to_num(data1[0,0,:,:]),origin="lower",cmap=cm.get_cmap('Purples', 500),norm=colors.PowerNorm(gamma=1),vmin=0.0,vmax=0.1)#65,55,65,0.9,0.9,0.9,0.9/x,x,x,0.1,0.03,0.025,0.025,0.025,0.1
#cbar=plt.colorbar(im, orientation='vertical',fraction=0.04,pad=0)
#cbar.set_label('Primary Beam Response')
ax1.tick_params(axis='both', which='major', labelsize=15,width=3,length=7,color='k')
ax1.tick_params(axis='both', which='minor', labelsize=15,width=1,length=7,color='k')
ax1.coords['ra'].set_axislabel('Right Ascension')
ax1.coords['dec'].set_axislabel('Declination',minpad=-0.1)
ax1.coords['ra'].set_major_formatter('hh:mm:ss.s')
ax1.set_ylim(150, 700)
ax1.set_xlim(100, 650)
plt.contour(X,Y,Z,levels,colors='w',transform=ax1.get_transform(wmap),linewidths=2)
plt.savefig(datadir+'for_paper/ALMA_field_map.pdf',bbox_inches='tight')
plt.show()
