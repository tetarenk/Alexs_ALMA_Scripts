from spectral_cube import SpectralCube
import numpy as np
from astropy.io import fits
from astropy import units as u
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
import glob

datadir = '/mnt/bigdata/tetarenk/ALMA_GRS1915_105/'

'''imfiles=[]
for name in glob.glob(datadir+'other_data/*_regrid.fits'):
	if 'hspire' in name:
		imfiles.append(name.strip('.fits')+'_fix.fits')
	elif 'ukids' in name:
		imfiles.append(name.strip('.fits')+'_fix.fits')
	elif 'tay2' in name:
		continue
	else:
    	imfiles.append(name)'''
#VLA contours
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
x=np.arange(0,len(data[0,:]))
y=np.arange(0,len(data[:,0]))
X, Y = np.meshgrid(x, y)
Z=data#[0,0,490:550,470:550]
levels=np.array([1,2,3,4,5,6,7])*0.000345

'''imfiles=['/mnt/bigdata/tetarenk/ALMA_GRS1915_105/other_data/residual_MG0450n005_024_all_mipsgal_regrid.fits',
 '/mnt/bigdata/tetarenk/ALMA_GRS1915_105/other_data/hspireplw529_herschel_regrid_fix.fits',
 '/mnt/bigdata/tetarenk/ALMA_GRS1915_105/other_data/hspirepmw529_herschel_regrid_fix.fits',
 '/mnt/bigdata/tetarenk/ALMA_GRS1915_105/other_data/hspirepsw529_herschel_regrid_fix.fits',
 '/mnt/bigdata/tetarenk/ALMA_GRS1915_105/other_data/v2.1_ds2_l045_13pca_medmap20_crop_bolocam_regrid.fits',
 '/mnt/bigdata/tetarenk/ALMA_GRS1915_105/other_data/GLM_04550-0075_mosaic_I3_glimpse_regrid.fits',
 '/mnt/bigdata/tetarenk/ALMA_GRS1915_105/other_data/GLM_04550-0075_mosaic_I4_glimpse_regrid.fits',
 '/mnt/bigdata/tetarenk/ALMA_GRS1915_105/other_data/GLM_04550-0075_mosaic_I1_glimpse_regrid.fits',
 '/mnt/bigdata/tetarenk/ALMA_GRS1915_105/other_data/GLM_04550-0075_mosaic_I2_glimpse_regrid.fits',
 '/mnt/bigdata/tetarenk/ALMA_GRS1915_105/other_data/236_471_1_1424_2_ukidssh_regrid_fix.fits',
 '/mnt/bigdata/tetarenk/ALMA_GRS1915_105/other_data/236_471_1_1424_3_ukidssk_regrid_fix.fits',
 '/mnt/bigdata/tetarenk/ALMA_GRS1915_105/other_data/236_471_1_1424_1_ukidssj_regrid_fix.fits']'''

imfiles=['/mnt/bigdata/tetarenk/ALMA_GRS1915_105/other_data/236_471_1_1424_1_ukidssj_regrid_fix.fits',\
'/mnt/bigdata/tetarenk/ALMA_GRS1915_105/other_data/v2.1_ds2_l045_13pca_medmap20_crop_bolocam_regrid.fits',\
'/mnt/bigdata/tetarenk/ALMA_GRS1915_105/other_data/hspireplw529_herschel_regrid_fix.fits',\
 '/mnt/bigdata/tetarenk/ALMA_GRS1915_105/other_data/hspirepmw529_herschel_regrid_fix.fits',\
 '/mnt/bigdata/tetarenk/ALMA_GRS1915_105/other_data/hspirepsw529_herschel_regrid_fix.fits',\
'/mnt/bigdata/tetarenk/ALMA_GRS1915_105/other_data/residual_MG0450n005_024_all_mipsgal_regrid.fits',\
'/mnt/bigdata/tetarenk/ALMA_GRS1915_105/other_data/GLM_04550-0075_mosaic_I4_glimpse_regrid.fits',\
'/mnt/bigdata/tetarenk/ALMA_GRS1915_105/other_data/GLM_04550-0075_mosaic_I3_glimpse_regrid.fits',\
'/mnt/bigdata/tetarenk/ALMA_GRS1915_105/other_data/GLM_04550-0075_mosaic_I2_glimpse_regrid.fits',\
'/mnt/bigdata/tetarenk/ALMA_GRS1915_105/other_data/GLM_04550-0075_mosaic_I1_glimpse_regrid.fits',\
'/mnt/bigdata/tetarenk/ALMA_GRS1915_105/other_data/236_471_1_1424_3_ukidssk_regrid_fix.fits',\
'/mnt/bigdata/tetarenk/ALMA_GRS1915_105/other_data/236_471_1_1424_2_ukidssh_regrid_fix.fits']
label2=['UKIDSS','BOLOCAM','HSPIRE','HSPIRE','HSPIRE','MIPSGAL','GLIMPSE','GLIMPSE','GLIMPSE','GLIMPSE',\
'UKIDSS','UKIDSS']
label3=['UKIRT','CSO','Herschel','Herschel','Herschel','Spitzer','Spitzer','Spitzer','Spitzer','Spitzer',\
'UKIRT','UKIRT']
label=['$1.33\\mu m$','$1.11mm$', '$500\\mu m$','$350\\mu m$','$250\\mu m$','$24\\mu m$','$8\\mu m$',\
'$5.8\\mu m$','$4.5\\mu m$','$3.6\\mu m$','$2.37\\mu m$','$1.78\\mu m$']

'''label=['$24\\mu m$','$500\\mu m$','$350\\mu m$','$250\\mu m$','$1.11mm$','$5.8\\mu m$',\
'$8\\mu m$','$3.6\\mu m$','$4.5\\mu m$','$1.78\\mu m$','$2.37\\mu m$','$1.33\\mu m$']
label2=['MIPSGAL','HSPIRE','HSPIRE','HSPIRE','Bolocam','GLIMPSE','GLIMPSE','GLIMPSE','GLIMPSE',\
'UKIDSS','UKIDSS','UKIDSS']
label3=['Spitzer','Herschel','Herschel','Herschel','CSO','Spitzer','Spitzer','Spitzer','Spitzer',\
'UKIRT','UKIRT','UKIRT']'''
fits_file1='/mnt/bigdata/tetarenk/VLA_grs1915_images/GRSVLA.fits'
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
x=np.arange(0,len(data[0,:]))
y=np.arange(0,len(data[:,0]))
X, Y = np.meshgrid(x, y)
Z=data#[0,0,490:550,470:550]
levels=np.array([1,2,3,4,5,6,7])*0.000345

fig=plt.figure(figsize=(17,10))
#plt.rcdefaults()
plt.rc('xtick.major', size=4)
plt.rc('xtick', color='w', labelsize='small')
for i in range(0,len(imfiles)):
	fits_file1=imfiles[i]
	hdulist1 = fits.open(fits_file1)[0]
	if i == 2:
		data1=hdulist1.data*(9.382805e-7)*(35.2**2)/(0.2**2)
	elif i ==3:
		data1=hdulist1.data*(9.382805e-7)*(23.9**2)/(0.2**2)
	elif i ==4:
		data1=hdulist1.data*(9.382805e-7)*(17.6**2)/(0.2**2)
	else:
		data1=hdulist1.data
	wmap1=wcs.WCS(hdulist1.header)
	coord0=SkyCoord('19h15m41.8s','+10d40m38s',frame='icrs')
	coord1=SkyCoord('19h15m34.9s','+10d41m52s',frame='icrs')
	x1=float(wmap1.wcs_world2pix(coord0.ra.value,coord0.dec.value,1)[0])
	y1=float(wmap1.wcs_world2pix(coord0.ra.value,coord0.dec.value,1)[1])
	x2=float(wmap1.wcs_world2pix(coord1.ra.value,coord1.dec.value,1)[0])
	y2=float(wmap1.wcs_world2pix(coord1.ra.value,coord1.dec.value,1)[1])
	ax=plt.subplot(3,4,i, projection=wmap1.celestial)
	im=plt.imshow(np.nan_to_num(data1),origin="lower",cmap=cm.get_cmap('jet', 500),norm=colors.PowerNorm(gamma=1.),vmin=0.0)
	cbar=plt.colorbar(im, orientation='vertical',fraction=0.035,pad=0)
	cbar.ax.tick_params(labelsize=6)
	plt.contour(X,Y,Z,levels,colors='w',transform=ax.get_transform(wmap))
#cbar.set_label('Jy')
	plt.tick_params(axis='both', which='major', labelsize=1,width=3,length=7,color='k')
	plt.tick_params(axis='both', which='minor', labelsize=1,width=1,length=7,color='k')
	if i in [2,3,4,6,7,8]:
		ax.coords['ra'].set_axislabel('')
		ax.coords['dec'].set_axislabel('',minpad=-0.1)
		ax.coords['ra'].set_major_formatter('hh:mm:ss.s')
		ax.coords['ra'].set_ticklabel_visible(False)
		ax.coords['dec'].set_ticklabel_visible(False)
	elif i in [10,11,0]:
		ax.coords['ra'].set_axislabel('Right Ascension')
		ax.coords['dec'].set_axislabel('',minpad=-0.1)
		ax.coords['ra'].set_major_formatter('hh:mm:ss.s')
		ax.coords['ra'].set_ticklabel_visible(True)
		ax.coords['dec'].set_ticklabel_visible(False)
	elif i ==9:
		ax.coords['ra'].set_axislabel('Right Ascension')
		ax.coords['dec'].set_axislabel('Declination',minpad=-0.1)
		ax.coords['ra'].set_major_formatter('hh:mm:ss.s')
		ax.coords['ra'].set_ticklabel_visible(True)
		ax.coords['dec'].set_ticklabel_visible(True)
		
	elif i in [1,5]:
		ax.coords['ra'].set_axislabel('')
		ax.coords['dec'].set_axislabel('Declination',minpad=-0.1)
		ax.coords['ra'].set_major_formatter('hh:mm:ss.s')
		ax.coords['ra'].set_ticklabel_visible(False)
		ax.coords['dec'].set_ticklabel_visible(True)
	ax.set_ylim(y1,y2)
	ax.set_xlim(x1,x2)
	if i in [2,3]:
		ax.text(30,475,label[i],color='k')
		ax.text(350,475,label2[i],color='k')
		ax.text(350,445,label3[i],color='k')
	else:
		ax.text(30,475,label[i],color='w')
		ax.text(350,475,label2[i],color='w')
		ax.text(350,445,label3[i],color='w')
	#ax.set_aspect('equal')
	#ax.set_xticklabels([])
	#ax.set_yticklabels([])
	#plt.setp(ax.get_yticklabels(),visible=False)
plt.subplots_adjust(wspace=0.2,hspace=0.1,top=0.9,bottom=0.1)
plt.savefig(datadir+'other_data/'+'all'+'_contour.pdf',bbox_inches='tight')
plt.show()
raw_input('stop')
#vla cont
fits_file1='/mnt/bigdata/tetarenk/VLA_grs1915_images/grs1915_rob_ms_nopb_tay2_redo_dec16.image.tt0.fits'
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
coord0=SkyCoord('19h15m41.3s','+10d41m01s',frame='icrs')
coord1=SkyCoord('19h15m36.9s','+10d41m49s',frame='icrs')
x1=float(wmap1.wcs_world2pix(coord0.ra.value,coord0.dec.value,0,0,1)[0])
y1=float(wmap1.wcs_world2pix(coord0.ra.value,coord0.dec.value,0,0,1)[1])
x2=float(wmap1.wcs_world2pix(coord1.ra.value,coord1.dec.value,0,0,1)[0])
y2=float(wmap1.wcs_world2pix(coord1.ra.value,coord1.dec.value,0,0,1)[1])
coord3=SkyCoord('19h15m37.1s','+10d41m04s',frame='icrs')
x3=float(wmap1.wcs_world2pix(coord3.ra.value,coord3.dec.value,0,0,1)[0])
y3=float(wmap1.wcs_world2pix(coord3.ra.value,coord3.dec.value,0,0,1)[1])


x=np.arange(0,len(data[0,:]))
y=np.arange(0,len(data[:,0]))
X, Y = np.meshgrid(x, y)
Z=data#[0,0,490:550,470:550]
levels=np.array([1,2,3,4,5,6,7])*0.000345
#evels=np.array([1,2,3,4,5,6,7])*0.000345
#sns.set_style("dark")
#cmap1 = mpl.colors.ListedColormap(sns.color_palette("colorblind",10))
#cmap2=colors_maps()
fig=plt.figure()
#plt.rcdefaults()
plt.rc('xtick.major', size=4)
plt.rc('xtick', color='w', labelsize='large')
ax1 = fig.add_subplot(111, projection=wmap1.celestial)
im=plt.imshow(np.nan_to_num(data1[0,0,:,:]),origin="lower",cmap=cm.get_cmap('jet', 500),norm=colors.PowerNorm(gamma=0.7),vmin=0.0)#65,55,65,0.9,0.9,0.9,0.9/x,x,x,0.1,0.03,0.025,0.025,0.025,0.1
cbar=plt.colorbar(im, orientation='vertical',fraction=0.0358,pad=0)
cbar.set_label('mJy/beam')
ax1.tick_params(axis='both', which='major', labelsize=15,width=3,length=7,color='k')
ax1.tick_params(axis='both', which='minor', labelsize=15,width=1,length=7,color='k')
ax1.coords['ra'].set_axislabel('Right Ascension')
ax1.coords['dec'].set_axislabel('Declination',minpad=-0.1)
ax1.coords['ra'].set_major_formatter('hh:mm:ss.s')
ax1.set_ylim(y1,y2)
ax1.set_xlim(x1,x2)
ax1.text(470,550,'4-8 GHz',color='w')
e1 = patches.Ellipse((x3,y3), 5.31, 4.50,angle=40, linewidth=2, fill=False,color='m')
ax1.add_patch(e1)
plt.contour(X,Y,Z,levels,colors='k',transform=ax1.get_transform(wmap))
plt.savefig(datadir+'other_data/vla_contour.pdf',bbox_inches='tight')
plt.show()

#spectralindex map
fits_file1='/mnt/bigdata/tetarenk/VLA_grs1915_images/grs1915_rob_ms_nopb_tay2_redo_dec16.image.alpha.fits'
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
coord0=SkyCoord('19h15m41.3s','+10d41m01s',frame='icrs')
coord1=SkyCoord('19h15m36.9s','+10d41m49s',frame='icrs')
x1=float(wmap1.wcs_world2pix(coord0.ra.value,coord0.dec.value,0,0,1)[0])
y1=float(wmap1.wcs_world2pix(coord0.ra.value,coord0.dec.value,0,0,1)[1])
x2=float(wmap1.wcs_world2pix(coord1.ra.value,coord1.dec.value,0,0,1)[0])
y2=float(wmap1.wcs_world2pix(coord1.ra.value,coord1.dec.value,0,0,1)[1])
coord3=SkyCoord('19h15m37.1s','+10d41m04s',frame='icrs')
x3=float(wmap1.wcs_world2pix(coord3.ra.value,coord3.dec.value,0,0,1)[0])
y3=float(wmap1.wcs_world2pix(coord3.ra.value,coord3.dec.value,0,0,1)[1])


x=np.arange(0,len(data[0,:]))
y=np.arange(0,len(data[:,0]))
X, Y = np.meshgrid(x, y)
Z=data#[0,0,490:550,470:550]
levels=np.array([1,2,3,4,5,6,7])*0.000345
#evels=np.array([1,2,3,4,5,6,7])*0.000345
#sns.set_style("dark")
#cmap1 = mpl.colors.ListedColormap(sns.color_palette("colorblind",10))
#cmap2=colors_maps()
fig=plt.figure()
#plt.rcdefaults()
plt.rc('xtick.major', size=4)
plt.rc('xtick', color='w', labelsize='large')
ax1 = fig.add_subplot(111, projection=wmap1.celestial)
im=plt.imshow((data1[0,0,:,:]),origin="lower",cmap=cm.get_cmap('jet', 500),vmin=-1,vmax=0.8)#65,55,65,0.9,0.9,0.9,0.9/x,x,x,0.1,0.03,0.025,0.025,0.025,0.1
cbar=plt.colorbar(im, orientation='vertical',fraction=0.0357,pad=0)
cbar.set_label('Spectral Index')
ax1.tick_params(axis='both', which='major', labelsize=15,width=3,length=7,color='k')
ax1.tick_params(axis='both', which='minor', labelsize=15,width=1,length=7,color='k')
ax1.coords['ra'].set_axislabel('Right Ascension')
ax1.coords['dec'].set_axislabel('Declination',minpad=-0.1)
ax1.coords['ra'].set_major_formatter('hh:mm:ss.s')
ax1.set_ylim(y1,y2)
ax1.set_xlim(x1,x2)
ax1.text(470,550,'4-8 GHz',color='k')
e1 = patches.Ellipse((x3,y3), 5.31, 4.50,angle=40, linewidth=2, fill=False,color='m')
ax1.add_patch(e1)
plt.contour(X,Y,Z,levels,colors='k',transform=ax1.get_transform(wmap))
plt.savefig(datadir+'other_data/vla_alpha_contour.pdf',bbox_inches='tight')
plt.show()
