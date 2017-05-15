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
import os


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
	print Xkms.shape
	mom=Xkms[:,:,:].moment(order=mind)
	mom_array=np.array(mom)
	os.system('rm -rf '+ddir+'moment_maps/'+line+'_moment'+str(mind)+flag+'.fits')
	fits.writeto(filename=ddir+'moment_maps/'+line+'_moment'+str(mind)+flag+'.fits',output_verify='ignore',\
	clobber=True,data=mom_array,header=header)


#for i in lines:
	#make_moment_maps(datadir+'alex_imaging_'+i+'/GRS1915_modelimg_'+i+'.image.pbcor.fits',0,datadir,i)
	#make_moment_maps(datadir+'alex_imaging_'+i+'/GRS1915_modelimg_'+i+'.image.pbcor.fits',1,datadir,i)
	#make_moment_maps(datadir+'alex_imaging_'+i+'/GRS1915_modelimg_'+i+'.image.pbcor.fits',2,datadir,i)

i='18CO'
ii='CH3OH_5'
flag='K'
moment=0
make_moment_maps(datadir+'alex_imaging_'+i+'_fix_all/GRS1915_modelimg_'+i+'.image.pbcor.fits',moment,datadir,i,flag)
#make_moment_maps(datadir+'alex_imaging_'+i+'/GRS1915_modelimg_'+i+'.image.pbcor.fits',1,datadir,i)
#make_moment_maps(datadir+'alex_imaging_'+i+'/GRS1915_modelimg_'+i+'.image.pbcor.fits',2,datadir,i)
#imview(datadir+'moment_maps/'+i+'_moment'+str(0)+'.fits')

fits_file1=datadir+'moment_maps/'+i+'_moment'+str(moment)+flag+'.fits'
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

'''#Make mask
#get spectra, into units of K and GHz
tmax=pyfits.getdata(datadir+'T_max_maps/'+i+'_tmax.fits')

#create mask - makes a 2D map- first tmax>=0, then erode to get rid of noisy edges, then sigma cut it
badmask=tmax>=0
#pixels that are NaN in the original data
#erode the badmask edge by thismuch to get rid of edges
rad=10
datam=nd.morphology.binary_erosion(badmask, np.ones((rad, rad)))
#igma cut- 5*6mJy rms=30 mJy or 0.03Jy
keep=(tmax*datam)>1.0
fig=plt.figure()
im=plt.imshow(keep)
plt.colorbar(im)
plt.gca().invert_yaxis()
#plt.savefig(datadir+'T_max_maps/keep321.png')
plt.show()'''

import scipy.ndimage as nd
rad=100#100,125,125,125,125,100
#data_mask=nd.morphology.binary_erosion(np.nan_to_num(data1), np.ones((rad, rad)))
'''tmax=pyfits.getdata(datadir+'moment_maps/'+i+'_moment'+str(0)+'K'+'.fits')
badmask=tmax>=0
#pixels that are NaN in the original data
#erode the badmask edge by thismuch to get rid of edges
rad=75
datam=nd.morphology.binary_erosion(badmask, np.ones((rad, rad)))
#igma cut- 5*6mJy rms=30 mJy or 0.03Jy
keep=(tmax*datam)>60'''


x=np.arange(0,len(data[0,:]))
y=np.arange(0,len(data[:,0]))
X, Y = np.meshgrid(x, y)
Z=data#[0,0,490:550,470:550]
#levels=np.array([1,2,3,4,5,6,7])*0.000345
levels=np.array([4,6,8,10,15,20,40,60])*0.00005
#n=np.array([0.2,0.6,0.8,1,2,3,4,5,6,7,8,9,10])
#levels=2**(n/2.)*0.000345
fig=plt.figure()
#plt.rcdefaults()
plt.rc('xtick.major', size=4)
#plt.rc('xtick', color='w', labelsize='large')
ax1 = fig.add_subplot(111, projection=wmap1.celestial)
im=plt.imshow(np.nan_to_num(data1),origin="lower",cmap=cm.get_cmap('hot_r', 500),vmin=0,vmax=10)#,vmax=2.0)#,vmax=68,norm=colors.PowerNorm(gamma=1.0))#,vmax=10.0)#,vmax=7)#,norm=colors.PowerNorm(gamma=1.5))
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
from matplotlib.patches import Rectangle
r1=Rectangle((385, 560), 15, 15,alpha=1, facecolor='none',edgecolor='white',lw=2,ls='solid',zorder=5)
r2=Rectangle((400, 570), 20, 25,alpha=1, facecolor='none',edgecolor='white',lw=2,ls='solid',zorder=5)
r3=Rectangle((340, 470), 20, 15,alpha=1, facecolor='none',edgecolor='white',lw=2,ls='solid',zorder=5)
r4=Rectangle((507, 447), 73, 45,alpha=1, facecolor='none',edgecolor='white',lw=2,ls='solid',zorder=5)
'''ax1.add_patch(r1)
ax1.add_patch(r2)
ax1.add_patch(r3)
ax1.add_patch(r4)
ax1.text(375,535, 'A',fontsize=15,color='white')
ax1.text(428,575, 'B',fontsize=15,color='white')
ax1.text(310,460, 'C',fontsize=15,color='white')
ax1.text(585,492, 'D',fontsize=15,color='white')'''
plt.contour(X,Y,Z,levels,colors='white',transform=ax1.get_transform(wmap),zorder=4)
#plt.savefig(datadir+'moment_maps/'+i+'_moment'+str(moment)+'_contourK.pdf',bbox_inches='tight')
plt.savefig(datadir+'for_paper/'+i+'_moment'+str(moment)+'_contourK.pdf',bbox_inches='tight')
plt.show()
#raw_input('stop')
#import scipy.ndimage as nd
#rad=150#100,125,125,125,125,100
#data_mask=nd.morphology.binary_erosion(np.nan_to_num(data1), np.ones((rad, rad)))
#import os
#os.system('cp -r '+datadir+'moment_maps/'+i+'_moment'+str(0)+'_contour.pdf /home/ubuntu/Dropbox')

'''#all lines big plot
imfiles=['/mnt/bigdata/tetarenk/ALMA_GRS1915_105/moment_maps/12CO_moment0.fits',\
'/mnt/bigdata/tetarenk/ALMA_GRS1915_105/moment_maps/13CO_moment0.fits',\
'/mnt/bigdata/tetarenk/ALMA_GRS1915_105/moment_maps/18CO_moment0J.fits',\
'/mnt/bigdata/tetarenk/ALMA_GRS1915_105/moment_maps/H2CO_303_202_moment0.fits',\
'/mnt/bigdata/tetarenk/ALMA_GRS1915_105/moment_maps/H2CO_321_220_moment0.fits',\
'/mnt/bigdata/tetarenk/ALMA_GRS1915_105/moment_maps/H2CO_322_221_moment0.fits',\
'/mnt/bigdata/tetarenk/ALMA_GRS1915_105/moment_maps/H30alpha_moment0J.fits',\
'/mnt/bigdata/tetarenk/ALMA_GRS1915_105/moment_maps/N2Dp_moment0.fits',\
'/mnt/bigdata/tetarenk/ALMA_GRS1915_105/moment_maps/SiO_moment0.fits']

label=['$^{12}{\\rm CO}$','$^{13}{\\rm CO}$','$^{18}{\\rm CO}$','${\\rm H}_{2}{\\rm CO}(303-202)$',\
'${\\rm H}_{2}{\\rm CO}(321-220)$','${\\rm H}_{2}{\\rm CO}(322-221)$','${\\rm H}_{30}\\alpha$',\
'${\\rmN}_{2}{\\rm D}+$','$SiO$']

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
minv=[0,0,0,0,0,0,0,0,0]
maxv=[7.2,3.2,1.4,0.28,0.105,0.200,1.5,0.08,0.08]

for i in range(0,len(imfiles)):
	fits_file1=imfiles[i]
	hdulist1 = fits.open(fits_file1)[0]
	data1=hdulist1.data
	wmap1=wcs.WCS(hdulist1.header)
	coord0=SkyCoord('19h15m41.8s','+10d40m38s',frame='icrs')
	coord1=SkyCoord('19h15m34.9s','+10d41m52s',frame='icrs')
	ax=plt.subplot(3,3,i+1, projection=wmap1.celestial)
	if i==3:
		data_mask=nd.morphology.binary_erosion(np.nan_to_num(data1), np.ones((100, 100)))
		data_mask2=nd.morphology.binary_erosion(np.nan_to_num(data1), np.ones((150, 150)))
	if i in [3,4,5,7,8]:
		#data_mask=nd.morphology.binary_erosion(np.nan_to_num(data1), np.ones((100, 100)))
		im=plt.imshow(np.nan_to_num(data1)*data_mask,origin="lower",cmap=cm.get_cmap('jet', 500),norm=colors.PowerNorm(gamma=1.),vmin=minv[i],vmax=maxv[i])
	elif i==6:
		#data_mask=nd.morphology.binary_erosion(np.nan_to_num(data1), np.ones((150, 150)))
		im=plt.imshow(np.nan_to_num(data1)*data_mask2,origin="lower",cmap=cm.get_cmap('jet', 500),norm=colors.PowerNorm(gamma=1.),vmin=minv[i],vmax=maxv[i])	
	else:
		im=plt.imshow(np.nan_to_num(data1),origin="lower",cmap=cm.get_cmap('jet', 500),norm=colors.PowerNorm(gamma=1.),vmin=minv[i],vmax=maxv[i])
	cbar=plt.colorbar(im, orientation='vertical',fraction=0.035,pad=0)
	cbar.ax.tick_params(labelsize=8)
	plt.contour(X,Y,Z,levels,colors='k',transform=ax.get_transform(wmap))
	plt.tick_params(axis='both', which='major', labelsize=1,width=3,length=7,color='k')
	plt.tick_params(axis='both', which='minor', labelsize=1,width=1,length=7,color='k')
	if i in [1,2,4,5]:
		ax.coords['ra'].set_axislabel('')
		ax.coords['dec'].set_axislabel('',minpad=-0.1)
		ax.coords['ra'].set_major_formatter('hh:mm:ss')
		ax.coords['ra'].set_ticklabel_visible(False)
		ax.coords['dec'].set_ticklabel_visible(False)
	elif i in [7,8]:
		ax.coords['ra'].set_axislabel('Right Ascension')
		ax.coords['dec'].set_axislabel('',minpad=-0.1)
		ax.coords['ra'].set_major_formatter('hh:mm:ss')
		ax.coords['ra'].set_ticklabel_visible(True)
		ax.coords['dec'].set_ticklabel_visible(False)
	elif i ==6:
		ax.coords['ra'].set_axislabel('Right Ascension')
		ax.coords['dec'].set_axislabel('Declination',minpad=-0.1)
		ax.coords['ra'].set_major_formatter('hh:mm:ss')
		ax.coords['ra'].set_ticklabel_visible(True)
		ax.coords['dec'].set_ticklabel_visible(True)
		
	elif i in [0,3]:
		ax.coords['ra'].set_axislabel('')
		ax.coords['dec'].set_axislabel('Declination',minpad=-0.1)
		ax.coords['ra'].set_major_formatter('hh:mm:ss.s')
		ax.coords['ra'].set_ticklabel_visible(False)
		ax.coords['dec'].set_ticklabel_visible(True)
	ax.set_ylim(150, 700)
	ax.set_xlim(100, 650)
	if i in [2,3]:
		ax.text(150,170,label[i],color='w')
	else:
		ax.text(150,170,label[i],color='w')
plt.subplots_adjust(wspace=0.2,hspace=0.1)
plt.savefig(datadir+'moment_maps/'+'all'+'_moment0.pdf',bbox_inches='tight')
plt.show()'''

