'''Big plot of VLA continuum on a spw by spw basis '''


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
direc='/mnt/bigdata/tetarenk/VLA_grs1915_data/lustre/aoc/ftp/e2earchive/'
'''data_vla=[]
for i in np.arange(32):
	if os.path.isfile(direc+'grs1915_rob_ms_nopb_tay1_redo_spw'+str(i)+'.fits'):
		if i not in [24,26]:
			hdulist1 = fits.open(direc+'grs1915_rob_ms_nopb_tay1_redo_spw'+str(i)+'.fits')[0]
			data=hdulist1.data[0,0,:,:]
			data_vla.append(data)
			print i
		
data_cube_vla=np.array(data_vla)
print data_cube_vla.shape'''

headervla = fits.getdata(direc+'grs1915_rob_ms_nopb_tay1_redo_spw'+str(2)+'.fits', header=True)[1]
'''fits.writeto(filename=datadir+'other_data/grs1915_vla_cube.fits',output_verify='ignore',\
	clobber=True,data=data_cube_vla,header=headervla)'''



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

label=['4.30GHz','4.40GHz','4.51GHz','4.62GHz','4.73GHz','4.84GHz','4.94GHz','5.05GHz','5.16GHz','5.27GHz','5.38GHz',\
'5.48GHz','5.59GHz','5.70GHz','5.81GHz','5.92GHz','6.02GHz','6.13GHz','6.24GHz','6.35GHz','6.46GHz','6.56GHz',\
'6.67GHz','6.78GHz','6.89GHz']
fig=plt.figure(figsize=(20,20))
#plt.rcdefaults()
plt.rc('xtick.major', size=4)
plt.rc('xtick', color='w', labelsize='small')
for i in range(0,data_cube_vla.shape[0]):
	j=i+1
	data1=data_cube_vla[i]
	wmap1=wcs.WCS(headervla)
	coord0=SkyCoord('19h15m41.8s','+10d40m38s',frame='icrs')
	coord1=SkyCoord('19h15m34.9s','+10d41m52s',frame='icrs')
	x1=float(wmap1.wcs_world2pix(coord0.ra.value,coord0.dec.value,0,0,1)[0])
	y1=float(wmap1.wcs_world2pix(coord0.ra.value,coord0.dec.value,0,0,1)[1])
	x2=float(wmap1.wcs_world2pix(coord1.ra.value,coord1.dec.value,0,0,1)[0])
	y2=float(wmap1.wcs_world2pix(coord1.ra.value,coord1.dec.value,0,0,1)[1])
	x=np.arange(0,len(data1[0,:]))
	y=np.arange(0,len(data1[:,0]))
	X, Y = np.meshgrid(x, y)
	Z=data1#[0,0,490:550,470:550]
	levels=np.array([1,2,3,4,5,6,7])*0.000345
	ax=plt.subplot(5,5,j, projection=wmap1.celestial)
	im=plt.imshow(np.nan_to_num(data1)*1e3,origin="lower",cmap=cm.get_cmap('jet', 500),norm=colors.PowerNorm(gamma=1.),vmin=0.0,vmax=5.6)
	#cbar=plt.colorbar(im, orientation='vertical',fraction=0.035,pad=0)
	#cbar.ax.tick_params(labelsize=6)
	plt.contour(X,Y,Z,levels,colors='w',transform=ax.get_transform(wmap))
#cbar.set_label('Jy')
	plt.tick_params(axis='both', which='major', labelsize=1,width=3,length=7,color='k')
	plt.tick_params(axis='both', which='minor', labelsize=1,width=1,length=7,color='k')
	if j in [2,3,4,5,7,8,9,10,12,13,14,15,17,18,19,20]:
		ax.coords['ra'].set_axislabel('')
		ax.coords['dec'].set_axislabel('',minpad=-0.1)
		ax.coords['ra'].set_major_formatter('hh:mm:ss')
		ax.coords['ra'].set_ticklabel_visible(False)
		ax.coords['dec'].set_ticklabel_visible(False)
	elif j in [22,23,24,25]:
		ax.coords['ra'].set_axislabel('Right Ascension')
		ax.coords['dec'].set_axislabel('',minpad=-0.1)
		ax.coords['ra'].set_major_formatter('hh:mm:ss')
		ax.coords['ra'].set_ticklabel_visible(True)
		ax.coords['dec'].set_ticklabel_visible(False)
	elif j ==21:
		ax.coords['ra'].set_axislabel('Right Ascension')
		ax.coords['dec'].set_axislabel('Declination',minpad=-0.1)
		ax.coords['ra'].set_major_formatter('hh:mm:ss')
		ax.coords['ra'].set_ticklabel_visible(True)
		ax.coords['dec'].set_ticklabel_visible(True)
		
	elif j in [1,6,11,16]:
		ax.coords['ra'].set_axislabel('')
		ax.coords['dec'].set_axislabel('Declination',minpad=-0.1)
		ax.coords['ra'].set_major_formatter('hh:mm:ss')
		ax.coords['ra'].set_ticklabel_visible(False)
		ax.coords['dec'].set_ticklabel_visible(True)
	ax.set_ylim(y1,y2)
	ax.set_xlim(x1,x2)
	ax.text(550,460,label[i],color='w')
	#ax.set_aspect('equal')
	#ax.set_xticklabels([])
	#ax.set_yticklabels([])
	#plt.setp(ax.get_yticklabels(),visible=False)
	#plt.setp(ax.get_xticklabels(), rotation=45, horizontalalignment='right')
cbar_ax = fig.add_axes([0.91, 0.13, 0.02, 0.74])
cbar=plt.colorbar(im, orientation='vertical',cax=cbar_ax)
cbar.ax.tick_params(labelsize=14)
cbar.set_label('mJy/beam')
plt.subplots_adjust(right=0.9,top=0.9,bottom=0.1,wspace=0.1,hspace=0.1)
#plt.tight_layout()
plt.savefig(datadir+'other_data/'+'vla_cube'+'_contour.pdf',bbox_inches='tight')
plt.show()

