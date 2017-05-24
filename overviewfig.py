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
from astropy.coordinates import SkyCoord,FK5,FK4
import matplotlib.patches as patches
import matplotlib
import matplotlib.colors as colors
import glob
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes,inset_axes,InsetPosition
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from matplotlib import gridspec
datadir = '/mnt/bigdata/tetarenk/ALMA_GRS1915_105/'


fits_file1='/mnt/bigdata/tetarenk/ALMA_GRS1915_105/other_data/hpacs_25HPPUNIMAPB_blue_1917_p1149_00_v1.0_1471608553396.fits'
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
hdulist1 = fits.open(fits_file1)[1]
data1=hdulist1.data
wmap1=wcs.WCS(hdulist1.header)
coord0=SkyCoord('19h15m00s','10d25m00s',frame='fk4',equinox='B1950')
coord1=SkyCoord('19h11m30s','11d25m00s',frame='fk4',equinox='B1950')
C1=coord0.transform_to(FK5(equinox='J2000'))
C2=coord1.transform_to(FK5(equinox='J2000'))
x1=float(wmap1.wcs_world2pix(C1.ra.value,C1.dec.value,1)[0])
y1=float(wmap1.wcs_world2pix(C1.ra.value,C1.dec.value,1)[1])
x2=float(wmap1.wcs_world2pix(C2.ra.value,C2.dec.value,1)[0])
y2=float(wmap1.wcs_world2pix(C2.ra.value,C2.dec.value,1)[1])

def conv_1950_2000(ra,dec):
	coord0=SkyCoord(ra,dec,frame='fk4',equinox='B1950')
	C1=coord0.transform_to(FK5(equinox='J2000'))
	return(C1.ra.hms,C1.dec.dms)
def conv_2000_1950(ra,dec):
	coord0=SkyCoord(ra,dec,frame='icrs')
	C1=coord0.transform_to(FK4(equinox='B1950'))
	return(C1.ra.hms,C1.dec.dms)

x=np.arange(0,len(data[0,:]))
y=np.arange(0,len(data[:,0]))
X, Y = np.meshgrid(x, y)
Z=data
levels=np.array([4,6,8,10,15,20,40,60])*0.00005

coord1915=SkyCoord('19h15m11.6s','10d56m44s',frame='icrs')
xgrs=float(wmap1.wcs_world2pix(coord1915.ra.value,coord1915.dec.value,1)[0])
ygrs=float(wmap1.wcs_world2pix(coord1915.ra.value,coord1915.dec.value,1)[1])
coord1915p5mra=coord1915.ra+5*u.arcmin
coord1915p5mdec=coord1915.dec+5*u.arcmin
xgrs2=float(wmap1.wcs_world2pix(coord1915p5mra.value,coord1915p5mdec.value,1)[0])
ygrs2=float(wmap1.wcs_world2pix(coord1915p5mra.value,coord1915p5mdec.value,1)[1])
lengarr=ygrs2-ygrs
endxarr=xgrs+lengarr*np.cos(64*np.pi/180.)
endyarr=ygrs+lengarr*np.sin(64*np.pi/180.)
endxarr2=xgrs-lengarr*np.cos(64*np.pi/180.)
endyarr2=ygrs-lengarr*np.sin(64*np.pi/180.)



fig,axes=plt.subplots()
gs=gridspec.GridSpec(4,4,width_ratios=[1,0.5,0.05,1])
gs2=gridspec.GridSpec(4,4,width_ratios=[1,0.5,0.05,1])
plt.rcdefaults()
plt.rc('xtick.major', size=4)
#plt.rc('xtick', color='w', labelsize='large')
ax1 = plt.subplot(gs[0:4, 0:2], projection=wmap1.celestial)
im=plt.imshow(np.nan_to_num(data1[:,:]),origin="lower",cmap=cm.get_cmap('hot_r', 500),norm=colors.PowerNorm(gamma=0.4),vmin=0.0,vmax=5.)#65,55,65,0.9,0.9,0.9,0.9/x,x,x,0.1,0.03,0.025,0.025,0.025,0.1
#cbar=plt.colorbar(im, orientation='vertical',fraction=0.04,pad=0)
#cbar.set_label('mJy/beam')
ax1.tick_params(axis='both', which='major', labelsize=14,width=3,length=7,color='k')
ax1.tick_params(axis='both', which='minor', labelsize=14,width=1,length=7,color='k')
ax1.coords['ra'].set_axislabel('Right Ascension')
ax1.coords['dec'].set_axislabel('Declination',minpad=-0.1)
ax1.coords['ra'].set_major_formatter('hh:mm:ss')
ax1.set_ylim(y1,y2)
ax1.set_xlim(x1,x2)
ax1.arrow(xgrs, ygrs, endxarr-xgrs, endyarr-ygrs, head_width=20, head_length=20, fc='k', ec='k',lw=1)
ax1.arrow(xgrs, ygrs, endxarr2-xgrs, endyarr2-ygrs, head_width=20, head_length=20, fc='k', ec='k',lw=1)
ax1.plot(xgrs,ygrs,marker='+',ls='',ms=14,mew=3,color='m')
ax2 = plt.subplot(gs2[1:3, -1],projection=wmap.celestial)
plt.rc('xtick', color='w', labelsize='large')
plt.rc('xtick.major', size=4)
coord0a=SkyCoord('19h15m40.8s','+10d40m58s',frame='icrs')
coord1a=SkyCoord('19h15m37.1s','+10d41m46s',frame='icrs')
x1a=float(wmap.wcs_world2pix(coord0a.ra.value,coord0a.dec.value,1)[0])
y1a=float(wmap.wcs_world2pix(coord0a.ra.value,coord0a.dec.value,1)[1])
x2a=float(wmap.wcs_world2pix(coord1a.ra.value,coord1a.dec.value,1)[0])
y2a=float(wmap.wcs_world2pix(coord1a.ra.value,coord1a.dec.value,1)[1])
coordt1=SkyCoord('19h15m39.8s','+10d40m00s',frame='icrs')
coordt2=SkyCoord('19h15m39.8s','+10d40m15s',frame='icrs')
x3=float(wmap.wcs_world2pix(coordt1.ra.value,coordt1.dec.value,1)[0])
y3=float(wmap.wcs_world2pix(coordt1.ra.value,coordt1.dec.value,1)[1])
x4=float(wmap.wcs_world2pix(coordt2.ra.value,coordt2.dec.value,1)[0])
y4=float(wmap.wcs_world2pix(coordt2.ra.value,coordt2.dec.value,1)[1])
fifty5=y4-y3
#plt.imshow(np.nan_to_num(data[:,:])*1000.,origin="lower",cmap=cm.get_cmap('hot_r', 500),norm=colors.PowerNorm(gamma=0.7),vmin=0.0,vmax=3.)
plt.contour(X,Y,Z,levels,colors='k',transform=ax2.get_transform(wmap))
ax2.tick_params(axis='both', which='major', labelsize=15,width=3,length=7,color='w')
ax2.tick_params(axis='both', which='minor', labelsize=15,width=1,length=7,color='w')
#ax2.coords['ra'].set_axislabel('Right Ascension')
#ax2.coords['dec'].set_axislabel('Declination',minpad=-0.1)
ax2.coords['ra'].set_major_formatter('hh:mm:ss')
ax2.coords['dec'].set_ticklabel_visible(False)
ax2.coords['ra'].set_ticklabel_visible(False)
ax2.set_ylim(y1a,y2a)
ax2.set_xlim(x1a,x2a)
ax2.errorbar(495,485,xerr=fifty5/2,marker=None,ls='',color='k')
ax2.text(490,480,'15"')
ax2.set_aspect('equal', 'datalim')
plt.savefig(datadir+'for_paper/overv.pdf',bbox_inches='tight')
plt.show()
