'''plot alma mosaic fields'''

from astropy import units as u
import numpy as np
from astropy.io import fits,ascii
from astropy.table import Table,Column,MaskedColumn
import math as ma
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from astropy import wcs
from astropy.coordinates import SkyCoord
import matplotlib.patches as patches
import matplotlib
import matplotlib.colors as colors
import os

datadir = '/mnt/bigdata/tetarenk/ALMA_GRS1915_105/'

data12=ascii.read('/home/ubuntu/Dropbox/12m_fields.txt', delimiter=' ',data_start=0,names=['ID','Code','Name','RA',\
'Dec', 'Epoch','SrcID','Rows'],guess='False')
data7=ascii.read('/home/ubuntu/Dropbox/7m_fields.txt', delimiter=' ',data_start=0,names=['ID','Code','Name','RA',\
'Dec', 'Epoch','SrcID','Rows'],guess='False')

vals_ra12=data12['RA']
vals_dec12=data12['Dec']
vals_id12=data12['ID']
vals_ra7=data7['RA']
vals_dec7=data7['Dec']
vals_id7=data7['ID']

#12m
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


coord0=SkyCoord('19h15m43.4s','+10d40m12s',frame='icrs')
coord1=SkyCoord('19h15m33.8s','+10d42m00s',frame='icrs')
x1=float(wmap1.wcs_world2pix(coord0.ra.value,coord0.dec.value,0,0,1)[0])
y1=float(wmap1.wcs_world2pix(coord0.ra.value,coord0.dec.value,0,0,1)[1])
x2=float(wmap1.wcs_world2pix(coord1.ra.value,coord1.dec.value,0,0,1)[0])
y2=float(wmap1.wcs_world2pix(coord1.ra.value,coord1.dec.value,0,0,1)[1])
coord3=SkyCoord('19h15m37.2s','+10d41m00s',frame='icrs')
x3=float(wmap1.wcs_world2pix(coord3.ra.value,coord3.dec.value,0,0,1)[0])
y3=float(wmap1.wcs_world2pix(coord3.ra.value,coord3.dec.value,0,0,1)[1])
coord4=SkyCoord('19h15m38.64s','+10d41m07.9s',frame='icrs')
x4=float(wmap1.wcs_world2pix(coord4.ra.value,coord4.dec.value,0,0,1)[0])
y4=float(wmap1.wcs_world2pix(coord4.ra.value,coord4.dec.value,0,0,1)[1])

x=np.arange(0,len(data[0,:]))
y=np.arange(0,len(data[:,0]))
X, Y = np.meshgrid(x, y)
Z=data#[0,0,490:550,470:550]
levels=np.array([4,6,8,10,15,20,40,60])*0.00005
fig=plt.figure()
plt.rcdefaults()
plt.rc('xtick.major', size=4)
ax1 = fig.add_subplot(111, projection=wmap1.celestial)
#im=plt.imshow(np.nan_to_num(data1[0,0,:,:]),origin="lower",cmap=cm.get_cmap('Purples', 500),norm=colors.PowerNorm(gamma=1),vmin=0.0,vmax=0.1)#65,55,65,0.9,0.9,0.9,0.9/x,x,x,0.1,0.03,0.025,0.025,0.025,0.1
#cbar=plt.colorbar(im, orientation='vertical',fraction=0.04,pad=0)
#cbar.set_label('Primary Beam Response')
ax1.tick_params(axis='both', which='major', labelsize=15,width=3,length=7,color='k')
ax1.tick_params(axis='both', which='minor', labelsize=15,width=1,length=7,color='k')
ax1.coords['ra'].set_axislabel('Right Ascension')
ax1.coords['dec'].set_axislabel('Declination',minpad=-0.1)
ax1.coords['ra'].set_major_formatter('hh:mm:ss.s')
ax1.set_ylim(y1, y2)
ax1.set_xlim(x1, x2)
#plt.contour(X,Y,Z,levels,colors='w',transform=ax1.get_transform(wmap),linewidths=2)
for i in range(0,len(vals_ra12)):
	ra=vals_ra12[i]
	dec=vals_dec12[i]
	ids=float(vals_id12[i])-3.
	print ra,dec
	coord=SkyCoord(ra[0:2]+'h'+ra[3:5]+'m'+ra[6:15]+'s',dec[0:3]+'d'+dec[4:6]+'m'+dec[7:15]+'s',frame='icrs')
	x=float(wmap1.wcs_world2pix(coord.ra.value,coord.dec.value,0,0,1)[0])
	y=float(wmap1.wcs_world2pix(coord.ra.value,coord.dec.value,0,0,1)[1])
	c1 = patches.Circle((x,y), (26.2/2)/0.2, linewidth=2, fill=False,color='b',zorder=1,ls='dotted')
	ax1.add_patch(c1)
	ax1.plot(x,y,marker='x',ms=12,color='r',mew=1.5,zorder=3)
	ax1.text(x-40,y-50,str(int(ids)),color='k',zorder=5,fontsize=15)
ax1.plot(x4,y4,marker='o',ms=12,color='m')
plt.savefig(datadir+'for_paper/ALMA_12mpoint_map.pdf',bbox_inches='tight')
ax1.set_aspect('equal', 'datalim')
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


coord0=SkyCoord('19h15m43.4s','+10d40m12s',frame='icrs')
coord1=SkyCoord('19h15m33.8s','+10d42m00s',frame='icrs')
x1=float(wmap1.wcs_world2pix(coord0.ra.value,coord0.dec.value,0,0,1)[0])
y1=float(wmap1.wcs_world2pix(coord0.ra.value,coord0.dec.value,0,0,1)[1])
x2=float(wmap1.wcs_world2pix(coord1.ra.value,coord1.dec.value,0,0,1)[0])
y2=float(wmap1.wcs_world2pix(coord1.ra.value,coord1.dec.value,0,0,1)[1])
coord3=SkyCoord('19h15m37.2s','+10d41m00s',frame='icrs')
x3=float(wmap1.wcs_world2pix(coord3.ra.value,coord3.dec.value,0,0,1)[0])
y3=float(wmap1.wcs_world2pix(coord3.ra.value,coord3.dec.value,0,0,1)[1])
coord4=SkyCoord('19h15m38.64s','+10d41m07.9s',frame='icrs')
x4=float(wmap1.wcs_world2pix(coord4.ra.value,coord4.dec.value,0,0,1)[0])
y4=float(wmap1.wcs_world2pix(coord4.ra.value,coord4.dec.value,0,0,1)[1])

x=np.arange(0,len(data[0,:]))
y=np.arange(0,len(data[:,0]))
X, Y = np.meshgrid(x, y)
Z=data#[0,0,490:550,470:550]
levels=np.array([4,6,8,10,15,20,40,60])*0.00005
fig=plt.figure()
plt.rcdefaults()
plt.rc('xtick.major', size=4)
ax1 = fig.add_subplot(111, projection=wmap1.celestial)
#im=plt.imshow(np.nan_to_num(data1[0,0,:,:]),origin="lower",cmap=cm.get_cmap('Purples', 500),norm=colors.PowerNorm(gamma=1),vmin=0.0,vmax=0.1)#65,55,65,0.9,0.9,0.9,0.9/x,x,x,0.1,0.03,0.025,0.025,0.025,0.1
#cbar=plt.colorbar(im, orientation='vertical',fraction=0.04,pad=0)
#cbar.set_label('Primary Beam Response')
ax1.tick_params(axis='both', which='major', labelsize=15,width=3,length=7,color='k')
ax1.tick_params(axis='both', which='minor', labelsize=15,width=1,length=7,color='k')
ax1.coords['ra'].set_axislabel('Right Ascension')
ax1.coords['dec'].set_axislabel('Declination',minpad=-0.1)
ax1.coords['ra'].set_major_formatter('hh:mm:ss.s')
ax1.set_ylim(y1, y2)
ax1.set_xlim(x1, x2)
#plt.contour(X,Y,Z,levels,colors='w',transform=ax1.get_transform(wmap),linewidths=2)
for i in range(0,len(vals_ra7)):
	ra=vals_ra7[i]
	dec=vals_dec7[i]
	ids=float(vals_id7[i])-4.
	print ra,dec
	coord=SkyCoord(ra[0:2]+'h'+ra[3:5]+'m'+ra[6:15]+'s',dec[0:3]+'d'+dec[4:6]+'m'+dec[7:15]+'s',frame='icrs')
	x=float(wmap1.wcs_world2pix(coord.ra.value,coord.dec.value,0,0,1)[0])
	y=float(wmap1.wcs_world2pix(coord.ra.value,coord.dec.value,0,0,1)[1])
	c1 = patches.Circle((x,y), (44.4/2)/1., linewidth=2, fill=False,color='b',zorder=1,ls='dotted')
	ax1.add_patch(c1)
	ax1.plot(x,y,marker='x',ms=12,color='r',mew=1.5,zorder=3)
	ax1.text(x-10,y-15,str(int(ids)),color='k',zorder=5,fontsize=15)
ax1.plot(x4,y4,marker='o',ms=12,color='m')
plt.savefig(datadir+'for_paper/ALMA_7mpoint_map.pdf',bbox_inches='tight')
ax1.set_aspect('equal', 'datalim')
plt.show()

