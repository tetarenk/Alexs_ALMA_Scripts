'''create pv diagrams'''

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
from astropy import wcs,coordinates
from astropy import units as u
from astropy.coordinates import SkyCoord
import matplotlib.patches as patches
import matplotlib
import matplotlib.colors as colors
import os
from pvextractor.gui import PVSlicer
from pvextractor import Path
from pvextractor import extract_pv_slice
from matplotlib import gridspec
import matplotlib as mpl

def mjrFormatter(x, cr2,cd2):
    return (x-(cr2/1000.))/(cd2/1000.)
def change_coords(ra,dec,wm):
	x=float(wm.wcs_world2pix(ra,dec,0,0,1)[0])
	y=float(wm.wcs_world2pix(ra,dec,0,0,1)[1])
	return(x,y)

datadir = '/mnt/bigdata/tetarenk/ALMA_GRS1915_105/'
line='18CO'
wd=5.
ra_line=['19h15m40.1s','19h15m37.6s']
dec_line=['10d40m47s','10d41m53s']
#ra_line=['19h15m39.8s','19h15m37.3s']
#dec_line=['10d40m44.56s','10d41m50.56s']
ff=datadir+'alex_imaging_'+line+'_fix/GRS1915_modelimg_'+line+'.image.pbcor.fits'
############################
'''pv = PVSlicer(Xkms)
pv.show()'''

def create_PV(ff,g,path):
	X=SpectralCube.read(ff)
	XK=X.to(u.K, equivalencies=X.beam.jtok_equiv(X.header['RESTFRQ']*u.Hz))
	Xkms=XK.with_spectral_unit(u.km/u.s,velocity_convention='radio')
	path=Path(g,width=wd*u.arcsec)
	slices = extract_pv_slice(Xkms, path)
	os.system('rm -rf '+datadir+'PVdiag/pv_'+line+'.fits')
	slicestop.writeto(datadir+'PVdiag/pv_'+line+'.fits')


g=SkyCoord(ra_line,dec_line,frame='icrs')
path=Path(g,width=wd*u.arcsec)
create_PV(ff,g,path)

'''gtop=SkyCoord(['19h15m40.39s','19h15m37.89s'],['10d40m49.43s','10d41m55.43s'],frame='icrs')
gbot=SkyCoord(['19h15m39.8s','19h15m37.3s'],['10d40m44.56s','10d41m50.56s'],frame='icrs')
pathtop=Path(gtop,width=wd*u.arcsec)
pathbot=Path(gbot,width=wd*u.arcsec)
slicestop = extract_pv_slice(Xkms, pathtop)
slicesbot = extract_pv_slice(Xkms, pathbot)
os.system('rm -rf '+datadir+'PVdiag/pvtop_'+line+'.fits')
slicestop.writeto(datadir+'PVdiag/pvtop_'+line+'.fits')
os.system('rm -rf '+datadir+'PVdiag/pvbot_'+line+'.fits')
slicesbot.writeto(datadir+'PVdiag/pvbot_'+line+'.fits')'''
gtop=SkyCoord(['19h15m40.39s','19h15m37.89s'],['10d40m49.43s','10d41m55.43s'],frame='icrs')
gbot=SkyCoord(['19h15m39.8s','19h15m37.3s'],['10d40m44.56s','10d41m50.56s'],frame='icrs')

#test pv slice
fitsfile = fits.open(datadir+'PVdiag/pvtop_'+line+'.fits')
hdu = fitsfile[0]
datapv=hdu.data
cd1=hdu.header['CDELT1']
cd2=hdu.header['CDELT2']
cr2=hdu.header['CRVAL2']
hdu.header['CUNIT1']='arcsec'
hdu.header['CDELT1']=cd1*3600.
w = wcs.WCS(hdu.header)
NTsep1=56.00539935050752
NTsep2=66.20744343524173
fig=plt.figure(figsize=(10,7))
#plt.rc_context({'axes.edgecolor':'white', 'xtick.color':'white', 'ytick.color':'white', 'figure.facecolor':'white'})
ax=plt.gca()
#ax2a=ax.twinx()
#ax=fig.add_axes([0.1,0.1,0.8,0.8],projection=w)
ax.imshow(datapv,aspect='auto',interpolation='nearest',cmap=plt.get_cmap('binary'),vmin=0)#,vmax=1.5e16)
ax.get_yaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x,p :"%.0f"%((cr2/1000.)+x*(cd2/1000.))))
ax.set_ylim(mjrFormatter(63, cr2,cd2),mjrFormatter(71, cr2,cd2))
ax.set_yticks([mjrFormatter(63, cr2,cd2),mjrFormatter(66, cr2,cd2),mjrFormatter(69, cr2,cd2)])
ax.get_xaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x,p :"%.0f"%(x*cd1*3600)))
#ax.set_xlim(1,375)
#plt.setp(ax.get_yticklabels(), rotation=55, horizontalalignment='right')
ax.set_ylabel('Velocity\n(km/s)')
ax.set_xlabel('Offset (arcsec)')
x=np.arange(0,len(datapv[0,:]))
y=np.arange(0,len(datapv[:,0]))
X, Y = np.meshgrid(x, y)
Z=datapv#[0,0,490:550,470:550]
#levels=np.array([3,4,5,6,7,8,9])*1.5e15
levels=np.array([0.2,0.5,0.7,0.8])*5
plt.contour(X,Y,Z,levels,colors='k')
ax.axvspan(NTsep1/(3600.*cd1),NTsep2/(3600.*cd1),color='m',alpha=0.3)
plt.savefig(datadir+'for_paper/pv_'+line+'.pdf',bbox_inches='tight')
#plt.savefig(datadir+'for_paper/pvcol18_'+line+'.pdf',bbox_inches='tight')
plt.show()


fits_file1='/mnt/bigdata/tetarenk/ALMA_GRS1915_105/'+'moment_maps/'+line+'_moment'+str(0)+'.fits'
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
coord0=g[0]
coord1=g[1]
coord0t=gtop[0]
coord1t=gtop[1]
#val=np.sqrt(((wd/2.)**2)/2.)
x1=float(wmap1.wcs_world2pix(coord0.ra.value,coord0.dec.value,0,0,1)[0])
y1=float(wmap1.wcs_world2pix(coord0.ra.value,coord0.dec.value,0,0,1)[1])
x2=float(wmap1.wcs_world2pix(coord1.ra.value,coord1.dec.value,0,0,1)[0])
y2=float(wmap1.wcs_world2pix(coord1.ra.value,coord1.dec.value,0,0,1)[1])
x1t=float(wmap1.wcs_world2pix(coord0t.ra.value,coord0t.dec.value,0,0,1)[0])
y1t=float(wmap1.wcs_world2pix(coord0t.ra.value,coord0t.dec.value,0,0,1)[1])
x2t=float(wmap1.wcs_world2pix(coord1t.ra.value,coord1t.dec.value,0,0,1)[0])
y2t=float(wmap1.wcs_world2pix(coord1t.ra.value,coord1t.dec.value,0,0,1)[1])

thet2val=180-90-np.arctan((y2-y1)/(x2-x1))*(180./np.pi)
xval=(wd/2.)*np.cos(thet2val*np.pi/180.)
yval=(wd/2.)*np.sin(thet2val*np.pi/180.)
cra_low=coord0.ra-xval*u.arcsec
cdec_low=coord0.dec-yval*u.arcsec
cra_high=coord0.ra+xval*u.arcsec
cdec_high=coord0.dec+yval*u.arcsec
x3=float(wmap1.wcs_world2pix(cra_low.value,cdec_low.value,0,0,1)[0])
y3=float(wmap1.wcs_world2pix(cra_low.value,cdec_low.value,0,0,1)[1])
x4=float(wmap1.wcs_world2pix(cra_high.value,cdec_high.value,0,0,1)[0])
y4=float(wmap1.wcs_world2pix(cra_high.value,cdec_high.value,0,0,1)[1])
cra_low1=coord0.ra-(wd/2.)*u.arcsec
cdec_low1=coord0.dec-(wd/2.)*u.arcsec
cra_high1=coord0.ra+(wd/2.)*u.arcsec
cdec_high1=coord0.dec+(wd/2.)*u.arcsec
x31=float(wmap1.wcs_world2pix(cra_low1.value,cdec_low1.value,0,0,1)[0])
y31=float(wmap1.wcs_world2pix(cra_low1.value,cdec_low1.value,0,0,1)[1])
x41=float(wmap1.wcs_world2pix(cra_high1.value,cdec_high1.value,0,0,1)[0])
y41=float(wmap1.wcs_world2pix(cra_high1.value,cdec_high1.value,0,0,1)[1])

thet2valt=180-90-np.arctan((y2t-y1t)/(x2t-x1t))*(180./np.pi)
xvalt=(wd/2.)*np.cos(thet2valt*np.pi/180.)
yvalt=(wd/2.)*np.sin(thet2valt*np.pi/180.)
cra_lowt=coord0t.ra-xvalt*u.arcsec
cdec_lowt=coord0t.dec-yvalt*u.arcsec
cra_hight=coord0t.ra+xvalt*u.arcsec
cdec_hight=coord0t.dec+yvalt*u.arcsec
x3t=float(wmap1.wcs_world2pix(cra_lowt.value,cdec_lowt.value,0,0,1)[0])
y3t=float(wmap1.wcs_world2pix(cra_lowt.value,cdec_lowt.value,0,0,1)[1])
x4t=float(wmap1.wcs_world2pix(cra_hight.value,cdec_hight.value,0,0,1)[0])
y4t=float(wmap1.wcs_world2pix(cra_hight.value,cdec_hight.value,0,0,1)[1])

#pv on inset in moment map
'''fig,axes=plt.subplots()
gs=gridspec.GridSpec(1,3,width_ratios=[1,0.1,1])
gs2=gridspec.GridSpec(3,3,width_ratios=[1,0.1,1])
#ax1=plt.subplot2grid((6,6),(0,3),colspan=4,rowspan=6,projection=wmap1.celestial)
ax1=plt.subplot(gs[2],projection=wmap1.celestial)
im=plt.imshow(np.nan_to_num(data1),origin="lower",cmap=cm.get_cmap('jet', 500),vmin=0.0)
#cbar=plt.colorbar(im, orientation='vertical',fraction=0.04,pad=0)
#cbar.set_label('J km/s')
ax1.tick_params(axis='both', which='major', labelsize=15,width=3,length=7,color='k')
ax1.tick_params(axis='both', which='minor', labelsize=15,width=1,length=7,color='k')
ax1.coords['ra'].set_axislabel('Right Ascension')
ax1.coords['dec'].set_axislabel('Declination',minpad=-0.7)
ax1.coords['ra'].set_major_formatter('hh:mm:ss')
ax1.set_ylim(150, 700)
ax1.set_xlim(100, 650)
plt.plot([x1,x2],[y1,y2],marker='',lw=1,color='k',ls='--')
rect = patches.Rectangle((x3,y3),np.sqrt((y2-y1)**2+(x2-x1)**2),y41-y31,angle=np.arctan((y2-y1)/(x2-x1))*(180./np.pi),\
linewidth=1,edgecolor='k',facecolor='none')
# Add the patch to the Axes
ax1.add_patch(rect)
plt.contour(X,Y,Z,levels,colors='k',transform=ax1.get_transform(wmap))
#ax2=plt.subplot2grid((6,6),(0,0),colspan=2,rowspan=6)#2,projection=w)
ax2=plt.subplot(gs2[1,0])#2,projection=w)
#ax2a=ax2.twinx()
ax2.imshow(datapv,interpolation='nearest',aspect='auto')
ax2.get_yaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x,p :"%.0f"%((cr2/1000.)+x*(cd2/1000.))))
ax2.set_ylim(mjrFormatter(60, cr2,cd2),mjrFormatter(78, cr2,cd2))
ax2.set_yticks([mjrFormatter(60, cr2,cd2),mjrFormatter(65, cr2,cd2),mjrFormatter(70, cr2,cd2),mjrFormatter(75, cr2,cd2)])
ax2.set_xlim(1,375)
ax2.get_xaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x,p :"%.0f"%(x*cd1*3600)))
ytickloc = ax2.get_yticks()
#ax2a.set_yticks(ytickloc)
#ax2a.set_ylim(ax2.get_ylim())
plt.setp(ax2.get_yticklabels(), rotation=55, horizontalalignment='right')
ax2.set_ylabel('Velocity (km/s)')
ax2.set_xlabel('Offset (arcsec)')
#plt.setp(ax2.get_yticklabels(),visible=False)
#ax2.invert_xaxis()
#ax2a.invert_yaxis()
#ax2.invert_yaxis()
for tick in ax2.get_xticklines():
	tick.set_color('white')
for tick in ax2.get_yticklines():
	tick.set_color('white')
for tick in ax1.get_xticklines():
	tick.set_color('white')
for tick in ax1.get_yticklines():
	tick.set_color('white')
plt.savefig(datadir+'PVdiag/pvmap_'+line+'.pdf',bbox_inches='tight')
plt.show()'''

#three plot for paper
cra_high1a=coord0.ra+(wd/2.)*u.arcsec
cdec_high1a=coord0.dec+(wd/2.)*u.arcsec
x41a=float(wmap1.wcs_world2pix(cra_high1a.value,cdec_high1a.value,0,0,1)[0])
y41a=float(wmap1.wcs_world2pix(cra_high1a.value,cdec_high1a.value,0,0,1)[1])
cra_low1a=coord0.ra-(wd/2.)*u.arcsec
cdec_low1a=coord0.dec-(wd/2.)*u.arcsec
x41b=float(wmap1.wcs_world2pix(cra_low1a.value,cdec_low1a.value,0,0,1)[0])
y41b=float(wmap1.wcs_world2pix(cra_low1a.value,cdec_low1a.value,0,0,1)[1])
fifty=y41a-y41b
fitsfile12 = fits.open(datadir+'PVdiag/pv_12CO.fits')
fitsfile13 = fits.open(datadir+'PVdiag/pv_13CO.fits')
fitsfile18 = fits.open(datadir+'PVdiag/pv_18CO.fits')#datadir+'PVdiag/pvtop_'+line+'.fits'
fitsfile12top = fits.open(datadir+'PVdiag/pvtop_12CO.fits')
fitsfile13top = fits.open(datadir+'PVdiag/pvtop_13CO.fits')
fitsfile18top = fits.open(datadir+'PVdiag/pvtop_18CO.fits')
hdu12 = fitsfile12[0]
hdu13 = fitsfile13[0]
hdu18 = fitsfile18[0]
datapv12=hdu12.data
datapv13=hdu13.data
datapv18=hdu18.data
hdu12t = fitsfile12top[0]
hdu13t = fitsfile13top[0]
hdu18t = fitsfile18top[0]
datapv12t=hdu12t.data
datapv13t=hdu13t.data
datapv18t=hdu18t.data
cd1_12=hdu12.header['CDELT1']
cd2_12=hdu12.header['CDELT2']
cr2_12=hdu12.header['CRVAL2']
hdu12.header['CUNIT1']='arcsec'
hdu12.header['CDELT1']=cd1_12*3600.
w12 = wcs.WCS(hdu12.header)
cd1_18=hdu18.header['CDELT1']
cd2_18=hdu18.header['CDELT2']
cr2_18=hdu18.header['CRVAL2']
hdu18.header['CUNIT1']='arcsec'
hdu18.header['CDELT1']=cd1_18*3600.
w18 = wcs.WCS(hdu18.header)
cd1_13=hdu13.header['CDELT1']
cd2_13=hdu13.header['CDELT2']
cr2_13=hdu13.header['CRVAL2']
hdu13.header['CUNIT1']='arcsec'
hdu13.header['CDELT1']=cd1_13*3600.
w13 = wcs.WCS(hdu12.header)


raNT='19h15m38.178s'
decNT='10d41m36.829s'
gt0=SkyCoord(raNT,decNT,frame='icrs')
XNT,YNT=change_coords(gt0.ra.value,gt0.dec.value,wmap1)
raNT='19h15m38.084s'
decNT='10d41m34.472s'
gt1=SkyCoord(raNT,decNT,frame='icrs')
XNT1,YNT1=change_coords(gt1.ra.value,gt1.dec.value,wmap1)
raNT='19h15m37.784s'
decNT='10d41m43.728s'
gt2=SkyCoord(raNT,decNT,frame='icrs')
XNT2,YNT2=change_coords(gt2.ra.value,gt2.dec.value,wmap1)

c1=g[0]
c2=gt0
NTsep=c1.separation(c2).arcsec
c1=g[0]
c2=gt1
NTsep1=c1.separation(c2).arcsec
c1=g[0]
c2=gt2
NTsep2=c1.separation(c2).arcsec
c1=gt1
c2=gt2
NTsep3=c1.separation(c2).arcsec

cra_high1a2=coord0.ra+(NTsep3/2.)*u.arcsec
cdec_high1a2=coord0.dec+(NTsep3/2.)*u.arcsec
x41a2=float(wmap1.wcs_world2pix(cra_high1a2.value,cdec_high1a2.value,0,0,1)[0])
y41a2=float(wmap1.wcs_world2pix(cra_high1a2.value,cdec_high1a2.value,0,0,1)[1])
cra_low1a2=coord0.ra-(NTsep3/2.)*u.arcsec
cdec_low1a2=coord0.dec-(NTsep3/2.)*u.arcsec
x41b2=float(wmap1.wcs_world2pix(cra_low1a2.value,cdec_low1a2.value,0,0,1)[0])
y41b2=float(wmap1.wcs_world2pix(cra_low1a2.value,cdec_low1a2.value,0,0,1)[1])
fifty2=y41a2-y41b2

fig,axes=plt.subplots(figsize=(10,6.))
gs=gridspec.GridSpec(3,10,width_ratios=[1.7,1.7,0.7,0.5,1.5,1.5,1.5,0.2,1.7,1.7])
gs2=gridspec.GridSpec(3,10,width_ratios=[1.7,1.7,0.7,0.5,1.5,1.5,1.5,0.2,1.7,1.7])
gs3=gridspec.GridSpec(3,10,width_ratios=[1.7,1.7,0.7,0.5,1.5,1.5,1.5,0.2,1.7,1.7])
ax1 = plt.subplot(gs3[0, 8:10])
ax1b=ax1.twinx()
ax1.imshow(datapv12,aspect='auto',interpolation='nearest',cmap=plt.get_cmap('binary'))
ax1.get_yaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x,p :"%.0f"%((cr2_12/1000.)+x*(cd2_12/1000.))))
ax1.set_ylim(mjrFormatter(63, cr2_12,cd2_12),mjrFormatter(71, cr2_12,cd2_12))
ax1.set_yticks([mjrFormatter(63, cr2_12,cd2_12),mjrFormatter(66, cr2_12,cd2_12),mjrFormatter(69, cr2_12,cd2_12)])
ax1b.get_yaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x,p :"%.0f"%((cr2_12/1000.)+x*(cd2_12/1000.))))
ax1b.set_ylim(mjrFormatter(63, cr2_12,cd2_12),mjrFormatter(71, cr2_12,cd2_12))
ax1b.set_yticks([mjrFormatter(63, cr2_12,cd2_12),mjrFormatter(66, cr2_12,cd2_12),mjrFormatter(69, cr2_12,cd2_12)])
#ax1.set_xlim(1,375)
ax1.get_xaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x,p :"%.0f"%(x*cd1_12*3600)))
#ax1.axvline(x=NTsep/(3600.*cd1_12),color='m',lw=2,ls='--')
#ax1.axvline(x=NTsep1/(3600.*cd1_12),color='m',lw=2,ls='--')
#ax1.axvline(x=NTsep2/(3600.*cd1_12),color='m',lw=2,ls='--')
ax1.axvspan(NTsep1/(3600.*cd1_12),NTsep2/(3600.*cd1_12),color='m',alpha=0.3)
plt.setp(ax1.get_xticklabels(), rotation=55, horizontalalignment='right')
plt.setp(ax1.get_xticklabels(),visible=False)
plt.setp(ax1.get_yticklabels(),visible=False)
x=np.arange(0,len(datapv12[0,:]))
y=np.arange(0,len(datapv12[:,0]))
X, Y = np.meshgrid(x, y)
Z=datapv12#[0,0,490:550,470:550]
levels=np.array([1,2,3,3.5,4])*5.
plt.contour(X,Y,Z,levels,colors='k')
#ax1.set_ylabel('Velocity (km/s)')
#ax1.set_xlabel('Offset (arcsec)')
ax2 = plt.subplot(gs3[1, 8:10])
ax2b=ax2.twinx()
ax2.imshow(datapv13,aspect='auto',interpolation='nearest',cmap=plt.get_cmap('binary'))
ax2.get_yaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x,p :"%.0f"%((cr2_13/1000.)+x*(cd2_13/1000.))))
ax2.set_ylim(mjrFormatter(63, cr2_13,cd2_13),mjrFormatter(71, cr2_13,cd2_13))
ax2.set_yticks([mjrFormatter(63, cr2_13,cd2_13),mjrFormatter(66, cr2_13,cd2_13),mjrFormatter(69, cr2_13,cd2_13)])
ax2b.get_yaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x,p :"%.0f"%((cr2_13/1000.)+x*(cd2_13/1000.))))
ax2b.set_ylim(mjrFormatter(63, cr2_13,cd2_13),mjrFormatter(71, cr2_13,cd2_13))
ax2b.set_yticks([mjrFormatter(63, cr2_13,cd2_13),mjrFormatter(66, cr2_13,cd2_13),mjrFormatter(69, cr2_13,cd2_13)])
#ax2.set_xlim(1,375)
ax2.get_xaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x,p :"%.0f"%(x*cd1_13*3600)))
#ax2.axvline(x=NTsep/(3600.*cd1_13),color='m',lw=2,ls='--')
ax2.axvspan(NTsep1/(3600.*cd1_13),NTsep2/(3600.*cd1_13),color='m',alpha=0.3)
plt.setp(ax2.get_xticklabels(), rotation=55, horizontalalignment='right')
plt.setp(ax2.get_xticklabels(),visible=False)
plt.setp(ax2.get_yticklabels(),visible=False)
x=np.arange(0,len(datapv13[0,:]))
y=np.arange(0,len(datapv13[:,0]))
X, Y = np.meshgrid(x, y)
Z=datapv13#[0,0,490:550,470:550]
levels=np.array([1,1.3,1.6,2,2.2,2.5])*5
plt.contour(X,Y,Z,levels,colors='k')
#ax1.set_ylabel('Velocity (km/s)')
#ax2.set_xlabel('Offset (arcsec)')
ax3 = plt.subplot(gs3[2, 8:10])
ax3b=ax3.twinx()
ax3.imshow(datapv18,aspect='auto',interpolation='nearest',cmap=plt.get_cmap('binary'))
ax3.get_yaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x,p :"%.0f"%((cr2_18/1000.)+x*(cd2_18/1000.))))
ax3.set_ylim(mjrFormatter(63, cr2_18,cd2_18),mjrFormatter(71, cr2_18,cd2_18))
ax3.set_yticks([mjrFormatter(63, cr2_18,cd2_18),mjrFormatter(66, cr2_18,cd2_18),mjrFormatter(69, cr2_18,cd2_18)])
ax3b.get_yaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x,p :"%.0f"%((cr2_18/1000.)+x*(cd2_18/1000.))))
ax3b.set_ylim(mjrFormatter(63, cr2_18,cd2_18),mjrFormatter(71, cr2_18,cd2_18))
ax3b.set_yticks([mjrFormatter(63, cr2_18,cd2_18),mjrFormatter(66, cr2_18,cd2_18),mjrFormatter(69, cr2_18,cd2_18)])
#ax3.set_xlim(1,375)
ax3.get_xaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x,p :"%.0f"%(x*cd1_18*3600)))
#ax3.axvline(x=NTsep/(3600.*cd1_18),color='m',lw=2,ls='--')
ax3.axvspan(NTsep1/(3600.*cd1_18),NTsep2/(3600.*cd1_18),color='m',alpha=0.3)
plt.setp(ax3.get_xticklabels(), rotation=55, horizontalalignment='right')
plt.setp(ax3.get_yticklabels(),visible=False)
x=np.arange(0,len(datapv18[0,:]))
y=np.arange(0,len(datapv18[:,0]))
X, Y = np.meshgrid(x, y)
Z=datapv18#[0,0,490:550,470:550]
levels=np.array([0.2,0.5,0.7,0.8])*5
plt.contour(X,Y,Z,levels,colors='k')
ax2b.set_ylabel('Velocity (km/s)',labelpad=10)
#ax3.xaxis.set_label_coords(0.8,-0.05)
ax3.set_xlabel('Offset (arcsec)')
ax4 = plt.subplot(gs2[0:3, 4:7],projection=wmap1.celestial)
ax4.coords['ra'].set_axislabel('Right Ascension')
ax4.coords['dec'].set_axislabel('Declination',minpad=-0.7)
#ax4.coords['ra'].set_ticklabel_visible(False)
#ax4.coords['dec'].set_ticklabel_visible(False)
ax4.coords['ra'].set_major_formatter('hh:mm:ss')
ax4.set_ylim(250, 700)
ax4.set_xlim(200, 500)
#plt.plot([x1,x2],[y1,y2],marker='',lw=1,color='b',ls='--')
rect = patches.Rectangle((x3,y3),np.sqrt((y2-y1)**2+(x2-x1)**2),fifty,angle=np.arctan((y2-y1)/(x2-x1))*(180./np.pi),\
linewidth=1,edgecolor='b',facecolor='none')
# Add the patch to the Axes
ax4.add_patch(rect)
rectt = patches.Rectangle((x3t,y3t),np.sqrt((y2t-y1t)**2+(x2t-x1t)**2),fifty,angle=np.arctan((y2t-y1t)/(x2t-x1t))*(180./np.pi),\
linewidth=1,edgecolor='b',facecolor='none',linestyle='dashed')
# Add the patch to the Axes
ax4.add_patch(rectt)
rect2 = patches.Rectangle((XNT1,YNT1),np.sqrt((YNT2-YNT1)**2+(XNT1-XNT1)**2),fifty,angle=np.arctan((YNT2-YNT1)/(XNT2-XNT1))*(180./np.pi),\
linewidth=1,edgecolor='m',facecolor='m',alpha=0.3)
# Add the patch to the Axes
ax4.add_patch(rect2)
	

thet2=180-90-np.arctan((y2-y1)/(x2-x1))*(180./np.pi)
x_top=x3-fifty*np.cos(thet2*np.pi/180.)
y_top=y3+fifty*np.sin(thet2*np.pi/180.)
x_bot=x3
y_bot=y3
xs_bot=[]
ys_bot=[]
xs_top=[]
ys_top=[]
thet2t=180-90-np.arctan((y2t-y1t)/(x2t-x1t))*(180./np.pi)
x_topt=x3t-fifty*np.cos(thet2t*np.pi/180.)
y_topt=y3t+fifty*np.sin(thet2t*np.pi/180.)
x_bott=x3t
y_bott=y3t
xs_bott=[]
ys_bott=[]
xs_topt=[]
ys_topt=[]
for i in range(0,7):
	xs_top.append(x_top+(2.*fifty+2*fifty*i)*np.cos(np.arctan((y2-y1)/(x2-x1))))
	ys_top.append(y_top+(2.*fifty+2*fifty*i)*np.sin(np.arctan((y2-y1)/(x2-x1))))
	xs_bot.append(x_bot+(2.*fifty+2*fifty*i)*np.cos(np.arctan((y2-y1)/(x2-x1))))
	ys_bot.append(y_bot+(2.*fifty+2*fifty*i)*np.sin(np.arctan((y2-y1)/(x2-x1))))
	xs_topt.append(x_topt+(2.*fifty+2*fifty*i)*np.cos(np.arctan((y2t-y1t)/(x2t-x1t))))
	ys_topt.append(y_topt+(2.*fifty+2*fifty*i)*np.sin(np.arctan((y2t-y1t)/(x2t-x1t))))
	xs_bott.append(x_bott+(2.*fifty+2*fifty*i)*np.cos(np.arctan((y2t-y1t)/(x2t-x1t))))
	ys_bott.append(y_bott+(2.*fifty+2*fifty*i)*np.sin(np.arctan((y2t-y1t)/(x2t-x1t))))
for i in range(0,len(xs_top)):
	plt.plot([xs_top[i],xs_bot[i]],[ys_top[i],ys_bot[i]],marker='',lw=1,color='r',ls='-')
for i in range(0,len(xs_topt)):
	plt.plot([xs_topt[i],xs_bott[i]],[ys_topt[i],ys_bott[i]],marker='',lw=1,color='r',ls='--')
plt.plot([x_top,x_bot],[y_top,y_bot],marker='',lw=1,color='r',ls='-')
plt.plot([x_topt,x_bott],[y_topt,y_bott],marker='',lw=1,color='r',ls='--')
x=np.arange(0,len(data[0,:]))
y=np.arange(0,len(data[:,0]))
X, Y = np.meshgrid(x, y)
Z=data#[0,0,490:550,470:550]
levels=np.array([4,6,8,10,15,20,40,60])*0.00005
plt.contour(X,Y,Z,levels,colors='k',transform=ax4.get_transform(wmap))
plt.subplots_adjust(wspace=None,hspace=0.05)
ax4.spines['left'].set_visible(False)
ax4.spines['right'].set_visible(False)
ax4.spines['top'].set_visible(False)
ax4.spines['bottom'].set_visible(False)
ax4.set_axis_off()
ax1.tick_params(axis='y', which='major',labelsize=11)
ax2.tick_params(axis='y', which='major',labelsize=11)
ax3.tick_params(axis='y', which='major',labelsize=11)
ax4.errorbar(250,300,xerr=fifty,marker=None,ls='',color='k')
ax4.text(240,275,'10"')
'''for tick in ax2.get_xticklines():
	tick.set_color('white')
for tick in ax2.get_yticklines():
	tick.set_color('white')
for tick in ax1.get_xticklines():
	tick.set_color('white')
for tick in ax1.get_yticklines():
	tick.set_color('white')
for tick in ax3.get_xticklines():
	tick.set_color('white')
for tick in ax3.get_yticklines():
	tick.set_color('white')
#for tick in ax4.get_xticklines():
	#tick.set_color('white')
#for tick in ax4.get_yticklines():
	#tick.set_color('white')'''
ax4.set_aspect('equal', 'datalim')

ax5 = plt.subplot(gs[0, :2])
#ax5b = ax5.twinx()
ax5.imshow(datapv12t,aspect='auto',interpolation='nearest',cmap=plt.get_cmap('binary'))
ax5.get_yaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x,p :"%.0f"%((cr2_12/1000.)+x*(cd2_12/1000.))))
ax5.set_ylim(mjrFormatter(63, cr2_12,cd2_12),mjrFormatter(71, cr2_12,cd2_12))
ax5.set_yticks([mjrFormatter(63, cr2_12,cd2_12),mjrFormatter(66, cr2_12,cd2_12),mjrFormatter(69, cr2_12,cd2_12)])
ax5.get_xaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x,p :"%.0f"%(x*cd1_12*3600)))
#ax5b.get_yaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x,p :"%.0f"%((cr2_12/1000.)+x*(cd2_12/1000.))))
#ax5b.set_ylim(mjrFormatter(63, cr2_12,cd2_12),mjrFormatter(71, cr2_12,cd2_12))
#ax5b.set_yticks([mjrFormatter(63, cr2_12,cd2_12),mjrFormatter(66, cr2_12,cd2_12),mjrFormatter(69, cr2_12,cd2_12)])
#ax1.axvline(x=NTsep/(3600.*cd1_12),color='m',lw=2,ls='--')
#ax1.axvline(x=NTsep1/(3600.*cd1_12),color='m',lw=2,ls='--')
#ax1.axvline(x=NTsep2/(3600.*cd1_12),color='m',lw=2,ls='--')
#ax5.axvspan(NTsep1/(3600.*cd1_12),NTsep2/(3600.*cd1_12),color='m',alpha=0.3)
plt.setp(ax5.get_xticklabels(), rotation=55, horizontalalignment='right')
plt.setp(ax5.get_xticklabels(),visible=False)
#plt.setp(ax5.get_yticklabels(),visible=False)
x=np.arange(0,len(datapv12[0,:]))
y=np.arange(0,len(datapv12[:,0]))
X, Y = np.meshgrid(x, y)
Z=datapv12t#[0,0,490:550,470:550]
levels=np.array([1,2,3,3.5,4])*5.
plt.contour(X,Y,Z,levels,colors='k')
#ax1.set_ylabel('Velocity (km/s)')
#ax1.set_xlabel('Offset (arcsec)')
ax6 = plt.subplot(gs[1, :2])
#ax6b=ax6.twinx()
ax6.imshow(datapv13t,aspect='auto',interpolation='nearest',cmap=plt.get_cmap('binary'))
ax6.get_yaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x,p :"%.0f"%((cr2_13/1000.)+x*(cd2_13/1000.))))
ax6.set_ylim(mjrFormatter(63, cr2_13,cd2_13),mjrFormatter(71, cr2_13,cd2_13))
ax6.set_yticks([mjrFormatter(63, cr2_13,cd2_13),mjrFormatter(66, cr2_13,cd2_13),mjrFormatter(69, cr2_13,cd2_13)])
#ax6b.get_yaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x,p :"%.0f"%((cr2_13/1000.)+x*(cd2_13/1000.))))
#ax6b.set_ylim(mjrFormatter(63, cr2_13,cd2_13),mjrFormatter(71, cr2_13,cd2_13))
#ax6b.set_yticks([mjrFormatter(63, cr2_13,cd2_13),mjrFormatter(66, cr2_13,cd2_13),mjrFormatter(69, cr2_13,cd2_13)])
#ax2.set_xlim(1,375)
ax6.get_xaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x,p :"%.0f"%(x*cd1_13*3600)))
#ax2.axvline(x=NTsep/(3600.*cd1_13),color='m',lw=2,ls='--')
#ax6.axvspan(NTsep1/(3600.*cd1_13),NTsep2/(3600.*cd1_13),color='m',alpha=0.3)
plt.setp(ax6.get_xticklabels(), rotation=55, horizontalalignment='right')
plt.setp(ax6.get_xticklabels(),visible=False)
#plt.setp(ax6.get_yticklabels(),visible=False)
x=np.arange(0,len(datapv13[0,:]))
y=np.arange(0,len(datapv13[:,0]))
X, Y = np.meshgrid(x, y)
Z=datapv13t#[0,0,490:550,470:550]
levels=np.array([1,1.3,1.6,2,2.2,2.5])*5
plt.contour(X,Y,Z,levels,colors='k')
#ax1.set_ylabel('Velocity (km/s)')
#ax2.set_xlabel('Offset (arcsec)')
ax7 = plt.subplot(gs[2, :2])
#ax7b=ax7.twinx()
ax7.imshow(datapv18t,aspect='auto',interpolation='nearest',cmap=plt.get_cmap('binary'))
ax7.get_yaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x,p :"%.0f"%((cr2_18/1000.)+x*(cd2_18/1000.))))
ax7.set_ylim(mjrFormatter(63, cr2_18,cd2_18),mjrFormatter(71, cr2_18,cd2_18))
ax7.set_yticks([mjrFormatter(63, cr2_18,cd2_18),mjrFormatter(66, cr2_18,cd2_18),mjrFormatter(69, cr2_18,cd2_18)])
#ax7b.get_yaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x,p :"%.0f"%((cr2_18/1000.)+x*(cd2_18/1000.))))
#ax7b.set_ylim(mjrFormatter(63, cr2_18,cd2_18),mjrFormatter(71, cr2_18,cd2_18))
#ax7b.set_yticks([mjrFormatter(63, cr2_18,cd2_18),mjrFormatter(66, cr2_18,cd2_18),mjrFormatter(69, cr2_18,cd2_18)])
#ax3.set_xlim(1,375)
ax7.get_xaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x,p :"%.0f"%(x*cd1_18*3600)))
#ax3.axvline(x=NTsep/(3600.*cd1_18),color='m',lw=2,ls='--')
#ax7.axvspan(NTsep1/(3600.*cd1_18),NTsep2/(3600.*cd1_18),color='m',alpha=0.3)
plt.setp(ax7.get_xticklabels(), rotation=55, horizontalalignment='right')
x=np.arange(0,len(datapv18[0,:]))
y=np.arange(0,len(datapv18[:,0]))
X, Y = np.meshgrid(x, y)
Z=datapv18t#[0,0,490:550,470:550]
levels=np.array([0.2,0.5,0.7,0.8])*5
plt.contour(X,Y,Z,levels,colors='k')
ax6.set_ylabel('Velocity (km/s)')
#ax3.xaxis.set_label_coords(0.8,-0.05)
ax7.set_xlabel('Offset (arcsec)')
#plt.setp(ax7.get_yticklabels(),visible=False)
plt.savefig(datadir+'for_paper/pv_allCO_2off.pdf',bbox_inches='tight')
plt.show()

