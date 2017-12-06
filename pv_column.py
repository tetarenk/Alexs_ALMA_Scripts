'''PV slice of C18O column density'''


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
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes,inset_axes,InsetPosition
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

def mjrFormatter(x, cr2,cd2):
    return (x-(cr2/1000.))/(cd2/1000.)

datadir = '/mnt/bigdata/tetarenk/ALMA_GRS1915_105/'
line='18CO'
wd=5.
ra_line=['19h15m40.1s','19h15m37.6s']
dec_line=['10d40m47s','10d41m53s']
ff=datadir+'columns_plots/18CO_cube_column.fits'
############################
'''pv = PVSlicer(Xkms)
pv.show()'''

X=SpectralCube.read(ff)
Xkms=X.with_spectral_unit(u.km/u.s,velocity_convention='radio')
g=SkyCoord(ra_line,dec_line,frame='icrs')
path=Path(g,width=wd*u.arcsec)

'''gtop=SkyCoord(['19h15m40.39s','19h15m37.89s'],['10d40m49.43s','10d41m55.43s'],frame='icrs')
gbot=SkyCoord(['19h15m39.8s','19h15m37.3s'],['10d40m44.56s','10d41m50.56s'],frame='icrs')
pathtop=Path(gtop,width=wd*u.arcsec)
pathbot=Path(gbot,width=wd*u.arcsec)
slicestop = extract_pv_slice(Xkms, pathtop)
slicesbot = extract_pv_slice(Xkms, pathbot)
os.system('rm -rf '+datadir+'PVdiag/pvtop_'+line+'.fits')
slices.writeto(datadir+'PVdiag/pvtop_'+line+'.fits')
os.system('rm -rf '+datadir+'PVdiag/pvbot_'+line+'.fits')
slices.writeto(datadir+'PVdiag/pvbot_'+line+'.fits')'''


slices = extract_pv_slice(Xkms, path)
os.system('rm -rf '+datadir+'PVdiag/pvcol_'+line+'.fits')
slices.writeto(datadir+'PVdiag/pvcol_'+line+'.fits')


fitsfile = fits.open(datadir+'PVdiag/pvcol_'+line+'.fits')
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
ax.imshow(datapv,aspect='auto',interpolation='nearest',cmap=plt.get_cmap('binary'),vmin=0,vmax=1.5e16)
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
levels=np.array([3,4,5,6,7,8,9])*1.5e15
#levels=np.array([1,2,3,3.5,4,5])*0.08
plt.contour(X,Y,Z,levels,colors='k')
ax.axvspan(NTsep1/(3600.*cd1),NTsep2/(3600.*cd1),color='m',alpha=0.3)
#plt.savefig(datadir+'for_paper/pv_'+line+'.pdf',bbox_inches='tight')
#plt.savefig(datadir+'for_paper/pvcol18_'+line+'.pdf',bbox_inches='tight')
plt.show()

fits_file1=datadir+'columns_plots/18CO_column.fits'
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
fig=plt.figure()
plt.rcdefaults()
plt.rc('xtick.major', size=4)
#plt.rc('xtick', color='w', labelsize='large')
ax1 = fig.add_subplot(111, projection=wmap1.celestial)
im=plt.imshow(np.nan_to_num(data1)/1e16,origin="lower",cmap=cm.get_cmap('hot_r', 500),\
norm=colors.PowerNorm(gamma=1),vmin=0.0,vmax=0.65)
cbar=plt.colorbar(im, orientation='vertical',fraction=0.04,pad=0)
cbar.set_label('$N_{{\\rm C}^{18}{\\rm O}}\\,\\times 10^{16} \\,{\\rm cm}^{-2}$')
ax1.tick_params(axis='both', which='major', labelsize=15,width=3,length=7,color='k')
ax1.tick_params(axis='both', which='minor', labelsize=15,width=1,length=7,color='k')
ax1.coords['ra'].set_axislabel('Right Ascension')
ax1.coords['dec'].set_axislabel('Declination',minpad=-0.1)
ax1.coords['ra'].set_major_formatter('hh:mm:ss.s')
ax1.set_ylim(150, 700)
ax1.set_xlim(100, 650)
plt.contour(X,Y,Z,levels,colors='w',transform=ax1.get_transform(wmap))
ax2=inset_axes(ax1,width=3,height=1,loc=1)
plt.rcdefaults()
ax2.plot([300,400],[300,600],ls='',lw=0)
ax2.set_xticks([])
ax2.set_yticks([])
plt.setp(ax2.get_yticklabels(), visible=False)
plt.setp(ax2.get_xticklabels(), visible=False)
ax3=plt.gcf().add_axes([0.5,0.5,1,1])
ax3.imshow(datapv,aspect='auto',interpolation='nearest',cmap=plt.get_cmap('binary'),vmin=0,vmax=1.5e16)
ax3.get_yaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x,p :"%.0f"%((cr2/1000.)+x*(cd2/1000.))))
ax3.set_ylim(mjrFormatter(63, cr2,cd2),mjrFormatter(71, cr2,cd2))
ax3.set_yticks([mjrFormatter(63, cr2,cd2),mjrFormatter(66, cr2,cd2),mjrFormatter(69, cr2,cd2)])
ax3.get_xaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x,p :"%.0f"%(x*cd1*3600)))
x=np.arange(0,len(datapv[0,:]))
y=np.arange(0,len(datapv[:,0]))
X, Y = np.meshgrid(x, y)
Z=datapv#[0,0,490:550,470:550]
levels=np.array([3,4,5,6,7,8,9])*1.5e15
plt.contour(X,Y,Z,levels,colors='k')
ax3.axvspan(NTsep1/(3600.*cd1),NTsep2/(3600.*cd1),color='m',alpha=0.3)
ip = InsetPosition(ax1,[0.2,0.1,0.7,0.25])
ax2.set_axes_locator(ip)
ax3.set_axes_locator(ip)
plt.setp(ax3.get_xticklabels(), rotation=45, horizontalalignment='right')
ax3.tick_params(axis='both', which='major', length=5,width=2)
for tick in ax3.get_xticklines():
    tick.set_color('k')
for tick in ax3.get_yticklines():
    tick.set_color('k')
ax3.patch.set_facecolor('grey')
ax3.patch.set_alpha(0.3)
ax3.set_ylabel('Velocity (km/s)')
ax3.set_xlabel('Offset (arcsec)',labelpad=2)
#mark_inset(ax1, ax2, loc1=2, loc2=1, fc="none", ec="0.3",lw=2)
#ax3.set_yticks([0.0,0.01,0.02,0.03,0.04,0.05,0.06])
plt.savefig(datadir+'for_paper/18CO_column_contour_pv.pdf',bbox_inches='tight')
plt.show()

