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
import pyfits
import scipy.ndimage as nd
import pyspeckit
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes,inset_axes,InsetPosition
from mpl_toolkits.axes_grid1.inset_locator import mark_inset


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
	fits.writeto(filename=ddir+'moment_maps/'+line+'_moment'+str(mind)+flag+'.fits',output_verify='ignore',\
	clobber=True,data=mom_array,header=header)
def make_spec(fitsfile,ddir,limits,flag):
	x0,x1,y0,y1=limits[0],limits[1],limits[2],limits[3]
	a=SpectralCube.read(fitsfile)
	a.allow_huge_operations=True
	a1=a.to(u.K, equivalencies=a.beam.jtok_equiv(a.header['RESTFRQ']*u.Hz))
	a2=a1.with_spectral_unit(u.km/u.s, velocity_convention='radio')
	if flag == 'reg':
		suba=a2[:, y0:y1, x0:x1]
		veloc=suba.spectral_axis
		#print veloc
		amps = np.mean(suba.filled_data[:],axis=(1,2))
	elif flag == 'pix':
		suba=a2[:, y0, x0]
		veloc=suba.spectral_axis
		#print veloc
		amps=a2[:, y0, x0].value
	return(veloc.to(u.km/u.s),amps)


i='CS'#322_221 needs 0-80 only,sio 375-464 in full cube
ii='CH3OH_5'
flag='K'
moment=0
make_moment_maps(datadir+'alex_imaging_'+i+'_fix/GRS1915_modelimg_'+i+'.image.pbcor.fits',moment,datadir,i,flag)
#lim=[309,348,513,550]#methanol
#lim2=[344,376,471,489]#methanol2
#lim=[315,344,523,546]#sio
#lim=[328,372,434,486]#h2co1
#lim2=[310,355,506,545]#h2co2
lim=[332,381,446,496]
lim2=[306,357,518,556]
vel,amp=make_spec(datadir+'alex_imaging_'+i+'_fix/GRS1915_modelimg_'+i+'.image.pbcor.fits',datadir,lim,'reg')
vel2,amp2=make_spec(datadir+'alex_imaging_'+i+'_fix/GRS1915_modelimg_'+i+'.image.pbcor.fits',datadir,lim2,'reg')
plt.plot(vel,amp,color='b')
plt.plot(vel2,amp2,color='r')
plt.show()
raw_input('stop')


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
rad=150#100,125,125,125,125,100
data_mask=nd.morphology.binary_erosion(np.nan_to_num(data1), np.ones((rad, rad)))
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
#evels=np.array([1,2,3,4,5,6,7])*0.000345
levels=np.array([4,6,8,10,15,20,40,60])*0.00005
fig=plt.figure()
plt.rcdefaults()
plt.rc('xtick.major', size=4)
#plt.rc('xtick', color='w', labelsize='large')
ax1 = fig.add_subplot(111, projection=wmap1.celestial)
im=plt.imshow(np.nan_to_num(data1)*data_mask,origin="lower",cmap=cm.get_cmap('hot_r', 500),vmin=0,vmax=7,norm=colors.PowerNorm(gamma=0.5))#,vmax=2)#,vmax=2.0)#,vmax=2)#,vmax=68,norm=colors.PowerNorm(gamma=1.0))#,vmax=10.0)#,vmax=7)#,norm=colors.PowerNorm(gamma=1.5))
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
plt.contour(X,Y,Z,levels,colors='k',transform=ax1.get_transform(wmap),zorder=4)
'''ax1.add_patch(r1)
ax1.add_patch(r2)
ax1.add_patch(r3)
ax1.add_patch(r4)
ax1.text(365,560, 'A',fontsize=15)
ax1.text(428,575, 'B',fontsize=15)
ax1.text(368,485, 'C',fontsize=15)
ax1.text(585,492, 'D',fontsize=15)'''
ax2=inset_axes(ax1,width=3,height=1,loc=1)
plt.rcdefaults()
ax2.plot([lim[0],lim[1]],[lim[2],lim[3]],ls='',lw=0)
ax2.set_xticks([])
ax2.set_yticks([])
plt.setp(ax2.get_yticklabels(), visible=False)
plt.setp(ax2.get_xticklabels(), visible=False)
ax3=plt.gcf().add_axes([0.5,0.5,1,1])
ax3.plot(vel, amp,color='b', ls='-',lw=3)
ax3.set_ylabel('$T_B$ (K)',labelpad=3)
ax3.set_xlabel('$v$ (km/s)',labelpad=2)
ax3.set_ylim(0,0.12)
#ax3.set_yticks([0.0,0.01,0.02,0.03,0.04])
ax3.set_xlim(60,75)
ip = InsetPosition(ax1,[0.15,0.15,0.3,0.2])
ax2.set_axes_locator(ip)
ax3.set_axes_locator(ip)
plt.setp(ax3.get_xticklabels(), rotation=45, horizontalalignment='right')
#ax3.xaxis.label.set_color('white')
#ax3.tick_params(axis='x', colors='white')
#ax3.yaxis.label.set_color('white')
#ax3.tick_params(axis='y', colors='white')
ax3.tick_params(axis='both', which='major', length=5,width=2)
for tick in ax3.get_xticklines():
    tick.set_color('k')
for tick in ax3.get_yticklines():
    tick.set_color('k')
ax3.patch.set_facecolor('grey')
ax3.patch.set_alpha(0.3)
mark_inset(ax1, ax2, loc1=2, loc2=1, fc="none", ec="0.3",lw=2)
ax3.set_yticks([0.0,0.02,0.04,0.06,0.08,0.1])

ax2b=inset_axes(ax1,width=3,height=1,loc=3)
plt.rcdefaults()
ax2b.plot([lim2[0],lim2[1]],[lim2[2],lim2[3]],ls='',lw=0)
ax2b.set_xticks([])
ax2b.set_yticks([])
plt.setp(ax2b.get_yticklabels(), visible=False)
plt.setp(ax2b.get_xticklabels(), visible=False)
ax3b=plt.gcf().add_axes([0,0,0.4,0.4])
ax3c=ax2b.twinx()
ax3d=ax2b.twiny()
ax3b.plot(vel2, amp2,color='b', ls='-',lw=3)
ax3c.set_ylabel('$T_B$ (K)',labelpad=3)
ax3d.set_xlabel('$v$ (km/s)',labelpad=2)
ax3b.set_ylim(0.008,0.08)
ax3c.set_ylim(0.008,0.08)
#ax3.set_yticks([0.0,0.01,0.02,0.03,0.04])
ax3b.set_xlim(53,80)
ax3d.set_xlim(53,80)
ax3d.set_xticks([55,60,65,70,75,80])
ax3c.set_yticks([0.0,0.02,0.04,0.06,0.08])
#ax3.set_yticks([0.0,0.02,0.04,0.06,0.08,0.1])
ipb = InsetPosition(ax1,[0.57,0.7,0.3,0.2])
ax2b.set_axes_locator(ipb)
ax3b.set_axes_locator(ipb)
plt.setp(ax3b.get_xticklabels(), rotation=45, horizontalalignment='right')
plt.setp(ax3d.get_xticklabels(), rotation=-55, horizontalalignment='right')
#ax3b.xaxis.label.set_color('white')
#ax3b.tick_params(axis='x', colors='white')
#ax3b.yaxis.label.set_color('white')
#ax3b.tick_params(axis='y', colors='white')
ax3b.tick_params(axis='both', which='major', length=5,width=2)
plt.setp(ax3b.get_yticklabels(), visible=False)
#ax3c.yaxis.label.set_color('white')
#ax3c.tick_params(axis='y', colors='white')
ax3c.tick_params(axis='both', which='major', length=5,width=2)
plt.setp(ax3b.get_xticklabels(), visible=False)
#ax3d.xaxis.label.set_color('white')
#ax3d.tick_params(axis='x', colors='white')
ax3d.tick_params(axis='both', which='major', length=5,width=2)
for tick in ax3b.get_xticklines():
    tick.set_color('k')
for tick in ax3b.get_yticklines():
    tick.set_color('k')
ax3b.patch.set_facecolor('grey')
ax3b.patch.set_alpha(0.3)
mark_inset(ax1, ax2b, loc1=2, loc2=3, fc="none", ec="0.3",lw=2)
plt.savefig(datadir+'for_paper/'+i+'_moment'+str(moment)+'_contourK_inset.pdf',bbox_inches='tight')
plt.show()
#raw_input('stop')
#import scipy.ndimage as nd
#rad=150#100,125,125,125,125,100
#data_mask=nd.morphology.binary_erosion(np.nan_to_num(data1), np.ones((rad, rad)))
#lim=[328,372,434,486]#h2co1
#lim2=[310,355,506,545]

