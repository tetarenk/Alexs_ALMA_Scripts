'''Creates Moment 0 maps with inset spectra'''


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
from matplotlib.patches import Rectangle
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
		amps = np.mean(suba.filled_data[:],axis=(1,2))
	elif flag == 'pix':
		suba=a2[:, y0, x0]
		veloc=suba.spectral_axis
		amps=a2[:, y0, x0].value
	return(veloc.to(u.km/u.s),amps)
def fit_gauss_to_spec(fitsfile,ddir,line,limits,aguess):
	x0,x1,y0,y1=limits[0],limits[1],limits[2],limits[3]
	a=SpectralCube.read(fitsfile)
	a1=a.to(u.K, equivalencies=a.beam.jtok_equiv(a.header['RESTFRQ']*u.Hz))
	a2=a1.with_spectral_unit(u.km/u.s, velocity_convention='radio')
	suba=a2[:, y0:y1, x0:x1]
	data=np.mean(suba.filled_data[:],axis=(1,2))
	asp=pyspeckit.Spectrum(data=np.array(data), xarr=suba.spectral_axis, xarrkwargs={'unit':'km/s'}, unit="$T_B$")
	#width guess
	#aguess=[aamp, acenter, awidth]
	asp.specfit(fittype='gaussian', guesses=aguess)
	asp.plotter(errstyle='fill')
	asp.specfit.plot_fit()
	asp.plotter.savefig(ddir+'spectra_plots/'+line+'_fitted.eps')
	#asp.specfit.modelpars/modelerrs give parameters and errors of spectrum fit
	return(asp.specfit.modelpars,asp.specfit.modelerrs)


'''Line name; 
Reminders: (1) H2CO322_221 has methanol line in spw as well (chans 0-80),
		   (2) continuum has CS line in spw as well
		   (3) SiO has noise away from line (chans 375-464) '''
i='CS'
#ii='CH3OH_5'
flag='K'
moment=0
make_moment_maps(datadir+'alex_imaging_'+i+'_fix2/GRS1915_modelimg_'+i+'.image.pbcor.fits',moment,\
	datadir,i,flag)
lim=[332,381,446,496]#cs1
lim2=[306,357,518,556]#cs2
'''Regions for line emission
lim=[309,348,513,550]#methanol
lim2=[344,376,471,489]#methanol2
lim=[315,344,523,546]#sio
lim2=[375,387,562,577]#sio2
lim=[328,372,434,486]#h2co1
lim2=[310,355,506,545]#h2co2
lim=[332,381,446,496]#cs1
lim2=[306,357,518,556]#cs2
lim=[385,400,560,575]#ridge
lim2=[400,420,570,595]#nt
lim3=[340,360,470,485]#hmsr
lim=[518,597,328,383]#off'''

vel,amp=make_spec(datadir+'alex_imaging_'+i+'_fix/GRS1915_modelimg_'+i+'.image.pbcor.fits',datadir,lim,'reg')
vel2,amp2=make_spec(datadir+'alex_imaging_'+i+'_fix/GRS1915_modelimg_'+i+'.image.pbcor.fits',datadir,lim2,'reg')
#check channel range over which you are integrating
plt.figure()
plt.plot(vel,amp,color='b')
plt.plot(vel2,amp2,color='r')
plt.show()

#write spectra files for fitting later
v=np.array(vel)
d=np.array(amp)
v2=np.array(vel2)
d2=np.array(amp2)
fileo=open(datadir+'for_paper/'+i+'_spectra_file_1.txt','w')
for k in range(0,len(vel)):
	fileo.write('{0} {1}\n'.format(v[k],d[k]))
fileo.close()
fileo=open(datadir+'for_paper/'+i+'_spectra_file_2.txt','w')
for k in range(0,len(vel2)):
	fileo.write('{0} {1}\n'.format(v2[k],d2[k]))
fileo.close()


#built in gaussian fitting if you like, rather use own mcmc code to do it!
'''fitsfile=datadir+'alex_imaging_'+i+'_fix_all/GRS1915_modelimg_'+i+'.image.pbcor.fits'
pars,errs=fit_gauss_to_spec(fitsfile,datadir,i,lim,[0.06,65,20])
print pars,errs'''


#read in fits files
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

#Make mask to cut noisy edges of map if needed
rad=130#100,120,125
data_mask=nd.morphology.binary_erosion(np.nan_to_num(data1), np.ones((rad, rad)))

#radio contours
x=np.arange(0,len(data[0,:]))
y=np.arange(0,len(data[:,0]))
X, Y = np.meshgrid(x, y)
Z=data
levels=np.array([4,6,8,10,15,20,40,60])*0.00005

#start plot
fig=plt.figure(figsize=(12,10))
plt.rcdefaults()
plt.rc('xtick.major', size=4)
#plt.rc('xtick', color='w', labelsize='large')
ax1 = fig.add_subplot(111, projection=wmap1.celestial)
im=plt.imshow(np.nan_to_num(data1)*data_mask,origin="lower",cmap=cm.get_cmap('hot_r', 500),\
	vmin=0,vmax=5.5,norm=colors.PowerNorm(gamma=0.55))
if flag=='K':
	cbar.set_label('K km/s',size=22)
else:
	cbar.set_label('Jy km/s',size=22)
cbar.ax.tick_params(labelsize=22)
ax1.tick_params(axis='both', which='major', labelsize=22,width=3,length=7,color='k')
ax1.tick_params(axis='both', which='minor', labelsize=22,width=1,length=7,color='k')
ax1.coords['ra'].set_axislabel('Right Ascension',size=22)
ax1.coords['dec'].set_axislabel('Declination',minpad=-0.1,size=22)
ax1.coords['ra'].set_major_formatter('hh:mm:ss.s')
ax1.coords['ra'].set_ticklabel(size=25)
ax1.coords['dec'].set_ticklabel(size=25)
ax1.set_ylim(150, 700)
ax1.set_xlim(100, 650)
r1=Rectangle((385, 560), 15, 15,alpha=1, facecolor='none',edgecolor='white',lw=2,ls='solid',zorder=5)
r2=Rectangle((400, 570), 20, 25,alpha=1, facecolor='none',edgecolor='white',lw=2,ls='solid',zorder=5)
r3=Rectangle((340, 470), 20, 15,alpha=1, facecolor='none',edgecolor='white',lw=2,ls='solid',zorder=5)
r4=Rectangle((507, 447), 73, 45,alpha=1, facecolor='none',edgecolor='white',lw=2,ls='solid',zorder=5)
plt.contour(X,Y,Z,levels,colors='k',transform=ax1.get_transform(wmap),zorder=4,linewidths=2)

#first inset spectra
ax2=inset_axes(ax1,width=3,height=1,loc=1)
plt.rcdefaults()
ax2.plot([lim[0],lim[1]],[lim[2],lim[3]],ls='',lw=0)
ax2.set_xticks([])
ax2.set_yticks([])
plt.setp(ax2.get_yticklabels(), visible=False)
plt.setp(ax2.get_xticklabels(), visible=False)
ax3=plt.gcf().add_axes([0.5,0.5,1,1])
ax3.plot(vel, amp,color='b', ls='-',lw=3)
ax3.set_ylabel('$T_B$ (K)',labelpad=3,size=22)
ax3.set_xlabel('$v$ (km/s)',labelpad=2,size=22)
ax3.set_ylim(0,0.1)
#ax3.set_yticks([0.0,0.01,0.02,0.03,0.04])
ax3.set_xlim(60,76)
ip = InsetPosition(ax1,[0.18,0.15,0.3,0.2])
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
ax3.set_yticks([0.0,0.04,0.08])
ax3.set_xticks([60,64,68,72,76])
ax3.tick_params(axis='both', which='major', labelsize=22,color='k')

#second inset spectra
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
ax3c.set_ylabel('$T_B$ (K)',labelpad=3,size=22)
ax3b.set_xlabel('$v$ (km/s)',labelpad=2,size=22)
#ax2b.set_ylim(0.008,0.08)
ax3b.set_ylim(0,0.08)
ax3c.set_ylim(0,0.08)
#ax3.set_yticks([0.0,0.01,0.02,0.03,0.04])
ax3b.set_xlim(55,80)
ax3d.set_xlim(55,80)
#ax3d.set_xticks([])
#ax3c.set_yticks([])
ax3c.tick_params(axis='both', which='major', labelsize=22,color='k')
ax3b.tick_params(axis='both', which='major', labelsize=22,color='k')
ax3c.set_xticks([60,65,70,75,80])
ax3d.set_yticks([0.0,0.02,0.04,0.06,0.08])
ax3d.set_xticks([60,65,70,75,80])
ax3c.set_yticks([0.0,0.02,0.04,0.06,0.08])
ax3b.set_xticks([60,65,70,75,80])
ax3b.set_yticks([0.0,0.02,0.04,0.06,0.08])
#ax3.set_yticks([0.0,0.03,0.06,0.09,0.12])
ipb = InsetPosition(ax1,[0.56,0.55,0.3,0.2])

ax2b.set_axes_locator(ipb)
ax3b.set_axes_locator(ipb)
plt.setp(ax3b.get_xticklabels(), rotation=45, horizontalalignment='right')
plt.setp(ax3d.get_xticklabels(), rotation=-45, horizontalalignment='right')
#ax3b.xaxis.label.set_color('white')
#ax3b.tick_params(axis='x', colors='white')
#ax3b.yaxis.label.set_color('white')
#ax3b.tick_params(axis='y', colors='white')
ax3b.tick_params(axis='both', which='major', length=5,width=2)
plt.setp(ax3b.get_yticklabels(), visible=False)
#plt.setp(ax2b.get_yticklabels(), visible=False)
#ax3c.yaxis.label.set_color('white')
#ax3c.tick_params(axis='y', colors='white')
ax3c.tick_params(axis='both', which='major', length=5,width=2)
plt.setp(ax3d.get_xticklabels(), visible=False)
#ax3d.xaxis.label.set_color('white')
#ax3d.tick_params(axis='x', colors='white')
ax3d.tick_params(axis='both', which='major', length=5,width=2)
#for tick in ax3b.get_xticklines():
    #tick.set_color('k')
#for tick in ax3b.get_yticklines():
    #tick.set_color('k')
ax3b.patch.set_facecolor('grey')
ax3b.patch.set_alpha(0.3)
mark_inset(ax1, ax2b, loc1=2, loc2=3, fc="none", ec="0.3",lw=2)
plt.savefig(datadir+'for_paper/'+i+'_moment'+str(moment)+'_contourK_inset_fix2.pdf',bbox_inches='tight')
plt.show()


