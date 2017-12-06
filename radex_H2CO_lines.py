'''RADEX fitting of thermometry lines'''

import pyspeckit
import matplotlib
import pyfits
from pyspeckit.spectrum import models
from spectral_cube import SpectralCube
import astropy.units as u
import numpy as np
import scipy.ndimage as nd
import astropy.table
from astropy.table import Table
import matplotlib.pyplot as plt
from astropy.io import fits
import pylab as pl
import math as ma
from astropy.table import hstack
import matplotlib.cm as cm
from astropy import wcs
from astropy.coordinates import SkyCoord
import matplotlib.patches as patches
import matplotlib.colors as colors


datadir = '/mnt/bigdata/tetarenk/ALMA_GRS1915_105/'

#Step 1:
#read in model grids to create H2CO radex fitter
#grids use the radex program to calculate brightness of H2CO under different conditions
grids_dir='/mnt/bigdata/tetarenk/ALMA_GRS1915_105/radex_grids/'
texgrid1a=pyfits.getdata(grids_dir+'303-202_321-220_5kms_temperature_para_tex1.fits')
taugrid1a=pyfits.getdata(grids_dir+'303-202_321-220_5kms_temperature_para_tau1.fits')
texgrid2a=pyfits.getdata(grids_dir+'303-202_321-220_5kms_temperature_para_tex2.fits')
taugrid2a=pyfits.getdata(grids_dir+'303-202_321-220_5kms_temperature_para_tau2.fits')
hdra=pyfits.getheader(grids_dir+'303-202_321-220_5kms_temperature_para_tau2.fits')
texgrid1b=pyfits.getdata(grids_dir+'303-202_322-221_5kms_temperature_para_tex1.fits')
taugrid1b=pyfits.getdata(grids_dir+'303-202_322-221_5kms_temperature_para_tau1.fits')
texgrid2b=pyfits.getdata(grids_dir+'303-202_322-221_5kms_temperature_para_tex2.fits')
taugrid2b=pyfits.getdata(grids_dir+'303-202_322-221_5kms_temperature_para_tau2.fits')
hdrb=pyfits.getheader(grids_dir+'303-202_322-221_5kms_temperature_para_tau2.fits')


#Step 2: Fit lines
#NOTE: -models.formaldehyde.formaldehyde_radex is the model that we are going to fit 
#-models.model.SpectralModel is a wrapper to deal with parinfo, multiple peaks, and annotations
#-all parameters after the first are passed to the model function

#fit all three lines
#tex/tau grid sets frequency range (in GHz) over which frequency range is valid 
formaldehyde_radex_fitter_both=models.model.SpectralModel(
	models.formaldehyde_mm.formaldehyde_mm_radex, 5,
	parnames=['temperature', 'column', 'density', 'center', 'width'],
	parvalues=[50,12,4.5,0,1],
	parlimited=[(True, True), (True, True), (True, True), (False, False), (True, False)],
	parlimits=[(5,205), (10,17), (2,7), (0,0), (0,0)],
	parsteps=[0.01, 0.01, 0.1, 0, 0], fitunits='Hz',
	texgrid=((218.15, 218.25, texgrid1b), (218.4, 218.55, texgrid2b),(218.7, 218.8, texgrid2a)),
	taugrid=((218.15, 218.25, taugrid1b), (218.4, 218.55, taugrid2b), (218.7, 218.8, taugrid2a)),
	hdr=hdrb, shortvarnames=("T", "N", "n", "v", "\\sigma"), grid_vwidth=5.0)

formaldehyde_radex_fitter_303322 = models.model.SpectralModel(
    models.formaldehyde_mm.formaldehyde_mm_radex, 5,
    parnames=['temperature','column','density','center','width'],
    parvalues=[50,12,4.5,0,1],
    parlimited=[(True,True), (True,True), (True,True), (False,False), (True,False)],
    parlimits=[(5,205), (10,17), (2,7), (0,0), (0,0)],
    parsteps=[0.01,0.01,0.1,0,0],
    fitunits='Hz',
    texgrid=((218.15,218.25,texgrid1b),(218.4,218.55,texgrid2b)), # specify the frequency range over which the grid is valid (in GHz)
    taugrid=((218.15,218.25,taugrid1b),(218.4,218.55,taugrid2b)),
    hdr=hdrb,
    shortvarnames=("T","N","n","v","\\sigma"), # specify the parameter names (TeX is OK)
    grid_vwidth=5.0,
)

formaldehyde_radex_fitter_303321 = models.model.SpectralModel(
    models.formaldehyde_mm.formaldehyde_mm_radex, 5,
    parnames=['temperature','column','density','center','width'],
    parvalues=[50,12,4.5,0,1],
    parlimited=[(True,True), (True,True), (True,True), (False,False), (True,False)],
    parlimits=[(5,205), (10,17), (2,7), (0,0), (0,0)],
    parsteps=[0.01,0.01,0.1,0,0],
    fitunits='Hz',
    texgrid=((218.15,218.25,texgrid1a),(218.7,218.8,texgrid2a)), # specify the frequency range over which the grid is valid (in GHz)
    taugrid=((218.15,218.25,taugrid1a),(218.7,218.8,taugrid2a)),
    hdr=hdra,
    shortvarnames=("T","N","n","v","\\sigma"), # specify the parameter names (TeX is OK)
    grid_vwidth=5.0,
)

#Step 3: Make mask
#get spectra, into units of K and GHz
tmax=pyfits.getdata(datadir+'T_max_maps/H2CO_303_202_tmax.fits')
tmax2=pyfits.getdata(datadir+'T_max_maps/H2CO_322_221_tmax.fits')
tmax3=pyfits.getdata(datadir+'T_max_maps/H2CO_321_220_tmax.fits')
#create mask - makes a 2D map- first tmax>=0, then erode to get rid of noisy edges, then sigma cut it
badmask=tmax>=0
#pixels that are NaN in the original data
#erode the badmask edge by thismuch to get rid of edges
rad=150
data=nd.morphology.binary_erosion(badmask, np.ones((rad, rad)))
#sigma cut- 5*5mJy,rms=25 mJy or 0.025Jy or 0.25K
keep=(tmax*data)>0.25#0.025
#print tmax pixels
maxim=np.nanmax(data*tmax*keep)#finds max value
maxind=np.nanargmax(data*tmax*keep)#finds position of max - 1d
ymax,xmax=np.unravel_index(maxind, tmax.shape)
print(xmax, ymax)
#plot to make sure mask is good for all three lines
print 'show mask...'
fig=plt.figure()
im=plt.imshow(keep)
plt.colorbar(im)
plt.gca().invert_yaxis()
plt.show()

#Step 4: Maps and spectra
a=SpectralCube.read(datadir+'alex_imaging_H2CO_303_202_fix/GRS1915_modelimg_H2CO_303_202.image.pbcor.fits')
a.allow_huge_operations=True
a1=a.to(u.K, equivalencies=a.beam.jtok_equiv(a.header['RESTFRQ']*u.Hz))
a2=a1.with_spectral_unit(u.Hz)
#a3=a2.with_mask(keep)
a3=a1.with_spectral_unit(u.km/u.s, velocity_convention='radio')

b=SpectralCube.read(datadir+'alex_imaging_H2CO_322_221_fix/GRS1915_modelimg_H2CO_322_221.image.pbcor.fits')
b.allow_huge_operations=True
b1=b.to(u.K, equivalencies=b.beam.jtok_equiv(b.header['RESTFRQ']*u.Hz))
b2=b1.with_spectral_unit(u.Hz)
#b3=b2.with_mask(keep)
b3=b1.with_spectral_unit(u.km/u.s, velocity_convention='radio')

c=SpectralCube.read(datadir+'alex_imaging_H2CO_321_220_fix/GRS1915_modelimg_H2CO_321_220.image.pbcor.fits')
c.allow_huge_operations=True
c1=c.to(u.K, equivalencies=c.beam.jtok_equiv(c.header['RESTFRQ']*u.Hz))
c2=c1.with_spectral_unit(u.Hz)
#c3=c2.with_mask(keep)
c3=c1.with_spectral_unit(u.km/u.s, velocity_convention='radio')

#create a table to store the spectrum parameters
column_names=['x', 'y', 'temp', 'column', 'density', 'center', 'width', 'temp errors', 'column errors', 'density errors', 'center errors', 'width errors']
column_types=['f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4']
table=Table(names=column_names, dtype=column_types)

print 'testing at peak pixel...'
#test at max pixel
w=ymax
v=xmax
v1=np.concatenate((a2[:, w, v].value, b2[:, w, v].value, c2[:, w, v].value))
w1=np.concatenate((a2[:, w, v].spectral_axis.value, b2[:, w, v].spectral_axis.value,c2[:, w, v].spectral_axis.value))
#need to take only values in concatenation and add units with xarrkwargs below
w1 /= 1e9
print 'showing spectra in frequency space where it will be fit...'
plt.figure()
plt.plot(w1, v1, label='test', color='b', ls='-',lw=3)
plt.ylabel('$T_B$ (K)')
plt.xlabel('$v$ (km/s)')
plt.show()

sp=pyspeckit.Spectrum(data=v1, xarr=w1, unit="$T_B$ (K)",xarrkwargs={'unit':'GHz'})
#add (register) fitters to the spectrum
sp.Registry.add_fitter('formaldehyde_mm_radex', formaldehyde_radex_fitter_both, 5)
#plot fit for all 3 ('both')
sp.plotter(figure=1)
sp.specfit(fittype='formaldehyde_mm_radex',
	guesses=[95, 14.5, 4.0, 67, 0.8],
	limits=[(50,550), (11,17), (3,8.5), (65,70), (0.5,10)],
	limited=[(True, True)]*5,
	fixed=[False, False, True, False, False])
#sp.plotter.savefig(datadir+'T_max_maps/AtestH2CO_spec_'+str(v)+'_'+str(w)+'.png')
print 'fit at peak pixel (', v,w, ') saved.'
print 'temp', sp.specfit.modelpars[0], '+/-',sp.specfit.modelerrs[0]
print 'column', sp.specfit.modelpars[1], '+/-',sp.specfit.modelerrs[1]
print 'density', sp.specfit.modelpars[2], '+/-',sp.specfit.modelerrs[2]
print 'cen', sp.specfit.modelpars[3], '+/-',sp.specfit.modelerrs[3]
print 'width', sp.specfit.modelpars[4], '+/-',sp.specfit.modelerrs[4]
raw_input('Press enter to continue with full fit')


#vlo=237
#vup=462
#wlow=404
#wup=618
print 'Beginning fitting loop...'
for v in range(0, a2.shape[1]):
	for w in range(0, a2.shape[2]):
		#combine three lines into one spectrum by concatenating numpy arrays, y-y axis,x-x axis 
		if keep[w,v]==True:
			v1=np.concatenate((a2[:, w, v].value, b2[:, w, v].value, c2[:, w, v].value))
			w1=np.concatenate((a2[:, w, v].spectral_axis.value, b2[:, w, v].spectral_axis.value,c2[:, w, v].spectral_axis.value))
			#need to take only values in concatenation and add units with xarrkwargs below
			w1 /= 1e9
			sp=pyspeckit.Spectrum(data=v1, xarr=w1, unit="$T_B$ (K)",xarrkwargs={'unit':'GHz'})
			#add (register) fitters to the spectrum
			sp.Registry.add_fitter('formaldehyde_mm_radex', formaldehyde_radex_fitter_both, 5)
			#plot fit for all 3 ('both')
			sp.plotter(figure=1)
			sp.specfit(fittype='formaldehyde_mm_radex',	guesses=[95, 14.5, 4.0, 67, 0.8], limits=[(50,550), (11,17), (3,8.5), (65,70), (0.5,10)], limited=[(True, True)]*5, fixed=[False, False, True, False, False])
			sp.plotter.savefig(datadir+'T_max_maps/figs/testH2CO_spec_'+str(v)+'_'+str(w)+'.png')
			#print v,w
			table.add_row()
			table[-1]['x']=v
			table[-1]['y']=w
			table[-1]['temp']=sp.specfit.modelpars[0]
			table[-1]['column']=sp.specfit.modelpars[1]
			table[-1]['density']=sp.specfit.modelpars[2]
			table[-1]['center']=sp.specfit.modelpars[3]
			table[-1]['width']=sp.specfit.modelpars[4]
			table[-1]['temp errors']=sp.specfit.modelerrs[0]
			table[-1]['column errors']=sp.specfit.modelerrs[1]
			table[-1]['density errors']=sp.specfit.modelerrs[2]
			table[-1]['center errors']=sp.specfit.modelerrs[3]
			table[-1]['width errors']=sp.specfit.modelerrs[4]
		else:
			table.add_row()
			table[-1]['x']=v
			table[-1]['y']=w
			table[-1]['temp']=np.nan
			table[-1]['column']=np.nan
			table[-1]['density']=np.nan
			table[-1]['center']=np.nan
			table[-1]['width']=np.nan
			table[-1]['temp errors']=np.nan
			table[-1]['column errors']=np.nan
			table[-1]['density errors']=np.nan
			table[-1]['center errors']=np.nan
			table[-1]['width errors']=np.nan
table.write(datadir+'T_max_maps/grs1915H2COparameters_fixdc.fits', overwrite=True)
t=Table.read(datadir+'T_max_maps/grs1915H2COparameters_fixdc.fits')
#test image of temp result
plt.figure()
plt.rcdefaults()
plt.scatter(t['x'], t['y'], c=t['temp'], s=2*t['temp'],edgecolor='none')
plt.colorbar(label='Temperature(K)')
plt.xlim(0,a2.shape[1])
plt.ylim(0,a2.shape[2])
plt.savefig(datadir+'T_max_maps/grs1915H2COtempmap.pdf')
plt.show()


#ra/dec versions
#temp
i='H2CO_303_202'
fits_file1=datadir+'T_max_maps/'+i+'_tmax.fits'
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
#data1=hdulist1.data
wmap1=wcs.WCS(hdulist1.header)
coord0=SkyCoord('19h15m40.8s','+10d40m58s',frame='icrs')
coord1=SkyCoord('19h15m37.1s','+10d41m46s',frame='icrs')
x1=float(wmap1.wcs_world2pix(coord0.ra.value,coord0.dec.value,0,0,1)[0])
y1=float(wmap1.wcs_world2pix(coord0.ra.value,coord0.dec.value,0,0,1)[1])
x2=float(wmap1.wcs_world2pix(coord1.ra.value,coord1.dec.value,0,0,1)[0])
y2=float(wmap1.wcs_world2pix(coord1.ra.value,coord1.dec.value,0,0,1)[1])

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
im=plt.imshow(np.swapaxes(np.reshape(t['temp'],(a2.shape[1],a2.shape[2])),0,1),origin="lower",cmap=cm.get_cmap('hot_r', 500),norm=colors.PowerNorm(gamma=1.),vmin=0.0)#,vmin=0.0,vmax=2.)10,2
cbar=plt.colorbar(im, orientation='vertical',fraction=0.04,pad=0)
cbar.set_label('Temperature (K)')
ax1.tick_params(axis='both', which='major', labelsize=15,width=3,length=7,color='k')
ax1.tick_params(axis='both', which='minor', labelsize=15,width=1,length=7,color='k')
ax1.coords['ra'].set_axislabel('Right Ascension')
ax1.coords['dec'].set_axislabel('Declination',minpad=-0.1)
ax1.coords['ra'].set_major_formatter('hh:mm:ss.s')
ax1.set_ylim(y1, y2)
ax1.set_xlim(x1, x2)
plt.contour(X,Y,Z,levels,colors='k',transform=ax1.get_transform(wmap))
plt.savefig(datadir+'for_paper/H2CO_fit_temp_contour.pdf',bbox_inches='tight')
plt.show()

#column
fig=plt.figure()
plt.rcdefaults()
plt.rc('xtick.major', size=4)
#plt.rc('xtick', color='w', labelsize='large')
ax1 = fig.add_subplot(111, projection=wmap1.celestial)
im=plt.imshow(np.swapaxes(np.reshape(t['column'],(a2.shape[1],a2.shape[2])),0,1),origin="lower",cmap=cm.get_cmap('hot_r', 500),norm=colors.PowerNorm(gamma=1.),vmin=13,vmax=15)#,vmin=0.0,vmax=2.)10,2
cbar=plt.colorbar(im, orientation='vertical',fraction=0.04,pad=0)
cbar.set_label('Column')
ax1.tick_params(axis='both', which='major', labelsize=15,width=3,length=7,color='k')
ax1.tick_params(axis='both', which='minor', labelsize=15,width=1,length=7,color='k')
ax1.coords['ra'].set_axislabel('Right Ascension')
ax1.coords['dec'].set_axislabel('Declination',minpad=-0.1)
ax1.coords['ra'].set_major_formatter('hh:mm:ss.s')
ax1.set_ylim(y1, y2)
ax1.set_xlim(x1, x2)
plt.contour(X,Y,Z,levels,colors='k',transform=ax1.get_transform(wmap))
plt.savefig(datadir+'other_data/H2CO_fit_column_contour.pdf',bbox_inches='tight')
plt.show()

#density
fig=plt.figure()
plt.rcdefaults()
plt.rc('xtick.major', size=4)
#plt.rc('xtick', color='w', labelsize='large')
ax1 = fig.add_subplot(111, projection=wmap1.celestial)
im=plt.imshow(np.swapaxes(np.reshape(t['density'],(a2.shape[1],a2.shape[2])),0,1),origin="lower",cmap=cm.get_cmap('hot_r', 500),norm=colors.PowerNorm(gamma=1.),vmin=0.0)#,vmin=0.0,vmax=2.)10,2
cbar=plt.colorbar(im, orientation='vertical',fraction=0.04,pad=0)
cbar.set_label('Temperature (K)')
ax1.tick_params(axis='both', which='major', labelsize=15,width=3,length=7,color='k')
ax1.tick_params(axis='both', which='minor', labelsize=15,width=1,length=7,color='k')
ax1.coords['ra'].set_axislabel('Right Ascension')
ax1.coords['dec'].set_axislabel('Declination',minpad=-0.1)
ax1.coords['ra'].set_major_formatter('hh:mm:ss.s')
ax1.set_ylim(y1, y2)
ax1.set_xlim(x1, x2)
plt.contour(X,Y,Z,levels,colors='k',transform=ax1.get_transform(wmap))
plt.savefig(datadir+'other_data/H2CO_fit_density_contour.pdf',bbox_inches='tight')
plt.show()

#TEMP CONTOURS ON VLA
i='H2CO_303_202'
fits_file1=datadir+'T_max_maps/'+i+'_tmax.fits'
hdulist = fits.open('/mnt/bigdata/tetarenk/VLA_grs1915_images/GRSVLA.fits')[0]
data=hdulist.data
wmap=wcs.WCS(hdulist.header)
hdulist1 = fits.open(fits_file1)[0]
hdulist1.header.remove('CRPIX3')
hdulist1.header.remove('CRVAL3')
hdulist1.header.remove('CDELT3')
hdulist1.header.remove('CUNIT3')
hdulist1.header.remove('CTYPE3')
hdulist1.header.remove('CRPIX4')
hdulist1.header.remove('CRVAL4')
hdulist1.header.remove('CDELT4')
hdulist1.header.remove('CUNIT4')
hdulist1.header.remove('CTYPE4')
hdulist1.header['WCSAXES']=2
wmap1=wcs.WCS(hdulist1.header)
#coord0=SkyCoord('19h15m41.3s','+10d41m01s',frame='icrs')
#coord1=SkyCoord('19h15m36.9s','+10d41m49s',frame='icrs')
coord0=SkyCoord('19h15m40.8s','+10d40m58s',frame='icrs')
coord1=SkyCoord('19h15m37.1s','+10d41m46s',frame='icrs')
x1=float(wmap.wcs_world2pix(coord0.ra.value,coord0.dec.value,0,0,1)[0])
y1=float(wmap.wcs_world2pix(coord0.ra.value,coord0.dec.value,0,0,1)[1])
x2=float(wmap.wcs_world2pix(coord1.ra.value,coord1.dec.value,0,0,1)[0])
y2=float(wmap.wcs_world2pix(coord1.ra.value,coord1.dec.value,0,0,1)[1])
coord3=SkyCoord('19h15m37.2s','+10d41m00s',frame='icrs')
x3=float(wmap.wcs_world2pix(coord3.ra.value,coord3.dec.value,0,0,1)[0])
y3=float(wmap.wcs_world2pix(coord3.ra.value,coord3.dec.value,0,0,1)[1])


x=np.arange(0,a2.shape[1])
y=np.arange(0,a2.shape[2])
X, Y = np.meshgrid(x, y)
Z=np.swapaxes(np.reshape(t['temp'],(a2.shape[1],a2.shape[2])),0,1)#[0,0,490:550,470:550]
levels=np.array([15,25,50,100,125,150,175,200,300])
fig=plt.figure()
plt.rcdefaults()
plt.rc('xtick.major', size=4)
ax1 = fig.add_subplot(111, projection=wmap.celestial)
im=plt.imshow(np.nan_to_num(data),origin="lower",cmap=cm.get_cmap('jet', 500),norm=colors.PowerNorm(gamma=0.7),vmin=0.0)#65,55,65,0.9,0.9,0.9,0.9/x,x,x,0.1,0.03,0.025,0.025,0.025,0.1
cbar=plt.colorbar(im, orientation='vertical',fraction=0.04,pad=0)
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
plt.contour(X,Y,Z,levels,colors='k',transform=ax1.get_transform(wmap1))
plt.savefig(datadir+'other_data/H2CO_vla_contour.pdf',bbox_inches='tight')
plt.show()



