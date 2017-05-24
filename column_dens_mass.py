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


datadir = '/mnt/bigdata/tetarenk/ALMA_GRS1915_105/'
lines=['13CO','12CO','18CO','H2CO_303_202','H2CO_322_221','H2CO_321_220','SiO','N2Dp','H30a']

#mass estimate from H2 column
def mass_est(H2col_fits_file,cell,dist,lims):
	x1,x2,y1,y2=lims[0],lims[1],lims[2],lims[3]
	#dist in kpc cell in arcsec
	col=np.nan_to_num(pyfits.getdata(H2col_fits_file))
	col_data0=col[x1:x2,y1:y2]
	col_data=col_data0[col_data0>0]
	numpix=np.count_nonzero(col_data)
	print numpix
	numarea=numpix*(cell*(dist*1000*3e18)/206265)**2
	mass=np.nansum(col_data)*(1.00794*2*1.6e-24)*(numarea)/(1.989e33)
	print np.nansum(col_data),numarea
	return(mass,numpix,numarea)

def make_moment_maps(fitsfile,mind,ddir,line):
	header = fits.getdata(fitsfile, header=True)[1]
	X=SpectralCube.read(fitsfile)
	X.allow_huge_operations=True
	X1=X.to(u.K, equivalencies=X.beam.jtok_equiv(X.header['RESTFRQ']*u.Hz))
	Xkms=X1.with_spectral_unit(u.km/u.s,velocity_convention='radio')
	#Xs=Xkms.spectral_slab(0*u.km/u.s,67*u.km/u.s)
	#print Xs.spectral_axis
	mom=Xkms.moment(order=mind)
	mom_array=np.array(mom)
	fits.writeto(filename=ddir+'columns_plots/'+line+'_moment'+str(mind)+'.fits',output_verify='ignore',\
	clobber=True,data=mom_array,header=header)
def make_TMAX(fitsfile,ddir,line):
	header = fits.getdata(fitsfile, header=True)[1]
	readfile=SpectralCube.read(fitsfile)
	readfile.allow_huge_operations=True
	readfile1=readfile.to(u.K, equivalencies=readfile.beam.jtok_equiv(readfile.header['RESTFRQ']*u.Hz))
	filekms=readfile1.with_spectral_unit(u.km/u.s,velocity_convention='radio')
	tmax=filekms.apply_numpy_function(np.nanmax,axis=0)
	fits.writeto(filename=ddir+'columns_plots/'+line+'_tmax.fits',output_verify='ignore',\
	clobber=True,data=tmax,header=header)

def Tex_2to1(nu0,Tmax0):
	#nu0 is reference frequency molecule in ghz,Tmax0 is reference molecule max T in K
	unitdict={'h':6.626068e-27,'k':1.3806503e-16,'c':2.99792458e10}
	h,k,c=unitdict['h'],unitdict['k'],unitdict['c']
	nu=nu0*1e9
	Tbg=2.7
	C0=h*nu/k
	Cbg=1./(np.exp((h*nu)/(k*Tbg))-1)
	Tex=C0/(np.log(1+(C0/(Tmax0+C0*Cbg))))
	return(Tex)
def Tau_2to1(nu1,Tmax1,Tex):
	#nu1 is target frequency molecule in ghz,Tmax1 is target molecule max T in K
	unitdict={'h':6.626068e-27,'k':1.3806503e-16,'c':2.99792458e10}
	h,k,c=unitdict['h'],unitdict['k'],unitdict['c']
	nu=nu1*1e9
	Tbg=2.7
	C1=h*nu/k
	A=1./(np.exp(C1/Tex)-1)
	B=1./(np.exp(C1/Tbg)-1)
	tau=-np.log(1-((Tmax1/C1)*((A-B)**(-1))))
	return(tau)
def Ntot_factor_2to1(nu1,tau,Tex):
	unitdict={'h':6.626068e-27,'k':1.3806503e-16,'c':2.99792458e10,'Be':5.521711e10}
	h,k,c,Be=unitdict['h'],unitdict['k'],unitdict['c'],unitdict['Be']
	nu=nu1*1e9
	C1=h*nu/k
	C2=2*h*Be/k
	A=tau/(1-np.exp(-tau))
	B=(1.5e14)/(1-np.exp(-C1/Tex))
	cons=A*B*np.exp(C2/Tex)
	return(cons)
Tex_2to1vec=np.vectorize(Tex_2to1)
Tau_2to1vec=np.vectorize(Tau_2to1)
Ntot_factor_2to1vec=np.vectorize(Ntot_factor_2to1)


#per spectra version

fits_file_cube18=datadir+'alex_imaging_'+'18CO'+'_fix/GRS1915_modelimg_'+'18CO'+'.image.pbcor.fits'
fits_file_cube12=datadir+'alex_imaging_'+'12CO'+'_fix/GRS1915_modelimg_'+'12CO'+'.image.pbcor.fits'
fits_file_cube13=datadir+'alex_imaging_'+'13CO'+'_fix/GRS1915_modelimg_'+'13CO'+'.image.pbcor.fits'

make_TMAX(fits_file_cube12,datadir,'12CO')
make_TMAX(fits_file_cube13,datadir,'13CO')
make_TMAX(fits_file_cube18,datadir,'18CO')
tmax12_fits=datadir+'columns_plots/12CO_tmax.fits'
tmax13_fits=datadir+'columns_plots/13CO_tmax.fits'
tmax18_fits=datadir+'columns_plots/18CO_tmax.fits'
make_moment_maps(fits_file_cube18,0,datadir,'18CO')
integ18_fits=datadir+'columns_plots/18CO_moment0.fits'

hdtmax12=fits.open(tmax12_fits)[0]
hdtmax13=fits.open(tmax13_fits)[0]
hdtmax18=fits.open(tmax18_fits)[0]
hdmom18=fits.open(integ18_fits)[0]
dat12=np.nan_to_num(hdtmax12.data)
dat13=np.nan_to_num(hdtmax13.data)
dat18t=np.nan_to_num(hdtmax18.data)
dat18=np.nan_to_num(hdmom18.data)
freq12=230.538
freq13=220.39868
fr18=np.empty((dat12.shape[0],dat12.shape[1]))
fr18.fill(219.56035)
tm=np.zeros((dat12.shape[0],dat12.shape[1]))
fr=np.zeros((dat12.shape[0],dat12.shape[1]))
for kk in range(0,(dat12.shape[0])):
	for ll in range(0,(dat12.shape[1])):
		val12=dat12[kk,ll]
		val13=dat13[kk,ll]
		if val12 >val13:
			tm[kk,ll]=val12
			fr[kk,ll]=freq12
		else:
			tm[kk,ll]=val13
			fr[kk,ll]=freq13
Tex_array=Tex_2to1vec(fr,tm)
Tau_array=Tau_2to1vec(fr18,dat18t,Tex_array)
Nfact_array=Ntot_factor_2to1vec(fr18,Tau_array,Tex_array)

colum=dat18*Nfact_array
header = fits.getdata(tmax12_fits, header=True)[1]
fits.writeto(filename=datadir+'columns_plots/18CO_column.fits',output_verify='ignore',\
	clobber=True,data=colum,header=header)
facH2=np.empty((dat12.shape[0],dat12.shape[1]))
facH2.fill(2.65e21)
columH2=dat18*facH2
fits.writeto(filename=datadir+'columns_plots/H2_column.fits',output_verify='ignore',\
	clobber=True,data=columH2,header=header)
Mass=mass_est(datadir+'columns_plots/H2_column.fits',0.2,8.6,[346,398,510,569])[0]

print 'The mass in the selected region is ', "%.2f" % Mass,'solar masses'

#ridge-[379,411,544,577]-602 solar masses

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
#plt.rcdefaults()
plt.rc('xtick.major', size=4)
plt.rc('xtick', color='w', labelsize='large')
ax1 = fig.add_subplot(111, projection=wmap1.celestial)
im=plt.imshow(np.nan_to_num(data1)/1e16,origin="lower",cmap=cm.get_cmap('hot_r', 500),\
norm=colors.PowerNorm(gamma=1),vmin=0.0,vmax=1.2)
cbar=plt.colorbar(im, orientation='vertical',fraction=0.04,pad=0)
cbar.set_label('$N_{18{\\rm CO}}\\,\\times 10^{16} \\,{\\rm cm}^{-2}$')
ax1.tick_params(axis='both', which='major', labelsize=15,width=3,length=7,color='k')
ax1.tick_params(axis='both', which='minor', labelsize=15,width=1,length=7,color='k')
ax1.coords['ra'].set_axislabel('Right Ascension')
ax1.coords['dec'].set_axislabel('Declination',minpad=-0.1)
ax1.coords['ra'].set_major_formatter('hh:mm:ss.s')
ax1.set_ylim(150, 700)
ax1.set_xlim(100, 650)
plt.contour(X,Y,Z,levels,colors='k',transform=ax1.get_transform(wmap))
plt.savefig(datadir+'for_paper/18CO_column_contour.pdf',bbox_inches='tight')
plt.show()

#per spectra version
fits_file1=datadir+'columns_plots/H2_column.fits'
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
im=plt.imshow(np.nan_to_num(data1)/1e22,origin="lower",cmap=cm.get_cmap('hot_r', 500),norm=colors.PowerNorm(gamma=1),vmax=2,vmin=0.0)
cbar=plt.colorbar(im, orientation='vertical',fraction=0.04,pad=0)
cbar.set_label('$N_{{\\rm H}2}\\,\\times 10^{22} \\,{\\rm cm}^{-2}$')
ax1.tick_params(axis='both', which='major', labelsize=15,width=3,length=7,color='k')
ax1.tick_params(axis='both', which='minor', labelsize=15,width=1,length=7,color='k')
ax1.coords['ra'].set_axislabel('Right Ascension')
ax1.coords['dec'].set_axislabel('Declination',minpad=-0.1)
ax1.coords['ra'].set_major_formatter('hh:mm:ss.s')
ax1.set_ylim(150, 700)
ax1.set_xlim(100, 650)
plt.contour(X,Y,Z,levels,colors='k',transform=ax1.get_transform(wmap))
plt.savefig(datadir+'for_paper/H2_column_contour_lobe.pdf',bbox_inches='tight')
plt.show()

#cube version
'''hdulist118 = fits.open(fits_file_cube18)[0]
data118=hdulist118.data
#data118H=hdulist118.data
a=SpectralCube.read(fits_file_cube18)
a.allow_huge_operations=True
a1=a.to(u.K, equivalencies=a.beam.jtok_equiv(a.header['RESTFRQ']*u.Hz))
a2=a1.with_spectral_unit(u.km/u.s, velocity_convention='radio')
suba=a2[:,200,200]
veloc=np.array(suba.spectral_axis.to(u.km/u.s))

for kk in range(0,data118.shape[1]):
	data118[0,kk,:,:]=data118[0,kk,:,:]*Nfact_array*veloc[kk]
	data118H[0,kk,:,:]=data118H[0,kk,:,:]*facH2*veloc[kk]
col18_cube=np.nan_to_num(data118)
header = fits.getdata(fits_file_cube18, header=True)[1]
fits.writeto(filename=datadir+'columns_plots/18CO_cube_column.fits',output_verify='ignore',\
	clobber=True,data=data118,header=header)'''






