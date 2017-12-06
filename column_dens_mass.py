'''Creates column density maps'''


from spectral_cube import SpectralCube
from astropy import units as u
import astropy.constants as con
import numpy as np
from astropy.io import fits
import pyfits
import pylab as pl
import math as ma
import matplotlib.pyplot as plt
from astropy.table import hstack
from mpl_toolkits.mplot3d.axes3d import Axes3D
from mpl_toolkits.axes_grid1.inset_locator import mark_inset,zoomed_inset_axes
#import aplpy, atpy
import matplotlib.cm as cm
from astropy import wcs
from astropy.coordinates import SkyCoord
import matplotlib.patches as patches
import matplotlib
import matplotlib.colors as colors
import mpmath


datadir = '/mnt/bigdata/tetarenk/ALMA_GRS1915_105/'
lines=['13CO','12CO','18CO','H2CO_303_202','H2CO_322_221','H2CO_321_220','SiO','N2Dp','H30a']

#mass estimate from H2 column
def mass_est(H2col_fits_file,cell,dist,lims):
	x1,x2,y1,y2=lims[0],lims[1],lims[2],lims[3]
	#dist in kpc cell in arcsec
	col=np.nan_to_num(pyfits.getdata(H2col_fits_file))
	col_data=col[x1:x2,y1:y2]
	col_data=col_data0[col_data0>0]
	numpix=1.
	print col_data.shape
	numarea=numpix*(cell*(dist*1000*3e18)/206265)**2
	mass=np.nansum(col_data)*(2*1.00794*1.6e-24)*(numarea)/(1.989e33)
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
def makeVarr(fitsfile):
	header = fits.getdata(fitsfile, header=True)[1]
	X=SpectralCube.read(fitsfile)
	X.allow_huge_operations=True
	X1=X.to(u.K, equivalencies=X.beam.jtok_equiv(X.header['RESTFRQ']*u.Hz))
	Xkms=X1.with_spectral_unit(u.km/u.s,velocity_convention='radio')
	Vel=np.array(Xkms.spectral_axis)
	Ts=np.array(Xkms)
	return(Ts,Vel)
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
	Tbg=2.73
	C12=h*nu/k
	Cbg=1./(np.exp((h*nu)/(k*Tbg))-1)
	Tex=C12/(np.log(1+(C12/(Tmax0+C12*Cbg))))
	return(Tex)
def Tau_2to1(nu1,TV,Tex):
	#nu1 is target frequency molecule in ghz,TV is array of target molecule T in K for a bunch of dV
	unitdict={'h':6.626068e-27,'k':1.3806503e-16,'c':2.99792458e10}
	h,k,c=unitdict['h'],unitdict['k'],unitdict['c']
	nu=nu1*1e9
	Tbg=2.73
	C18=h*nu/k
	A=1./(np.exp(C18/Tex)-1)
	B=1./(np.exp(C18/Tbg)-1)
	tau=-np.log(1-((TV/C18)*((A-B)**(-1))))
	return(tau)
def NJ(nu1,Tex,Tau_array,dV):
	#TV is array of target molecule T in K for a bunch of dV (array also in km/s)
	unitdict={'h':6.626068e-27,'k':1.3806503e-16,'c':2.99792458e10}
	h,k,c=unitdict['h'],unitdict['k'],unitdict['c']
	nu=nu1*1e9
	prefac = (8 * np.pi * 1*u.GHz**3 / (con.c**3)) * u.km/u.s * u.s
	C=(prefac.to(u.cm**-2)).value
	col=C*(3./2.)*(1.165e11)*(0.122)**(-2)*(1-np.exp(-h*nu/(k*Tex)))**(-1)
	inte=np.trapz(Tau_array,dV,axis=0)
	return(col*inte)
def Zpart(Tex_array):
	ZZ=np.zeros((Tex_array.shape[0],Tex_array.shape[1]))
	for kk in range(0,(Tex_array.shape[0])):
		for ll in range(0,(Tex_array.shape[1])):
			val=Tex_array[kk,ll]
			ZZ[kk,ll]=mpmath.nsum(lambda x: (2*x+1)*mpmath.exp(-2.65*x*(x+1)/val), [0,6])
	return(ZZ)
def Ntot_2to1(Tex,Nj):
	ZZ=Zpart(Tex)/(3.)
	NN=Nj*ZZ*np.exp(2.65*2./Tex)
	return(NN)
Tex_2to1vec=np.vectorize(Tex_2to1)
Tau_2to1vec=np.vectorize(Tau_2to1)



#per spectra version

fits_file_cube18=datadir+'alex_imaging_'+'18CO'+'_fix/GRS1915_modelimg_'+'18CO'+'.image.pbcor.fits'
fits_file_cube12=datadir+'alex_imaging_'+'12CO'+'_fix/GRS1915_modelimg_'+'12CO'+'.image.pbcor.fits'
fits_file_cube13=datadir+'alex_imaging_'+'13CO'+'_fix/GRS1915_modelimg_'+'13CO'+'.image.pbcor.fits'

make_TMAX(fits_file_cube12,datadir,'12CO')
make_TMAX(fits_file_cube13,datadir,'13CO')
tmax12_fits=datadir+'columns_plots/12CO_tmax.fits'
tmax13_fits=datadir+'columns_plots/13CO_tmax.fits'
make_moment_maps(fits_file_cube18,0,datadir,'18CO')
integ18_fits=datadir+'columns_plots/18CO_moment0.fits'

hdtmax12=fits.open(tmax12_fits)[0]
hdtmax13=fits.open(tmax13_fits)[0]
hdmom18=fits.open(integ18_fits)[0]
dat12=np.nan_to_num(hdtmax12.data)
dat13=np.nan_to_num(hdtmax13.data)
dat18m=np.nan_to_num(hdmom18.data)
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
raw_input('s')
Tex_array=Tex_2to1vec(fr,tm)
TV,dV=makeVarr(fits_file_cube18)
Tau_array=Tau_2to1vec(fr18,TV,Tex_array)
NJ_array=NJ(fr18,Tex_array,Tau_array,dV)
Nfact_array=Ntot_2to1(Tex_array,NJ_array)


colum=Nfact_array
header = fits.getdata(tmax12_fits, header=True)[1]
fits.writeto(filename=datadir+'columns_plots/18CO_column2.fits',output_verify='ignore',\
	clobber=True,data=colum,header=header)
facH2=np.empty((dat12.shape[0],dat12.shape[1]))
facH2.fill(1.5e6)
columH2=colum*facH2
fits.writeto(filename=datadir+'columns_plots/H2_column2.fits',output_verify='ignore',\
	clobber=True,data=columH2,header=header)
Mass=mass_est(datadir+'columns_plots/H2_column2.fits',0.2,6.,[250.75,417,433.25,599.5])[0]

print 'The mass in the selected region is ', "%.2f" % Mass,'solar masses'
raw_input('stop')
#ridge-[379,411,544,577]-602 solar masses

#only VLA contours
fits_file1=datadir+'columns_plots/18CO_column2.fits'
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
xx=np.arange(0,len(data1[0,:]))
yy=np.arange(0,len(data1[:,0]))
XX, YY = np.meshgrid(xx, yy)
ZZ=data1#[0,0,490:550,470:550]
levelss=np.array([5.5,6,7,8])*1e15
#evels=np.array([1,2,3,4,5,6,7])*0.000345
fig=plt.figure()
plt.rcdefaults()
plt.rc('xtick.major', size=4)
#plt.rc('xtick', color='w', labelsize='large')
ax1 = fig.add_subplot(111, projection=wmap1.celestial)
im=plt.imshow(np.nan_to_num(data1)/1e16,origin="lower",cmap=cm.get_cmap('hot_r', 500),\
norm=colors.PowerNorm(gamma=1),vmin=0.0,vmax=0.64)
cbar=plt.colorbar(im, orientation='vertical',fraction=0.04,pad=0)
cbar.set_label('$N_{{\\rm C}^{18}{\\rm O}}\\,\\times 10^{16} \\,{\\rm cm}^{-2}$')
ax1.tick_params(axis='both', which='major', labelsize=15,width=3,length=7,color='k')
ax1.tick_params(axis='both', which='minor', labelsize=15,width=1,length=7,color='k')
ax1.coords['ra'].set_axislabel('Right Ascension')
ax1.coords['dec'].set_axislabel('Declination',minpad=-0.1)
ax1.coords['ra'].set_major_formatter('hh:mm:ss.s')
ax1.set_xlim(200, 500)
ax1.set_ylim(380, 650)
plt.contour(X,Y,Z,levels,colors='w',transform=ax1.get_transform(wmap),linewidths=1.5)
plt.contour(XX,YY,ZZ,levelss,colors='#00FF00',linewidths=1.5)
plt.savefig(datadir+'for_paper/18CO_column_contour2.pdf',bbox_inches='tight')
plt.show()



#VLA and column density contours
'''Get rid of extra WCS axes for contours'''
#filename=datadir+'columns_plots/18CO_column.fits'
#fh = fits.open(filename)
#data = fh[0].data.squeeze() # drops the size-1 axes
#header = fh[0].header
#mywcs = wcs.WCS(header)
#new_header = mywcs.to_header()
#new_fh = fits.PrimaryHDU(data=data, header=new_header)
#new_fh.writeto(datadir+'columns_plots/18CO_column_squeeze.fits')
fits_file1=datadir+'columns_plots/18CO_column_squeeze.fits'
hdulist1 = fits.open(fits_file1)[0]
hdulist1.header.remove('CRPIX3')
hdulist1.header.remove('CRVAL3')
hdulist1.header.remove('CDELT3')
hdulist1.header.remove('CUNIT3')
hdulist1.header.remove('CTYPE3')
hdulist1.header.remove('CRPIX4')
hdulist1.header.remove('CRVAL4')
hdulist1.header.remove('CDELT4')
#hdulist.header.remove('CUNIT4')
hdulist1.header.remove('CTYPE4')
hdulist1.header['WCSAXES']=2
data1=hdulist1.data
wmap1=wcs.WCS(hdulist1.header)
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
x1=float(wmap.wcs_world2pix(coord0.ra.value,coord0.dec.value,1)[0])
y1=float(wmap.wcs_world2pix(coord0.ra.value,coord0.dec.value,1)[1])
x2=float(wmap.wcs_world2pix(coord1.ra.value,coord1.dec.value,1)[0])
y2=float(wmap.wcs_world2pix(coord1.ra.value,coord1.dec.value,1)[1])

levelss=np.array([5.25,5.5,6,6.5,7,8,9])*1e15
fig=plt.figure()
plt.rcdefaults()
plt.rc('xtick.major', size=4)
#plt.rc('xtick', color='w', labelsize='large')
ax1 = fig.add_subplot(111, projection=wmap.celestial)
#im=plt.imshow(np.nan_to_num(data[:,:])*10000.,origin="lower",cmap=cm.get_cmap('hot_r', 500),\
#norm=colors.PowerNorm(gamma=0.7),vmin=0,vmax=3.)
#cbar=plt.colorbar(im, orientation='vertical',fraction=0.04,pad=0)
#cbar.set_label('${\\rm mJy\, bm}^{-1}$')
ax1.tick_params(axis='both', which='major', labelsize=15,width=3,length=7,color='k')
ax1.tick_params(axis='both', which='minor', labelsize=15,width=1,length=7,color='k')
ax1.coords['ra'].set_axislabel('Right Ascension')
ax1.coords['dec'].set_axislabel('Declination',minpad=-0.1)
ax1.coords['ra'].set_major_formatter('hh:mm:ss.s')
ax1.set_xlim(x1,x2)
ax1.set_ylim(y1,y2)
plt.contour(X,Y,Z,levels,colors='k',linewidths=1.5)
plt.contour(XX,YY,ZZ,levelss,colors='m',transform=ax1.get_transform(wmap1),linewidths=2)
ax1.set_aspect('equal', 'datalim')
plt.savefig(datadir+'for_paper/18CO_column_contour_vla.pdf',bbox_inches='tight')
plt.show()


#H2 version, vla contours only
fits_file1=datadir+'columns_plots/H2_column2.fits'
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
im=plt.imshow(np.nan_to_num(data1)/1e22,origin="lower",cmap=cm.get_cmap('hot_r', 500),norm=colors.PowerNorm(gamma=1),vmax=3,vmin=0.0)
cbar=plt.colorbar(im, orientation='vertical',fraction=0.04,pad=0)
cbar.set_label('$N_{{\\rm H}2}\\,\\times 10^{22} \\,{\\rm cm}^{-2}$')
ax1.tick_params(axis='both', which='major', labelsize=15,width=3,length=7,color='k')
ax1.tick_params(axis='both', which='minor', labelsize=15,width=1,length=7,color='k')
ax1.coords['ra'].set_axislabel('Right Ascension')
ax1.coords['dec'].set_axislabel('Declination',minpad=-0.1)
ax1.coords['ra'].set_major_formatter('hh:mm:ss.s')
ax1.set_ylim(150, 700)
ax1.set_xlim(100, 650)
plt.contour(X,Y,Z,levels,colors='w',transform=ax1.get_transform(wmap))
axins=zoomed_inset_axes(ax1,2,loc=4)
im2=axins.imshow(np.nan_to_num(data1)/1e22,origin="lower",cmap=cm.get_cmap('hot_r', 500),norm=colors.PowerNorm(gamma=1),vmax=3,vmin=0.0)
#plt.contour(X,Y,Z,levels,colors='w',transform=axins.get_transform(wmap))
axins.set_ylim(550,655)
axins.set_xlim(342,497)
plt.yticks(visible=False)
plt.xticks(visible=False)
plt.setp(axins,xticks=[],yticks=[])
mark_inset(ax1,axins,loc1=2,loc2=1,fc="none",ec="0.4",lw=2)
#plt.savefig(datadir+'for_paper/H2_column_contour_lobe.pdf',bbox_inches='tight')
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


