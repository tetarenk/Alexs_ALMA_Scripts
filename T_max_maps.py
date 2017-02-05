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
#from mpl_toolkits.mplot3d.axes3d import Axes3D
#import aplpy, atpy
import matplotlib.cm as cm
from astropy import wcs
from astropy.coordinates import SkyCoord
import matplotlib.patches as patches
import matplotlib
import matplotlib.colors as colors
import seaborn as sns
import matplotlib as mpl


def colors_maps():
	CBcdict={
		'Bl':(0,0,0),
		'Or':(.9,.6,0),
		'SB':(.35,.7,.9),
		'bG':(0,.6,.5),
		'Ye':(.95,.9,.25),
		'Bu':(0,.45,.7),
		'Ve':(.8,.4,0),
		'rP':(.8,.6,.7),
	}

	##Single color gradient maps
	def lighter(colors):
		li=lambda x: x+.5*(1-x)
		return (li(colors[0]),li(colors[1]),li(colors[2]))

	def darker(colors):
		return (.5*colors[0],.5*colors[1],.5*colors[2])

	CBLDcm={}
	for key in CBcdict:
		    CBLDcm[key]=matplotlib.colors.LinearSegmentedColormap.from_list('CMcm'+key,[lighter(CBcdict[key]),darker(CBcdict[key])])


	##Two color gradient maps
	CB2cm={}
	for key in CBcdict:
		for key2 in CBcdict:
		    if key!=key2: CB2cm[key+key2]=matplotlib.colors.LinearSegmentedColormap.from_list('CMcm'+key+key2,[CBcdict[key],CBcdict[key2]])

	##Two color gradient maps with white in the middle
	CBWcm={}
	for key in CBcdict:
		for key2 in CBcdict:
		    if key!=key2: CBWcm[key+key2]=matplotlib.colors.LinearSegmentedColormap.from_list('CMcm'+key+key2,[CBcdict[key],(1,1,1),CBcdict[key2]])

	##Two color gradient maps with Black in the middle
	CBBcm={}
	for key in CBcdict:
		for key2 in CBcdict:
		    if key!=key2: CBBcm[key+key2]=matplotlib.colors.LinearSegmentedColormap.from_list('CMcm'+key+key2,[CBcdict[key],(0,0,0),CBcdict[key2]])
	return(CB2cm)
datadir = '/mnt/bigdata/tetarenk/ALMA_GRS1915_105/'
lines=['13CO','12CO','18CO','H2CO_303_202','H2CO_322_221','H2CO_321_220','SiO','N2Dp','H30a']

def make_TMAX(fitsfile,ddir,line):
	header = fits.getdata(fitsfile, header=True)[1]
	readfile=SpectralCube.read(fitsfile)
	filekms=readfile.with_spectral_unit(u.km/u.s,velocity_convention='radio')
	tmax=filekms.apply_numpy_function(np.nanmax,axis=0)
	print tmax
	fits.writeto(filename=ddir+'T_max_maps/'+line+'_tmax.fits',output_verify='ignore',\
	clobber=True,data=tmax,header=header)


#for i in lines:
	#make_TMAX(datadir+'alex_imaging_'+i+'/GRS1915_modelimg_'+i+'.image.pbcor.fits',datadir,i)
#lines=['13CO','12CO','18CO','H2CO_303_202','H2CO_322_221','H2CO_321_220','SiO','N2Dp','H30a']
i='H30alpha'
make_TMAX(datadir+'alex_imaging_'+i+'_fix/GRS1915_modelimg_'+i+'.image.pbcor.fits',datadir,i)
#imview(datadir+'T_max_maps/'+i+'_tmax.fits')

#convert vla image so only 2 dimensional
'''filename='/mnt/bigdata/tetarenk/VLA_grs1915_images/grs1915_rob_ms_pb_tay2.fits'
fh = fits.open(filename)
data = fh[0].data.squeeze() # drops the size-1 axes
header = fh[0].header
mywcs = wcs.WCS(header)
new_header = mywcs.to_header()
new_fh = fits.PrimaryHDU(data=data, header=new_header)
new_fh.writeto('GRSVLA.fits')'''

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
data1=hdulist1.data
wmap1=wcs.WCS(hdulist1.header)
coord0=SkyCoord('19h15m40.8s','+10d40m58s',frame='icrs')
coord1=SkyCoord('19h15m37.1s','+10d41m44s',frame='icrs')
x1=float(wmap1.wcs_world2pix(coord0.ra.value,coord0.dec.value,0,0,1)[0])
y1=float(wmap1.wcs_world2pix(coord0.ra.value,coord0.dec.value,0,0,1)[1])
x2=float(wmap1.wcs_world2pix(coord1.ra.value,coord1.dec.value,0,0,1)[0])
y2=float(wmap1.wcs_world2pix(coord1.ra.value,coord1.dec.value,0,0,1)[1])


x=np.arange(0,len(data[0,:]))
y=np.arange(0,len(data[:,0]))
X, Y = np.meshgrid(x, y)
Z=data#[0,0,490:550,470:550]
levels=np.array([1,2,3,4,5,6,7])*0.000345
#evels=np.array([1,2,3,4,5,6,7])*0.000345
#sns.set_style("dark")
#cmap1 = mpl.colors.ListedColormap(sns.color_palette("colorblind",10))
#cmap2=colors_maps()
fig=plt.figure()
plt.rcdefaults()
plt.rc('xtick.major', size=4)
plt.rc('xtick', color='w', labelsize='large')
ax1 = fig.add_subplot(111, projection=wmap1.celestial)
im=plt.imshow(np.nan_to_num(data1),origin="lower",cmap=cm.get_cmap('jet', 500),norm=colors.PowerNorm(gamma=1.),vmin=0.,vmax=0.025)#65,55,65,0.9,0.9,0.9,0.9/x,x,x,0.1,0.045,0.065
cbar=plt.colorbar(im, orientation='vertical',fraction=0.04,pad=0)
cbar.set_label('Jy/beam')
ax1.tick_params(axis='both', which='major', labelsize=15,width=3,length=7,color='k')
ax1.tick_params(axis='both', which='minor', labelsize=15,width=1,length=7,color='k')
ax1.coords['ra'].set_axislabel('Right Ascension')
ax1.coords['dec'].set_axislabel('Declination',minpad=-0.1)
ax1.coords['ra'].set_major_formatter('hh:mm:ss.s')
ax1.set_ylim(150, 700)
ax1.set_xlim(100, 650)
plt.contour(X,Y,Z,levels,colors='k',transform=ax1.get_transform(wmap))
plt.savefig(datadir+'T_max_maps/'+i+'_tmax_contour.pdf',bbox_inches='tight')
plt.show()

import scipy.ndimage as nd
rad=125#100,125,125,150,150,100
#data_mask=nd.morphology.binary_erosion(np.nan_to_num(data1), np.ones((rad, rad)))
#grsmask=data>1*0.00035
'''import os
os.system('cp -r '+datadir+'T_max_maps/'+i+'_tmax_contour.pdf /home/ubuntu/Dropbox')'''




