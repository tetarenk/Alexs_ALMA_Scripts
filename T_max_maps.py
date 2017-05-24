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
	fileKkms=filekms.to(u.K, equivalencies=filekms.beam.jtok_equiv(filekms.header['RESTFRQ']*u.Hz))
	tmax=fileKkms.apply_numpy_function(np.nanmax,axis=0)
	print tmax
	fits.writeto(filename=ddir+'T_max_maps/'+line+'_tmax.fits',output_verify='ignore',\
	clobber=True,data=tmax,header=header)

i='H2CO_303_202'
#make all Tmax maps
#for i in lines:
	#make_TMAX(datadir+'alex_imaging_'+i+'/GRS1915_modelimg_'+i+'.image.pbcor.fits',datadir,i)
#lines=['13CO','12CO','18CO','H2CO_303_202','H2CO_322_221','H2CO_321_220','SiO','N2Dp','H30a']
#i='H30alpha'
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

i='12CO'
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
#coord0=SkyCoord('19h15m40.8s','+10d40m58s',frame='icrs')
#coord1=SkyCoord('19h15m37.1s','+10d41m44s',frame='icrs')
coord0=SkyCoord('19h15m40.8s','+10d40m58s',frame='icrs')
coord1=SkyCoord('19h15m37.1s','+10d41m46s',frame='icrs')

x1=float(wmap1.wcs_world2pix(coord0.ra.value,coord0.dec.value,0,0,1)[0])
y1=float(wmap1.wcs_world2pix(coord0.ra.value,coord0.dec.value,0,0,1)[1])
x2=float(wmap1.wcs_world2pix(coord1.ra.value,coord1.dec.value,0,0,1)[0])
y2=float(wmap1.wcs_world2pix(coord1.ra.value,coord1.dec.value,0,0,1)[1])

import scipy.ndimage as nd
rad=80#100,125,125,150,150,100
#data_mask=nd.morphology.binary_erosion(np.nan_to_num(data1), np.ones((rad, rad)))

x=np.arange(0,len(data[0,:]))
y=np.arange(0,len(data[:,0]))
X, Y = np.meshgrid(x, y)
Z=data#[0,0,490:550,470:550]
levels=np.array([4,6,8,10,15,20,40,60])*0.00005
#evels=np.array([1,2,3,4,5,6,7])*0.000345
#sns.set_style("dark")
#cmap1 = mpl.colors.ListedColormap(sns.color_palette("colorblind",10))
#cmap2=colors_maps()
fig=plt.figure()
plt.rcdefaults()
plt.rc('xtick.major', size=4)
#plt.rc('xtick', color='w', labelsize='large')
ax1 = fig.add_subplot(111, projection=wmap1.celestial)
im=plt.imshow(np.nan_to_num(data1),origin="lower",cmap=cm.get_cmap('hot_r', 500),norm=colors.PowerNorm(gamma=1),vmin=0.,vmax=25)#,vmax=1.2)#65,55,65,0.9,0.9,0.9,0.9/x,x,x,0.1,0.045,0.065
cbar=plt.colorbar(im, orientation='vertical',fraction=0.04,pad=0)
cbar.set_label('K')
ax1.tick_params(axis='both', which='major', labelsize=15,width=3,length=7,color='k')
ax1.tick_params(axis='both', which='minor', labelsize=15,width=1,length=7,color='k')
ax1.coords['ra'].set_axislabel('Right Ascension')
ax1.coords['dec'].set_axislabel('Declination',minpad=-0.1)
ax1.coords['ra'].set_major_formatter('hh:mm:ss.s')
ax1.set_ylim(150, 700)
ax1.set_xlim(100, 650)
ax1.set_ylim(y1, y2)
ax1.set_xlim(x1, x2)
from matplotlib.patches import Rectangle
r1=Rectangle((385, 560), 15, 15,alpha=1, facecolor='none',edgecolor='white',lw=2,ls='solid',zorder=5)
r2=Rectangle((400, 570), 20, 25,alpha=1, facecolor='none',edgecolor='white',lw=2,ls='solid',zorder=5)
r3=Rectangle((340, 470), 20, 15,alpha=1, facecolor='none',edgecolor='white',lw=2,ls='solid',zorder=5)
r4=Rectangle((507, 447), 73, 45,alpha=1, facecolor='none',edgecolor='white',lw=2,ls='solid',zorder=5)
'''ax1.add_patch(r1)
ax1.add_patch(r2)
ax1.add_patch(r3)
ax1.add_patch(r4)
ax1.text(365,560, 'A',fontsize=15)
ax1.text(428,575, 'B',fontsize=15)
ax1.text(368,485, 'C',fontsize=15)
ax1.text(585,492, 'D',fontsize=15)'''
plt.contour(X,Y,Z,levels,colors='k',transform=ax1.get_transform(wmap))
plt.savefig(datadir+'for_paper/'+i+'_tmax_contour.pdf',bbox_inches='tight')
plt.show()

import scipy.ndimage as nd
rad=125#100,125,125,150,150,100
#data_mask=nd.morphology.binary_erosion(np.nan_to_num(data1), np.ones((rad, rad)))
#grsmask=data>1*0.00035
'''import os
os.system('cp -r '+datadir+'T_max_maps/'+i+'_tmax_contour.pdf /home/ubuntu/Dropbox')'''

#all lines big plot
'''imfiles=['/mnt/bigdata/tetarenk/ALMA_GRS1915_105/T_max_maps/12CO_tmax.fits',\
'/mnt/bigdata/tetarenk/ALMA_GRS1915_105/T_max_maps/13CO_tmax.fits',\
'/mnt/bigdata/tetarenk/ALMA_GRS1915_105/T_max_maps/18CO_tmax.fits',\
'/mnt/bigdata/tetarenk/ALMA_GRS1915_105/T_max_maps/H2CO_303_202_tmax.fits',\
'/mnt/bigdata/tetarenk/ALMA_GRS1915_105/T_max_maps/H2CO_321_220_tmax.fits',\
'/mnt/bigdata/tetarenk/ALMA_GRS1915_105/T_max_maps/H2CO_322_221_tmax.fits',\
'/mnt/bigdata/tetarenk/ALMA_GRS1915_105/T_max_maps/H30alpha_tmax.fits',\
'/mnt/bigdata/tetarenk/ALMA_GRS1915_105/T_max_maps/N2Dp_tmax.fits',\
'/mnt/bigdata/tetarenk/ALMA_GRS1915_105/T_max_maps/SiO_tmax.fits']

label=['$^{12}{\\rm CO}$','$^{13}{\\rm CO}$','$^{18}{\\rm CO}$','${\\rm H}_{2}{\\rm CO}(303-202)$',\
'${\\rm H}_{2}{\\rm CO}(321-220)$','${\\rm H}_{2}{\\rm CO}(322-221)$','${\\rm H}_{30}\\alpha$',\
'${\\rmN}_{2}{\\rm D}+$','$SiO$']

fits_file1='/mnt/bigdata/tetarenk/VLA_grs1915_images/GRSVLA.fits'
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
x=np.arange(0,len(data[0,:]))
y=np.arange(0,len(data[:,0]))
X, Y = np.meshgrid(x, y)
Z=data#[0,0,490:550,470:550]
levels=np.array([1,2,3,4,5,6,7])*0.000345

fig=plt.figure(figsize=(17,10))
#plt.rcdefaults()
plt.rc('xtick.major', size=4)
plt.rc('xtick', color='w', labelsize='small')
minv=[0,0,0,0,0,0,0,0,0]
maxv=[2.4,1.2,0.54,0.1,0.045,0.045,0.064,0.064,0.064]

for i in range(0,len(imfiles)):
	fits_file1=imfiles[i]
	hdulist1 = fits.open(fits_file1)[0]
	data1=hdulist1.data
	wmap1=wcs.WCS(hdulist1.header)
	coord0=SkyCoord('19h15m41.8s','+10d40m38s',frame='icrs')
	coord1=SkyCoord('19h15m34.9s','+10d41m52s',frame='icrs')
	ax=plt.subplot(3,3,i+1, projection=wmap1.celestial)
	if i==3:
		data_mask=nd.morphology.binary_erosion(np.nan_to_num(data1), np.ones((100, 100)))
		data_mask2=nd.morphology.binary_erosion(np.nan_to_num(data1), np.ones((150, 150)))
	if i in [3,4,5,7,8]:
		#data_mask=nd.morphology.binary_erosion(np.nan_to_num(data1), np.ones((100, 100)))
		im=plt.imshow(np.nan_to_num(data1)*data_mask,origin="lower",cmap=cm.get_cmap('jet', 500),norm=colors.PowerNorm(gamma=1.),vmin=minv[i],vmax=maxv[i])
	elif i==6:
		#data_mask=nd.morphology.binary_erosion(np.nan_to_num(data1), np.ones((150, 150)))
		im=plt.imshow(np.nan_to_num(data1)*data_mask2,origin="lower",cmap=cm.get_cmap('jet', 500),norm=colors.PowerNorm(gamma=1.),vmin=minv[i],vmax=maxv[i])	
	else:
		im=plt.imshow(np.nan_to_num(data1),origin="lower",cmap=cm.get_cmap('jet', 500),norm=colors.PowerNorm(gamma=1.),vmin=minv[i],vmax=maxv[i])
	cbar=plt.colorbar(im, orientation='vertical',fraction=0.035,pad=0)
	cbar.ax.tick_params(labelsize=8)
	plt.contour(X,Y,Z,levels,colors='k',transform=ax.get_transform(wmap))
	plt.tick_params(axis='both', which='major', labelsize=1,width=3,length=7,color='k')
	plt.tick_params(axis='both', which='minor', labelsize=1,width=1,length=7,color='k')
	if i in [1,2,4,5]:
		ax.coords['ra'].set_axislabel('')
		ax.coords['dec'].set_axislabel('',minpad=-0.1)
		ax.coords['ra'].set_major_formatter('hh:mm:ss')
		ax.coords['ra'].set_ticklabel_visible(False)
		ax.coords['dec'].set_ticklabel_visible(False)
	elif i in [7,8]:
		ax.coords['ra'].set_axislabel('Right Ascension')
		ax.coords['dec'].set_axislabel('',minpad=-0.1)
		ax.coords['ra'].set_major_formatter('hh:mm:ss')
		ax.coords['ra'].set_ticklabel_visible(True)
		ax.coords['dec'].set_ticklabel_visible(False)
	elif i ==6:
		ax.coords['ra'].set_axislabel('Right Ascension')
		ax.coords['dec'].set_axislabel('Declination',minpad=-0.1)
		ax.coords['ra'].set_major_formatter('hh:mm:ss')
		ax.coords['ra'].set_ticklabel_visible(True)
		ax.coords['dec'].set_ticklabel_visible(True)
		
	elif i in [0,3]:
		ax.coords['ra'].set_axislabel('')
		ax.coords['dec'].set_axislabel('Declination',minpad=-0.1)
		ax.coords['ra'].set_major_formatter('hh:mm:ss.s')
		ax.coords['ra'].set_ticklabel_visible(False)
		ax.coords['dec'].set_ticklabel_visible(True)
	ax.set_ylim(150, 700)
	ax.set_xlim(100, 650)
	if i in [2,3]:
		ax.text(150,170,label[i],color='w')
	else:
		ax.text(150,170,label[i],color='w')
plt.subplots_adjust(wspace=0.2,hspace=0.1)
plt.savefig(datadir+'T_max_maps/'+'all'+'_tmax.pdf',bbox_inches='tight')
plt.show()'''


