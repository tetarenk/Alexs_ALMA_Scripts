'''Other continuum maps'''

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
from mpl_toolkits.mplot3d.axes3d import Axes3D
#import aplpy, atpy
import matplotlib.cm as cm
from astropy import wcs
from astropy.coordinates import SkyCoord
import matplotlib.patches as patches
import matplotlib
import matplotlib.colors as colors
import glob
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredEllipse

datadir = '/mnt/bigdata/tetarenk/ALMA_GRS1915_105/'

'''imfiles=[]
for name in glob.glob(datadir+'other_data/*_regrid.fits'):
	if 'hspire' in name:
		imfiles.append(name.strip('.fits')+'_fix.fits')
	elif 'ukids' in name:
		imfiles.append(name.strip('.fits')+'_fix.fits')
	elif 'tay2' in name:
		continue
	else:
    	imfiles.append(name)'''
#VLA contours
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
#levels=np.array([1,2,3,4,5,6,7])*0.000345
levels=np.array([4,6,8,10,15,20,40,60])*0.00005



imfiles=['/mnt/bigdata/tetarenk/ALMA_GRS1915_105/other_data/236_471_1_1424_1_ukidssj_regrid_fix.fits',\
'/mnt/bigdata/tetarenk/ALMA_GRS1915_105/other_data/v2.1_ds2_l045_13pca_medmap20_crop_bolocam_regrid.fits',\
'/mnt/bigdata/tetarenk/ALMA_GRS1915_105/other_data/hspireplw529_herschel_regrid_fix.fits',\
 '/mnt/bigdata/tetarenk/ALMA_GRS1915_105/other_data/hspirepmw529_herschel_regrid_fix.fits',\
 '/mnt/bigdata/tetarenk/ALMA_GRS1915_105/other_data/hspirepsw529_herschel_regrid_fix.fits',\
'/mnt/bigdata/tetarenk/ALMA_GRS1915_105/other_data/residual_MG0450n005_024_all_mipsgal_regrid.fits',\
'/mnt/bigdata/tetarenk/ALMA_GRS1915_105/other_data/2892p106_ab41-w4-int-3_ra288.9125_dec10.6875_asec300.000.fits',\
'/mnt/bigdata/tetarenk/ALMA_GRS1915_105/other_data/2892p106_ab41-w3-int-3_ra288.9125_dec10.6875_asec300.000.fits',\
'/mnt/bigdata/tetarenk/ALMA_GRS1915_105/other_data/GLM_04550-0075_mosaic_I4_glimpse_regrid.fits',\
'/mnt/bigdata/tetarenk/ALMA_GRS1915_105/other_data/GLM_04550-0075_mosaic_I3_glimpse_regrid.fits',\
'/mnt/bigdata/tetarenk/ALMA_GRS1915_105/other_data/2892p106_ab41-w2-int-3_ra288.9125_dec10.6875_asec300.000.fits',\
'/mnt/bigdata/tetarenk/ALMA_GRS1915_105/other_data/GLM_04550-0075_mosaic_I2_glimpse_regrid.fits',\
'/mnt/bigdata/tetarenk/ALMA_GRS1915_105/other_data/GLM_04550-0075_mosaic_I1_glimpse_regrid.fits',\
'/mnt/bigdata/tetarenk/ALMA_GRS1915_105/other_data/2892p106_ab41-w1-int-3_ra288.9125_dec10.6875_asec300.000.fits',\
'/mnt/bigdata/tetarenk/ALMA_GRS1915_105/other_data/236_471_1_1424_3_ukidssk_regrid_fix.fits',\
'/mnt/bigdata/tetarenk/ALMA_GRS1915_105/other_data/236_471_1_1424_2_ukidssh_regrid_fix.fits']
label2=['UKIDSS','BOLOCAM','HSPIRE','HSPIRE','HSPIRE','MIPSGAL','WISE','WISE','GLIMPSE','GLIMPSE','WISE','GLIMPSE','GLIMPSE',\
'WISE','UKIDSS','UKIDSS']
label3=['UKIRT','CSO','Herschel','Herschel','Herschel','Spitzer','NASA','NASA','Spitzer','Spitzer','NASA','Spitzer','Spitzer',\
'NASA','UKIRT','UKIRT']
label=['$1.33\\mu m$','$1.11mm$', '$500\\mu m$','$350\\mu m$','$250\\mu m$','$24\\mu m$','$22\\mu m$','$12\\mu m$','$8\\mu m$',\
'$5.8\\mu m$','$4.6\\mu m$','$4.5\\mu m$','$3.6\\mu m$','$3.4\\mu m$','$2.37\\mu m$','$1.78\\mu m$']

#all images
'''fits_file1='/mnt/bigdata/tetarenk/VLA_grs1915_images/GRSVLA.fits'
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
minv=[920,0,100,260,521,0,180,1000,27,4,20,0,0.2,20,2485,6013]
maxv=[23500,0.6,360,1300,4000,2000,900,6000,5000,5000,2600,3000,2000,2400,32000,54000]
pl=[0.35,1,1,1,1,0.45,0.8,0.5,0.4,0.32,0.4,0.24,0.25,0.5,0.35,0.35]
for i in range(0,len(imfiles)):
	fits_file1=imfiles[i]
	hdulist1 = fits.open(fits_file1)[0]
	if i == 2:
		data1=hdulist1.data*(9.382805e-7)*(35.2**2)/(0.2**2)
		minv0=0#minv[i]*(9.382805e-7)*(35.2**2)/(0.2**2)
		maxv0=10.5#maxv[i]*(9.382805e-7)*(35.2**2)/(0.2**2)
	elif i ==3:
		data1=hdulist1.data*(9.382805e-7)*(23.9**2)/(0.2**2)
		minv0=0#minv[i]*(9.382805e-7)*(23.9**2)/(0.2**2)
		maxv0=18.5#maxv[i]*(9.382805e-7)*(23.9**2)/(0.2**2)
	elif i ==4:
		data1=hdulist1.data*(9.382805e-7)*(17.6**2)/(0.2**2)
		minv0=0#minv[i]*(9.382805e-7)*(17.6**2)/(0.2**2)
		maxv0=maxv[i]*(9.382805e-7)*(17.6**2)/(0.2**2)
	else:
		data1=hdulist1.data
	wmap1=wcs.WCS(hdulist1.header)
	coord0=SkyCoord('19h15m41.8s','+10d40m38s',frame='icrs')
	coord1=SkyCoord('19h15m34.9s','+10d41m52s',frame='icrs')
	x1=float(wmap1.wcs_world2pix(coord0.ra.value,coord0.dec.value,1)[0])
	y1=float(wmap1.wcs_world2pix(coord0.ra.value,coord0.dec.value,1)[1])
	x2=float(wmap1.wcs_world2pix(coord1.ra.value,coord1.dec.value,1)[0])
	y2=float(wmap1.wcs_world2pix(coord1.ra.value,coord1.dec.value,1)[1])
	ax=plt.subplot(4,4,i, projection=wmap1.celestial)
	if i in [2,3,4]:
		im=plt.imshow(np.nan_to_num(data1),origin="lower",cmap=cm.get_cmap('jet', 500),norm=colors.PowerNorm(gamma=pl[i]),vmin=minv0,vmax=maxv0)
	else:
		im=plt.imshow(np.nan_to_num(data1),origin="lower",cmap=cm.get_cmap('jet', 500),norm=colors.PowerNorm(gamma=pl[i]),vmin=minv[i],vmax=maxv[i])
	cbar=plt.colorbar(im, orientation='vertical',fraction=0.035,pad=0)
	cbar.ax.tick_params(labelsize=8)

	if i==5:
		cbar.set_ticks([10,100,300,600,900,1200,1500])
		cbar.set_ticklabels([10,100,300,600,900,1200,1500])
	elif i==1:
		cbar.set_ticks([0,0.12,0.24,0.36,0.48,0.6])
		cbar.set_ticklabels([0,0.12,0.24,0.36,0.48,0.6])
	elif i==8:
		cbar.set_ticks([100,300,900,1800,2700,3600,4500])
		cbar.set_ticklabels([100,300,900,1800,2700,3600,4500])
	elif i==9:
		cbar.set_ticks([100,300,900,1800,2700,3600,4500])
		cbar.set_ticklabels([100,300,900,1800,2700,3600,4500])
	elif i==11:
		cbar.set_ticks([1,10,50,100,300,600,1200,2200])
		cbar.set_ticklabels([1,10,50,100,300,600,1200,2200])
	elif i==12:
		cbar.set_ticks([0.5,5,25,100,200,500,1000,1700])
		cbar.set_ticklabels([0.5,5,25,100,200,500,1000,1700])
	elif i==14:
		cbar.set_ticks([1000,4000,8000,12000,16000,20000,25000])
		cbar.set_ticklabels([1000,4000,8000,12000,16000,20000,25000])
	elif i==15:
		cbar.set_ticks([1000,10000,15000,20000,30000,40000,50000])
		cbar.set_ticklabels([1000,10000,15000,20000,30000,40000,50000])
	elif i==0:
		cbar.set_ticks([500,1500,3000,6000,10000,15000,20000])
		cbar.set_ticklabels([500,1500,3000,6000,10000,15000,20000])
	elif i==6:
		cbar.set_ticks([200,400,600,800])
		cbar.set_ticklabels([200,400,600,800])
	elif i==7:
		cbar.set_ticks([1200,2000,3000,4000,5000,5500])
		cbar.set_ticklabels([1200,2000,3000,4000,5000,5500])
	elif i==10:
		cbar.set_ticks([20,100,300,600,900,1300,1800,2200,2500])
		cbar.set_ticklabels([20,100,300,600,900,1300,1800,2200,2500])
	elif i==13:
		cbar.set_ticks([20,100,300,600,900,1300,1800,2200])
		cbar.set_ticklabels([20,100,300,600,900,1300,1800,2200])
	plt.contour(X,Y,Z,levels,colors='w',transform=ax.get_transform(wmap))
#cbar.set_label('Jy')
	plt.tick_params(axis='both', which='major', labelsize=1,width=3,length=7,color='k')
	plt.tick_params(axis='both', which='minor', labelsize=1,width=1,length=7,color='k')
	if i in [2,3,4,6,7,8,10,11,12]:
		ax.coords['ra'].set_axislabel('')
		ax.coords['dec'].set_axislabel('',minpad=-0.1)
		ax.coords['ra'].set_major_formatter('hh:mm:ss')
		ax.coords['ra'].set_ticklabel_visible(False)
		ax.coords['dec'].set_ticklabel_visible(False)
	elif i in [14,15,0]:
		ax.coords['ra'].set_axislabel('Right Ascension')
		ax.coords['dec'].set_axislabel('',minpad=-0.1)
		ax.coords['ra'].set_major_formatter('hh:mm:ss')
		ax.coords['ra'].set_ticklabel_visible(True)
		ax.coords['dec'].set_ticklabel_visible(False)
	elif i ==13:
		ax.coords['ra'].set_axislabel('Right Ascension')
		ax.coords['dec'].set_axislabel('Declination',minpad=-0.1)
		ax.coords['ra'].set_major_formatter('hh:mm:ss')
		ax.coords['ra'].set_ticklabel_visible(True)
		ax.coords['dec'].set_ticklabel_visible(True)
		
	elif i in [1,5,9]:
		ax.coords['ra'].set_axislabel('')
		ax.coords['dec'].set_axislabel('Declination',minpad=-0.1)
		ax.coords['ra'].set_major_formatter('hh:mm:ss')
		ax.coords['ra'].set_ticklabel_visible(False)
		ax.coords['dec'].set_ticklabel_visible(True)
	ax.set_ylim(y1,y2)
	ax.set_xlim(x1,x2)
	if i in [2,3]:
		ax.text(30,475,label[i],color='k')
		ax.text(350,475,label2[i],color='k')
		ax.text(350,445,label3[i],color='k')
	elif i in [6,7,10,13]:
		ax.text(85,130,label[i],color='w')
		ax.text(130,130,label2[i],color='w')
		ax.text(130,125,label3[i],color='w')
	else:
		ax.text(30,475,label[i],color='w')
		ax.text(350,475,label2[i],color='w')
		ax.text(350,445,label3[i],color='w')
	#ax.set_aspect('equal')
	#ax.set_xticklabels([])
	#ax.set_yticklabels([])
	#plt.setp(ax.get_yticklabels(),visible=False)
plt.subplots_adjust(wspace=0.2,hspace=0.1,top=0.9,bottom=0.1)
plt.savefig(datadir+'other_data/'+'all'+'_contour.pdf',bbox_inches='tight')
plt.show()'''



#vla continuum
fits_file1='/mnt/bigdata/tetarenk/VLA_grs1915_images/grs1915_rob_ms_nopb_tay2_redo_dec16.image.tt0.fits'
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
coord3=SkyCoord('19h15m37.3s','+10d41m01s',frame='icrs')
x3=float(wmap1.wcs_world2pix(coord3.ra.value,coord3.dec.value,0,0,1)[0])
y3=float(wmap1.wcs_world2pix(coord3.ra.value,coord3.dec.value,0,0,1)[1])


x=np.arange(0,len(data[0,:]))
y=np.arange(0,len(data[:,0]))
X, Y = np.meshgrid(x, y)
Z=data
levels=np.array([4,6,8,10,15,20,40,60])*0.00005
fig=plt.figure()
plt.rcdefaults()
plt.rc('xtick.major', size=4)
#plt.rc('xtick', color='w', labelsize='large')
ax1 = fig.add_subplot(111, projection=wmap1.celestial)
im=plt.imshow(np.nan_to_num(data1[0,0,:,:])*1000.,origin="lower",cmap=cm.get_cmap('hot_r', 500),norm=colors.PowerNorm(gamma=0.7),vmin=0.0,vmax=3.)#65,55,65,0.9,0.9,0.9,0.9/x,x,x,0.1,0.03,0.025,0.025,0.025,0.1
cbar=plt.colorbar(im, orientation='vertical',fraction=0.04,pad=0)
cbar.set_label('mJy/beam')
ax1.tick_params(axis='both', which='major', labelsize=15,width=3,length=7,color='k')
ax1.tick_params(axis='both', which='minor', labelsize=15,width=1,length=7,color='k')
ax1.coords['ra'].set_axislabel('Right Ascension')
ax1.coords['dec'].set_axislabel('Declination',minpad=-0.1)
ax1.coords['ra'].set_major_formatter('hh:mm:ss.s')
ax1.set_ylim(y1,y2)
ax1.set_xlim(x1,x2)
e1 = patches.Ellipse((x3,y3), 5.31, 4.50,angle=-50.1963, linewidth=2, fill=False,color='m')
ax1.add_patch(e1)
ae = AnchoredEllipse(ax1.transData, width=5.31, height=4.5, angle=-50.1963,loc=4, pad=0.5, borderpad=0.4, frameon=True)
ax1.add_artist(ae)
plt.contour(X,Y,Z,levels,colors='k',transform=ax1.get_transform(wmap))
plt.savefig(datadir+'for_paper/vla_contour.pdf',bbox_inches='tight')
plt.show()

# VLA spectral index map
fits_file1='/mnt/bigdata/tetarenk/VLA_grs1915_images/grs1915_rob_ms_nopb_tay2_redo_dec16.image.alpha.fits'
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
coord3=SkyCoord('19h15m37.3s','+10d41m01s',frame='icrs')
x3=float(wmap1.wcs_world2pix(coord3.ra.value,coord3.dec.value,0,0,1)[0])
y3=float(wmap1.wcs_world2pix(coord3.ra.value,coord3.dec.value,0,0,1)[1])


x=np.arange(0,len(data[0,:]))
y=np.arange(0,len(data[:,0]))
X, Y = np.meshgrid(x, y)
Z=data
levels=np.array([4,6,8,10,15,20,40,60])*0.00005
fig=plt.figure()
plt.rcdefaults()
plt.rc('xtick.major', size=4)
#plt.rc('xtick', color='w', labelsize='large')
ax1 = fig.add_subplot(111, projection=wmap1.celestial)
im=plt.imshow((data1[0,0,:,:]),origin="lower",cmap=cm.get_cmap('hot_r', 500),vmin=-1.,vmax=0.6)#,vmax=0.8)#65,55,65,0.9,0.9,0.9,0.9/x,x,x,0.1,0.03,0.025,0.025,0.025,0.1
cbar=plt.colorbar(im, orientation='vertical',fraction=0.04,pad=0)
cbar.set_label('Spectral Index')
ax1.tick_params(axis='both', which='major', labelsize=15,width=3,length=7,color='k')
ax1.tick_params(axis='both', which='minor', labelsize=15,width=1,length=7,color='k')
ax1.coords['ra'].set_axislabel('Right Ascension')
ax1.coords['dec'].set_axislabel('Declination',minpad=-0.1)
ax1.coords['ra'].set_major_formatter('hh:mm:ss.s')
ax1.set_ylim(y1,y2)
ax1.set_xlim(x1,x2)
e1 = patches.Ellipse((x3,y3), 5.31, 4.50,angle=-50.1963, linewidth=2, fill=False,color='m')
ax1.add_patch(e1)
ae = AnchoredEllipse(ax1.transData, width=5.31, height=4.5, angle=-50.1963,loc=4, pad=0.5, borderpad=0.4, frameon=True)
ax1.add_artist(ae)
plt.contour(X,Y,Z,levels,colors='k',transform=ax1.get_transform(wmap))
plt.savefig(datadir+'for_paper/vla_alpha_contour.pdf',bbox_inches='tight')
plt.show()

#new full plot for paper (4 ims only)

imfiles=['/mnt/bigdata/tetarenk/ALMA_GRS1915_105/other_data/hpacs_25HPPUNIMAPB_blue_1917_p1149_00_v1.0_1471608553396.fits',\
'/mnt/bigdata/tetarenk/ALMA_GRS1915_105/other_data/2892p106_ab41-w4-int-3_ra288.9125_dec10.6875_asec300.000.fits',\
'/mnt/bigdata/tetarenk/ALMA_GRS1915_105/other_data/GLM_04550-0075_mosaic_I4_glimpse_regrid.fits',\
'/mnt/bigdata/tetarenk/ALMA_GRS1915_105/other_data/236_471_1_1424_3_ukidssk_regrid_fix.fits']
label=['$70\\mu m$','$22\\mu m$','$8\\mu m$','$2.37\\mu m$']
label2=['PACS','WISE','GLIMPSE','UKIDSS']
label3=['Hershel','NASA','Spitzer','UKIRT']


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
levels=np.array([4,6,8,10,15,20,40,60])*0.00005

fig=plt.figure()
plt.rcdefaults()
plt.rc('xtick.major', size=4)
#plt.rc('xtick', color='w', labelsize='small')
minv=[0,180,27,2485]
maxv=[5,1100,5000,32000]
pl=[1,1,0.8,0.35]
for i in range(0,len(imfiles)):
	fits_file1=imfiles[i]
	if i==0:
		hdulist1 = fits.open(fits_file1)[1]
	else:
		hdulist1 = fits.open(fits_file1)[0]
	if i == 1:
		data1=hdulist1.data*(8.287e-9)*1e3*(12**2)/(1.375**2)
		minvs=minv[i]*(8.287e-9)*1e3*(12**2)/(1.375**2)
		maxvs=maxv[i]*(8.287e-9)*1e3*(12**2)/(1.375**2)
	elif i==0:
		data1=hdulist1.data*(6.6**2)/(3.2**2)#*1e3*(6.6**2)/(3.2**2)*(2.401998072e-4)
		minvs=minv[i]*(6.6**2)/(3.2**2)#*1e3*(6.6**2)/(3.2**2)*(2.401998072e-4)
		maxvs=maxv[i]*(6.6**2)/(3.2**2)#*1e3*(6.6**2)/(3.2**2)*(2.401998072e-4)
	elif i==2:
		data1=hdulist1.data*(9.382805e-7)*1e3*(2.5**2)/(1.2**2)
		minvs=minv[i]*(9.382805e-7)*1e3*(2.5**2)/(1.2**2)
		maxvs=20#maxv[i]*(9.382805e-7)*1e3*(2.5**2)/(1.2**2)
	elif i==3:
		data1=hdulist1.data*(9.382805e-7)*1e3*(0.4**2)/(0.2**2)
		minvs=minv[i]*(9.382805e-7)*1e3*(0.4**2)/(0.2**2)
		maxvs=100#maxv[i]*(9.382805e-7)*1e3*(0.4**2)/(0.2**2)
	else:
		data1=hdulist1.data
		minvs=minv[i]
		maxvs=maxv[i]
	wmap1=wcs.WCS(hdulist1.header)
	coord0=SkyCoord('19h15m41.8s','+10d40m38s',frame='icrs')
	coord1=SkyCoord('19h15m34.9s','+10d41m52s',frame='icrs')
	x1=float(wmap1.wcs_world2pix(coord0.ra.value,coord0.dec.value,1)[0])
	y1=float(wmap1.wcs_world2pix(coord0.ra.value,coord0.dec.value,1)[1])
	x2=float(wmap1.wcs_world2pix(coord1.ra.value,coord1.dec.value,1)[0])
	y2=float(wmap1.wcs_world2pix(coord1.ra.value,coord1.dec.value,1)[1])
	ax=plt.subplot(2,2,i+1, projection=wmap1.celestial)
	im=plt.imshow(np.nan_to_num(data1),origin="lower",cmap=cm.get_cmap('hot_r', 500),norm=colors.PowerNorm(gamma=pl[i]),vmin=minvs,vmax=maxvs)
	cbar=plt.colorbar(im, orientation='vertical',fraction=0.0353,pad=0)
	cbar.ax.tick_params(labelsize=8)
	if i==3:
		cbar.set_ticks([10,15,30,50,60,80])
		cbar.set_ticklabels([10,15,20,50,60,80])
	plt.contour(X,Y,Z,levels,colors='k',transform=ax.get_transform(wmap))
	plt.tick_params(axis='both', which='major', labelsize=1,width=3,length=7,color='k')
	plt.tick_params(axis='both', which='minor', labelsize=1,width=1,length=7,color='k')
	if i==1:
		ax.coords['ra'].set_axislabel('')
		ax.coords['dec'].set_axislabel('',minpad=-0.1)
		ax.coords['ra'].set_major_formatter('hh:mm:ss')
		ax.coords['ra'].set_ticklabel_visible(False)
		ax.coords['dec'].set_ticklabel_visible(False)
		#cbar.set_label('mJy/beam')
	elif i==3:
		ax.coords['ra'].set_axislabel('Right Ascension')
		ax.coords['dec'].set_axislabel('',minpad=-0.1)
		ax.coords['ra'].set_major_formatter('hh:mm:ss')
		ax.coords['ra'].set_ticklabel_visible(True)
		ax.coords['dec'].set_ticklabel_visible(False)
		#cbar.set_label('mJy/beam')
	elif i ==2:
		ax.coords['ra'].set_axislabel('Right Ascension')
		ax.coords['dec'].set_axislabel('Declination',minpad=-0.1)
		ax.coords['ra'].set_major_formatter('hh:mm:ss')
		ax.coords['ra'].set_ticklabel_visible(True)
		ax.coords['dec'].set_ticklabel_visible(True)
		#cbar.set_label('mJy/beam')
		
	elif i==0:
		ax.coords['ra'].set_axislabel('')
		ax.coords['dec'].set_axislabel('Declination',minpad=-0.1)
		ax.coords['ra'].set_major_formatter('hh:mm:ss')
		ax.coords['ra'].set_ticklabel_visible(False)
		ax.coords['dec'].set_ticklabel_visible(True)
		#cbar.set_label('Jy/beam')
	ax.set_ylim(y1,y2)
	ax.set_xlim(x1,x2)
	if i==0:
		ax.text(2195,700,label[i],color='k')
		ax.text(2215,700,label2[i],color='k')
		ax.text(2215,698,label3[i],color='k')
	if i==1:
		ax.text(85,132,label[i],color='k')
		ax.text(140,132,label2[i],color='k')
		ax.text(140,127,label3[i],color='k')
	else:
		ax.text(30,475,label[i],color='k')
		ax.text(350,475,label2[i],color='k')
		ax.text(350,445,label3[i],color='k')
	#ax.set_aspect('equal')
	#ax.set_xticklabels([])
	#ax.set_yticklabels([])
	#plt.setp(ax.get_yticklabels(),visible=False)
#plt.subplots_adjust(wspace=0.2,hspace=0.1,top=0.9,bottom=0.1)
plt.savefig(datadir+'for_paper/'+'all'+'_contour.pdf',bbox_inches='tight')
plt.show()

#pacs map
fits_file1='/mnt/bigdata/tetarenk/ALMA_GRS1915_105/other_data/hpacs_25HPPUNIMAPB_blue_1917_p1149_00_v1.0_1471608553396.fits'
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
hdulist1 = fits.open(fits_file1)[1]
data1=hdulist1.data*(6.6**2)/(3.2**2)
wmap1=wcs.WCS(hdulist1.header)
coord0=SkyCoord('19h15m45.907s','+10d40m07.796s',frame='icrs')
coord1=SkyCoord('19h15m34.491s','+10d42m59.3s',frame='icrs')
x1=float(wmap1.wcs_world2pix(coord0.ra.value,coord0.dec.value,1)[0])
y1=float(wmap1.wcs_world2pix(coord0.ra.value,coord0.dec.value,1)[1])
x2=float(wmap1.wcs_world2pix(coord1.ra.value,coord1.dec.value,1)[0])
y2=float(wmap1.wcs_world2pix(coord1.ra.value,coord1.dec.value,1)[1])
coord3=SkyCoord('19h15m37.3s','+10d41m01s',frame='icrs')
x3=float(wmap1.wcs_world2pix(coord3.ra.value,coord3.dec.value,1)[0])
y3=float(wmap1.wcs_world2pix(coord3.ra.value,coord3.dec.value,1)[1])


x=np.arange(0,len(data[0,:]))
y=np.arange(0,len(data[:,0]))
X, Y = np.meshgrid(x, y)
Z=data
levels=np.array([4,6,8,10,15,20,40,60])*0.00005
xp2=np.arange(0,len(data1[0,:]))
yp2=np.arange(0,len(data1[:,0]))
X2, Y2 = np.meshgrid(xp2, yp2)
Z2=data1
levels2=np.array([0.35,0.5,0.6,0.8,1,2,4,6])*1.
fig=plt.figure()
plt.rcdefaults()
plt.rc('xtick.major', size=4)
#plt.rc('xtick', color='w', labelsize='large')
ax1 = fig.add_subplot(111, projection=wmap1.celestial)
im=plt.imshow(np.nan_to_num(data1[:,:]),origin="lower",cmap=cm.get_cmap('hot', 500),norm=colors.PowerNorm(gamma=0.45),vmin=0*(6.6**2)/(3.2**2),vmax=5*(6.6**2)/(3.2**2))#,vmax=0.8)#65,55,65,0.9,0.9,0.9,0.9/x,x,x,0.1,0.03,0.025,0.025,0.025,0.1
cbar=plt.colorbar(im, orientation='vertical',fraction=0.04,pad=0)
cbar.set_label('Jy/beam')
ax1.tick_params(axis='both', which='major', labelsize=15,width=3,length=7,color='k')
ax1.tick_params(axis='both', which='minor', labelsize=15,width=1,length=7,color='k')
ax1.coords['ra'].set_axislabel('Right Ascension')
ax1.coords['dec'].set_axislabel('Declination',minpad=-0.1)
ax1.coords['ra'].set_major_formatter('hh:mm:ss.s')
ax1.set_ylim(y1,y2)
ax1.set_xlim(x1,x2)
#ax1.text(470,550,'4-8 GHz',color='k')
#e1 = patches.Ellipse((x3,y3), 5.31, 4.50,angle=-50.1963, linewidth=2, fill=False,color='m')
#ax1.add_patch(e1)
#ae = AnchoredEllipse(ax1.transData, width=5.31, height=4.5, angle=-50.1963,loc=4, pad=0.5, borderpad=0.4, frameon=True)
#ax1.add_artist(ae)
plt.contour(X2,Y2,Z2,levels2,colors='w')
plt.contour(X,Y,Z,levels,colors='k',transform=ax1.get_transform(wmap))
plt.savefig(datadir+'for_paper/pacs_contour.pdf',bbox_inches='tight')
plt.show()

