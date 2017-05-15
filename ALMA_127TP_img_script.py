###################################
#ALMA Spectral Line Imaging Script
###################################
'''CASA script to be used for imaging 12m + 7m + TP ALMA Spectral Line Data
INPUT: Parameter file detailing all data and imaging parameters (param_dir_file set below)
OUTPUT: (1) 12m + 7m + TP image (modelimage version) -- GRS1915_modelimg_[LINE].image.pbcor
		(2) 12m + 7m + TP image (feather version)-- GRS1915_Featherall7_[LINE].image.pbcor
		(3) 12m + 7m image -- GRS1915_12m7m_[LINE].image.pbcor
		(4) 12m + 7m + TP image (combined Kauffman method) -- GRS1915_modelimgnew_[LINE].image.pbcor
NOTES: - All output images & intermediate data products are put in my_dir directory set below.
	   - All output images are also converted to fits format (just append .fits to end of images 1,2,3,4 above)
	   - This script utilizes methods developed by a bunch of different people.
	     -More info on the Kauffman method here--> http://tinyurl.com/zero-spacing
	     -More info on feathering in the CASA Guide here-->https://casaguides.nrao.edu/index.php/M100_Band3_Combine_4.5
       -It is recommended to install the CASA-Python executable wrapper in order to install the 
        python packages needed in this script within CASA --> https://github.com/radio-astro-tools/casa-python
       -If using CASA 4.7 uncomment extra header fixes below after any call of imtrans, as apparently in CASA 4.7
        the imtrans task deletes the beam and units in the image headers.
Written by: Alex J. Tetarenko
Last Updated: Dec 15 2016'''

print '##################################################'
print 'Welcome to Alexs ALMA Spectral Line Imaging Script'
print '##################################################\n'


#packages you may need
from astropy.io import fits
from spectral_cube import SpectralCube
from astropy.wcs.utils import add_stokes_axis_to_wcs
import radio_beam
import re
import sys
import imp
import os
import linecache
import find
import warnings
warnings.filterwarnings('ignore')

#define output directory
my_dir='/mnt/bigdata/tetarenk/ALMA_GRS1915_105/alex_imaging_CS_fix/'
if not os.path.isdir(my_dir):
	os.system('sudo mkdir '+my_dir)
	os.system('sudo chown ubuntu '+my_dir)
	os.system('sudo chmod -R u+r '+my_dir) 
	os.system('sudo chmod -R u+w '+my_dir)
	os.system('sudo chmod -R u+x '+my_dir)
print 'You have set your output directory to ', my_dir
print 'All output images & intermediate data products are put in this directory.\n'

#param file location
param_dir_file='/mnt/bigdata/tetarenk/ALMA_GRS1915_105/alexs_scripts/params_CS.txt'
print 'You have set your param file to ', param_dir_file
print 'Please make sure all parameters are correct, they will change for each line!\n'

###################################
#Defining Parameters Section
###################################
print 'Reading in parameters...'
#read in param file
def getVar(filename):
	'''Easy way to read in variables from parameter file'''
	f=open(filename)
	global data_params
	data_params=imp.load_source('data_params','',f)
	f.close()
getVar(param_dir_file)
#data set params
data_12m=data_params.data_12m
data_7m=data_params.data_7m
data_TP=data_params.data_TP
field_12m=data_params.field_12m
field_7m=data_params.field_7m
numf12=str(int(field_12m.split('~')[1])-int(field_12m.split('~')[0]))
numf7=str(int(field_7m.split('~')[1])-int(field_7m.split('~')[0]))
remakems=data_params.remakems
remakemsconts=data_params.remakemsconts
#line params
line=data_params.line
spw_sd=data_params.spw_sd
spw_arr=data_params.spw_arr
spw_7m=data_params.spw_7m
#template image params
myniter_temp=data_params.myniter_temp
mynchan_temp=data_params.mynchan_temp
mystart_temp=data_params.mystart_temp
mywidth_temp=data_params.mywidth_temp
#general image params
mythreshold=data_params.mythreshold
myimsize=data_params.myimsize
mycell=data_params.mycell
myrestfreq=data_params.myrestfreq
#deep clean params
myniter=data_params.myniter
mynchan=data_params.mynchan
mystart=data_params.mystart
mywidth=data_params.mywidth
whichmethod=data_params.whichmethod
#optional params
if mynchan==-1:
	chans_select=''
else:
	chans_select=str(mystart)+'~'+str(mystart+mynchan-1)
###################################
print '**********************************************'
print 'You have chosen to image the ', line, ' line.'
print '**********************************************'
###################################
#Examining and Splitting Data Section
###################################
#listobs first
if not os.path.isfile(my_dir+'12m_listfile.txt'):
	print 'Making 12m listobs file...'
	listobs(data_params.data_12m,listfile=my_dir+'12m_listfile.txt')
else:
	print '12m listobs file already exists'
if not os.path.isfile(my_dir+'7m_listfile.txt'):
	print 'Making 7m listobs file...'
	listobs(data_params.data_7m,listfile=my_dir+'7m_listfile.txt')
else:
	print '7m listobs file already exists'
seelo=raw_input('Do you want to see the 12m and 7m listobs? y or n-->')
if seelo=='y':
	os.system('pluma '+my_dir+'12m_listfile.txt &')
	os.system('pluma '+my_dir+'7m_listfile.txt &')
	#raw_input('')
	raw_input('Please press enter to continue when you are done.')
	print 'Moving on to making split data sets.'
elif seelo=='n':
	print 'Okay. Moving on to making split data sets.'

#make split data sets if they dont exist
if not os.path.isdir(my_dir+'data_12m_'+line+'.ms'):
	print 'Making split 12m MS...'
	split(vis=data_12m,outputvis=my_dir+'data_12m_'+line+'.ms',spw=spw_arr,datacolumn='data',field=field_12m)
elif remakems=='T':
	print 'Remaking 12m MS...'
	os.system('rm -rf '+my_dir+'data_12m_'+line+'.ms')
	split(vis=data_12m,outputvis=my_dir+'data_12m_'+line+'.ms',spw=spw_arr,datacolumn='data',field=field_12m)
else:
	print 'Split 12m MS already exists'
if not os.path.isdir(my_dir+'data_7m_'+line+'.ms'):
	print 'Making split 7m MS...'
	split(vis=data_7m,outputvis=my_dir+'data_7m_'+line+'.ms',spw=spw_7m,datacolumn='data',field=field_7m)
elif remakems=='T':
	print 'Remaking 7m MS...'
	os.system('rm -rf '+my_dir+'data_7m_'+line+'.ms')
	split(vis=data_7m,outputvis=my_dir+'data_7m_'+line+'.ms',spw=spw_7m,datacolumn='data',field=field_7m)
else:
	print 'Split 7m MS already exists'
if not os.path.isdir(my_dir+'data_7m_'+line+'_comb.ms'):
	print 'Combining 7m spws...'
	mstransform(vis=my_dir+'data_7m_'+line+'.ms',outputvis=my_dir+'data_7m_'+line+'_comb.ms',combinespws=True,spw='',\
		datacolumn='data')
elif remakems=='T':
	print 'Remaking 7m combined MS...'
	os.system('rm -rf '+my_dir+'data_7m_'+line+'_comb.ms')
	mstransform(vis=my_dir+'data_7m_'+line+'.ms',outputvis=my_dir+'data_7m_'+line+'_comb.ms',combinespws=True,spw='',\
		datacolumn='data')
else:
	print 'Combined Split 7m MS already exists'
docheck=raw_input('Do you want to do checks on the 12m and combined 7m MSs? y or n-->')
if docheck=='y':
#checks
	print 'Displaying 7m combined listobs...'
	listobs(my_dir+'data_7m_'+line+'_comb.ms',listfile=my_dir+'7m_listfile_'+line+'_comb.txt',overwrite=True)
	os.system('pluma '+my_dir+'7m_listfile_'+line+'_comb.txt &')
	raw_input('Please press enter to continue with plots.')
	print 'First the amp/vel plot for 12m.'
	plotms(vis=my_dir+'data_12m_'+line+'.ms',yaxis='amp',xaxis='velocity',spw='', \
		avgtime='1e8',avgscan=True,coloraxis='spw',avgchannel='5',transform=True,\
		freqframe='LSRK',restfreq=myrestfreq, plotfile=my_dir+'A12m_'+line+'_vel.png',overwrite=True)
	raw_input('Please press enter to continue when you are done.')
	print 'Next the amp/vel plot for 7m.'
	plotms(vis=my_dir+'data_7m_'+line+'_comb.ms',yaxis='amp',xaxis='velocity',spw='', \
		avgtime='1e8',avgscan=True,coloraxis='spw',avgchannel='5',transform=True,\
		freqframe='LSRK',restfreq=myrestfreq, plotfile=my_dir+'A7m_'+line+'_vel.png',overwrite=True)
	raw_input('Please press enter to continue when you are done.')
	print 'Next the wt/uvdist plot for 12m.'
	plotms(vis=my_dir+'data_12m_'+line+'.ms',yaxis='wt',xaxis='uvdist',spw='',coloraxis='spw',\
		plotfile=my_dir+'A12m_'+line+'_WT.png',overwrite=True)
	raw_input('Please press enter to continue when you are done.')
	print 'Next the wt/uvdist plot for 7m.'
	plotms(vis=my_dir+'data_7m_'+line+'_comb.ms',yaxis='wt',xaxis='uvdist',spw='',coloraxis='spw',\
		plotfile=my_dir+'A7m_'+line+'_WT.png',overwrite=True)
	raw_input('Please press enter to continue when you are done.')
	print 'Next the amp/uvdist plot for 12m.'
	plotms(vis=my_dir+'data_12m_'+line+'.ms',yaxis='amp',xaxis='uvdist',spw='', avgscan=True,\
		avgchannel='5', coloraxis='spw',plotfile=my_dir+'A12m_'+line+'_UVdist.png',overwrite=True)
	raw_input('Please press enter to continue when you are done.')
	print 'Next the amp/uvdist plot for 7m.'
	plotms(vis=my_dir+'data_7m_'+line+'_comb.ms',yaxis='amp',xaxis='uvdist',spw='', avgscan=True,\
		avgchannel='5', coloraxis='spw',plotfile=my_dir+'A7m_'+line+'_UVdist.png',overwrite=True)
	raw_input('Please press enter to continue when you are done.')
	print 'Moving on.'
elif docheck=='n':
	print 'Okay. Moving on.'

#continuum subtraction if needed
docontsub=raw_input('Do you need to do continuum subtraction before imaging?y or n-->')
if docontsub=='y':
	print 'Plotting 12m spectra. Please use tool to identify channel range of line.'
	plotms(vis=my_dir+'data_12m_'+line+'.ms',yaxis='amp',xaxis='channel',spw='', \
		avgtime='1e8',avgscan=True,coloraxis='spw',plotfile=my_dir+'cs12m_'+line+'_chan.png',overwrite=True)
	raw_input('Please press enter to continue when you are done.')
	line_chans_12m=raw_input('Please enter line channel range for 12m; e.g.., 2~10.-->')
	orde12=raw_input('Please enter order for 12m cont fit (recomended 0)-->')
	print 'Next plot 7m spectra. Please use tool to identify channel range of line.'
	plotms(vis=my_dir+'data_7m_'+line+'_comb.ms',yaxis='amp',xaxis='channel',spw='', \
		avgtime='1e8',avgscan=True,coloraxis='spw',plotfile=my_dir+'cs7m_'+line+'_chan.png',overwrite=True)
	raw_input('Please press enter to continue when you are done.')
	line_chans_7m=raw_input('Please enter line channel range for 7m; e.g.., 2~10.-->')
	orde7=raw_input('Please enter order for 7m cont fit (recomended 0)-->')
	if not os.path.isdir(my_dir+'data_12m_'+line+'.ms.contsub'):
		print 'Making 12m continuum subtracted data set...'
		uvcontsub(vis = my_dir+'data_12m_'+line+'.ms',fitspw='0:'+line_chans_12m,excludechans = True,fitorder = int(orde12),solint='int')
	elif remakemsconts=='T':
		print 'Remaking 12m continuum subtracted data set...'
		os.system('rm -rf '+my_dir+'data_12m_'+line+'.ms.contsub')
		uvcontsub(vis = my_dir+'data_12m_'+line+'.ms',fitspw='0:'+line_chans_12m,excludechans = True,fitorder = int(orde12),solint='int')
	if not os.path.isdir(my_dir+'data_7m_'+line+'_comb.ms.contsub'):
		print 'Making 7m continuum subtracted data set...'
		uvcontsub(vis = my_dir+'data_7m_'+line+'_comb.ms',fitspw='0:'+line_chans_7m,excludechans = True,fitorder = int(orde7),solint='int')
	elif remakemsconts=='T':
		print 'Remaking 7m continuum subtracted data set...'
		os.system('rm -rf '+my_dir+'data_7m_'+line+'_cont.ms.contsub')
		uvcontsub(vis = my_dir+'data_7m_'+line+'_cont.ms',fitspw='0:'+line_chans_7m,excludechans = True,fitorder = int(orde7),solint='int')
	MS_12m=my_dir+'data_12m_'+line+'.ms.contsub'
	MS_7m=my_dir+'data_7m_'+line+'_comb.ms.contsub'
	print 'Plot continuum subtracted data...'
	print 'First 12m...'
	plotms(vis=MS_12m,yaxis='amp',xaxis='channel',spw='', \
		avgtime='1e8',avgscan=True,coloraxis='spw',plotfile=my_dir+'cs12m_'+line+'_chan.png',overwrite=True)
	raw_input('Please press enter to continue.')
	print 'Next 7m ...'
	plotms(vis=MS_7m,yaxis='amp',xaxis='channel',spw='', \
		avgtime='1e8',avgscan=True,coloraxis='spw',plotfile=my_dir+'cs12m_'+line+'_chan.png',overwrite=True)
	raw_input('Please press enter to continue.')
	print 'Okay. Moving on to imaging.'
	print 'The rest of the script is not interactive so please go do something else while'
	print 'it runs.'
	print '*******************************************************************************'
elif docontsub=='n':
	MS_12m=my_dir+'data_12m_'+line+'.ms'
	MS_7m=my_dir+'data_7m_'+line+'_comb.ms'
	print 'Okay. Moving on to imaging.'
	print 'The rest of the script is not interactive so please go do something else while'
	print 'it runs.'
	print '*******************************************************************************'
###################################

###################################
#Imaging Section
###################################

###template clean [(1)]
print 'Step 1: Model Image CLEAN'
os.system('rm -rf '+my_dir+'Template_'+line+'_fullchannel_nopb_all7m*')
print 'Making template 12m + 7m image...'
clean(vis=[MS_12m,MS_7m],\
	imagename=my_dir+'Template_'+line+'_fullchannel_nopb_all7m',\
	outlierfile="",field=['0~'+numf12, '0~'+numf7],spw=['0','0'], selectdata=True,timerange="",uvrange="",\
	antenna="",scan="",observation="",intent="", mode="channel",resmooth=False,gridmode="",\
	wprojplanes=-1,facets=1,cfcache="cfcache.dir", rotpainc=5.0,painc=360.0,aterm=True,\
	psterm=False,mterm=True,wbawp=False,conjbeams=True, epjtable="",interpolation="linear",\
	niter=myniter_temp,gain=0.1,threshold=mythreshold,psfmode="clark",imagermode="csclean", \
	ftmachine="mosaic",mosweight=False,scaletype="SAULT",multiscale=[], negcomponent=-1,\
	smallscalebias=0.6, interactive=False,mask=[], nchan=mynchan_temp,start=mystart_temp,width=mywidth_temp,\
	outframe="lsrk",veltype="radio", imsize=[myimsize, myimsize],cell=mycell, \
	phasecenter="J2000 19h15m38.305s 10d41m01.018s", restfreq=myrestfreq, stokes="I",\
	weighting="natural",robust=0.0,uvtaper=False,outertaper=[''],innertaper=['1.0'],modelimage='',\
	restoringbeam=[''], pbcor=False,minpb=0.2,usescratch=False,noise="1.0Jy", npixels=0,\
	npercycle=100,cyclefactor=1.5,cyclespeedup=-1,nterms=1,reffreq="",chaniter=False, \
	flatnoise=True,allowchunk=False)

#regrid the sd image to the template image of the line
os.system('rm -rf '+my_dir+'Radio_Peak_of_IRAS_19132+1035.spw'+spw_sd+'.I.sd.image')
os.system('rm -rf '+my_dir+'ReGrid'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image')
print 'Regriding SD image to 12m + 7m template...'
importfits(fitsimage=data_TP+'.spw'+spw_sd+'.I.sd.fits', \
	imagename=my_dir+'Radio_Peak_of_IRAS_19132+1035.spw'+spw_sd+'.I.sd.image')
imregrid(imagename=my_dir+'Radio_Peak_of_IRAS_19132+1035.spw'+spw_sd+'.I.sd.image',\
 template=my_dir+'Template_'+line+'_fullchannel_nopb_all7m.image', \
 output=my_dir+'ReGrid'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image', asvelocity=False,overwrite=True)
imhead(my_dir+'ReGrid'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image')

#either method works the same, was just there for testing!
if whichmethod=='C':
#Make sure TP image has same non-uniform PB response (lower in outskirts) as 12_7m and fix stokes axis/headers
	## (a) CASA task method
	#cut out clean masked stuff/noisy edges to make things same size for immath, esentially use TP regridded image size
	print 'Fixing Stokes axis and garbled headers...'
	os.system('rm -rf '+my_dir+'Template_'+line+'_fullchannel_nopb_all7m.subim')
	os.system('rm -rf '+my_dir+'ReGrid'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.subim')
	os.system('rm -rf '+my_dir+'Template_'+line+'_fullchannel_nopb_all7m.flux.subim')
	os.system('rm -rf '+my_dir+'ReGrid'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.subim.trans')
	os.system('rm -rf '+my_dir+'ReGrid'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.subim.trans.debpb')
	imsubimage(imagename=my_dir+'Template_'+line+'_fullchannel_nopb_all7m.image',\
		outfile=my_dir+'Template_'+line+'_fullchannel_nopb_all7m.subim',\
		region='box [[19:15:42.85546, +010.40.10.6834], [19:15:34.57909, +010.42.06.5645]]')
	imsubimage(imagename=my_dir+'ReGrid'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image',\
		outfile=my_dir+'ReGrid'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.subim',\
		region='box [[19:15:42.85546, +010.40.10.6834], [19:15:34.57909, +010.42.06.5645]]')
	#make cutout of PB flux image too
	imsubimage(imagename=my_dir+'Template_'+line+'_fullchannel_nopb_all7m.flux',\
		outfile=my_dir+'Template_'+line+'_fullchannel_nopb_all7m.flux.subim',\
		region='box [[19:15:42.85546, +010.40.10.6834], [19:15:34.57909, +010.42.06.5645]]')
	#fix stokes axis
	imtrans(imagename=my_dir+'ReGrid'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.subim',\
		outfile=my_dir+'ReGrid'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.subim.trans', order='0132')
	#fix headers: only needed in CASA 4.7 and above!!!
	#tt=imhead(data_TP+'.spw'+spw_sd+'.I.sd.fits')
	#U=tt['unit']
	#RMAJ=str(tt['restoringbeam']['major']['value'])+'arcsec'
	#RMIN=str(tt['restoringbeam']['minor']['value'])+'arcsec'
	#RPA=str(tt['restoringbeam']['positionangle']['value'])+'deg'
	#imhead(my_dir+'ReGrid'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.subim.trans',mode='put',hdkey='bunit',hdvalue=U)
	#imhead(my_dir+'ReGrid'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.subim.trans',mode='add',hdkey='bmaj',hdvalue=RMAJ)
	#imhead(my_dir+'ReGrid'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.subim.trans',mode='put',hdkey='bmin',hdvalue=RMIN)
	#imhead(my_dir+'ReGrid'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.subim.trans',mode='put',hdkey='bpa',hdvalue=RPA)
	#multiply by PB
	print 'Make sure TP image has same non-uniform PB response (lower in outskirts) as 12m + 7m...'
	immath(imagename=[my_dir+'ReGrid'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.subim.trans',\
		my_dir+'Template_'+line+'_fullchannel_nopb_all7m.flux.subim'],expr='IM0*IM1',\
		outfile=my_dir+'ReGrid'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.subim.trans.debpb')
	modelimg=my_dir+'ReGrid'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.subim.trans.debpb'
elif whichmethod=='S':
	## (b) Spectral cube method
	#export the template .flux file and the regridded sd file to fits
	os.system('rm -rf '+my_dir+'Template_'+line+'_fullchannelflux_nopb_all7m.fits')
	os.system('rm -rf '+my_dir+'ReGrid'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.fits')
	os.system('rm -rf '+my_dir+'ReGrid_fluxmultregrid0.spw'+spw_sd+'.I.sd1.fits')
	os.system('rm -rf '+my_dir+'ReGrid_fluxmultregrid0.spw'+spw_sd+'.I.sd1.image')
	os.system('rm -rf '+my_dir+'ReGrid_fluxmultregrid1.spw'+spw_sd+'.I.sd1.image')
	exportfits(imagename=my_dir+'Template_'+line+'_fullchannel_nopb_all7m.flux',\
		fitsimage=my_dir+'Template_'+line+'_fullchannelflux_nopb_all7m.fits', overwrite=True)
	exportfits(imagename=my_dir+'ReGrid'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image',\
		fitsimage=my_dir+'ReGrid'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.fits', overwrite=True)
	#with spectralcube multiply the flux and sd
	print 'Make sure TP image has same non-uniform PB response (lower in outskirts) as 12m + 7m...'
	flux=SpectralCube.read(my_dir+'Template_'+line+'_fullchannelflux_nopb_all7m.fits')
	flux.allow_huge_operations=True
	sd=SpectralCube.read(my_dir+'ReGrid'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.fits')
	sd.allow_huge_operations=True
	end=flux*sd
	#need to add a 4th axis - Stokes and modify the headers to match the ones of the original sd file
	print 'Fixing Stokes axis and garbled headers...'
	hdu = fits.PrimaryHDU(end.filled_data[:].value.reshape((1,)+end.shape),\
		header=add_stokes_axis_to_wcs(end.wcs, 3).to_header())
	#fix garbled headers part 1
	a=SpectralCube.read(data_TP+'.spw'+spw_sd+'.I.sd.fits')
	hdu.header['BUNIT']='Jy/beam'
	hdu.header['BMAJ']=a.header['BMAJ']
	hdu.header['BMIN']=a.header['BMIN']
	hdu.header['BPA']=a.header['BPA']
	hdu.header['BZERO']=a.header['BZERO']
	hdu.header['BSCALE']=a.header['BSCALE']
	hdu.header['BTYPE']=a.header['BTYPE']
	hdu.header['TELESCOP']=a.header['TELESCOP']
	hdu.header['INSTRUME']=a.header['INSTRUME']
	hdu.header['OBSERVER']=a.header['OBSERVER']
	hdu.header['TIMESYS']=a.header['TIMESYS']
	hdu.header['OBSDEC']=a.header['OBSDEC']
	hdu.header['OBSRA']=a.header['OBSRA']
	hdu.header['VELREF']=a.header['VELREF']
	hdu.writeto(my_dir+'ReGrid_fluxmultregrid0.spw'+spw_sd+'.I.sd1.fits', clobber=True)
	importfits(fitsimage=my_dir+'ReGrid_fluxmultregrid0.spw'+spw_sd+'.I.sd1.fits',\
		imagename=my_dir+'ReGrid_fluxmultregrid0.spw'+spw_sd+'.I.sd1.image')
	#get back to a good casa image by imtrans 
	imtrans(imagename=my_dir+'ReGrid_fluxmultregrid0.spw'+spw_sd+'.I.sd1.image',\
		outfile=my_dir+'ReGrid_fluxmultregrid1.spw'+spw_sd+'.I.sd1.image', order='0132')
	#fix garbled headers part 2, only needed in CASA 4.7 and above!!!
	#tt=imhead(data_TP+'.spw'+spw_sd+'.I.sd.fits')
	#U=tt['unit']
	#RMAJ=str(tt['restoringbeam']['major']['value'])+'arcsec'
	#RMIN=str(tt['restoringbeam']['minor']['value'])+'arcsec'
	#RPA=str(tt['restoringbeam']['positionangle']['value'])+'deg'
	#imhead(my_dir+'ReGrid_fluxmultregrid1.spw'+spw_sd+'.I.sd1.image',mode='put',hdkey='bunit',hdvalue=U)
	#imhead(my_dir+'ReGrid_fluxmultregrid1.spw'+spw_sd+'.I.sd1.image',mode='add',hdkey='bmaj',hdvalue=RMAJ)
	#imhead(my_dir+'ReGrid_fluxmultregrid1.spw'+spw_sd+'.I.sd1.image',mode='put',hdkey='bmin',hdvalue=RMIN)
	#imhead(my_dir+'ReGrid_fluxmultregrid1.spw'+spw_sd+'.I.sd1.image',mode='put',hdkey='bpa',hdvalue=RPA)
	modelimg=my_dir+'ReGrid_fluxmultregrid1.spw'+spw_sd+'.I.sd1.image'
#modelimage CLEAN
print 'Making 12m + 7m + TP image using regridded TP as model image...'
os.system('rm -rf '+my_dir+'CLEAN'+line+"_combinewithtp*")

clean(vis=[MS_12m,MS_7m],\
	imagename=my_dir+'CLEAN'+line+"_combinewithtp",\
	outlierfile="",field=['0~'+numf12, '0~'+numf7],spw=['0','0'], selectdata=True,timerange="",uvrange="",\
	antenna="",scan="",observation="",intent="", mode="channel",resmooth=False,gridmode="",\
	wprojplanes=-1,facets=1,cfcache="cfcache.dir", rotpainc=5.0,painc=360.0,aterm=True,\
	psterm=False,mterm=True,wbawp=False,conjbeams=True, epjtable="",interpolation="linear",\
	niter=myniter,gain=0.1,threshold=mythreshold,psfmode="clark",imagermode="csclean", \
	ftmachine="mosaic",mosweight=False,scaletype="SAULT",multiscale=[0, 1, 5, 10, 15], negcomponent=-1,\
	smallscalebias=0.6, interactive=False,mask=[], nchan=mynchan,start=mystart,width=mywidth,\
	outframe="lsrk",veltype="radio", imsize=[myimsize, myimsize],cell=mycell, \
	phasecenter="J2000 19h15m38.305s 10d41m01.018s", restfreq=myrestfreq, stokes="I",\
	weighting="natural",robust=0.0,uvtaper=False,outertaper=[''],innertaper=['1.0'],modelimage=modelimg,\
	restoringbeam=[''], pbcor=False,minpb=0.2,usescratch=False,noise="1.0Jy", npixels=0,\
	npercycle=100,cyclefactor=1.5,cyclespeedup=-1,nterms=1,reffreq="",chaniter=False, \
	flatnoise=True,allowchunk=False)
#PB correction
print "Doing PB correction..."
os.system('rm -rf '+my_dir+'GRS1915_modelimg_'+line+'.image.pbcor')
os.system('rm -rf '+my_dir+'GRS1915_modelimg_'+line+'.image.pbcor.fits')
immath(imagename=[my_dir+'CLEAN'+line+"_combinewithtp.image",
                  my_dir+'CLEAN'+line+"_combinewithtp.flux"],
       expr='IM0/IM1',
       outfile=my_dir+'GRS1915_modelimg_'+line+'.image.pbcor')
#make fits for radex fits and Tmax maps
exportfits(imagename=my_dir+'GRS1915_modelimg_'+line+'.image.pbcor',\
		fitsimage=my_dir+'GRS1915_modelimg_'+line+'.image.pbcor.fits', overwrite=True)

###Kauffman method [(2)]
print 'Step 2: Performing Kauffman Method...'
print 'Getting positive interferometer-only clean components...'
# get the positive interferometer-only clean components
os.system('rm -rf '+my_dir+'CLEAN'+line+"_combinewithtp.model.subim")
os.system('rm -rf '+my_dir+'CLEAN'+line+"_combinewithtp.flux.subim")
os.system('rm -rf '+my_dir+'ReGrid'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.subim.trans.debpb.chancut')
os.system('rm -rf '+my_dir+'ReGrid'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.subim.trans.debpb.chancut.pix')
imsubimage(imagename=my_dir+'CLEAN'+line+"_combinewithtp.model",
           outfile=my_dir+'CLEAN'+line+"_combinewithtp.model.subim",
           region='box [[19:15:42.85546, +010.40.10.6834], [19:15:34.57909, +010.42.06.5645]]')
imsubimage(imagename=my_dir+'CLEAN'+line+"_combinewithtp.flux",
           outfile=my_dir+'CLEAN'+line+"_combinewithtp.flux.subim",
           region='box [[19:15:42.85546, +010.40.10.6834], [19:15:34.57909, +010.42.06.5645]]')
imsubimage(imagename=my_dir+'ReGrid'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.subim.trans.debpb',\
 outfile=my_dir+'ReGrid'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.subim.trans.debpb.chancut',chans=chans_select)
rest_ghz=float(myrestfreq.replace('GHz',''))
sd_res=74.88/(rest_ghz/100.)/(12./10.)
cell_arc=float(mycell.replace('arcsec',''))
Conv=(cell_arc**2)/(1.133*sd_res**2)
#oldfile=my_dir+'ReGrid'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.subim.trans.debpb.chancut'
#im1=ia.imagecalc(outfile=my_dir+'ReGrid'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.subim.trans.debpb.chancut.pix',\
 #pixels=oldfile+'*'+str(Conv),overwrite=True)
#im1.done()
#ia.close()
immath(imagename=[my_dir+'ReGrid'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.subim.trans.debpb.chancut'],
       expr='IM0*'+str(Conv),
       outfile=my_dir+'ReGrid'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.subim.trans.debpb.chancut.pix')

imhead(my_dir+'ReGrid'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.subim.trans.debpb.chancut.pix',\
	mode='put',hdkey='bunit',hdvalue='Jy/pixel')
os.system('rm -rf '+my_dir+'deconvolved-sdinput.intmodel')
os.system('rm -rf '+my_dir+'deconvolved-sdinput-masked.intmodel')
immath(outfile=my_dir+'deconvolved-sdinput.intmodel',
   imagename=[my_dir+'CLEAN'+line+"_combinewithtp.model.subim",
              my_dir+'ReGrid'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.subim.trans.debpb.chancut.pix'],
   mode='evalexpr',
   expr='iif((IM0-IM1) >= 0.00, IM0-IM1, 0.0)')
# remove those components if they are at the map edge
print 'Remove components if they are at the map edge...'
immath(outfile=my_dir+'deconvolved-sdinput-masked.intmodel',
   imagename=[my_dir+'deconvolved-sdinput.intmodel',
               my_dir+'CLEAN'+line+"_combinewithtp.flux.subim"],
   mode='evalexpr',
   expr='iif((IM1) >= 0.25, IM0, 0.0)')
imhead(my_dir+'deconvolved-sdinput-masked.intmodel',
       mode='put',
       hdkey='bunit', hdvalue='Jy/pixel')
# smooth the interferometer-only components to the synthesized beam
print 'Smooth the interferometer-only components to the synthesized beam...'
test_imhead=imhead(my_dir+'CLEAN'+line+"_combinewithtp.image")
numchans=test_imhead['shape'][3]
#if multiple beams (per channel) they dont change much so you middle channel one
if 'perplanebeams' in test_imhead:
	SynthBeamMaj=test_imhead['perplanebeams']['beams']['*'+str(numchans/2)]['*0']['major']['value']
	SynthBeamMin=test_imhead['perplanebeams']['beams']['*'+str(numchans/2)]['*0']['minor']['value']
	SynthBeamPA=test_imhead['perplanebeams']['beams']['*'+str(numchans/2)]['*0']['positionangle']['value']
else:
	SynthBeamMaj = imhead(my_dir+'CLEAN'+line+"_combinewithtp.image", mode='get', hdkey='beammajor')['value']
	SynthBeamMin = imhead(my_dir+'CLEAN'+line+"_combinewithtp.image", mode='get', hdkey='bmin')['value']
	SynthBeamPA = imhead(my_dir+'CLEAN'+line+"_combinewithtp.image", mode='get', hdkey='bpa')['value']
os.system('rm -rf '+my_dir+'deconvolved-sdinput.intimage')
imsmooth(imagename=my_dir+'deconvolved-sdinput-masked.intmodel',
         outfile=my_dir+'deconvolved-sdinput.intimage',
         kernel='gauss',
         major=str(SynthBeamMaj)+'arcsec',
         minor=str(SynthBeamMin)+'arcsec',
         pa=str(SynthBeamPA)+'deg')

# feather the data
print 'Feather smoothed interferometer component image with regridded SD image...'
os.system('rm -rf '+my_dir+'deconvolved-combi.model')
feather(imagename = my_dir+'deconvolved-combi.model',
        highres = my_dir+'deconvolved-sdinput.intimage',
        lowres = modelimg)
# clean up and keep the final model
os.system('rm -rf '+my_dir+'combimodel.image')
os.system('mv '+my_dir+'deconvolved-combi.model '+my_dir+'combimodel.image')
#CLEAN round 2
print 'CLEANing with feathered image as modelimage...'
os.system('rm -rf '+my_dir+'CLEAN'+line+"_combinewithtpnew*")

clean(vis=[MS_12m,MS_7m],\
	imagename=my_dir+'CLEAN'+line+"_combinewithtpnew",\
	outlierfile="",field=['0~'+numf12, '0~'+numf7],spw=['0','0'], selectdata=True,timerange="",uvrange="",\
	antenna="",scan="",observation="",intent="", mode="channel",resmooth=False,gridmode="",\
	wprojplanes=-1,facets=1,cfcache="cfcache.dir", rotpainc=5.0,painc=360.0,aterm=True,\
	psterm=False,mterm=True,wbawp=False,conjbeams=True, epjtable="",interpolation="linear",\
	niter=myniter,gain=0.1,threshold=mythreshold,psfmode="clark",imagermode="csclean", \
	ftmachine="mosaic",mosweight=False,scaletype="SAULT",multiscale=[0, 1, 5, 10, 15], negcomponent=-1,\
	smallscalebias=0.6, interactive=False,mask=[], nchan=mynchan,start=mystart,width=mywidth,\
	outframe="lsrk",veltype="radio", imsize=[myimsize, myimsize],cell=mycell, \
	phasecenter="J2000 19h15m38.305s 10d41m01.018s", restfreq=myrestfreq, stokes="I",\
	weighting="natural",robust=0.0,uvtaper=False,outertaper=[''],innertaper=['1.0'],modelimage=my_dir+'combimodel.image',\
	restoringbeam=[''], pbcor=False,minpb=0.2,usescratch=False,noise="1.0Jy", npixels=0,\
	npercycle=100,cyclefactor=1.5,cyclespeedup=-1,nterms=1,reffreq="",chaniter=False, \
	flatnoise=True,allowchunk=False)
os.system('rm -rf '+my_dir+'GRS1915_modelimgnew_'+line+'.image.pbcor')
os.system('rm -rf '+my_dir+'GRS1915_modelimgnew_'+line+'.image.pbcor.fits')
print 'Doing PB correction...'
immath(imagename=[my_dir+'CLEAN'+line+"_combinewithtpnew.image",
                  my_dir+'CLEAN'+line+"_combinewithtpnew.flux"],
       expr='IM0/IM1',
       outfile=my_dir+'GRS1915_modelimgnew_'+line+'.image.pbcor')
#make fits for radex fits and Tmax maps
exportfits(imagename=my_dir+'GRS1915_modelimgnew_'+line+'.image.pbcor',\
		fitsimage=my_dir+'GRS1915_modelimgnew_'+line+'.image.pbcor.fits', overwrite=True)

###Feathering method [(3) and (4)]
print 'Step 3: Performing Feathering Method...'
#first clean 12m+7m
print 'Making 12m + 7m image...'
os.system('rm -rf '+my_dir+'CLEAN'+line+"_combine12m7m_nopb_all7*")

clean(vis=[MS_12m,MS_7m],\
	imagename=my_dir+'CLEAN'+line+"_combine12m7m_nopb_all7",\
	outlierfile="",field=['0~'+numf12, '0~'+numf7],spw=['0','0'], selectdata=True,timerange="",uvrange="",\
	antenna="",scan="",observation="",intent="", mode="channel",resmooth=False,gridmode="",\
	wprojplanes=-1,facets=1,cfcache="cfcache.dir", rotpainc=5.0,painc=360.0,aterm=True,\
	psterm=False,mterm=True,wbawp=False,conjbeams=True, epjtable="",interpolation="linear",\
	niter=myniter,gain=0.1,threshold=mythreshold,psfmode="clark",imagermode="csclean", \
	ftmachine="mosaic",mosweight=False,scaletype="SAULT",multiscale=[0, 1, 5, 10, 15], negcomponent=-1,\
	smallscalebias=0.6, interactive=False,mask=[], nchan=mynchan,start=mystart,width=mywidth,\
	outframe="lsrk",veltype="radio", imsize=[myimsize, myimsize],cell=mycell, \
	phasecenter="J2000 19h15m38.305s 10d41m01.018s", restfreq=myrestfreq, stokes="I",\
	weighting="natural",robust=0.0,uvtaper=False,outertaper=[''],innertaper=['1.0'],modelimage='',\
	restoringbeam=[''], pbcor=False,minpb=0.2,usescratch=False,noise="1.0Jy", npixels=0,\
	npercycle=100,cyclefactor=1.5,cyclespeedup=-1,nterms=1,reffreq="",chaniter=False, \
	flatnoise=True,allowchunk=False)
os.system('rm -rf '+my_dir+'CLEAN'+line+'_combine12m7m_nopb_all7.image.subim')
os.system('rm -rf '+my_dir+'CLEAN'+line+'_combine12m7m_nopb_all7.flux.subim')
os.system('rm -rf '+my_dir+'GRS1915_12m7m_'+line+'.image.pbcor')
os.system('rm -rf '+my_dir+'GRS1915_12m7m_'+line+'.image.pbcor.fits')
print 'Doing PB correction on 12m + 7m image...'
imsubimage(imagename=my_dir+'CLEAN'+line+'_combine12m7m_nopb_all7.image',
           outfile=my_dir+'CLEAN'+line+'_combine12m7m_nopb_all7.image.subim',
           region='box [[19:15:42.85546, +010.40.10.6834], [19:15:34.57909, +010.42.06.5645]]')
imsubimage(imagename=my_dir+'CLEAN'+line+'_combine12m7m_nopb_all7.flux',
           outfile=my_dir+'CLEAN'+line+'_combine12m7m_nopb_all7.flux.subim',
           region='box [[19:15:42.85546, +010.40.10.6834], [19:15:34.57909, +010.42.06.5645]]')
immath(imagename=[my_dir+'CLEAN'+line+'_combine12m7m_nopb_all7.image.subim',
                 my_dir+'CLEAN'+line+'_combine12m7m_nopb_all7.flux.subim'],
       expr='IM0/IM1',
       outfile=my_dir+'GRS1915_12m7m_'+line+'.image.pbcor')
#make fits for radex fits and Tmax maps
exportfits(imagename=my_dir+'GRS1915_12m7m_'+line+'.image.pbcor',\
		fitsimage=my_dir+'GRS1915_12m7m_'+line+'.image.pbcor.fits', overwrite=True)


os.system('rm -rf '+my_dir+'GRS1915_Featherall_'+line+'.image')
os.system('rm -rf '+my_dir+'GRS1915_Featherall7_'+line+'.image.pbcor')
os.system('rm -rf '+my_dir+'GRS1915_Featherall7_'+line+'.image.pbcor.fits')
print 'Feathering 12m + 7m image and regridded TP image...'
feather(imagename=my_dir+'GRS1915_Featherall_'+line+'.image',
        highres=my_dir+'CLEAN'+line+"_combine12m7m_nopb_all7.image.subim",
        lowres=modelimg)
print 'Doing PB correction...'
immath(imagename=[my_dir+'GRS1915_Featherall_'+line+'.image',
                 my_dir+'CLEAN'+line+'_combine12m7m_nopb_all7.flux.subim'],
       expr='IM0/IM1',
       outfile=my_dir+'GRS1915_Featherall7_'+line+'.image.pbcor')
#make fits for radex fits and Tmax maps
exportfits(imagename=my_dir+'GRS1915_Featherall7_'+line+'.image.pbcor',\
		fitsimage=my_dir+'GRS1915_Featherall7_'+line+'.image.pbcor.fits', overwrite=True)

print 'Cleaning up...'
os.system('rm -rf *.log')
os.system('rm -rf *.last')

print '********************************************************************'
print 'The script is finished. Please inspect the resulting data products.'
print '********************************************************************'
###################################
