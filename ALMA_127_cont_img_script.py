#################################################
#ALMA Continuum Imaging Script
#################################################
'''CASA script to be used for imaging 12m and 7m ALMA Continuum SPWs
INPUT: Parameter file detailing all data and imaging parameters (param_dir_file set below)
OUTPUT: (1) 12m continuum image -- GRS1915_12m_cont.image.pbcor
		(2) 7m continuum image -- GRS1915_7m_cont.image.pbcor
NOTES: - All output images & intermediate data products are put in my_dir directory set below.
	   - All output images are also converted to fits format (just append .fits to end of images 1,2 above)
	   - This script flags any lines in the continuum spw, then does an mfs clean with nterms
		 for the wide bandwidth.
Written by: Alex J. Tetarenko
Last Updated: Dec 15 2016'''

print '##################################################'
print 'Welcome to Alexs ALMA Continuum Imaging Script'
print '##################################################\n'


#packages you may need
from astropy.io import fits
from astropy.wcs.utils import add_stokes_axis_to_wcs
import re
import sys
import imp
import os
import linecache
import find
import warnings
warnings.filterwarnings('ignore')

#define output directory
my_dir='/mnt/bigdata/tetarenk/ALMA_GRS1915_105/alex_imaging_contspw_fix_redo/'
if not os.path.isdir(my_dir):
	os.system('sudo mkdir '+my_dir)
	os.system('sudo chown ubuntu '+my_dir)
	os.system('sudo chmod -R u+r '+my_dir) 
	os.system('sudo chmod -R u+w '+my_dir)
	os.system('sudo chmod -R u+x '+my_dir)
print 'You have set your output directory to ', my_dir
print 'All output images & intermediate data products are put in this directory.\n'

#param file location
param_dir_file='/mnt/bigdata/tetarenk/ALMA_GRS1915_105/alexs_scripts/params_cont.txt'
print 'You have set your param file to ', param_dir_file
print 'Please make sure all parameters are correct, they will change for each data set!\n'

#################################################
#Defining Parameters Section
#################################################
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
field_12m=data_params.field_12m
field_7m=data_params.field_7m
numf12=str(int(field_12m.split('~')[1])-int(field_12m.split('~')[0]))
numf7=str(int(field_7m.split('~')[1])-int(field_7m.split('~')[0]))
remakems=data_params.remakems
#cont spw params
spw_arr=data_params.spw_arr
spw_7m=data_params.spw_7m
#general image params
mythreshold12=data_params.mythreshold12
myimsize12=data_params.myimsize12
mycell12=data_params.mycell12
myniter12=data_params.myniter12
mythreshold7=data_params.mythreshold7
myimsize7=data_params.myimsize7
mycell7=data_params.mycell7
myniter7=data_params.myniter7
mythreshold127=data_params.mythreshold127
myimsize127=data_params.myimsize127
mycell127=data_params.mycell127
myniter127=data_params.myniter127
mynterms=data_params.mynterms
mymask12=data_params.mymask12
mymask7=data_params.mymask7
mymask127=data_params.mymask127
line='cont'
#################################################

#################################################
#Examining and Splitting Data Section
#################################################
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
	os.system('rm -rf '+my_dir+'data_12m_'+line+'.ms.flagversions')
	split(vis=data_12m,outputvis=my_dir+'data_12m_'+line+'.ms',spw=spw_arr,datacolumn='data',field=field_12m)
else:
	print 'Split 12m MS already exists'
if not os.path.isdir(my_dir+'data_7m_'+line+'.ms'):
	print 'Making split 7m MS...'
	split(vis=data_7m,outputvis=my_dir+'data_7m_'+line+'.ms',spw=spw_7m,datacolumn='data',field=field_7m)
elif remakems=='T':
	print 'Remaking 7m MS...'
	os.system('rm -rf '+my_dir+'data_7m_'+line+'.ms')
	os.system('rm -rf '+my_dir+'data_7m_'+line+'.ms.flagversions')
	split(vis=data_7m,outputvis=my_dir+'data_7m_'+line+'.ms',spw=spw_7m,datacolumn='data',field=field_7m)
else:
	print 'Split 7m MS already exists'
if not os.path.isdir(my_dir+'data_7m_'+line+'_comb.ms'):
	print 'Combining 7m spws...'
	mstransform(vis=my_dir+'data_7m_'+line+'.ms',outputvis=my_dir+'data_7m_'+line+'_comb.ms',combinespws=True,spw='',\
		datacolumn='data')
	listobs(my_dir+'data_7m_'+line+'_comb.ms',listfile=my_dir+'7m_listfile_'+line+'_comb.txt',overwrite=True)
elif remakems=='T':
	print 'Remaking 7m combined MS...'
	os.system('rm -rf '+my_dir+'data_7m_'+line+'_comb.ms')
	os.system('rm -rf '+my_dir+'data_7m_'+line+'_comb.ms.flagversions')
	mstransform(vis=my_dir+'data_7m_'+line+'.ms',outputvis=my_dir+'data_7m_'+line+'_comb.ms',combinespws=True,spw='',\
		datacolumn='data')
	listobs(my_dir+'data_7m_'+line+'_comb.ms',listfile=my_dir+'7m_listfile_'+line+'_comb.txt',overwrite=True)
else:
	print 'Combined Split 7m MS already exists'
#################################################

#################################################
#Flagging line emission in continuum spws section
#################################################
print 'Now flagging the lines in the continuum spw...'
flagorno=raw_input('Do you need to look at the plots (y) or do you know the channel range to flag (n)?-->')
if flagorno =='y':
	print 'First plot 12m spectra. Please use tool to identify channel range of line to flag.'
	plotms(vis=my_dir+'data_12m_'+line+'.ms',yaxis='amp',xaxis='channel',spw='', \
		avgtime='1e8',avgscan=True,coloraxis='spw',plotfile=my_dir+'A12m_'+line+'_chan.png',overwrite=True)
	raw_input('Please press enter to continue when you are done.')
	flag_chans_12m=raw_input('Please enter channel range to flag for 12m; e.g.., 2~10.-->')
	print 'Next plot 7m spectra. Please use tool to identify channel range of line to flag.'
	plotms(vis=my_dir+'data_7m_'+line+'_comb.ms',yaxis='amp',xaxis='channel',spw='', \
		avgtime='1e8',avgscan=True,coloraxis='spw',plotfile=my_dir+'A7m_'+line+'_chan.png',overwrite=True)
	raw_input('Please press enter to continue when you are done.')
	flag_chans_7m=raw_input('Please enter channel range to flag for 7m; e.g.., 2~10.-->')
	print 'Flagging selected line emission in both 12m and 7m data...'
elif flagorno =='n':
	flag_chans_12m=raw_input('Please enter channel range to flag for 12m; e.g.., 2~10.-->')
	flag_chans_7m=raw_input('Please enter channel range to flag for 7m; e.g.., 2~10.-->')
	print 'Flagging selected line emission in both 12m and 7m data...'
flagdata(vis=my_dir+'data_12m_'+line+'.ms', mode='manual', spw='0:'+flag_chans_12m)
flagdata(vis=my_dir+'data_7m_'+line+'_comb.ms', mode='manual', spw='0:'+flag_chans_7m)

seeagain=raw_input('Do you want to examine the flagged data plots? y or n-->')
if seeagain=='y':
	plotms(vis=my_dir+'data_12m_'+line+'.ms',yaxis='amp',xaxis='channel',spw='', \
	avgtime='1e8',avgscan=True,coloraxis='spw',plotfile=my_dir+'A12m_'+line+'_chan_flagged.png',overwrite=True)
	raw_input('Please press enter to continue when you are done.')
	plotms(vis=my_dir+'data_7m_'+line+'_comb.ms',yaxis='amp',xaxis='channel',spw='', \
	avgtime='1e8',avgscan=True,coloraxis='spw',plotfile=my_dir+'A7m_'+line+'_chan_flagged.png',overwrite=True)
	raw_input('Please press enter to continue when you are done.')
	print 'Continuing to MFS imaging...'
else:
	print 'Okay, Continuing to MFS imaging...'
#################################################

#################################################
#Imaging Section
#################################################
print 'Starting with 12m...'
if mymask12=='':
	print 'Using interactive mode so you can make a mask...'
	print 'Cleaning...'
	os.system('rm -rf '+my_dir+'GRS1915_cont12m*')
	clean(vis=my_dir+'data_12m_'+line+'.ms',imagename=my_dir+'GRS1915_cont12m',mode='mfs',\
	imagermode='csclean',imsize=[myimsize12,myimsize12],cell=mycell12,spw='',gain=0.1,\
	weighting='natural',nterms=mynterms, mask=mymask12,usescratch=False\
	,interactive=T,threshold=mythreshold12,niter=0,pbcor=False,minpb=0.2,\
	phasecenter="J2000 19h15m38.305s 10d41m01.018s")
	raw_input('Please press enter to continue when you are done.')
else:
	print 'Cleaning...'
	print 'This make take awhile, please go do something else for a while.'
	os.system('rm -rf '+my_dir+'GRS1915_cont12m*')
	clean(vis=my_dir+'data_12m_'+line+'.ms',imagename=my_dir+'GRS1915_cont12m',mode='mfs',\
	imagermode='csclean',imsize=[myimsize12,myimsize12],cell=mycell12,spw='',gain=0.1,\
	weighting='natural',nterms=mynterms, mask=mymask12,usescratch=False\
	,interactive=F,threshold=mythreshold12,niter=myniter12,pbcor=False,minpb=0.2,\
	phasecenter="J2000 19h15m38.305s 10d41m01.018s")
	raw_input('Please press enter to continue when you are done.')
os.system('rm -rf '+my_dir+'GRS1915_12m_cont.image')
os.system('rm -rf '+my_dir+'GRS1915_12m_cont.image.pbcor')
os.system('rm -rf '+my_dir+'GRS1915_12m_cont.image.pbcor.fits')
os.system('sudo mv '+my_dir+'GRS1915_cont12m.image.tt0 '+my_dir+'GRS1915_12m_cont.image')
immath(imagename=[my_dir+'GRS1915_12m_cont.image',
                 my_dir+'GRS1915_cont12m.flux'],
       expr='IM0/IM1',
       outfile=my_dir+'GRS1915_12m_cont.image.pbcor')
exportfits(imagename=my_dir+'GRS1915_12m_cont.image.pbcor',fitsimage=my_dir+'GRS1915_12m_cont.image.pbcor.fits')

print 'Next the 7m...'
if mymask7=='':
	print 'Using interactive mode so you can make a mask...'
	print 'Cleaning...'
	os.system('rm -rf '+my_dir+'GRS1915_cont7m*')
	clean(vis=my_dir+'data_7m_'+line+'_comb.ms',imagename=my_dir+'GRS1915_cont7m',mode='mfs',\
	imagermode='csclean',imsize=[myimsize7,myimsize7],cell=mycell7,spw='',gain=0.1,\
	weighting='natural',nterms=mynterms, mask=mymask7,usescratch=False\
	,interactive=T,threshold=mythreshold7,niter=0,pbcor=False,minpb=0.2,\
	phasecenter="J2000 19h15m38.305s 10d41m01.018s")
else:
	print 'Cleaning...'
	print 'This make take awhile, please go do something else for a while.'
	os.system('rm -rf '+my_dir+'GRS1915_cont7m*')
	clean(vis=my_dir+'data_7m_'+line+'_comb.ms',imagename=my_dir+'GRS1915_cont7m',mode='mfs',\
	imagermode='csclean',imsize=[myimsize7,myimsize7],cell=mycell7,spw='',gain=0.1,\
	weighting='natural',nterms=mynterms, mask=mymask7,usescratch=False\
	,interactive=F,threshold=mythreshold7,niter=myniter7,pbcor=False,minpb=0.2,\
	phasecenter="J2000 19h15m38.305s 10d41m01.018s")
os.system('rm -rf '+my_dir+'GRS1915_7m_cont.image')
os.system('rm -rf '+my_dir+'GRS1915_7m_cont.image.pbcor')
os.system('rm -rf '+my_dir+'GRS1915_7m_cont.image.pbcor.fits')
os.system('sudo mv '+my_dir+'GRS1915_cont7m.image.tt0 '+my_dir+'GRS1915_7m_cont.image')
immath(imagename=[my_dir+'GRS1915_7m_cont.image',
                 my_dir+'GRS1915_cont7m.flux'],
       expr='IM0/IM1',
       outfile=my_dir+'GRS1915_7m_cont.image.pbcor')
exportfits(imagename=my_dir+'GRS1915_7m_cont.image.pbcor',fitsimage=my_dir+'GRS1915_7m_cont.image.pbcor.fits')

print 'Next the 12m + 7m...'
if mymask12=='':
	print 'Using interactive mode so you can make a mask...'
	print 'Cleaning...'
	os.system('rm -rf '+my_dir+'GRS1915_cont12m7m*')
	clean(vis=[my_dir+'data_12m_'+line+'.ms',my_dir+'data_7m_'+line+'_comb.ms'],\
	imagename=my_dir+'GRS1915_cont12m7m',mode='mfs',psfmode="clark",ftmachine="mosaic",\
	imagermode='csclean',imsize=[myimsize127,myimsize127],cell=mycell127,gain=0.1,\
	weighting='natural',nterms=mynterms,usescratch=False\
	,interactive=T,threshold=mythreshold127,niter=0,pbcor=False,minpb=0.2,\
	phasecenter="J2000 19h15m38.305s 10d41m01.018s",field=['0~'+numf12, '0~'+numf7],spw=['0','0'])
	raw_input('Please press enter to continue when you are done.')

else:
	print 'Cleaning...'
	print 'This make take awhile, please go do something else for a while.'
	os.system('rm -rf '+my_dir+'GRS1915_cont12m7m*')
	clean(vis=[my_dir+'data_12m_'+line+'.ms',my_dir+'data_7m_'+line+'_comb.ms'],\
	imagename=my_dir+'GRS1915_cont12m7m',mode='mfs',psfmode="clark",ftmachine="mosaic",\
	imagermode='csclean',imsize=[myimsize127,myimsize127],cell=mycell127,gain=0.1,\
	weighting='natural',nterms=mynterms, mask=mymask127,usescratch=False\
	,interactive=F,threshold=mythreshold127,niter=myniter127,pbcor=False,minpb=0.2,\
	phasecenter="J2000 19h15m38.305s 10d41m01.018s",field=['0~'+numf12, '0~'+numf7],spw=['0','0'])
	raw_input('Please press enter to continue when you are done.')
os.system('rm -rf '+my_dir+'GRS1915_12m7m_cont.image')
os.system('rm -rf '+my_dir+'GRS1915_12m7m_cont.image.pbcor')
os.system('rm -rf '+my_dir+'GRS1915_12m7m_cont.image.pbcor.fits')
os.system('sudo mv '+my_dir+'GRS1915_cont12m7m.image.tt0 '+my_dir+'GRS1915_12m7m_cont.image')
immath(imagename=[my_dir+'GRS1915_12m7m_cont.image',
                 my_dir+'GRS1915_cont12m7m.flux'],
       expr='IM0/IM1',
       outfile=my_dir+'GRS1915_12m7m_cont.image.pbcor')
exportfits(imagename=my_dir+'GRS1915_12m7m_cont.image.pbcor',fitsimage=my_dir+'GRS1915_12m7m_cont.image.pbcor.fits')


print 'Cleaning up...'
os.system('rm -rf *.log')
os.system('rm -rf *.last')

print '********************************************************************'
print 'The script is finished. Please inspect the resulting data products.'
print '********************************************************************'
#################################################
