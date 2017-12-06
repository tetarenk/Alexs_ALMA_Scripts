#CASA Re-run Pams imaging of ALMA data, shows process

from astropy.io import fits
from spectral_cube import SpectralCube
from astropy.wcs.utils import add_stokes_axis_to_wcs
import radio_beam


#####
#step 1- combine and image the 7m + 12m by creating all channel template image
#####
data_12m='/mnt/bigdata/tetarenk/ALMA_GRS1915_105/calib_MS/calibrated_final.ms'
data_7m='/mnt/bigdata/tetarenk/ALMA_GRS1915_105/calib_MS_7m/calibrated_final.ms'
data_TP='/mnt/bigdata/tetarenk/ALMA_GRS1915_105/2015.1.00976.S/science_goal.uid___A001_X2d6_X61/group.uid___A001_X2d6_X62/member.uid___A001_X2d6_X69/product/Radio_Peak_of_IRAS_19132+1035'
my_dir='/mnt/bigdata/tetarenk/ALMA_GRS1915_105/alex_redo/'
#listobs
listobs(data_12m,listfile='12m_listfile.txt')
listobs(data_7m,listfile='7m_listfile.txt')
os.system('pluma 12m_listfile.txt &')
os.system('pluma 7m_listfile.txt &')


##define variables
line='13CO'
spw_sd='33'#for single dish
spw_arr='8'#for array
spw_7m='8,26'

myniter=1
mythreshold='8mJy'
mynchan=-1
mystart=''
mywidth=1
myimsize=800
mycell='0.2arcsec'
myrestfreq='220.39868GHz'

split(vis=data_12m,outputvis='data_12m_CO13.ms',spw=spw_arr,datacolumn='data',field='4~36')
split(vis=data_7m,outputvis='data_7m_CO13.ms',spw=spw_7m,datacolumn='data',field='5~18')
mstransform(vis='data_7m_CO13.ms',outputvis='data_7m_CO13_comb.ms',combinespws=True,spw='',datacolumn='data')
listobs('data_7m_CO13_comb.ms')
plotms(vis='data_7m_CO13_comb.ms',yaxis='amp',xaxis='velocity',spw='', avgtime='1e8',avgscan=True,coloraxis='spw',avgchannel='5',
       transform=True,freqframe='LSRK',restfreq=myrestfreq, plotfile='All127_13CO_vel.png',overwrite=True)
#plotms(vis='CO13_12m_7m_concat.ms',yaxis='wt',xaxis='uvdist',spw='',
       #coloraxis='spw',plotfile='combine_CO_WT.png',overwrite=True)
#plotms(vis='CO13_12m_7m_concat.ms',yaxis='amp',xaxis='uvdist',spw='0~1:100~500', avgscan=True,
       #avgchannel='200', coloraxis='spw',plotfile='combine_uvdist.png',overwrite=True)


clean(vis=['data_12m_CO13.ms','data_7m_CO13_comb.ms'],imagename=line+'_fullchannel_nopb_all7m',outlierfile="",field=['0~32', '0~13'],spw=['0','0'], selectdata=True,timerange="",uvrange="",antenna="",scan="",observation="",intent="", mode="channel",resmooth=False,gridmode="",wprojplanes=-1,facets=1,cfcache="cfcache.dir", rotpainc=5.0,painc=360.0,aterm=True,psterm=False,mterm=True,wbawp=False,conjbeams=True, epjtable="",interpolation="linear", niter=myniter,gain=0.1,threshold=mythreshold,psfmode="clark",imagermode="csclean", ftmachine="mosaic",mosweight=False,scaletype="SAULT",multiscale=[], negcomponent=-1,smallscalebias=0.6, interactive=False,mask=[], nchan=mynchan,start=mystart,width=mywidth,outframe="lsrk",veltype="radio", imsize=[myimsize, myimsize],cell=mycell, phasecenter="J2000 19h15m38.305s 10d41m01.018s", restfreq=myrestfreq, stokes="I",weighting="natural",robust=0.0,uvtaper=False,outertaper=[''],innertaper=['1.0'],modelimage='', restoringbeam=[''], pbcor=False,minpb=0.2,usescratch=False,noise="1.0Jy", npixels=0,npercycle=100,cyclefactor=1.5,cyclespeedup=-1,nterms=1,reffreq="",chaniter=False, flatnoise=True,allowchunk=False)

#####
#step 2- regrid the sd image to the template image of the line
#Note: use the already made [800,800] 0.2 arcsec _fullchannel.* templates
#####
importfits(fitsimage=data_TP+'.spw'+spw_sd+'.I.sd.fits', imagename=my_dir+'Radio_Peak_of_IRAS_19132+1035.spw'+spw_sd+'.I.sd.image')
imregrid(imagename=my_dir+'Radio_Peak_of_IRAS_19132+1035.spw'+spw_sd+'.I.sd.image', template=line+'_fullchannel_nopb_all7m.image', output='B17m'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image', asvelocity=False,overwrite=True)
imhead('B17m'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image')
#convert to Jy/pixel
#rest_ghz=float(myrestfreq.replace('GHz',''))
#sd_res=74.88/(rest_ghz/100.)/(12./10.)
#cell_arc=float(mycell.replace('arcsec',''))
#Conv=(cell_arc**2)/(1.133*sd_res**2)

#im1=ia.imagecalc(outfile=line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd_scaled.image', pixels=str(Conv)+'*'+'B'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image',overwrite=True)
#im1.done()
#ia.close()
#imhead(line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd_scaled.image',mode='put',hdkey='bunit',hdvalue='Jy/pixel')
#alternate method
#cut out clean masked stuff
imsubimage(imagename=line+'_fullchannel_nopb_all7m.image',
           outfile=line+'_fullchannel_nopb_all7m.subim',
           region='box [[19:15:42.85546, +010.40.10.6834], [19:15:34.57909, +010.42.06.5645]]')
imsubimage(imagename='B17m'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image',
           outfile='B17m'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.subim',
           region='box [[19:15:42.85546, +010.40.10.6834], [19:15:34.57909, +010.42.06.5645]]')

#multiply by PB
imsubimage(imagename=line+'_fullchannel_nopb_all7m.flux',
           outfile='B17m'+line+'_fullchannel_nopb_all7m.flux.subim',
           region='box [[19:15:42.85546, +010.40.10.6834], [19:15:34.57909, +010.42.06.5645]]')
imtrans(imagename='B17m'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.subim', outfile='B17m'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.subim.trans', order='0132')
tt=imhead(data_TP+'.spw'+spw_sd+'.I.sd.fits')
U=tt['unit']
RMAJ=str(tt['restoringbeam']['major']['value'])+'arcsec'
RMIN=str(tt['restoringbeam']['minor']['value'])+'arcsec'
RPA=str(tt['restoringbeam']['positionangle']['value'])+'deg'

imhead('B17m'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.subim.trans',mode='put',hdkey='bunit',hdvalue=U)
imhead('B17m'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.subim.trans',mode='add',hdkey='bmaj',hdvalue=RMAJ)
imhead('B17m'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.subim.trans',mode='put',hdkey='bmin',hdvalue=RMIN)
imhead('B17m'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.subim.trans',mode='put',hdkey='bpa',hdvalue=RPA)

immath(imagename=['B17m'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.subim.trans',
                  'B17m'+line+'_fullchannel_nopb_all7m.flux.subim'],
       expr='IM0*IM1',
       outfile='B17m'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.subim.trans.debpb')
#tryjy.pix
rest_ghz=float(myrestfreq.replace('GHz',''))
sd_res=74.88/(rest_ghz/100.)/(12./10.)
cell_arc=float(mycell.replace('arcsec',''))
Conv=(cell_arc**2)/(1.133*sd_res**2)
im1=ia.imagecalc(outfile='B17m'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.subim.trans.debpb.pix', pixels=str(Conv)+'*B17m'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.subim.trans.debpb',overwrite=True)
im1.done()
ia.close()
imhead('B17m'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.subim.trans.debpb.pix',mode='put',hdkey='bunit',hdvalue='Jy/pixel')



#####
#step 3- make modelimage
#####

#export the template .flux file and the regridded sd file to fits
exportfits(imagename=line+'_fullchannel_nopb.flux', fitsimage=line+'_fullchannelflux_nopb.fits', overwrite=True)
exportfits(imagename='B1'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image', fitsimage=line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd1.fits', overwrite=True)

#with spectralcube multiply the flux and sd
flux=SpectralCube.read(line+'_fullchannelflux_nopb.fits')
flux.allow_huge_operations=True
sd=SpectralCube.read(line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd1.fits')
sd.allow_huge_operations=True
end=flux*sd

#need to add a 4th axis - Stokes and modify the headers to match the ones of the original sd file
hdu = fits.PrimaryHDU(end.filled_data[:].value.reshape((1,)+end.shape), header=add_stokes_axis_to_wcs(end.wcs, 3).to_header())

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

hdu.writeto(my_dir+'_fluxmultregrid0.spw'+spw_sd+'.I.sd1.fits', clobber=True)

importfits(fitsimage=my_dir+'_fluxmultregrid0.spw'+spw_sd+'.I.sd1.fits',imagename=my_dir+'_fluxmultregrid0.spw'+spw_sd+'.I.sd1.image')

#get back to a good casa image by imtrans 
imtrans(imagename=my_dir+'_fluxmultregrid0.spw'+spw_sd+'.I.sd1.image', outfile=my_dir+'fluxmultregrid1.spw'+spw_sd+'.I.sd1.image', order='0132')

tt=imhead(data_TP+'.spw'+spw_sd+'.I.sd.fits')
U=tt['unit']
RMAJ=str(tt['restoringbeam']['major']['value'])+'arcsec'
RMIN=str(tt['restoringbeam']['minor']['value'])+'arcsec'
RPA=str(tt['restoringbeam']['positionangle']['value'])+'deg'

imhead(my_dir+'fluxmultregrid1.spw'+spw_sd+'.I.sd1.image',mode='put',hdkey='bunit',hdvalue=U)
imhead(my_dir+'fluxmultregrid1.spw'+spw_sd+'.I.sd1.image',mode='add',hdkey='bmaj',hdvalue=RMAJ)
imhead(my_dir+'fluxmultregrid1.spw'+spw_sd+'.I.sd1.image',mode='put',hdkey='bmin',hdvalue=RMIN)
imhead(my_dir+'fluxmultregrid1.spw'+spw_sd+'.I.sd1.image',mode='put',hdkey='bpa',hdvalue=RPA)


#convert to Jy/pixel
##rest_ghz=float(myrestfreq.replace('GHz',''))
#sd_res=74.88/(rest_ghz/100.)/(12./10.)
#cell_arc=float(mycell.replace('arcsec',''))
#Conv=(cell_arc**2)/(1.133*sd_res**2)

#im1=ia.imagecalc(outfile=my_dir+'fluxmultregrid2.spw'+spw_sd+'.I.sd.image', pixels=str(Conv)+'*fluxmultregrid1.spw'+spw_sd+'.I.sd.image',overwrite=True)
#im1.done()
#ia.close()
#imhead(my_dir+'fluxmultregrid2.spw'+spw_sd+'.I.sd.image',mode='put',hdkey='bunit',hdvalue='Jy/pixel')


#imtrans(imagename=my_dir+'fluxmultregrid0.spw'+spw_sd+'.I.sd.image', outfile=my_dir+'fluxmultregrid1.spw'+spw_sd+'.I.sd.image', order='0132')

#####
#step 4- clean using this new file as a modelimage
#####
#13CO 2 is with Jy/beam, 3 is with Jy/pixel in model image, 21 is with Jy/beam and no pb cor, 22 is with JY/beam, casa method PB, all spws from 7m and no pbcor, 22pix same as 22 but jy/pixwl, tp2new is same as 22 but with added alternate method

spw_arr='8' #spw for array
spw_7m='8,17,26'
myniter=20000
mythreshold='8mJy'
mynchan=60
mystart=40
mywidth=1
myimsize=800
mycell='0.2arcsec'
myrestfreq='220.39868GHz'

#clean(vis=[data_12m,data_7m],imagename=line+"_combinewithtp21",outlierfile="",field=['4~36', '5~18'],spw=spw_arr, selectdata=True,timerange="",uvrange="",antenna="",scan="",observation="",intent="",  mode="channel",resmooth=False,gridmode="",wprojplanes=-1,facets=1,cfcache="cfcache.dir", rotpainc=5.0,painc=360.0,aterm=True,psterm=False,mterm=True,wbawp=False,conjbeams=True, epjtable="",interpolation="linear", niter=myniter,gain=0.1,threshold=mythreshold,psfmode="clark",imagermode="csclean", ftmachine="mosaic",mosweight=False,scaletype="SAULT",multiscale=[0, 1, 5, 10, 15], negcomponent=-1,smallscalebias=0.6, interactive=False,mask=[], nchan=mynchan,start=mystart,width=mywidth,outframe="lsrk",veltype="radio", imsize=[myimsize, myimsize],cell=mycell, phasecenter="J2000 19h15m38.305s 10d41m01.018s", restfreq=myrestfreq, stokes="I",weighting="natural",robust=0.0,uvtaper=False,outertaper=[''],innertaper=['1.0'],modelimage=my_dir+'fluxmultregrid1.spw'+spw_sd+'.I.sd1.image', restoringbeam=[''], pbcor=False,minpb=0.2,usescratch=False,noise="1.0Jy", npixels=0,npercycle=100,cyclefactor=2.0,cyclespeedup=-1,nterms=1,reffreq="",chaniter=False, flatnoise=True,allowchunk=False)
clean(vis=['data_12m_CO13.ms','data_7m_CO13_comb.ms'],imagename=line+"_combinewithtp22pix",outlierfile="",field=['0~32', '0~13'],spw=['0','0'], selectdata=True,timerange="",uvrange="",antenna="",scan="",observation="",intent="",  mode="channel",resmooth=False,gridmode="",wprojplanes=-1,facets=1,cfcache="cfcache.dir", rotpainc=5.0,painc=360.0,aterm=True,psterm=False,mterm=True,wbawp=False,conjbeams=True, epjtable="",interpolation="linear", niter=myniter,gain=0.1,threshold=mythreshold,psfmode="clark",imagermode="csclean", ftmachine="mosaic",mosweight=False,scaletype="SAULT",multiscale=[0, 1, 5, 10, 15], negcomponent=-1,smallscalebias=0.6, interactive=False,mask=[], nchan=mynchan,start=mystart,width=mywidth,outframe="lsrk",veltype="radio", imsize=[myimsize, myimsize],cell=mycell, phasecenter="J2000 19h15m38.305s 10d41m01.018s", restfreq=myrestfreq, stokes="I",weighting="natural",robust=0.0,uvtaper=False,outertaper=[''],innertaper=['1.0'],modelimage='B17m'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.subim.trans.debpb.pix', restoringbeam=[''], pbcor=False,minpb=0.2,usescratch=False,noise="1.0Jy", npixels=0,npercycle=100,cyclefactor=2.0,cyclespeedup=-1,nterms=1,reffreq="",chaniter=False, flatnoise=True,allowchunk=False)

immath(imagename=['B'+line+"_combinewithtp22pix.image",
                  'B'+line+"_combinewithtp22pix.flux"],
       expr='IM0/IM1',
       outfile='GRS1915_modelimg7mpix_CO13.image.pbcor')

###other method
# get the positive interferometer-only clean components
imsubimage(imagename='B'+line+"_combinewithtp22.model",
           outfile='B'+line+"_combinewithtp22.model.subim",
           region='box [[19:15:42.85546, +010.40.10.6834], [19:15:34.57909, +010.42.06.5645]]')
imsubimage(imagename='B'+line+"_combinewithtp22.flux",
           outfile='B'+line+"_combinewithtp22.flux0.subim",
           region='box [[19:15:42.85546, +010.40.10.6834], [19:15:34.57909, +010.42.06.5645]]')
imsubimage(imagename='B17m'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.subim.trans.debpb', outfile='B17m'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.subim.trans.debpb.chancut',chans='40~99')
rest_ghz=float(myrestfreq.replace('GHz',''))
sd_res=74.88/(rest_ghz/100.)/(12./10.)
cell_arc=float(mycell.replace('arcsec',''))
Conv=(cell_arc**2)/(1.133*sd_res**2)
im1=ia.imagecalc(outfile='B17m'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.subim.trans.debpb.chancut.pix', pixels=str(Conv)+'*B17m'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.subim.trans.debpb.chancut',overwrite=True)
im1.done()
ia.close()
imhead('B17m'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.subim.trans.debpb.chancut.pix',mode='put',hdkey='bunit',hdvalue='Jy/pixel')

immath(outfile='deconvolved-sdinput.intmodel',
   imagename=['B'+line+"_combinewithtp22.model.subim",
              'B17m'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.subim.trans.debpb.chancut.pix'],
   mode='evalexpr',
   expr='iif((IM0-IM1) >= 0.00, IM0-IM1, 0.0)')

# remove those components if they are at the map edge
immath(outfile='deconvolved-sdinput-masked.intmodel',
   imagename=['deconvolved-sdinput.intmodel',
               'B'+line+"_combinewithtp22.flux0.subim"],
   mode='evalexpr',
   expr='iif((IM1) >= 0.25, IM0, 0.0)')

imhead('deconvolved-sdinput-masked.intmodel',
       mode='put',
       hdkey='bunit', hdvalue='Jy/pixel')


# smooth the interferometer-only components to the synthesized beam
SynthBeamMaj = imhead('B'+line+"_combinewithtp22.image", mode='get', hdkey='bmaj')['value']
SynthBeamMin = imhead('B'+line+"_combinewithtp22.image", mode='get', hdkey='bmin')['value']
SynthBeamPA = imhead('B'+line+"_combinewithtp22.image", mode='get', hdkey='bpa')['value']


imsmooth(imagename='deconvolved-sdinput-masked.intmodel',
         outfile='deconvolved-sdinput.intimage',
         kernel='gauss',
         major=str(SynthBeamMaj)+'arcsec',
         minor=str(SynthBeamMin)+'arcsec',
         pa=str(SynthBeamPA)+'deg')

# feather the data
os.system('rm -rf deconvolved-combi.model')
feather(imagename = 'deconvolved-combi.model',
        highres = 'deconvolved-sdinput.intimage',
        lowres = 'B17m'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.subim.trans.debpb')

# clean up and keep the final model
os.system('rm -rf combimodel.image')
os.system('mv deconvolved-combi.model combimodel.image')


clean(vis=['data_12m_CO13.ms','data_7m_CO13_comb.ms'],imagename='B'+line+"_combinewithtp2new",outlierfile="",field=['0~32', '0~13'],spw=['0','0'], selectdata=True,timerange="",uvrange="",antenna="",scan="",observation="",intent="",  mode="channel",resmooth=False,gridmode="",wprojplanes=-1,facets=1,cfcache="cfcache.dir", rotpainc=5.0,painc=360.0,aterm=True,psterm=False,mterm=True,wbawp=False,conjbeams=True, epjtable="",interpolation="linear", niter=myniter,gain=0.1,threshold=mythreshold,psfmode="clark",imagermode="csclean", ftmachine="mosaic",mosweight=False,scaletype="SAULT",multiscale=[0, 1, 5, 10, 15], negcomponent=-1,smallscalebias=0.6, interactive=False,mask=[], nchan=mynchan,start=mystart,width=mywidth,outframe="lsrk",veltype="radio", imsize=[myimsize, myimsize],cell=mycell, phasecenter="J2000 19h15m38.305s 10d41m01.018s", restfreq=myrestfreq, stokes="I",weighting="natural",robust=0.0,uvtaper=False,outertaper=[''],innertaper=['1.0'],modelimage='combimodel.image', restoringbeam=[''], pbcor=False,minpb=0.2,usescratch=False,noise="1.0Jy", npixels=0,npercycle=100,cyclefactor=2.0,cyclespeedup=-1,nterms=1,reffreq="",chaniter=False, flatnoise=True,allowchunk=False)

immath(imagename=['B'+line+"_combinewithtp2new.image",
                  'B'+line+"_combinewithtp2new.flux"],
       expr='IM0/IM1',
       outfile='GRS1915_modelimg7mnew_CO13.image.pbcor')

########


im.open('data_12m_CO13.ms')  
im.selectvis(spw='8',field='4~32')  
im.defineimage(nx=800, cellx='0.2arcsec', phasecenter="J2000 19h15m38.305s 10d41m01.018s")  
im.setvp(dovp=T)  
im.setoptions(ftmachine='mosaic')  
#im.setscales(nscales=3)  
#im.setsdoptions(scale=0.9);  
im.makemodelfromsd(sdimage=my_dir+'_fluxmultregrid0.spw'+spw_sd+'.I.sd1.image', modelimage='co13_model')  
  
 


