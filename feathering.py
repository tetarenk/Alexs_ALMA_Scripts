# CASA Feathering ALMA data, shows process

#####step 1- combine and image the 7m + 12m
data_12m='/mnt/bigdata/tetarenk/ALMA_GRS1915_105/calib_MS/calibrated_final.ms'
data_7m='/mnt/bigdata/tetarenk/ALMA_GRS1915_105/calib_MS_7m/calibrated_final.ms'
data_TP='/mnt/bigdata/tetarenk/ALMA_GRS1915_105/2015.1.00976.S/science_goal.uid___A001_X2d6_X61/group.uid___A001_X2d6_X62/member.uid___A001_X2d6_X69/product/Radio_Peak_of_IRAS_19132+1035'
my_dir='/mnt/bigdata/tetarenk/ALMA_GRS1915_105/alex_redo/'

#listobs
listobs(data_12m,listfile='12m_listfile.txt')
listobs(data_7m,listfile='7m_listfile.txt')
os.system('sudo pluma 12m_listfile.txt &')
os.system('sudo pluma 7m_listfile.txt &')

#show pointings
au.plotmosaic(data_12m,sourceid='3',figfile='12m_mosaic.png')
au.plotmosaic(data_7m,sourceid='4',figfile='7m_mosaic.png')

####13CO
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

plotms(vis=data_12m,yaxis='wt',xaxis='uvdist',spw='0',
       coloraxis='spw',plotfile='12m_WT_13CO.png',overwrite=True)
#
plotms(vis=data_7m,yaxis='wt',xaxis='uvdist',spw='8,17,26',
       coloraxis='spw',plotfile='7m_WT_13CO.png',overwrite=True)

plotms(vis=data_12m,yaxis='amp',xaxis='velocity',spw='8', avgtime='1e8',avgscan=True,coloraxis='spw',avgchannel='5',
       transform=True,freqframe='LSRK',restfreq=myrestfreq, plotfile='12m_13CO_vel.png',overwrite=True)
plotms(vis=data_7m,yaxis='amp',xaxis='velocity',spw='8,17,26', avgtime='1e8',avgscan=True,coloraxis='spw',avgchannel='5',
       transform=True,freqframe='LSRK',restfreq=myrestfreq, plotfile='7m_13CO_vel.png',overwrite=True)
#SiO myrestfreq='217.105GHz'
plotms(vis=data_12m,yaxis='amp',xaxis='velocity',spw='3', avgtime='1e8',avgscan=True,coloraxis='spw',avgchannel='5',
       transform=True,freqframe='LSRK',restfreq=myrestfreq, plotfile='12m_SiO_vel.png',overwrite=True)
plotms(vis=data_7m,yaxis='amp',xaxis='velocity',spw='3,12,21', avgtime='1e8',avgscan=True,coloraxis='spw',avgchannel='5',
       transform=True,freqframe='LSRK',restfreq=myrestfreq, plotfile='7m_SiO_vel.png',overwrite=True)
#12CO myrestfreq='230.538GHz'
plotms(vis=data_12m,yaxis='amp',xaxis='velocity',spw='0', avgtime='1e8',avgscan=True,coloraxis='spw',avgchannel='5',
       transform=True,freqframe='LSRK',restfreq=myrestfreq, plotfile='12m_12CO_vel.png',overwrite=True)
plotms(vis=data_7m,yaxis='amp',xaxis='velocity',spw='0,9,18', avgtime='1e8',avgscan=True,coloraxis='spw',avgchannel='5',
       transform=True,freqframe='LSRK',restfreq=myrestfreq, plotfile='7m_12CO_vel.png',overwrite=True)
#18CO myrestfreq='219.56035GHz'
plotms(vis=data_12m,yaxis='amp',xaxis='velocity',spw='7', avgtime='1e8',avgscan=True,coloraxis='spw',avgchannel='5',
       transform=True,freqframe='LSRK',restfreq=myrestfreq, plotfile='12m_18CO_vel.png',overwrite=True)
plotms(vis=data_7m,yaxis='amp',xaxis='velocity',spw='7,16,25', avgtime='1e8',avgscan=True,coloraxis='spw',avgchannel='5',
       transform=True,freqframe='LSRK',restfreq=myrestfreq, plotfile='7m_18CO_vel.png',overwrite=True)
#N2D+ myrestfreq='231.32183GHz'
plotms(vis=data_12m,yaxis='amp',xaxis='velocity',spw='1', avgtime='1e8',avgscan=True,coloraxis='spw',avgchannel='5',
       transform=True,freqframe='LSRK',restfreq=myrestfreq, plotfile='12m_N2D_vel.png',overwrite=True)
plotms(vis=data_7m,yaxis='amp',xaxis='velocity',spw='1,10,19', avgtime='1e8',avgscan=True,coloraxis='spw',avgchannel='5',
       transform=True,freqframe='LSRK',restfreq=myrestfreq, plotfile='7m_N2D_vel.png',overwrite=True)
#H30alpha myrestfreq='231.90093GHz'
plotms(vis=data_12m,yaxis='amp',xaxis='velocity',spw='2', avgtime='1e8',avgscan=True,coloraxis='spw',avgchannel='5',
       transform=True,freqframe='LSRK',restfreq=myrestfreq, plotfile='12m_Hal_vel.png',overwrite=True)
plotms(vis=data_7m,yaxis='amp',xaxis='velocity',spw='2,11,20', avgtime='1e8',avgscan=True,coloraxis='spw',avgchannel='5',
       transform=True,freqframe='LSRK',restfreq=myrestfreq, plotfile='7m_Hal_vel.png',overwrite=True)
#303-202 myrestfreq='218.22219GHz'
plotms(vis=data_12m,yaxis='amp',xaxis='velocity',spw='4', avgtime='1e8',avgscan=True,coloraxis='spw',avgchannel='5',
       transform=True,freqframe='LSRK',restfreq=myrestfreq, plotfile='12m_303_vel.png',overwrite=True)
plotms(vis=data_7m,yaxis='amp',xaxis='velocity',spw='4,13,22', avgtime='1e8',avgscan=True,coloraxis='spw',avgchannel='5',
       transform=True,freqframe='LSRK',restfreq=myrestfreq, plotfile='7m_303_vel.png',overwrite=True)
#322-221 myrestfreq='218.47563GHz'
plotms(vis=data_12m,yaxis='amp',xaxis='velocity',spw='5', avgtime='1e8',avgscan=True,coloraxis='spw',avgchannel='5',
       transform=True,freqframe='LSRK',restfreq=myrestfreq, plotfile='12m_322_vel.png',overwrite=True)
plotms(vis=data_7m,yaxis='amp',xaxis='velocity',spw='5,14,23', avgtime='1e8',avgscan=True,coloraxis='spw',avgchannel='5',
       transform=True,freqframe='LSRK',restfreq=myrestfreq, plotfile='7m_322_vel.png',overwrite=True)
#321-220 myrestfreq='218.76007GHz'
plotms(vis=data_12m,yaxis='amp',xaxis='velocity',spw='6', avgtime='1e8',avgscan=True,coloraxis='spw',avgchannel='5',
       transform=True,freqframe='LSRK',restfreq=myrestfreq, plotfile='12m_321_vel.png',overwrite=True)
plotms(vis=data_7m,yaxis='amp',xaxis='velocity',spw='6,15,24', avgtime='1e8',avgscan=True,coloraxis='spw',avgchannel='5',
       transform=True,freqframe='LSRK',restfreq=myrestfreq, plotfile='7m_321_vel.png',overwrite=True)


#make combined 12m and 7m image, dont do PB 
#clean(vis=[data_12m,data_7m],imagename=line+"_combine12m7m_nppb",outlierfile="",field=['4~36', '5~18'],spw=spw_arr, selectdata=True,timerange="",uvrange="",antenna="",scan="",observation="",intent="",  mode="channel",resmooth=False,gridmode="",wprojplanes=-1,facets=1,cfcache="cfcache.dir", rotpainc=5.0,painc=360.0,aterm=True,psterm=False,mterm=True,wbawp=False,conjbeams=True, epjtable="",interpolation="linear", niter=myniter,gain=0.1,threshold=mythreshold,psfmode="clark",imagermode="csclean", ftmachine="mosaic",mosweight=False,scaletype="SAULT",multiscale=[0, 1, 5, 10, 15], negcomponent=-1,smallscalebias=0.6, interactive=False,mask=[], nchan=mynchan,start=mystart,width=mywidth,outframe="lsrk",veltype="radio", imsize=[myimsize, myimsize],cell=mycell, phasecenter="J2000 19h15m38.305s 10d41m01.018s", restfreq=myrestfreq, stokes="I",weighting="natural",robust=0.0,uvtaper=False,outertaper=[''],innertaper=['1.0'],modelimage='', restoringbeam=[''], pbcor=False,minpb=0.2,usescratch=False,noise="1.0Jy", npixels=0,npercycle=100,cyclefactor=2.0,cyclespeedup=-1,nterms=1,reffreq="",chaniter=False, flatnoise=True,allowchunk=False)
clean(vis=['data_12m_CO13.ms','data_7m_CO13_comb.ms'],imagename=line+"_combine12m7m_nopb_all7",outlierfile="",field=['0~32', '0~13'],spw=['0','0'], selectdata=True,timerange="",uvrange="",antenna="",scan="",observation="",intent="",  mode="channel",resmooth=False,gridmode="",wprojplanes=-1,facets=1,cfcache="cfcache.dir", rotpainc=5.0,painc=360.0,aterm=True,psterm=False,mterm=True,wbawp=False,conjbeams=True, epjtable="",interpolation="linear", niter=myniter,gain=0.1,threshold=mythreshold,psfmode="clark",imagermode="csclean", ftmachine="mosaic",mosweight=False,scaletype="SAULT",multiscale=[0, 1, 5, 10, 15], negcomponent=-1,smallscalebias=0.6, interactive=False,mask=[], nchan=mynchan,start=mystart,width=mywidth,outframe="lsrk",veltype="radio", imsize=[myimsize, myimsize],cell=mycell, phasecenter="J2000 19h15m38.305s 10d41m01.018s", restfreq=myrestfreq, stokes="I",weighting="natural",robust=0.0,uvtaper=False,outertaper=[''],innertaper=['1.0'],modelimage='', restoringbeam=[''], pbcor=False,minpb=0.2,usescratch=False,noise="1.0Jy", npixels=0,npercycle=100,cyclefactor=2.0,cyclespeedup=-1,nterms=1,reffreq="",chaniter=False, flatnoise=True,allowchunk=False)
#regrid TP
importfits(fitsimage=data_TP+'.spw'+spw_sd+'.I.sd.fits', imagename=my_dir+'Radio_Peak_of_IRAS_19132+1035.spw'+spw_sd+'.I.sd.image')
imregrid(imagename=my_dir+'Radio_Peak_of_IRAS_19132+1035.spw'+spw_sd+'.I.sd.image', template=line+"_combine12m7m_nppb.image", output='F'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.regrid', asvelocity=False,overwrite=True)
imhead('F'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.regrid')

#cut out clean masked stuff
imsubimage(imagename=line+"_combine12m7m_nppb.image",
           outfile=line+"_combine12m7m_nppb.image.subim",
           region='box [[19:15:42.85546, +010.40.10.6834], [19:15:34.57909, +010.42.06.5645]]')
imsubimage(imagename='F'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.regrid',
           outfile='F'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.regrid.subim',
           region='box [[19:15:42.85546, +010.40.10.6834], [19:15:34.57909, +010.42.06.5645]]')

#multiply by PB
imsubimage(imagename=line+"_combine12m7m_nppb.flux",
           outfile=line+"_combine12m7m_nppb.flux.subim",
           region='box [[19:15:42.85546, +010.40.10.6834], [19:15:34.57909, +010.42.06.5645]]')
imtrans(imagename='F'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.regrid.subim', outfile='F'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.regrid.subim.trans', order='0132')


#imtrans(imagename='F'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.regrid', outfile='F'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.regrid.t2', order='0132')
#imhead('F'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.regrid.t2',mode='put',hdkey='bunit',hdvalue=U)
#imhead('F'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.regrid.t2',mode='add',hdkey='bmaj',hdvalue=RMAJ)
#imhead('F'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.regrid.t2',mode='put',hdkey='bmin',hdvalue=RMIN)
#imhead('F'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.regrid.t2',mode='put',hdkey='bpa',hdvalue=RPA)
#imsubimage(imagename='F'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.regrid.t2',
           #outfile='F'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.regrid.t2.subim',
           #region='ellipse[[19:15:38.60125, +010.40.59.1847], [45.8646arcsec, 44.1345arcsec], 90.00000000deg]')



tt=imhead(data_TP+'.spw'+spw_sd+'.I.sd.fits')
U=tt['unit']
RMAJ=str(tt['restoringbeam']['major']['value'])+'arcsec'
RMIN=str(tt['restoringbeam']['minor']['value'])+'arcsec'
RPA=str(tt['restoringbeam']['positionangle']['value'])+'deg'

imhead('F'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.regrid.subim.trans',mode='put',hdkey='bunit',hdvalue=U)
imhead('F'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.regrid.subim.trans',mode='add',hdkey='bmaj',hdvalue=RMAJ)
imhead('F'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.regrid.subim.trans',mode='put',hdkey='bmin',hdvalue=RMIN)
imhead('F'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.regrid.subim.trans',mode='put',hdkey='bpa',hdvalue=RPA)

immath(imagename=['F'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.regrid.subim.trans',
                  'F'+line+"_combine12m7m_nppb.flux.subim"],
       expr='IM0*IM1',
       outfile='F'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.regrid.subim.trans.depb')

#feathering task
#feather(imagename='GRS_Feather_CO13.image',
        #highres=line+"_combine12m7m.image.subim",
        #lowres='F'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.regrid.depb')
#feather(imagename='GRS_Feather_CO133.image',
        #highres=line+"_combine12m7m.image",
        #lowres='F'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.regrid.t2.subim')
feather(imagename='GRS1915_Feather_CO13.image',
        highres=line+"_combine12m7m_nppb.image.subim",
        lowres='F'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.regrid.subim.trans.depb')

immath(imagename=['GRS1915_Feather_CO13.image',
                  'F'+line+"_combine12m7m_nppb.flux.subim"],
       expr='IM0/IM1',
       outfile='GRS1915_Feather_CO13.image.pbcor')

#######including all 7m and using same total power as model image
imsubimage(imagename=line+'_combine12m7m_nopb_all7.image',
           outfile=line+'_combine12m7m_nopb_all7.image.subim',
           region='box [[19:15:42.85546, +010.40.10.6834], [19:15:34.57909, +010.42.06.5645]]')
imsubimage(imagename=line+'_combine12m7m_nopb_all7.flux',
           outfile='B17'+line+'_combine12m7m_nopb_all7.flux.subim',
           region='box [[19:15:42.85546, +010.40.10.6834], [19:15:34.57909, +010.42.06.5645]]')
feather(imagename='GRS1915_Featherall7_CO13.image',
        highres=line+"_combine12m7m_nopb_all7.image.subim",
        lowres='B17m'+line+'_TP'+'_regrid.spw'+spw_sd+'.I.sd.image.subim.trans.debpb')

immath(imagename=['GRS1915_Featherall7_CO13.image',
                 'B17'+line+'_combine12m7m_nopb_all7.flux.subim'],
       expr='IM0/IM1',
       outfile='GRS1915_Featherall7_CO13.image.pbcor')
