from spectral_cube import SpectralCube
import matplotlib.pyplot as plt
import astropy.units as u
import numpy as np
import pyfits
import scipy.ndimage as nd
import pyspeckit

datadir = '/mnt/bigdata/tetarenk/ALMA_GRS1915_105/'
lines=['13CO','12CO','18CO','H2CO_303_202','H2CO_322_221','H2CO_321_220','SiO','N2Dp','H30a']


def make_spec(fitsfile,ddir,limits,flag):
	x0,x1,y0,y1=limits[0],limits[1],limits[2],limits[3]
	a=SpectralCube.read(fitsfile)
	a.allow_huge_operations=True
	a1=a.to(u.K, equivalencies=a.beam.jtok_equiv(a.header['RESTFRQ']*u.Hz))
	a2=a1.with_spectral_unit(u.km/u.s, velocity_convention='radio')
	if flag == 'reg':
		suba=a2[:, y0:y1, x0:x1]
		veloc=suba.spectral_axis
		#print veloc
		amps = np.mean(suba.filled_data[:],axis=(1,2))
	elif flag == 'pix':
		suba=a2[:, y0, x0]
		veloc=suba.spectral_axis
		#print veloc
		amps=a2[:, y0, x0].value
	return(veloc.to(u.km/u.s),amps)

def find_peak_pix(tmaxfits,vals=100):
	#get tmax map to find brightest pixel
	tmax=pyfits.getdata(tmaxfits)
	#mask where tmax is 0
	mask=(tmax>=0.0)
	#ndimage package - multidimensional imaging. 
	#Non-zero (true) elements of the mask are dilated, np.ones are used to dilate
	#essentially got rid of noise at edges which would interfere with finding max of the signal.
	# ~ reverses the true and false elements
	data=nd.morphology.binary_erosion(mask, np.ones((vals, vals)))
	#plt.imshow(data*tmax) #shows tmax map with edges removed. Need to adjust size of np.ones based off image)
	maxim=np.nanmax(data*tmax)#finds max value
	maxind=np.nanargmax(data*tmax)#finds position of max - 1d
	#get back pixels in (x,y) of max, then assign ymax and xmax to those values
	ymax,xmax=np.unravel_index(maxind, tmax.shape)
	print xmax,ymax
	#fig=plt.figure()
	#a=plt.imshow(tmax)
	#cbar=fig.colorbar(a)
	#plt.gca().invert_yaxis()
	#plt.show()
	#raw_input('stop')
	return(xmax,ymax)

'''fig=plt.figure()
a=plt.imshow(tmax)
cbar=fig.colorbar(a)
plt.gca().invert_yaxis()
plt.show()'''


def fit_gauss_to_spec(fitsfile,ddir,line,tmaxfits):
	a=SpectralCube.read(fitsfile)
	a1=a.to(u.K, equivalencies=a.beam.jtok_equiv(a.header['RESTFRQ']*u.Hz))
	a2=a1.with_spectral_unit(u.km/u.s, velocity_convention='radio')
	xmax,ymax=find_peak_pix(tmaxfits)
	asp=pyspeckit.Spectrum(data=a2[:, ymax, xmax].value, xarr=a2.spectral_axis, xarrkwargs={'unit':'km/s'}, unit="$T_B$")
	aamp=a2[:, ymax, xmax].value.max()
	#amplitude guess
	acenter=(a2[:, ymax, xmax].value*a2.spectral_axis.value).sum()/a2[:, ymax, xmax].value.sum()
	#center guess
	awidth=a2[:, ymax, xmax].value.sum()/aamp/np.sqrt(2*np.pi)
	#width guess
	aguess=[aamp, acenter, awidth]
	asp.specfit(fittype='gaussian', guesses=aguess)
	asp.plotter(errstyle='fill')
	asp.specfit.plot_fit()
	asp.plotter.savefig(ddir+'spectra_plots/'+line+'_fitted.eps')
	#asp.specfit.modelpars/modelerrs give parameters and errors of spectrum fit
	return(asp.specfit.modelpars,asp.specfit.modelerrs)
	

'''
#single spectra for a specific line and a specific pixel or region
line='13CO'
lim=[385,400,560,575]
vel,amp=make_spec(datadir+'alex_imaging_'+line+'_fix/GRS1915_modelimg_'+line+'.image.pbcor.fits',datadir,lim,'reg')
fig=plt.figure()
plt.rcdefaults()
plt.plot(vel, amp, label=line, color='b', ls='-.',lw=3)
plt.ylabel('$T_B$ (K)')
plt.xlabel('$v$ (km/s)')
plt.legend(numpoints=1,loc='upper left')
plt.ylim(0,20)
#fig.savefig(datadir+'spectra_plots/'+line+'.png',format='png')
plt.show()
raw_input('enter')

#single spectra over HMSR for a specific line, do only for full spw cleans
line='18CO'
lim=[340,360,470,485]
vel,amp=make_spec(datadir+'alex_imaging_'+line+'_fix_all/GRS1915_modelimg_'+line+'.image.pbcor.fits',datadir,lim,'reg')
fig=plt.figure()
plt.plot(vel, amp, label=line, color='b', ls='-.',lw=3)
plt.ylabel('$T_B$ (K)')
plt.xlabel('$v$ (km/s)')
plt.legend(numpoints=1,loc='upper left')
plt.ylim(0,20)
plt.savefig(datadir+'spectra_plots/'+line+'_hmsr.png',format='png')
plt.show()
print np.trapz(y=amp,x=vel)
raw_input('enter')

#single spectra over HMSR for a specific line at peak pixel
line='H30alpha'
lim=[361,361,469,469]
vel,amp=make_spec(datadir+'alex_imaging_'+line+'_fix/GRS1915_modelimg_'+line+'.image.pbcor.fits',datadir,lim,'pix')
fig=plt.figure()
plt.plot(vel, amp, label=line, color='b', ls='-.',lw=3)
plt.ylabel('$T_B$ (K)')
plt.xlabel('$v$ (km/s)')
plt.legend(numpoints=1,loc='upper left')
plt.ylim(0,0.2)
plt.savefig(datadir+'spectra_plots/'+line+'_hmsr_peak.png',format='png')
plt.show()
raw_input('enter')

#getting formaldhyde triplet spectra at max pixel
tmaxfits=datadir+'T_max_maps/H2CO_303_202_tmax.fits'
xmax,ymax=find_peak_pix(tmaxfits)
line1='H2CO_303_202'
line2='H2CO_322_221'
line3='H2CO_321_220'
lim=[xmax,xmax,ymax,ymax]
vel1,amp1=make_spec(datadir+'alex_imaging_'+line1+'_fix/GRS1915_modelimg_'+line1+'.image.pbcor.fits',datadir,lim,'pix')
vel2,amp2=make_spec(datadir+'alex_imaging_'+line2+'_fix/GRS1915_modelimg_'+line2+'.image.pbcor.fits',datadir,lim,'pix')
vel3,amp3=make_spec(datadir+'alex_imaging_'+line3+'_fix/GRS1915_modelimg_'+line3+'.image.pbcor.fits',datadir,lim,'pix')
fig=plt.figure()
plt.plot(vel1, amp1, label=line1, color='m', ls='-',lw=3)
plt.plot(vel2, amp2, label=line2, color='r', ls='-.',lw=3)
plt.plot(vel3, amp3, label=line3, color='b', ls='--',lw=3)
plt.ylabel('$T_B$ (K)')
plt.xlabel('$v$ (km/s)')
plt.xlim(64,71)
plt.ylim(0,1.5)
plt.legend(numpoints=1,loc='upper left')
#plt.ylim(0,2)
plt.savefig(datadir+'spectra_plots/H2CO_max.pdf',format='pdf')
plt.show()
raw_input('enter')'''
#guassian fitting of a line at max pixel
'''tmaxfits=datadir+'T_max_maps/H2CO_303_202_tmax.fits'
line='13CO'
fitsfile=datadir+'alex_imaging_'+line+'_fix/GRS1915_modelimg_'+line+'.image.pbcor.fits'
pars,errs=fit_gauss_to_spec(fitsfile,datadir,line,tmaxfits)'''


#regions CO big figure
fig=plt.figure()
ax=fig.add_subplot(111)
ax.set_xticks([])
ax.set_yticks([])
plt.setp(ax.get_yticklabels(), visible=False)
plt.setp(ax.get_xticklabels(), visible=False)
ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['right'].set_color('none')
ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
ax.set_ylabel('$T_B$ (K)',labelpad=20)
ax.set_xlabel('$v$ (km/s)',labelpad=20)
line13='13CO'
line12='12CO'
line18='18CO'
ax1=fig.add_subplot(221)#ridge
lim=[385,400,560,575]
v13,a13=make_spec(datadir+'alex_imaging_'+line13+'_fix/GRS1915_modelimg_'+line13+'.image.pbcor.fits',datadir,lim,'reg')
v12,a12=make_spec(datadir+'alex_imaging_'+line12+'_fix/GRS1915_modelimg_'+line12+'.image.pbcor.fits',datadir,lim,'reg')
v18,a18=make_spec(datadir+'alex_imaging_'+line18+'_fix/GRS1915_modelimg_'+line18+'.image.pbcor.fits',datadir,lim,'reg')
ax1.plot(v12, a12, label='$^{12}{\\rm CO}$', color='#848484', ls='--',lw=3)
ax1.plot(v13, a13, label='$^{13}{\\rm CO}$', color='#424242', ls=':',lw=3)
ax1.plot(v18, a18, label='${\\rm C}^{18}{\\rm O}$', color='k', ls='-',lw=3)
plt.setp(ax1.get_xticklabels(), visible=False)
ax1.text(58, 22, 'A')
ax1.set_ylim(0,25)
ax2=fig.add_subplot(222,sharey=ax1)#nonthermal
lim=[400,420,570,595]
v13,a13=make_spec(datadir+'alex_imaging_'+line13+'_fix/GRS1915_modelimg_'+line13+'.image.pbcor.fits',datadir,lim,'reg')
v12,a12=make_spec(datadir+'alex_imaging_'+line12+'_fix/GRS1915_modelimg_'+line12+'.image.pbcor.fits',datadir,lim,'reg')
v18,a18=make_spec(datadir+'alex_imaging_'+line18+'_fix/GRS1915_modelimg_'+line18+'.image.pbcor.fits',datadir,lim,'reg')
ax2.plot(v12, a12, label='$^{12}{\\rm CO}$', color='#848484', ls='--',lw=3)
ax2.plot(v13, a13, label='$^{13}{\\rm CO}$', color='#424242', ls=':',lw=3)
ax2.plot(v18, a18, label='${\\rm C}^{18}{\\rm O}$', color='k', ls='-',lw=3)
plt.setp(ax2.get_xticklabels(), visible=False)
plt.setp(ax2.get_yticklabels(), visible=False)
ax2.text(58, 22, 'B')
ax2.set_ylim(0,25)
ax3=fig.add_subplot(223,sharex=ax1)#hmsr
lim=[340,360,470,485]
v13,a13=make_spec(datadir+'alex_imaging_'+line13+'_fix/GRS1915_modelimg_'+line13+'.image.pbcor.fits',datadir,lim,'reg')
v12,a12=make_spec(datadir+'alex_imaging_'+line12+'_fix/GRS1915_modelimg_'+line12+'.image.pbcor.fits',datadir,lim,'reg')
v18,a18=make_spec(datadir+'alex_imaging_'+line18+'_fix/GRS1915_modelimg_'+line18+'.image.pbcor.fits',datadir,lim,'reg')
ax3.plot(v12, a12, label='$^{12}{\\rm CO}$', color='#848484', ls='--',lw=3)
ax3.plot(v13, a13, label='$^{13}{\\rm CO}$', color='#424242', ls=':',lw=3)
ax3.plot(v18, a18, label='${\\rm C}^{18}{\\rm O}$', color='k', ls='-',lw=3)
ax3.text(58, 22, 'C')
ax3.set_ylim(0,25)
ax4=fig.add_subplot(224)#off
#lim=[490,505,490,505]
lim=[507,580,447,492]
v13,a13=make_spec(datadir+'alex_imaging_'+line13+'_fix/GRS1915_modelimg_'+line13+'.image.pbcor.fits',datadir,lim,'reg')
v12,a12=make_spec(datadir+'alex_imaging_'+line12+'_fix/GRS1915_modelimg_'+line12+'.image.pbcor.fits',datadir,lim,'reg')
v18,a18=make_spec(datadir+'alex_imaging_'+line18+'_fix/GRS1915_modelimg_'+line18+'.image.pbcor.fits',datadir,lim,'reg')
ax4.plot(v12, a12, label='$^{12}{\\rm CO}$', color='#848484', ls='--',lw=3)
ax4.plot(v13, a13, label='$^{13}{\\rm CO}$', color='#424242', ls=':',lw=3)
ax4.plot(v18, a18, label='${\\rm C}^{18}{\\rm O}$', color='k', ls='-',lw=3)
ax4.text(58, 22, 'D')
ax4.set_ylim(0,25)
plt.setp(ax4.get_yticklabels(), visible=False)
ax4.set_xticks([60,65,70,75,80,85])
ax3.set_xticks([55,60,65,70,75,80])
ax3.set_yticks([0,5,10,15,20])
ax1.set_yticks([0,5,10,15,20,25])
ax1.tick_params(axis='both', which='major', length=5,width=2)
ax2.tick_params(axis='both', which='major', length=5,width=2)
ax3.tick_params(axis='both', which='major', length=5,width=2)
ax4.tick_params(axis='both', which='major', length=5,width=2)
#plt.tight_layout()
#plt.legend((l1, l2,l3), (line13, line12,line18), loc='upper center',ncol=3)
fig.subplots_adjust(top=0.9,wspace=0.,hspace=0.)
ax2.legend(loc='upper right',ncol=1)
plt.savefig(datadir+'for_paper/COspectra.pdf')
plt.show()



