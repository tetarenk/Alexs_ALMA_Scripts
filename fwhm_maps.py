'''creates FWHM maps'''

from spectral_cube import SpectralCube
from astropy import units as u
import numpy as np
from astropy.io import fits
import math as ma
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from astropy import wcs
from astropy.coordinates import SkyCoord
import matplotlib.patches as patches
import matplotlib
import matplotlib.colors as colors
import pyspeckit
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes,inset_axes,InsetPosition
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import emcee
from astropy.modeling import models, fitting
from astropy.stats import mad_std
import scipy.stats as ss
import os

datadir = '/mnt/bigdata/tetarenk/ALMA_GRS1915_105/'

def error_estimate(x,data,vel_range):
    vel_low,vel_up=vel_range[0],vel_range[1]
    inds=np.where(np.logical_and(x>=vel_low,x<=vel_up))
    return(mad_std(data[inds]))
def oneDGauss(x,amp,mean,std):
    g = models.Gaussian1D(amplitude=amp, mean=mean, stddev=std)
    return(g(x))
def confidenceInterval(y,sig):
    median=np.median(y)
    pct15=np.percentile(y,15)
    pct85=np.percentile(y,85)
    list1=np.array([median,median-pct15,pct85-median])
    return list1
def lp(p,vel,data,error,guess):
    amp,cen,sig=p[0],p[1],p[2]
    mod=oneDGauss(vel,amp,cen,sig)
    if sig <0.:
        return(-np.inf)
    re=(mod-data)**2/(2*error**2)
    chi2_tot=np.nansum(re)
    return(-chi2_tot)
def chisq(p,vel,data,error):
    amp,cen,sig=p[0],p[1],p[2]
    mod=oneDGauss(vel,amp,cen,sig)
    chi2=(mod-data)**2/(2*error**2)
    return(np.nansum(chi2)/len(vel))
def fit_spec(X,Y,E,guess,nBurn,nSample):
	ndim=3
	nwalkers = ndim*2
	p0 = np.zeros((nwalkers,ndim))
	for i in np.arange(ndim):
		p0[:,i]=(((np.random.randn(nwalkers))*0.01)+guess[i])
	sampler = emcee.EnsembleSampler(nwalkers,ndim,lp,args=[X,Y,E,guess],threads=1)
	pos,prob,state = sampler.run_mcmc(p0,nBurn)
	sampler.reset()
	pos,prob,state = sampler.run_mcmc(pos,nSample)
	pars=[]
	for k in range(0,ndim):
		a=confidenceInterval(sampler.flatchain[:,k],1)
		pars.append(a[0])
	chi=chisq(pars,X,Y,E)
	fwhm=2*np.sqrt(2*np.log(2))*confidenceInterval(sampler.flatchain[:,2],1)[0]
	return(fwhm,chi)
def make_spec(fitsfile,ddir,offs,guess,line):
	ff=fits.open(fitsfile)[0]
	data=ff.data
	shapes=data.shape
	header = fits.getdata(fitsfile, header=True)[1]
	a=SpectralCube.read(fitsfile)
	a.allow_huge_operations=True
	a1=a.to(u.K, equivalencies=a.beam.jtok_equiv(a.header['RESTFRQ']*u.Hz))
	a2=a1.with_spectral_unit(u.km/u.s, velocity_convention='radio')
	fwhm_lst=np.empty((shapes[3],shapes[2]))
	chi_lst=np.empty((shapes[3],shapes[2]))
	for i in range(0,shapes[3]):
		for j in range(0,shapes[2]):
			if i>=216 and i<=451:
				if j>=401 and j<=626:
					#print i,j
					suba=a2[:, j, i]
					veloc=suba.spectral_axis
					amps=np.array(a2[:, j, i].value)
					Vels=np.array(veloc.to(u.km/u.s))
					e=np.empty(len(amps))
					e.fill(error_estimate(Vels,amps,offs))
					fwhm,chi=fit_spec(Vels,amps,e,[np.nanmax(amps),guess[1],guess[2]],50,500)
					fwhm_lst[j][i]=fwhm
					chi_lst[j][i]=chi
					#print fwhm,chi
				else:
					fwhm_lst[j][i]=np.nan
					chi_lst[j][i]=np.nan
			else:
				fwhm_lst[j][i]=np.nan
				chi_lst[j][i]=np.nan
	os.system('rm -rf '+ddir+'spectra_plots/'+line+'_chi.fits')
	fits.writeto(filename=ddir+'spectra_plots/'+line+'_chi.fits',output_verify='ignore',\
	clobber=True,data=chi_lst,header=header)
	os.system('rm -rf '+ddir+'spectra_plots/'+line+'_fwhm.fits')
	fits.writeto(filename=ddir+'spectra_plots/'+line+'_fwhm.fits',output_verify='ignore',\
	clobber=True,data=fwhm_lst,header=header)


#216,451,401,626--CO limits of cont emission
line='18CO'
fitsfile=datadir+'alex_imaging_'+line+'_fix/GRS1915_modelimg_'+line+'.image.pbcor.fits'
offs=[73,85]
guess=[10.,67.5,2.]


make_spec(fitsfile,datadir,offs,guess,line)	
			
