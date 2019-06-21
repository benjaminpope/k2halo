import numpy as np
import matplotlib.pyplot as plt

import matplotlib 
import matplotlib as mpl
from matplotlib import rc
from matplotlib.gridspec import GridSpec
from matplotlib.pyplot import figure, subplots, subplot
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size

import lightkurve 
from lightkurve import KeplerLightCurve, KeplerTargetPixelFile
from lightkurve.convenience import estimate_cdpp
from k2sc.standalone import k2sc_lc
import k2sc

import warnings ## eeek!
warnings.filterwarnings("ignore")

import halophot
from halophot.halo_tools import *

from astropy.table import Table
from astropy.io import fits

import fitsio

from argparse import ArgumentParser

mpl.style.use('seaborn-colorblind')

#To make sure we have always the same matplotlib settings
#(the ones in comments are the ipython notebook settings)

mpl.rcParams['figure.figsize']=(8.0,6.0)    #(6.0,4.0)
mpl.rcParams['font.size']=18               #10 
mpl.rcParams['savefig.dpi']= 200             #72 
mpl.rcParams['axes.labelsize'] = 16
mpl.rcParams['axes.labelsize'] = 16
mpl.rcParams['xtick.labelsize'] = 12
mpl.rcParams['ytick.labelsize'] = 12
mpl.rcParams["font.family"] = "Times New Roman"


'''-----------------------------------------------------------------
correct_halo.py

This executable Python script allows you to detrend any single object 
made using halophot with k2sc

-----------------------------------------------------------------'''

def match_cadences(halocads,lccads):
    indices =np.array([1 if j in lccads else 0 for j in halocads])
    return np.where(indices==1)[0]

def plot_weightmap_overlay(ax1,weightmap,name,title=False,mask=None):
        norm = np.size(weightmap)

        cmap = mpl.cm.seismic
        cmap.set_bad('k',1.)

        im = np.log10(weightmap.T*norm)
        pic = ax1.imshow(im.T,cmap=cmap, vmin=-2*np.nanmax(im),vmax=2*np.nanmax(im),
            interpolation='nearest',origin='lower')        
        if [mask] != None:
            mask_show = np.ma.masked_array(np.zeros_like(mask),mask=(mask))
            ax1.imshow(mask_show,cmap=mpl.cm.gray,alpha=0.3)

        if title:
            plt.title(r'TV-min Weightmap %s' % name)

        # cbaraxes, kw = mpl.colorbar.make_axes(ax1,location='right',pad=0.01)
        # plt.colorbar(pic,cax=cbaraxes)

        # cbaraxes.yaxis.set_ticks_position('right')
        aspect = 20
        pad_fraction = 0.5

        divider = make_axes_locatable(ax1)
        width = axes_size.AxesY(ax1, aspect=1./aspect)
        pad = axes_size.Fraction(pad_fraction, width)
        cax = divider.append_axes("right", size=width, pad=pad)
        plt.colorbar(pic, cax=cax)

        ax1.yaxis.set_ticks_position('left')        

def plot_lc(ax1,time,lc,name,trend=None,title=False):
        m = (lc>0.) * (np.isfinite(lc))

        ax1.plot(time[m],lc[m]/np.nanmedian(lc[m]),'.')
        dt = np.nanmedian(time[m][1:]-time[m][:-1])
        ax1.set_xlim(time[m].min()-dt,time[m].max()+dt)
        if trend is not None:
            ax1.plot(time[m],trend[m]/np.nanmedian(trend[m]),'-',color=colours[2])
            plt.legend(labels=['Flux','Trend'])
        ax1.set_xlabel('Time')
        ax1.set_ylabel('Relative Flux')
        if title:
            plt.title(r'%s' % name)

def plot_k2sc(lc,image,weightmap,save_file=None,formal_name='test',mask=None):
    min_p,max_p=1./24.,20.

    PW,PH = 8.27, 11.69
    
    frequency, power, spower = get_pgram(lc.time,lc.corr_flux-lc.tr_time,min_p=min_p,max_p=max_p)
    
    rc('axes', labelsize=7, titlesize=8)
    rc('font', size=6)
    rc('xtick', labelsize=7)
    rc('ytick', labelsize=7)
    rc('lines', linewidth=1)
    fig = plt.figure(figsize=(PW,PH))
    gs1 = GridSpec(3,2)
    gs1.update(top=0.95, bottom = 2/3.*1.05,hspace=0.0,left=0.09,right=0.96)
    gs2 = GridSpec(1,2)
    gs2.update(top=2/3.*0.97,bottom=1/3.*1.07,hspace=0.35,left=0.09,right=0.96)
    gs3 = GridSpec(2,2)
    gs3.update(top=1/3.*0.96,bottom=0.04,hspace=0.07,left=0.09,right=0.96)

    ax_lctime = subplot(gs1[0,:])
    ax_lcpos = subplot(gs1[1,:],sharex=ax_lctime)
    ax_lcwhite = subplot(gs1[2,:],sharex=ax_lctime)
    ax_fluxmap = subplot(gs2[0,0])
    ax_weightmap = subplot(gs2[0,1])
    ax_periodogram   = subplot(gs3[0,:])
    ax_logpgram    = subplot(gs3[1,:])

    plot_lc(ax_lctime,lc.time,lc.flux-lc.tr_time+np.nanmedian(lc.tr_time),formal_name,trend=lc.tr_position)
    plot_lc(ax_lcpos,lc.time,lc.flux-lc.tr_position+np.nanmedian(lc.tr_position),formal_name,trend=lc.tr_time)
    plot_lc(ax_lcwhite,lc.time,lc.corr_flux-lc.tr_time,formal_name+': Whitened')
    plot_weightmap_overlay(ax_weightmap,weightmap,formal_name,mask=mask)
    plot_fluxmap(ax_fluxmap,(image),formal_name)
    plot_pgram(ax_periodogram,frequency,power,spower,formal_name)        
    plot_log_pgram(ax_logpgram,frequency,power,spower,formal_name)  

    fig.suptitle(formal_name+' Detrended',y=0.99,fontsize=20)
    ax_periodogram.set_title('Periodograms')
    ax_fluxmap.set_title('Flux Map')
    ax_weightmap.set_title('TV-Min Weight Map')

    if save_file is not None:
        plt.savefig(save_file)

'''-----------------------------------------------------------------

An example call is 

python correct_halo.py -name Ascella -c 7 --do-plot
-----------------------------------------------------------------'''

if __name__ == '__main__':
    ap = ArgumentParser(description='halophot: K2 halo photometry with total variation.')
    ap.add_argument('epic',default=200000000,type=str,help='EPIC Number')
    ap.add_argument('--do-plot', action = 'store_true', default = True, \
                    help = 'produce plots')

    args = ap.parse_args()

    campaign = 6
    epic = args.epic

    tpf = lightkurve.open('../data/normal/ktwo%s-c06_lpd-targ.fits.gz' % (epic))
    lc = tpf.to_lightcurve('aperture')
    lc.pos_corr1 = tpf.pos_corr1
    lc.pos_corr2 = tpf.pos_corr2
    lc.primary_header = tpf.hdu[0].header
    lc.data_header = tpf.hdu[1].header


    lc_pipeline = lightkurve.open('../data/normal/ktwo%s-c06_llc.fits' % (epic))

    savedir='../release/c%d'

    cdpp_pdc = lc_pipeline.get_lightcurve('PDCSAP_FLUX').flatten().estimate_cdpp()
    cdpp_sap = lc_pipeline.get_lightcurve('SAP_FLUX').flatten().estimate_cdpp()
    print('CDPP: %.2f (PDC), %.2f (SAP)' % (cdpp_pdc,cdpp_sap))

    #load a lightkurve object with all the desired metadata

    tpf.__class__ = halophot.halo_tools.halo_tpf
    ws, lc_out = tpf.halo(aperture_mask=None,thresh=1.1,sigclip=True,objective='tv',lag=10,bitmask='hardest')
    weightmap = ws['weightmap']
    cdpp_halo = lc_out.flatten().estimate_cdpp()

    print('Halo CDPP:',cdpp_halo)

    enclosed_weight = np.sum(weightmap*tpf.pipeline_mask.T)
    print('enclosed_weight:',enclosed_weight)

    lc_out.__class__ = k2sc_lc
    lc_out.campaign = 6

    lc_out.k2sc()

    lc_dummy = lc_out.copy()
    lc_dummy.flux = lc_dummy.corr_flux

    cdpp_k2sc_halo = lc_dummy.flatten().estimate_cdpp()
    print('K2SC+Halo CDPP:',cdpp_k2sc_halo)

    # back to pdc
    print('Doing PDC')
    lc_pdc = lc_pipeline.get_lightcurve('PDCSAP_FLUX').remove_nans()
    lc_pdc.__class__ = k2sc_lc
    lc_pdc.campaign = 6
    lc_pdc.primary_header = tpf.hdu[0].header
    lc_pdc.data_header = tpf.hdu[1].header

    lc_pdc.k2sc()

    lc_dummy2 = lc_pdc.copy()
    lc_dummy2.flux = lc_dummy2.corr_flux

    cdpp_k2sc_pdc = lc_dummy2.flatten().estimate_cdpp()
    print('K2SC PDC CDPP:',cdpp_k2sc_pdc)

    with open('epic_%s.txt' % epic,'w') as f:
        f.write('%s %f %f %f %f %f %f %f' %  (epic, lc_pdc.primary_header['kepmag'], cdpp_sap, cdpp_pdc, cdpp_k2sc_pdc, cdpp_halo, cdpp_k2sc_halo, enclosed_weight))
    print('Written data to epic_%s.txt' % epic)

    to_save = ['time', 'flux', 'flux_err','centroid_col', 'centroid_row', 'quality', 'cadenceno','pos_corr1', 'pos_corr2','tr_position', 'tr_time','corr_flux']

    dummy = fits.getheader('../data/normal/ktwo%s-c06_lpd-targ.fits.gz' % (epic)) # get the old header from the TPF
    dummy['NAXIS']=1
    dummy['ORIGIN']=('Pope/NYU','institution responsible for creating this file ')
    dummy['CREATOR'] = ('halophot + k2halo/scripts/correct_unsaturated.py', 'pipeline')

    dummy['halo'] =(halophot.__version__,'halophot version')
    dummy['obj']=('tv','halophot objective')
    dummy['sub']=(1,'halophot subsampling')
    dummy['starname']=(epic,'Star Identifier')

    hdu = fits.PrimaryHDU(weightmap,dummy) # can't save a masked array yet so just using pixelmap
    cols = [fits.Column(name=key,format="D",array=lc_dummy.__dict__[key]) for key in to_save]
    tab = fits.BinTableHDU.from_columns(cols)

    hdul = fits.HDUList([hdu, tab])
    hdul.writeto('../data/normal/hlsp_halo_k2_llc_%s_-c%s_kepler_v1_lc.fits' % (epic,campaign),overwrite=True)
    print('Written Halo K2SC light curve to  ../data/normal/hlsp_halo_k2_llc_%s_-c%s_kepler_v1_lc.fits' % (epic,campaign))


    if args.do_plot:
        plot_k2sc(lc_out,np.nanmean(tpf.flux,axis=0),weightmap.T,formal_name='(EPIC %s)' % epic,mask=tpf.pipeline_mask,
            save_file='../data/normal/epic_%s.png' % (epic))
        print('Saved figure to ../data/normal/epic_%s.png' % (epic))





