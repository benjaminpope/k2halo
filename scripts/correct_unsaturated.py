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

def plot_k2sc(lc,image,weightmap,save_file=None,formal_name='test'):
    min_p,max_p=1./24.,20.

    PW,PH = 8.27, 11.69
    m = np.isfinite(lc.corr_flux)
    poly = np.poly1d(np.polyfit(lc.time[m],lc.corr_flux[m],15))
    trend = poly(lc.time)
    frequency, power, spower = get_pgram(lc.time,lc.corr_flux-trend+np.nanmedian(trend),min_p=min_p,max_p=max_p)
    font = {'fontname':'Times New Roman'}

    rc('axes', labelsize=7, titlesize=8)
    rc('font', size=6, family="Times New Roman")
    rc('xtick', labelsize=7)
    rc('ytick', labelsize=7)
    rc('lines', linewidth=1)
    fig = plt.figure(figsize=(PW,PH))
    gs1 = GridSpec(3,2)
    gs1.update(top=0.95, bottom = 2/3.*1.05,hspace=0.0,left=0.09,right=0.96)
    gs2 = GridSpec(1,2)
    gs2.update(top=2/3.*1.01,bottom=1/3.*1.1,hspace=0.35,left=0.09,right=0.96)
    gs3 = GridSpec(2,2)
    gs3.update(top=1/3.*1.01,bottom=0.04,hspace=0.23,left=0.09,right=0.96)

    ax_lctime = subplot(gs1[0,:])
    ax_lcpos = subplot(gs1[1,:],sharex=ax_lctime)
    ax_lcwhite = subplot(gs1[2,:],sharex=ax_lctime)
    ax_fluxmap = subplot(gs2[0,0])
    ax_weightmap = subplot(gs2[0,1])
    ax_periodogram   = subplot(gs3[0,:])
    ax_logpgram    = subplot(gs3[1,:])

    plot_lc(ax_lctime,lc.time,lc.flux-lc.tr_time+np.nanmedian(lc.tr_time),formal_name,trends=[lc.tr_position])
    plot_lc(ax_lcpos,lc.time,lc.flux-lc.tr_position+np.nanmedian(lc.tr_position),formal_name,trends=[lc.tr_time,trend])
    plot_lc(ax_lcwhite,lc.time,(lc.corr_flux-lc.tr_time)+np.nanmedian(lc.tr_time),formal_name+': Whitened')
    plot_weightmap(ax_weightmap,weightmap,formal_name)
    plot_fluxmap(ax_fluxmap,image,formal_name)
    plot_pgram(ax_periodogram,frequency,power,spower,formal_name)        
    plot_log_pgram(ax_logpgram,frequency,power,spower,formal_name)  

    fig.suptitle(formal_name,y=0.99,fontsize=20,**font)
    ax_periodogram.set_title('Periodograms')
    ax_fluxmap.set_title('Flux Map')
    ax_weightmap.set_title('TV-Min Weight Map')

    if save_file is not None:
        try:
            for fname in save_file:
                plt.savefig(fname)
        except:
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

    # read in our halo work


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
        f.write('%s %f %f %f %f %f %f' %  (epic, cdpp_sap, cdpp_pdc, cdpp_k2sc_pdc, cdpp_halo, cdpp_k2sc_halo, enclosed_weight))
    print('Written data to epic_%s.txt' % epic)

    if args.do_plot:
        plot_k2sc(lc_out,np.nanmean(tpf.flux,axis=0),weightmap.T,formal_name='(EPIC %s) Detrended' % epic,
            save_file=['../data/normal/epic_%s.png' % (epic)])
        print('Saved figure to ../data/normal/epic_%s.png' % (epic))





