import numpy as np
import matplotlib 
import matplotlib as mpl

import lightkurve 
from lightkurve import KeplerLightCurve, KeplerTargetPixelFile
from k2sc.standalone import k2sc_lc
import k2sc

import halophot
from halophot.halo_tools import translate_greek

from astropy.table import Table
from astropy.io import fits

import fitsio

from argparse import ArgumentParser

import matplotlib as mpl

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
    plot_weightmap(ax_weightmap,weightmap,formal_name)
    plot_fluxmap(ax_fluxmap,image,formal_name)
    plot_pgram(ax_periodogram,frequency,power,spower,formal_name)        
    plot_log_pgram(ax_logpgram,frequency,power,spower,formal_name)  

    fig.suptitle(formal_name,y=0.99,fontsize=20)
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
    ap.add_argument('-name', default='test',type=str,help='Target name')
    ap.add_argument('-c', '--campaign', metavar='C',default=13, type=int, 
        help='Campaign number')
    ap.add_argument('--do-plot', action = 'store_true', default = True, \
                    help = 'produce plots')

    args = ap.parse_args()

    campaign = args.campaign
    ddir = '../reduced/c%d/' % campaign
    starname = args.name
    fname = ddir+'%s_halo_lc_o1.fits' % starname

    f = fitsio.FITS(fname)
    hdr = fitsio.read_header(fname)

    # read in our halo work
    all_stars = Table.read('../data/haloC%d.csv' % campaign,format='ascii')
    star = all_stars[all_stars['Name']==starname]
    epic = star['EPIC ID'].data.data[0]

    #load a lightkurve object with all the desired metadata
    tpf = lightkurve.open('../data/ktwo%d-c%02d_lpd-targ.fits.gz' % (epic,campaign))
    lc = tpf.to_lightcurve('aperture')
    lc.pos_corr1 = tpf.pos_corr1
    lc.pos_corr2 = tpf.pos_corr2
    lc.primary_header = tpf.hdu[0].header
    lc.data_header = tpf.hdu[1].header

    # match cadences
    inds_lc = match_cadences(lc.cadenceno,f[1]['cadence'][:])
    lc = lc[inds_lc]
    lc.pos_corr1 = lc.pos_corr1[inds_lc]
    lc.pos_corr2 = lc.pos_corr2[inds_lc]

    inds_flux = match_cadences(f[1]['cadence'][:],lc.cadenceno)
    lc.flux = f[1]['corr_flux'][inds_flux]

    # now the magic happens
    lc.__class__ = k2sc_lc
    lc.k2sc()

    # save data
    to_save = ['time', 'flux', 'flux_err','centroid_col', 'centroid_row', 'quality', 'cadenceno','pos_corr1', 'pos_corr2','tr_position', 'tr_time','corr_flux']

    dummy = fits.getheader('../data/ktwo%d-c%02d_lpd-targ.fits.gz' % (epic,campaign)) # get the old header from the TPF
    dummy['NAXIS']=1
    dummy['halo'] =(halophot.__version__,'halophot version')
    dummy['order']=(1,'halophot TV order')
    dummy['sub']=(1,'halophot subsampling')
    dummy['starname']=(starname,'Star Identifier')

    hdu = fits.PrimaryHDU(f[0][:,:],dummy) # can't save a masked array yet so just using pixelmap
    cols = [fits.Column(name=key,format="D",array=lc.__dict__[key]) for key in to_save]
    tab = fits.BinTableHDU.from_columns(cols)

    hdul = fits.HDUList([hdu, tab])
    hdul.writeto('hlsp_halo_k2_llc_%s_-c%d.fits' % (epic,campaign),overwrite=True)

    if args.do_plot:
        plot_k2sc(lc,np.nanmean(tpf.flux,axis=0),f[0][:,:].T,formal_name=translate_greek(args.name).replace('_',' ')+' Detrended')