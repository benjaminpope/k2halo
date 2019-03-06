# k2halo

##The K2 Halo Campaign Data Release and Paper

We are delighted to present light curves of the very brightest stars observed by K2 - including the first magnitude stars Aldebaran and Spica, the Seven Sisters of the Pleiades, the Hyades giants and the blue supergiant ρ Leonis! These light curves were obtained through a special 'halo photometry' method, which uses the unsaturated 'halo' around a saturated star and applies a Total Variation minimization (TV-min) algorithm to obtain nearly normal quality photometry. These are then additionally corrected with the [`k2sc`](https://github.com/OxES/k2sc) Gaussian Process systematics-correction code to further correct pointing residuals.

We hope you are as excited as we are to look at the largest ever sample of the brightest ever stars to be observed with high-quality space photometry! 

### Contents

In this repository we include the entire set of scripts and metadata that went into producing the K2 Halo Campaign data release.

* In notebooks/ are the Jupyter notebooks used to generate all the scripts and explore the data. 

* In data/ are the catalogues and scripts to download the TPFs.

* In reduced/ are the `halophot` reduction scripts, together with subdirectories for each campaiggn including the halo light curve fits files and diagnostic plots.

* In scripts/ are the `k2sc` call script and HPC sbatch scripts to deploy this to the cluster.

* In release/ are the final `k2sc` light curves in MAST High Level Science Product format, together with `k2sc` + `halophot` diagnostic plots. __Use these for your science!__

### Citation

The halo apertures were kindly provided by the K2 team as part of the Guest Observer programs GO19051, GO18051, GO17051, GO16003, GO15003, GO14003, GO13047, GO12047, GO11047, GO10025, GO9923, GO8025, GO7081, GO6081, and as a DDT program in Campaign 4. 

If you use these data pre-publication, please contact Pope, White, Huber to include them as co-authors.

A data release paper is in preparation.

Please also cite the original White et al. halo paper
	
	@ARTICLE{White2017,
	   author = {{White}, T.~R. and {Pope}, B.~J.~S. and {Antoci}, V. and {P{\'a}pics}, P.~I. and
	  {Aerts}, C. and {Gies}, D.~R. and {Gordon}, K. and {Huber}, D. and
	  {Schaefer}, G.~H. and {Aigrain}, S. and {Albrecht}, S. and {Barclay}, T. and
	  {Barentsen}, G. and {Beck}, P.~G. and {Bedding}, T.~R. and {Fredslund Andersen}, M. and
	  {Grundahl}, F. and {Howell}, S.~B. and {Ireland}, M.~J. and
	  {Murphy}, S.~J. and {Nielsen}, M.~B. and {Silva Aguirre}, V. and
	  {Tuthill}, P.~G.},
	    title = "{Beyond the Kepler/K2 bright limit: variability in the seven brightest members of the Pleiades}",
	  journal = {\mnras},
	archivePrefix = "arXiv",
	   eprint = {1708.07462},
	 primaryClass = "astro-ph.SR",
	 keywords = {asteroseismology, techniques: photometric, stars: early type, stars: variables: general, open clusters and associations: individual: Pleiades},
	     year = 2017,
	    month = nov,
	   volume = 471,
	    pages = {2882-2901},
	      doi = {10.1093/mnras/stx1050},
	   adsurl = {http://adsabs.harvard.edu/abs/2017MNRAS.471.2882W},
	  adsnote = {Provided by the SAO/NASA Astrophysics Data System}
	}