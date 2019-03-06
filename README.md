# k2halo

The K2 Halo Campaign Data Release and Paper

### Contents

In this repository we include the entire set of scripts and metadata that went into producing the K2 Halo Campaign data release.

* In notebooks/ are the Jupyter notebooks used to generate all the scripts and explore the data. 

* In data/ are the catalogues and scripts to download the TPFs.

* In reduced/ are the `halophot` reduction scripts, together with subdirectories for each campaiggn including the halo light curve fits files and diagnostic plots.

* In scripts/ are the `k2sc` call script and HPC sbatch scripts to deploy this to the cluster.

* In release/ are the final `k2sc` light curves in MAST High Level Science Product format, together with `k2sc` + `halophot` diagnostic plots.