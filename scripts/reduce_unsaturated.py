import glob
from astropy.table import Table 
from tqdm import tqdm

files = glob.glob('epic_*.txt')

with open('all_epics.txt','w') as f:

	f.write('epic, Kp, cdpp_sap, cdpp_pdc, cdpp_k2sc_pdc, cdpp_halo, cdpp_k2sc_halo, enclosed_weight\n')

	for fname in tqdm(files):
		s = open(fname,'r').read().replace(' ',', ')
		f.write(s+'\n')

print('Done!')