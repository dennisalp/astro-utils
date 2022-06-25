import sys
import numpy as np
import pdb
import pandas as pd

if len(sys.argv) == 3:
    aa = pd.read_hdf(sys.argv[1])
    bb = pd.read_hdf(sys.argv[2])
else:
    aa = pd.read_hdf(sys.argv[1], sys.argv[3])
    bb = pd.read_hdf(sys.argv[2], sys.argv[3])
    
for ii, key in enumerate(aa.keys()):
    if np.any((aa[key] == bb[key]) | (np.isnan(aa[key]) & np.isnan(bb[key]))):
        continue
    pdb.set_trace()
    print('Difference.')
    sys.exit(1)
print('No difference.')
