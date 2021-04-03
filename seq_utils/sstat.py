#!/usr/bin/env python

#from Bio import SeqIO
import sys
import numpy as np
import pandas as pd
e=[float(x) for x in sys.stdin]

ee=pd.DataFrame.from_dict([{'len' : len(e), \
'min' : min(e), \
'max' : max(e), \
'median' : np.median(e), \
'avg' : np.average(e), \
'std' : np.std(e), \
'sum' : sum(e) \
}])

print(ee)
