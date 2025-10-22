"""Copies filt/dat/reprs.json into cluster_reprs.txt as a text list
with one representative per line"""

import json

IN = "../reprs.json"
OUT = "cluster_reprs.txt"

reprs = json.load(open(IN))
with open(OUT, 'w') as f:
	# needs newline at end of file or last entry gets discounted
	f.writelines([r + '\n' for r in reprs])

