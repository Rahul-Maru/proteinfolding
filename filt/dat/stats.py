"""Statistical analysis of the representative binding sites"""

import json
import os
from typing import Counter


INF = "filt/dat/reprs.json"

reprs = json.load(open(INF))
print(len(reprs))

# only 400 1-2 char ligands? some of them should've been removed
ligs = [b[15:-4] for b in reprs]
ligset = set(ligs)

counter = Counter(ligs).most_common()
print(len(counter), "unique ligands")
json.dump({c[0]: c[1] for c in counter}, open("filt/dat/ligcounts.json", 'w'))

with open("filt/dat/struct/unwanted-ligs.txt") as f:
	bad_ligs = [lig.strip() for lig in f.readlines()]

# print(bad_ligs)

interlopers = {l for l in ligset if l in bad_ligs}
print(interlopers)