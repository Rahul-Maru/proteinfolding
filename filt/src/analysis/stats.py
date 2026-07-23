"""Statistical analysis of the representative binding sites"""

import json
from collections import Counter
import matplotlib.pyplot as plt
import numpy as np


INF = "dat/reprs.json"

reprs = json.load(open(INF))
print(len(reprs))

# NOTE: only 400 1-2 char ligands?
ligs = [b[15:-4] for b in reprs]
ligset = set(ligs)

# ligand-frequency pairs
counter = Counter(ligs).most_common()
print(len(counter), "unique ligands")
json.dump({c[0]: c[1] for c in counter}, open("dat/ligcounts.json", 'w'))

with open("dat/struct/unwanted-ligs") as f:
	bad_ligs = [lig.strip() for lig in f.readlines()]

# print(bad_ligs)

# list of all unwanted ligands in the dataset
interlopers = {l for l in ligset if l in bad_ligs}
print(interlopers)

x_dat = np.array(range(1, len(counter) + 1))
y_dat = np.array([c[1] for c in counter])

plt.plot(x_dat, y_dat)
plt.xlabel("rank")
plt.ylabel("frequency")
plt.title("frequency of ligands")

plt.yscale('log')

plt.show()

