import json
import os
import re
import time

def fwrite(fp, x):
	with open(f"filt/dat/{fp}", 'w') as f:
		f.write(f"{x}")


with open("filt/dat/struct/f_bsites.txt") as f:
	lins = f.readlines()
	print(len(lins))
	f_bsites = {b.strip() for b in lins}
	print(len(f_bsites))

clust = json.load(open("filt/dat/f-clusters-by-bsite-70.json"))
print(len(clust))
clust_list = [site for c in clust for site in c]
print(len(clust_list))
clust_set = set(clust_list)
print(len(clust_set))

fwrite("bsites_missing_f.txt", (sites_m := f_bsites-clust_set))
print(len(sites_m))
print(len(sites_m) - 66968 - 1167)

with open("filt/dat/nonentities.txt") as f:
	nonentities = set(f.readlines()[0][2:-2].split("', '"))
	print(len(nonentities))

with open("filt/dat/removed.txt") as f:
	removed = set(f.readlines()[0][2:-2].split("', '"))
	print(len(removed))

with open("filt/dat/pdbs_missing_sites.txt") as f:
	unclustered_pdbs_sites = set(f.readlines()[0][2:-2].split("', '"))
	print(len(unclustered_pdbs_sites))

c =0
unclust_fsites = set()
for i, site in enumerate(unclustered_pdbs_sites):
	if site in f_bsites:
		unclust_fsites.add(site)

print(len(unclust_fsites))

print(len(nonentities) + len(removed) + len(unclust_fsites))
print(len(nonentities|removed|unclust_fsites))
print(len(sites_m - (nonentities|removed|unclust_fsites)))

fwrite("fbsites_missing.txt", sites_m - (nonentities|removed|unclust_fsites))