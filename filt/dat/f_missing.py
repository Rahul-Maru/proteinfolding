import json
import os
import time

def fwrite(fp, x):
	with open(f"filt/dat/{fp}", 'w') as f:
		f.write(f"{x}")

with open("filt/dat/struct/f_bsites.txt") as f:
	lins = f.readlines()
	print(len(lins))
	bsites = {b.strip() for b in lins}
	print(len(bsites))

clust = json.load(open("filt/dat/f-clusters-by-bsite-70.json"))
print(len(clust))
clust_list = [site for c in clust for site in c]
print(len(clust_list))
clust_set = set(clust_list)
print(len(clust_set))

fwrite("bsites_missing_f.txt", (sites_m:=bsites - clust_set))

seen = set()
duplicates = set()
found = set()
for item in clust_list:
    if item in found:
        duplicates.add(item)
    else:
        found.add(item)   

print(len(duplicates))

fwrite("dupes_f.txt", duplicates)
