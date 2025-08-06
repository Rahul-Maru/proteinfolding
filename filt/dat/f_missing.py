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
	bsites = {b.strip() for b in lins}
	print(len(bsites))

clust = json.load(open("filt/dat/f-clusters-by-bsite-70.json"))
print(len(clust))
clust_list = [site for c in clust for site in c]
print(len(clust_list))
clust_set = set(clust_list)
print(len(clust_set))

fwrite("bsites_missing_f.txt", (sites_m := bsites-clust_set))
print(len(sites_m))
print(len(sites_m) - 66968 - 1167)

# CHAIN_PATTERN = re.compile(rf"COMPND\s+\d*\s*CHAIN:")
# SOURCE_PATTERN = re.compile(r"SOURCE")

# removed = []
# nonentities = []

# c = 0
# for bsite in sites_m:
# 	c+=1
# 	if c%1000 == 0:
# 		print(c)
# 	tokens = bsite.strip().split("_")
# 	pdb_id = f"{tokens[0]}.ent"
# 	chains = []

# 	try:
# 		with open(f"filt/dat/struct/pdb/{pdb_id[4:6]}/{pdb_id}", "r") as f:
# 			for line in f:
# 				if SOURCE_PATTERN.match(line):
# 					break
# 				if CHAIN_PATTERN.match(line):
# 					chains.extend(line.replace(";", "").split("CHAIN:")[1].strip().split(', '))
# 	except FileNotFoundError as e:
# 		removed.append(bsite)
# 		continue

# 	if tokens[1] not in chains:
# 		nonentities.append(bsite)

# fwrite("nonentities.txt", nonentities)
# fwrite("removed.txt", removed)
# print(len(nonentities))
# print(len(removed))

with open("filt/dat/nonentities.txt") as f:
	nonentities = f.readlines()[0][2:-2].split("', '")
	print(len(nonentities))

with open("filt/dat/removed.txt") as f:
	removed = f.readlines()[0][2:-2].split("', '")
	print(len(removed))

with open("filt/dat/pdbs_missing_sites.txt") as f:
	unclustered_pdbs_sites = f.readlines()[0][2:-2].split("', '")
	print(len(unclustered_pdbs_sites))
