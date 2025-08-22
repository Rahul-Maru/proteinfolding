import json
import os
import re
import time

def fwrite(fp, x):
	"""Easier way to write files"""
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

print() #-----------------

fwrite("bsites_missing_f.txt", (sites_m := f_bsites-clust_set))
print(len(sites_m))

with open("filt/dat/nonentities.txt") as f:
	nonentities = set(f.readlines()[0][2:-2].split("', '"))
	print(len(nonentities))

with open("filt/dat/pdbs_missing_sites.txt") as f:
	unclustered_pdbs_sites = set(f.readlines()[0][2:-2].split("', '"))
	print(len(unclustered_pdbs_sites))

c =0
unclust_fsites = set()
for i, site in enumerate(unclustered_pdbs_sites):
	if site in f_bsites:
		unclust_fsites.add(site)

print(len(unclust_fsites))
fwrite("unclustered_fsites.txt", unclust_fsites)

print(len(nonentities) + len(unclust_fsites))
print(len(nonentities|unclust_fsites))
print(len(nonentities&unclust_fsites))
fwrite("nonentity_noncluster.txt", nonentities&unclust_fsites)
print(len(sites_m - (nonentities|unclust_fsites)))

fsites_m = sites_m - (nonentities|unclust_fsites)

fwrite("fbsites_missing.txt", fsites_m)

SOURCE_PATTERN = re.compile(r"SOURCE")
ENTY_PATTERN = re.compile(r"COMPND\s+\d*\s*MOL_ID:\s*(\d+)")
MOL_PATTERN = re.compile(r"COMPND\s+\d*\s*MOLECULE:\s*(.*)")
CHAIN_PATTERN = re.compile(r"COMPND\s+\d*\s*CHAIN:(.*)")

enty_names = []
for bsite in fsites_m:
	tokens = bsite.strip().split("_")
	pdb_id = f"{tokens[0]}.ent"
	enty_n = -1
	enty_nam = ""

	with open(f"filt/dat/struct/pdb/{pdb_id[4:6]}/{pdb_id}", "r") as f:
		for line in f:
			if SOURCE_PATTERN.match(line):
				print("FAILED - ", bsite, enty_n)
				break

			if (ent := ENTY_PATTERN.match(line)):
				enty_n = int(ent.group(1))
			elif (mol := MOL_PATTERN.match(line)):
				enty_nam = mol.group(1)
			elif (ch := CHAIN_PATTERN.match(line)):
				chains = ch.group(1).replace(';', '').strip().split(', ')
				if tokens[1] in chains:
					enty_names.append((f"{tokens[0][3:].upper()}_{enty_n}", bsite, enty_nam))
					break
		else:
			print("UNFOUND ", bsite)

# for nam in enty_names:
# 	print(f"\n---{nam[0]} –– {nam[1]}---")
# 	print(nam[2])

print(len(enty_names))
enties, _, names = map(list, zip(*enty_names))
fwrite("entities_missing_names.txt", '\n'.join(names))


CLUSTF = "filt/dat/clusters-by-entity-70.txt"
with open(CLUSTF, "r") as f:
	lines = f.readlines()
	# ignore non-PDB entries
	clusters = [[enty.strip() for enty in clust.split(" ") if len(enty) <= 8] for clust in lines]
	clusters = [clust for clust in clusters if len(clust) > 0]
	print(len(clusters))

clusters = [enty for c in clusters for enty in c]
print(len(clusters), clusters[0])
clusters = set(clusters)
print(len(clusters))

missing_enties = []
for ent in enties:
	if ent not in clusters:
		missing_enties.append(ent)

print(len(missing_enties))


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