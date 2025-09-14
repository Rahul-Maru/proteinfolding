import json
import os
import re
import time

def fwrite(fp, x, text):
	"""Easier way to write files"""
	with open(f"filt/dat/{fp}", 'w') as f:
		print(f"\nwriting to {fp}: {text}")
		f.write(f"{x}")

print("––f_bsites.txt––")
print(" - filtered list of all binding sites")

with open("filt/dat/struct/f_bsites.txt") as f:
	lins = f.readlines()
	print("number of lines in f_bsites.txt: ", len(lins))
	f_bsites = {b.strip() for b in lins}
	print("^ without whitespace (number of binding sites): ", len(f_bsites))

print("\n––clusters-by-bsite-70.json––")
print("- clustered list of filtered binding sites")
clust = json.load(open("filt/dat/f-clusters-by-bsite-70.json"))
print("number of ligand-clusters: ", len(clust))
clust_list = [site for c in clust for site in c]
print("number of clustered sites: ", len(clust_list))
clust_set = set(clust_list)
print("^ without duplicates: ", len(clust_set))

print() #-----------------

sites_m = f_bsites - clust_set

print("filtered binding sites that were not clustered: ", len(sites_m))
fwrite("bsites_missing_f.txt", sites_m, "^")

print("\n––nonentities.txt––")
print(" - list of binding sites whose chains do not belong to a named entity")
with open("filt/dat/nonentities.txt") as f:
	nonentities = set(f.readlines()[0][2:-2].split("', '"))
	print("number of binding sites not in an entity: ", len(nonentities))

print("\n––pdbs_missing_sites.txt––")
print(" - list of binding sites whose parent PDBs are not in the clusterfile")
with open("filt/dat/pdbs_missing_sites.txt") as f:
	unclustered_pdbs_sites = set(f.readlines()[0][2:-2].split("', '"))
	print("number of such binding sites: ", len(unclustered_pdbs_sites))

unclust_fsites = set()
for site in unclustered_pdbs_sites:
	if site in f_bsites:
		unclust_fsites.add(site)

print("number filtered binding sites with unclustered parent PDBs: ", len(unclust_fsites))
fwrite("unclustered_fsites.txt", unclust_fsites, "^")

print(f"{len(nonentities) + len(unclust_fsites)=}")
print("union of nonentities and unclustered f binding sites: ", len(nonentities|unclust_fsites))
print("intersection of nonentities and unclustered filtered binding sites: ", len(nonentities&unclust_fsites))
fwrite("nonentity_noncluster.txt", nonentities&unclust_fsites, "^")

fsites_m = sites_m - (nonentities|unclust_fsites)
print("number of missing binding sites except nonentities or those with unclustered parent PDBs: ", len(fsites_m))
fwrite("fbsites_missing.txt", fsites_m, "^ (fsites_m)")

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

print("number of entities represented by ^: ", len(enty_names))
enties, _, names = map(list, zip(*enty_names))
fwrite("entities_missing_names.txt", '\n'.join(names), "list of names of aforementioned entities")


CLUSTF = "filt/dat/clusters-by-entity-70.txt"
with open(CLUSTF, "r") as f:
	lines = f.readlines()
	# ignore non-PDB entries
	clusters = [[enty.strip() for enty in clust.split(" ") if len(enty) <= 8] for clust in lines]
	clusters = [clust for clust in clusters if len(clust) > 0]
	print("number of clusters: ", len(clusters))

clusters = [enty[:4] for c in clusters for enty in c]
print("number of entities in cluster-file (only PDB): ", len(clusters))
print("sample item:", clusters[0])
clusters = set(clusters)
print("number of PDBs in cluster-file: ", len(clusters))

print()

missing_enties = []
for ent in enties:
	if ent not in clusters:
		missing_enties.append(ent)

print("number of unaccounted-for entities not in the clusterfile: ", len(missing_enties))


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