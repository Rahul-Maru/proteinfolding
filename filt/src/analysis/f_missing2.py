""" SECOND VER (supersedes f_mising.py)
Script that investigates why some sites present in the filtered binding sites directory are
not in the final clustered list of binding sites.
"""

import json
import re
from bio import fwrite


print("––f_bsites.txt––")
print(" - filtered list of all binding sites")

with open("filt/dat/struct/f_bsites.txt") as f:
	lins = f.readlines()
	f_bsites = {b.strip() for b in lins}
	print("^ without whitespace (number of binding sites): ", len(f_bsites))

print("\n––clusters-by-bsite-70.json––")
print("- clustered list of filtered binding sites")
clust = json.load(open("filt/dat/f-clusters-by-bsite-70.json"))
clust_set = {site for c in clust for site in c}
print("number of clustered sites: ", len(clust_set))

print() #-----------------

sites_m = f_bsites - clust_set

print("filtered binding sites that were not clustered: ", len(sites_m))

print("\n––nonentities.txt––")
print(" - list of binding sites whose chains do not belong to a named entity")

# reads said list from file as computing it is time-consuming
with open("filt/dat/nonentities.txt") as f:
	nonentities = set(f.readlines()[0][2:-2].split("', '"))
	print("number of binding sites not in an entity: ", len(nonentities))

# unclustered binding sites whose chains ARE part of an entity
# there must exist a different reason as to why these sites were not clustered
sites_m2 = sites_m - nonentities
print("still-unaccounted for bsites: ", len(sites_m2))



SOURCE_PATTERN = re.compile(r"SOURCE")
CHAIN_PATTERN = re.compile(r"COMPND\s+\d*\s*CHAIN:(.*)")
ENTY_PATTERN = re.compile(r"COMPND\s+\d*\s*MOL_ID:\s*(\d+)")
MOL_PATTERN = re.compile(r"COMPND\s+\d*\s*MOLECULE:\s*(.*)")

enty_names = []
# gets all entities of the second set of missing binding sites
for bsite in sites_m2:
	tokens = bsite.strip().split("_")
	pdb_id = f"{tokens[0]}.ent"
	enty_n = -1
	enty_nam = ""

	with open(f"filt/dat/struct/pdb/{pdb_id[4:6]}/{pdb_id}", "r") as f:
		for line in f:
			if SOURCE_PATTERN.match(line):
				print("FAILED - ", bsite)
				break

			if (ent := ENTY_PATTERN.match(line)):
				enty_n = int(ent.group(1))
			elif (mol := MOL_PATTERN.match(line)):
				enty_nam = mol.group(1)
			elif (ch := CHAIN_PATTERN.match(line)):
				l = ch.group(1).strip()
				chains = []
				if l[-1] == ';':
					chains += l[:-1].split(', ')
				else:
					for line2 in f:
						l2 = line2[11:].strip()
						chains += l2[:-1].split(', ')
						if l2[-1] == ';':
							break

				if tokens[1] in chains:
					enty_names.append((f"{tokens[0][3:].upper()}_{enty_n}", bsite, enty_nam))
					break

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
