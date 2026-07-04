"""OBSOLETE (see f_missing2.py)
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
fwrite("bsites_missing_f.txt", sites_m)

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
fwrite("unclustered_fsites.txt", unclust_fsites)

print(f"{len(nonentities) + len(unclust_fsites)=}")
print("union of nonentities and unclustered f binding sites: ", len(nonentities|unclust_fsites))
print("intersection of nonentities and unclustered filtered binding sites: ", len(nonentities&unclust_fsites))
fwrite("nonentity_noncluster.txt", nonentities&unclust_fsites)

fsites_m = sites_m - (nonentities|unclust_fsites)
print("number of missing binding sites except nonentities or those with unclustered parent PDBs: ", len(fsites_m))
fwrite("fbsites_missing.txt", fsites_m, "^ (fsites_m)")

SOURCE_PATTERN = re.compile(r"SOURCE")
CHAIN_PATTERN = re.compile(r"COMPND\s+\d*\s*CHAIN:(.*)")
ENTY_PATTERN = re.compile(r"COMPND\s+\d*\s*MOL_ID:\s*(\d+)")
MOL_PATTERN = re.compile(r"COMPND\s+\d*\s*MOLECULE:\s*(.*)")

enty_names = []
for bsite in fsites_m:
	tokens = bsite.strip().split("_")
	pdb_id = f"{tokens[0]}.ent"
	enty_n = -1
	enty_nam = ""

	with open(f"filt/dat/struct/pdb/{pdb_id[4:6]}/{pdb_id}", "r") as f:
		for line in f:
			# if a SOURCE record is found then we have gone past the relevent section
			if SOURCE_PATTERN.match(line):
				print("FAILED - ", bsite)
				break

			if (ent := ENTY_PATTERN.match(line)):
				# we are now under entity x
				enty_n = int(ent.group(1))
			elif (mol := MOL_PATTERN.match(line)):
				# stores the name of the entity
				enty_nam = mol.group(1)
			elif (ch := CHAIN_PATTERN.match(line)):
				# gets the chains associated with the entity
				l = ch.group(1).strip()
				chains = []
				if l[-1] == ';' or l[-1] == ',':
					chains += l[:-1].split(', ')
				else:
					chains += l.split(', ')
				if l[-1] == ',':
					for line2 in f:
						l2 = line2[11:].strip()
						if l2[-1] == ';' or l2[-1] == ',':
							chains += l2[:-1].split(', ')
						else:
							chains += l2.split(', ')
						if l2[-1] != ',':
							break

					else:
						print("NO BREAK -", bsite)

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

print("\n––clusters-by-entity-70.txt––")
print("- original cluster file")
CLUSTF = "filt/dat/clusters-by-entity-70.txt"
with open(CLUSTF, "r") as f:
	lines = f.readlines()
	# ignore non-PDB entries
	clusters = [[e for enty in clust.split() if len(e:=enty.strip()) <= 8] for clust in lines]
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


removed = []
nonentities = []

c = 0
# finds all binding sites whose chains are not part of an entity
for bsite in sites_m:
	c+=1
	if c%10000 == 0:
		print(c)
	tokens = bsite.strip().split("_")
	pdb_id = f"{tokens[0]}.ent"
	chains = []

	try:
		with open(f"filt/dat/struct/pdb/{pdb_id[4:6]}/{pdb_id}", "r") as f:
			for line in f:
				# if a SOURCE record is found then we have gone past the relevent section
				if SOURCE_PATTERN.match(line):
					break

				if CHAIN_PATTERN.match(line):
					# find all the chains in the entity section
					l = line[18:].strip()
					if l[-1] == ';' or l[-1] == ',':
						chains += l[:-1].split(', ')
					else:
						chains += l.split(', ')
					if l[-1] == ',':
						for line2 in f:
							l2 = line2[11:].strip()
							if l2[-1] == ';' or l2[-1] == ',':
								chains += l2[:-1].split(', ')
							else:
								chains += l2.split(', ')
							if l2[-1] != ',':
								break

						else:
							print("NO BREAK -", bsite)
	except FileNotFoundError as e:
		removed.append(bsite)
		continue

	# check if the chain is not in the list
	if tokens[1] not in chains:
		nonentities.append(bsite)

fwrite("nonentities.txt", nonentities)
fwrite("removed.txt", removed)
print("number of non-entities: ", len(nonentities))
print("number of sites that should have been removed: ", len(removed))
