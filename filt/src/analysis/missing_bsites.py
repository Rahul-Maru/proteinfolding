"""OBSOLETE (see f_missing2.py)
Script that investigates why some sites present in the binding sites directory are
not in the final clustered list of binding sites.
"""

import json
import os


with open("dat/struct/bsites.txt") as f:
	lins = f.readlines()
	print("––bsites.txt––")
	print(" - list of all binding sites")
	print("number of lines in bsites.txt: ", len(lins))
	bsites = {b.strip() for b in lins}
	print("# of lines w/o whitespace (= number of binding sites)", len(bsites))

print("\n––clusters-by-bsite-70.json––")
print("- clustered list of binding sites")
clust = json.load(open("dat/clusters-by-bsite-70.json"))
print("number of (ligand-based) clusters: ", len(clust))
clust = [site for c in clust for site in c]
print("number of clustered sites: ", len(clust))
clust = set(clust)
print("# sites w/o duplicates: ", len(clust))

with open("dat/bsites_missing_b.txt", 'w') as f:
	print("\nwrote to bsites_missing_b: binding sites that were not clustered")
	f.write(f"{(sites_m:=bsites - clust)}")

with open("dat/bsites_missing_pdbs.txt", "w") as f:
	site_pdbs = set()
	for site in sites_m:
		site_pdbs.add(site[3:7])

	print("wrote to bsites_missing_pdbs: parent pdbs of ^")
	f.write(f"{site_pdbs}")

# ———————————
print('———————————')

with open("dat/struct/pdb.txt") as f:
	print("––pdb.txt––")
	print(" - list of all PDB entries")
	lins = f.readlines()
	print("number of lines in pdb.txt: ", len(lins))
	pdbs = {b.strip()[3:-4].upper() for b in lins}
	print("# of lines w/o whitespace (number of PDBs)", len(pdbs))


print("\n––clusters-by-entity-70.txt––")
print("- original cluster-file")
print("- clustered list of molecular entities (non-PDB items removed)")

CLUSTF = "dat/clusters-by-entity-70.txt"
with open(CLUSTF, "r") as f:
	lines = f.readlines()
	# ignore non-PDB entries
	clusters = [[e for enty in clust.split() if len(e:=enty.strip()) <= 8] for clust in lines]
	clusters = [clust for clust in clusters if clust]
	print("number of clusters: ", len(clusters))

entys = [enty for c in clusters for enty in c]
print("number of entities in cluster-file (only PDB ID): ", len(entys))
print("sample item:", entys[0])
# TODO rename this file to something that isn't misleading
print("writing entys to dat/clusters-by-entity-70.json")
json.dump(entys, open("dat/clusters-by-entity-70.json", "w"))
pdbs_clust = {enty[:4] for enty in entys}
print("number of PDBs in cluster-file: ", len(pdbs_clust))


with open("dat/pdbs_missing_p.txt", 'w') as f:
	print("\nwriting to pdbs_missing_p.txt: pdbs absent from the clusterfile")
	f.write(f"{(s2:=pdbs - pdbs_clust)}")

with open("dat/pdbs_missing_c.txt", 'w') as f:
	print("writing to pdbs_missing_c.txt: pdbs in the clusterfile but not in the directory")
	f.write(f"{pdbs_clust - pdbs}")

with open('dat/notfound.txt') as f:
	print("\n––notfound.txt––")
	lins = f.readlines()[0].split(', ')
	print("number of PDBs not found while clustering: ", len(lins))
	# time.sleep(0.7)
	# print(lins)

sites = []
# s2 = pdbs absent from the clusterfile
for i, pdb in enumerate(s2):
	pdb = pdb.lower()
	dir = os.listdir(f"dat/struct/binding-sites/{pdb[1:3]}")

	# get all binding sites of this pdb
	for f in dir:
		fsplit = f.split("_")
		# if the binding site belongs to this PDB, add it
		if pdb == fsplit[0][3:]:
			sites.append(f)

print("number of binding sites whose parent PDBs are not in the cluster-file: ", len(sites))

with open("dat/pdbs_missing_sites.txt", 'w') as f:
	print("writing to pdbs_missing_sites.txt: ^ the above binding sites")
	f.write(f"{sites}")

SITE_PDBS = {k.upper() for k in site_pdbs}
print("\npdbs whose binding sites were not clustered because \
they are themselves not clustered: ", len(s2 & SITE_PDBS))
