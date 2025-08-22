import json
import os
import time


with open("filt/dat/struct/bsites.txt") as f:
	lins = f.readlines()
	print("––bsites.txt––")
	print(" - list of all binding sites")
	print("number of lines in bsites.txt: ", len(lins))
	bsites = {b.strip() for b in lins}
	print("^ without whitespace (number of binding sites)", len(bsites))

print("\n––clusters-by-bsite-70.json––")
print("- clustered list of binding sites")
clust = json.load(open("filt/dat/clusters-by-bsite-70.json"))
print("number of clusters: ", len(clust))
clust = [site for c in clust for site in c]
print("number of clustered sites: ", len(clust))
clust = set(clust)
print("^ without duplicates: ", len(clust))

with open("filt/dat/bsites_missing_b.txt", 'w') as f:
	print("\nwrote to bsites_missing_b: binding sites that were not clustered")
	f.write(f"{(sites_m:=bsites - clust)}")

with open("filt/dat/bsites_missing_pdbs.txt", "w") as f:
	site_pdbs = set()
	for site in sites_m:
		site_pdbs.add(site[3:7])
	
	print("wrote to bsites_missing_pdbs: parent pdbs of ^")
	f.write(f"{site_pdbs}")

# ———————————
print('———————————')

with open("filt/dat/struct/pdbs.txt") as f:
	print("––pdbs.txt––")
	print(" - list of all pdb entries")
	lins = f.readlines()
	print("number of lines in pdbs.txt: ", len(lins))
	pdbs = {b.strip()[3:-4].upper() for b in lins}
	print("^ without whitespace (number of pdbs)", len(pdbs))


print("\n––clusters-by-entity-70.json––")
print("- clustered list of molecular entities (non-PDB items removed)")

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
print("^ without duplicates: ", len(clusters))


with open("filt/dat/pdbs_missing_p.txt", 'w') as f:
	print("\nwriting to pdbs_missing_p.txt: pdbs absent from the clusterfile")
	f.write(f"{(s2:=pdbs - clusters)}")

with open("filt/dat/pdbs_missing_c.txt", 'w') as f:
	print("writing to pdbs_missing_c.txt: pdbs in the clusterfile but not in the directory")
	f.write(f"{clusters - pdbs}")

with open('filt/dat/notfound.txt') as f:
	print("\n––notfound.txt––")
	lins = f.readlines()[0].split(', ')
	print("number of PDBs not found while clustering: ", len(lins))
	# time.sleep(0.7)
	# print(lins)

sites = []
# pdbs absent from the clusterfile
for i, pdb in enumerate(s2):
	pdb = pdb.lower()
	dir = os.listdir(f"filt/dat/struct/binding-sites/{pdb[1:3]}")

	# get all binding sites of the pdb
	for f in dir:
		fsplit = f.split("_")
		# if the binding site belongs to this PDB, add it
		if pdb == fsplit[0][3:]:
			sites.append(f)

print("number of binding sites whose parent PDBs are not in the clusterfile: ", len(sites))

with open("filt/dat/pdbs_missing_sites.txt", 'w') as f:
	print("writing to pdbs_missing_sites.txt: ^ the above binding sites")
	f.write(f"{sites}")

SITE_PDBS = {k.upper() for k in site_pdbs}
print("\npdbs whose binding sites were not clustered because they are themselves not clustered: ", len(s2 & SITE_PDBS))
