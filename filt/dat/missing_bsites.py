import json
import os
import time


with open("filt/dat/struct/bsites.txt") as f:
	lins = f.readlines()
	print(len(lins))
	bsites = {b.strip() for b in lins}
	print(len(bsites))

clust = json.load(open("filt/dat/clusters-by-bsite-70.json"))
print(len(clust))
clust = [site for c in clust for site in c]
print(len(clust))
clust = set(clust)
print(len(clust))

with open("filt/dat/bsites_missing_b.txt", 'w') as f:
	f.write(f"{(sites_m:=bsites - clust)}")

with open("filt/dat/bsites_missing_pdbs.txt", "w") as f:
	site_pdbs = set()
	for site in sites_m:
		site_pdbs.add(site[3:7])
	f.write(f"{site_pdbs}")

# ———————————
print('———————————')

with open("filt/dat/struct/pdbs.txt") as f:
	lins = f.readlines()
	print(len(lins))
	pdbs = {b.strip()[3:-4].upper() for b in lins}
	print(len(pdbs))

CLUSTF = "filt/dat/clusters-by-entity-70.txt"
with open(CLUSTF, "r") as f:
	lines = f.readlines()
	# ignore non-PDB entries
	clusters = [[enty.strip() for enty in clust.split(" ") if len(enty) <= 8] for clust in lines]
	clusters = [clust for clust in clusters if len(clust) > 0]
	print(len(clusters))

clusters = [enty[:4] for c in clusters for enty in c]
print(len(clusters), clusters[0])
clusters = set(clusters)
print(len(clusters))


with open("filt/dat/pdbs_missing_p.txt", 'w') as f:
	f.write(f"{(s2:=pdbs - clusters)}")

with open("filt/dat/pdbs_missing_c.txt", 'w') as f:
	f.write(f"{clusters - pdbs}")

with open('filt/dat/notfound.txt') as f:
	lins = f.readlines()[0].split(', ')
	print(len(lins))
	# time.sleep(0.7)
	# print(lins)

sites = []
for i, pdb in enumerate(pdbs - clusters):
	pdb = pdb.lower()
	mid = f"{pdb[1:3]}"
	dir = os.listdir(f"filt/dat/struct/binding-sites/{mid}")
	# print(dir)
	for f in dir:
		fsplit = f.split("_")
		# print(fsplit)
		# get the corresponding binding site
		if pdb == fsplit[0][3:]:
			sites.append(f)

print(len(sites))

with open("filt/dat/pdbs_missing_sites.txt", 'w') as f:
	f.write(f"{sites}")

print(len(s2-{k.upper() for k in site_pdbs}))
