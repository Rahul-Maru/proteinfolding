import json
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

fwrite("bsites_missing_f.txt", (sites_m:=bsites - clust_set))

CHAIN_PATTERN = re.compile(rf"COMPND\s+\d*\s*CHAIN:")


for bsite in sites_m:
	tokens = bsite.strip().split("_")
	pdb_id = f"{tokens[0]}.ent"
	with open(f"filt/dat/struct/pdb/{pdb_id[4:6]}/{pdb_id}", "r") as f:
		for line in f:
			if CHAIN_PATTERN.match(line):
				chains = line.replace(";", "").split("CHAIN:")[1].strip()
				break