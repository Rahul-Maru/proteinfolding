import json
import os


cluster_ents = set(json.load(open("filt/dat/clusters-by-entity-70.json")))

with open("filt/dat/struct/ents.txt") as f:
	pdb_ents = set(f.readlines())

pdb_ents = {ent.strip().upper() for ent in pdb_ents}


missing_entys = cluster_ents-pdb_ents
missing_entys2 = pdb_ents-cluster_ents
print(len(cluster_ents), len(pdb_ents), len(missing_entys), len(missing_entys2), len(pdb_ents&cluster_ents))

with open("filt/dat/missing_entys.txt", "w") as f_out:
	f_out.write(str(missing_entys))

missing_entys_pdbs = {ent[:4] for ent in missing_entys}
missing_entys_pdbs2 = {ent[:4] for ent in missing_entys2}

INDIR = "pdb"
missing_pdbs = set()
for ent in missing_entys_pdbs:
	if not os.path.exists(f"{INDIR}/{ent}"):
		missing_pdbs.add(ent)

print(len(missing_entys), len(missing_entys_pdbs), len(missing_pdbs))

print(len(missing_entys2), len(missing_entys_pdbs2))
