"""Script that investigates the differences between the molecular entites
of the pdb entries and the entity clusters."""

import json
import os

cluster_ents = set(json.load(open("filt/dat/clusters-by-entity-70.json")))
print("number of entities in cluster file: ", len(cluster_ents))

with open("filt/dat/struct/ents.txt") as f:
	pdb_ents = set(f.readlines())

pdb_ents = {ent.strip().upper() for ent in pdb_ents}
print("number of entites in the pdb: ", len(pdb_ents))

missing_entys = cluster_ents-pdb_ents
missing_entys2 = pdb_ents-cluster_ents
print("1. number of entities in cluster but not pdb: ", len(missing_entys))
print("2. number of entities in pdb but not cluster: ", len(missing_entys2))
print("3. number of entities in both pdb and cluster: ", len(pdb_ents&cluster_ents))

with open("filt/dat/missing_entys.txt", "w") as f_out:
	f_out.write(str(missing_entys))
	print("writing 1. to missing_entys.txt\n")

missing_entys_pdbs = {ent[:4] for ent in missing_entys}
missing_entys_pdbs2 = {ent[:4] for ent in missing_entys2}

INDIR = "filt/dat/struct/pdb"
missing_pdbs = set()
for pdb in missing_entys_pdbs:
    if not os.path.exists(f"{INDIR}/{pdb[1:3]}/{pdb}"):
        missing_pdbs.add(pdb)

print("1.: ", len(missing_entys))
print("number of pdbs represented by 1.: ", len(missing_entys_pdbs))
print("confirmation that all missing entities in 1. come from missing pdbs", len(missing_pdbs))

print("2.: ", len(missing_entys2))
print("number of pdbs represented by 2.: ", len(missing_entys_pdbs2))


