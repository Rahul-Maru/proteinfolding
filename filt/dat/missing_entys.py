import json


cluster_ents = set(json.load(open("filt/dat/clusters-by-entity-70.json")))

with open("filt/dat/struct/ents.txt", "r") as f:
	pdb_ents = set(f.readlines())

pdb_ents = {ent.strip() for ent in pdb_ents}


missing_entys = cluster_ents-pdb_ents

with open("filt/dat/missing_entys.txt", "w") as f_out:
	f_out.write(str(missing_entys))
