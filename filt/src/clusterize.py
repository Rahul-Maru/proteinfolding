from collections import defaultdict
import os
import re

def clusterize():
	"""
	map entities to binding sites
	for each cluster
	  go through each enty and add it to the lig-cluster
	"""

	with open("filt/dat/clusters-by-entity-70.txt", "r") as f:
		lines = f.readlines()

	clusters = [[enty.strip() for enty in clust.split(" ") if len(enty) <= 8] for clust in lines]
	clusters = [clust for clust in clusters if len(clust) > 0]

	new_clusters = []
	for i, clust in enumerate(clusters):
		lig_clusters = defaultdict(list)
		for entry in clust:
			pdb_id = f"pdb{entry[:4].lower()}.ent"
			enty_n = entry[5:]
			chains = enty_to_chains(pdb_id, enty_n)

			print('---')
			print(f"{pdb_id} - entity {enty_n}")

			all_ligs = []
			for chain in chains:
				ligs = get_lig(pdb_id, chain)
				all_ligs.extend(ligs)
				print(f"{chain}: {ligs}")

			for lig in all_ligs:
				lig_clusters[lig].append(entry)

		new_clusters.extend(lig_clusters.values())
		print(f"done with cluster {i}")

	return new_clusters

def enty_to_chains(pdb, enty):
	enty_pattern = re.compile(rf"COMPND\s+\d*\s*MOL_ID:\s*{enty};")
	chain_pattern = re.compile(rf"COMPND\s+\d*\s*CHAIN:")

	with open(f"filt/dat/struct/pdb/{pdb[4:6]}/{pdb}", "r") as f:
		lines = f.readlines()

	f = False
	for i, line in enumerate(lines):
		if enty_pattern.match(line):
			f = True
		if f:
			if chain_pattern.search(line):
				chains = line.replace(';', '').split("CHAIN:")[1].strip()
				break
	else:
		raise ValueError(f"{"Chains" if f else f"MOL_ID {enty}"} not found in {pdb}")

	return chains.split(', ')

def get_lig(pdb, chain):
	# TODO filter out unwanted ligands
	pdb = pdb[:-4]
	mid = f"{pdb[4:6]}"
	dir = os.listdir(f"filt/dat/struct/binding-sites/{mid}")
	ligs = []
	for f in dir:
		fsplit = f.split("_")
		if pdb == fsplit[0] and chain == fsplit[1]:
			print(f)
			ligs.append(fsplit[3][:-4])

	return ligs

if __name__ == "__main__":
	clusterize()
