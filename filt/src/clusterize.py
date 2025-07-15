from collections import defaultdict
import os
import re

CLUSTF = "filt/dat/clusters-by-entity-70.txt"
PDBDIR = "filt/dat/struct/pdb"

# CLUSTF = "filt/dat/clusters-by-entity-70.txt"
# PDBDIR = "filt/dat/struct/pdb"

def clusterize():
	"""
	map entities to binding sites
	for each cluster
	  go through each enty and add it to the lig-cluster
	"""

	# mirror: https://cdn.rcsb.org/resources/sequence/clusters/clusters-by-entity-70.txt
	with open(CLUSTF, "r") as f:
		lines = f.readlines()

	clusters = [[enty.strip() for enty in clust.split(" ") if len(enty) <= 8] for clust in lines]
	clusters = [clust for clust in clusters if len(clust) > 0]

	new_clusters = []
	elim = []
	for i, clust in enumerate(clusters):
		lig_clusters = defaultdict(list)
		for entry in clust:
			pdb_id = f"pdb{entry[:4].lower()}.ent"
			enty_n = entry[5:]

			print('---')
			try:
				chains = enty_to_chains(pdb_id, enty_n)

				print(f"{pdb_id} - entity {enty_n}")

				all_ligs = []
				for chain in chains:
					ligs = get_lig(pdb_id, chain)
					all_ligs.extend(ligs)
					print(f"{chain}: {ligs}")

				for lig in all_ligs:
					lig_clusters[lig].append(entry)

			except FileNotFoundError:
				# weeds out entries not available in the PDB format (or are otherwise absent from the dir) 
				elim.append(pdb_id)
				print(f"file not found: {pdb_id}")

		new_clusters.extend(lig_clusters.values())
		print(f"\n------done with cluster {i}------\n")

	with open('notfound', 'w'):
		f.write(', '.join(elim))

	return new_clusters

def enty_to_chains(pdb, enty):
	# match entity section in pdb file
	enty_pattern = re.compile(rf"COMPND\s+\d*\s*MOL_ID:\s*{enty};")
	# match chain lists under the entity section
	chain_pattern = re.compile(rf"COMPND\s+\d*\s*CHAIN:")

	with open(f"{PDBDIR}/{pdb[4:6]}/{pdb}", "r") as f:
		lines = f.readlines()

	f = False
	# go through the lines in the pdb file one by one
	for i, line in enumerate(lines):
		if enty_pattern.match(line):
			f = True
		if f:
			# once the right entity section is found, find the chain list and break
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
		# get the corresponding binding site
		if pdb == fsplit[0] and chain == fsplit[1]:
			print(f)
			# get the ligand name
			ligs.append(fsplit[3][:-4])

	return ligs

if __name__ == "__main__":
	clusterize()
