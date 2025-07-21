from collections import defaultdict
import json
import os
import re
import time

# CLUSTF = "test_hsm/clusters-by-entity-70.txt"
# PDBDIR = "hsm/pdbs"
# NF_FILE = 'test_hsm/notfound'
# OUT_FILE = "test_hsm/clusters-by-bsite-70.json"

CLUSTF = "filt/dat/clusters-by-entity-70.txt"
PDBDIR = "filt/dat/struct/pdb"
NF_FILE = 'filt/dat/notfound.txt'
OUT_FILE = "filt/dat/clusters-by-bsite-70.json"



def clusterize():
	"""
	map entities to binding sites
	for each cluster
	  go through each enty and add it to the lig-cluster
	"""

	# mirror: https://cdn.rcsb.org/resources/sequence/clusters/clusters-by-entity-70.txt
	with open(CLUSTF, "r") as f:
		lines = f.readlines()

	# ignore non-PDB entries
	clusters = [[enty.strip() for enty in clust.split(" ") if len(enty) <= 8] for clust in lines]
	clusters = [clust for clust in clusters if len(clust) > 0]

	new_clusters = []
	elim = set()

	# print(f"{len(clusters)} clusters")

	st = time.time()
	for i, clust in enumerate(clusters):
		print(f"beginning cluster {i}")
		# time.sleep(0.3)

		# current cluster broken into sub-clusters by ligand
		lig_clusters = defaultdict(list)

		st2 = time.time()
		for entry in clust:
			pdb_id = f"pdb{entry[:4].lower()}.ent"
			enty_n = entry[5:]

			# print('---')
			try:
				chains = enty_to_chains(pdb_id, enty_n)

				# print(f"{pdb_id} - entity {enty_n}")

				for chain in chains:
					sites = get_ligs_bsites(pdb_id, chain)
					for lig, site in sites:
						lig_clusters[lig].append(site)
					# print(f"{chain}: {sites}")

			except (FileNotFoundError, ValueError) as e:
				# weeds out entries not available in the PDB format (or are otherwise absent from the dir)
				elim.add(pdb_id)
				print(e)

		new_clusters.extend(lig_clusters.values())
		print(f"\n------done with cluster {i} in {time.time() - st2}s------\n")

	print(f"DONE IN {time.time() - st}")

	with open(NF_FILE, 'w') as f2:
		f2.write(', '.join(elim))

	json.dump(new_clusters, open(OUT_FILE, "w"))
	print(new_clusters)

	return new_clusters

def enty_to_chains(pdb, enty):
	# match entity section in pdb file
	enty_pattern = re.compile(rf"COMPND\s+\d*\s*MOL_ID:\s*{enty};")
	# match chain lists under the entity section
	chain_pattern = re.compile(rf"COMPND\s+\d*\s*CHAIN:")

	path = f"{PDBDIR}/{pdb[4:6]}/{pdb}"
	try:
		with open(path, "r") as f:
			lines = f.readlines()
	except FileNotFoundError:
		raise FileNotFoundError(f"file not found: {PDBDIR}/{pdb}")

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

def get_ligs_bsites(pdb, chain):
	# TODO filter out unwanted ligands
	pdb = pdb[:-4]
	mid = f"{pdb[4:6]}"
	dir = os.listdir(f"filt/dat/struct/binding-sites/{mid}")
	sites = []
	for f in dir:
		fsplit = f.split("_")
		# get the corresponding binding site
		if pdb == fsplit[0] and chain == fsplit[1]:
			# print("hi", f)
			# get the ligand name
			sites.append((fsplit[3][:-4], f))


	return sites

if __name__ == "__main__":
	clusterize()
