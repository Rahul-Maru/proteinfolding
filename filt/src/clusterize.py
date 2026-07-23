from collections import defaultdict
import json
import re
import time

# CLUSTF = "test_hsm/clusters-by-entity-70.txt"
# PDBDIR = "hsm/pdbs"

# NF_FILE = 'test_hsm/notfound'
# OUT_FILE = "test_hsm/clusters-by-bsite-70.json"

SITEF = "dat/struct/f_bsites.txt"

CLUSTF = "dat/clusters-by-entity-70.txt"
PDBDIR = "dat/struct/pdb"

NF_FILE = "dat/f-notfound.txt"
OUT_FILE = "dat/f-clusters-by-bsite-70.json"

bsites = defaultdict(list)
# pdbs = {}

# match entity section in pdb file, load for all possible entity numbers
ENTY_PATTERN = lambda enty : re.compile(rf"COMPND\s+\d*\s*MOL_ID:\s*{enty};")
# match chain lists under the entity section
CHAIN_PATTERN = re.compile(r"COMPND\s+\d*\s*CHAIN:")
# match the source section, implying that all relevant lines have been searched
SOURCE_PATTERN = re.compile(r"SOURCE")


def clusterize():
	"""
	map entities to binding sites
	for each cluster
	  go through each enty and add it to the lig-cluster
	"""

	st = time.time()

	load_bsites()

	# mirror: https://cdn.rcsb.org/resources/sequence/clusters/clusters-by-entity-70.txt
	with open(CLUSTF, "r") as f:
		lines = f.readlines()

	# ignore non-PDB entries by removing ones that are too long
	clusters = [
		[e for enty in clust.split() if len(e:=enty.strip()) <= 8] for clust in lines
	]
	clusters = [clust for clust in clusters if len(clust) > 0]

	new_clusters = []
	elim = set()

	st2 = time.time()

	for i, clust in enumerate(clusters):
		# current cluster broken into sub-clusters by ligand
		lig_clusters = defaultdict(list)

		for j, entry in enumerate(clust):
			nam = entry[:4].lower()
			pdb_f = f"pdb{nam}.ent"
			enty_n = int(entry[5:])
			try:
				# all the chains in this entity
				chains = enty_to_chains(pdb_f, enty_n)

				# put each site into the appropriate ligand-based cluster
				for chain in chains:
					sites = bsites[(nam, chain)]
					for lig, site in sites:
						lig_clusters[lig].append(site)

			except (FileNotFoundError, ValueError):
				# weeds out entries not available in the PDB format
				#   (or that are otherwise absent from the dir)
				elim.add(pdb_f)

			# if j % 200 == 0:
			#     print(f"{j}  {(t:=time.time() - st2)}s")

		new_clusters.extend(lig_clusters.values())

		if time.time() - st2 >= 15:
			print(f"\n------done up until cluster {i} in {time.time() - st2}s------\n")
			st2 = time.time()

	json.dump(new_clusters, open(OUT_FILE, "w"))
	with open(NF_FILE, "w") as f2:
		f2.write(", ".join(elim))

	print(f"DONE IN {time.time() - st}")

	return new_clusters


def load_bsites():
	t = time.time()
	with open(SITEF) as f:
		for line in f:
			line = line.strip()
			tokens = line.split("_")
			bsites[(tokens[0][3:], tokens[1])].append((tokens[3][:-4], line))

	print(f"loaded bsites in {time.time() - t}s")


def enty_to_chains(pdb, enty):
	"""Gets the chains from a given molecular entity"""

	path = f"{PDBDIR}/{pdb[4:6]}/{pdb}"
	try:
		with open(path, "r") as f:
			# go through the lines in the pdb file one by one
			enty_pattern = ENTY_PATTERN(enty)
			found = False
			for line in f:
				if not found and enty_pattern.match(line):
					found = True
					continue
				if found and CHAIN_PATTERN.match(line):
					# once the right entity section is found, find the chain list and return it
					l = line[18:].strip()
					chains = []
					if l[-1] == ';' or l[-1] == ',':
						chains += l[:-1].split(', ')
					else:
						chains += l.split(', ')
					if l[-1] == ',':
						for line2 in f:
							l2 = line2[11:].strip()
							if l2[-1] == ';' or l2[-1] == ',':
								chains += l2[:-1].split(', ')
							else:
								chains += l2.split(', ')
							if l2[-1] != ',':
								break
						else:
							print("NO BREAK -", pdb)

					return chains
			else:
				raise ValueError(
					f"{"Chains" if found else f"MOL_ID {enty}"} not found in {pdb}"
				)

	except FileNotFoundError:
		raise FileNotFoundError(f"file not found: {path}")

if __name__ == "__main__":
	clusterize()
