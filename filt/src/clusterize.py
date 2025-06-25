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

	for clust in clusters:
		for enty in clust:
			pdb_id = f"pdb{enty[:4].lower()}.ent"
			enty = enty[5:]
			print('---')
			print(pdb_id, enty)
			get_lig(pdb_id, enty_to_bsites(pdb_id, enty))


def enty_to_bsites(pdb, enty):
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
				print(f"{pdb}: {chains.split(', ')}")
				break
	else:
		raise ValueError(f"{"Chains" if f else f"MOL_ID {enty}"} not found in {pdb}")

	return chains.split(', ')
	
def get_lig(pdb, chain): pass


if __name__ == "__main__":
	clusterize()
