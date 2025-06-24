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
	print(clusters)

	for clust in clusters:
		for enty in clust:
			pdb_id = f"pdb{enty[:4].lower()}.ent"
			enty = enty[5:]
			get_lig(enty_to_bsites(pdb_id, enty))
			print(pdb_id, enty)


def enty_to_bsites(pdb, enty):
	with open(f"filt/dat/struct/pdb/{pdb[4:6]}/{pdb}", "r") as f:
		lines = f.readlines()

	i = lines.index(f"COMPND    MOL_ID: {enty}; ") + 2
	# TODO: can this wrap lines
	chain_info = lines[i][:-1].split(":")[1].split(', ')

def get_lig(*_): pass

if __name__ == "__main__":
	clusterize()
