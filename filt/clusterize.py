def main():
	"""
	map entities to binding sites
	for each cluster
	  go through each enty and add it to the lig-cluster
	"""

	with open("clusters-by-entity-70.txt", "r") as f:
		lines = f.readlines()

	clusters = [enty.strip() for clust in lines for enty in clust.split(" ") if enty[:5] not in ["AF_AF"]]
	print(clusters)

	for clust in clusters:
		for enty in clust:
			get_lig(enty_to_bsites(enty))


def enty_to_bsites(*_): pass
def get_lig(*_): pass

if __name__ == "__main__":
	main()