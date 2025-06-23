from collections import defaultdict
from res import get_res


def main():
	format_pdb = lambda x : x

	with open("clusters-by-entity-70.txt", "r") as f:
		lines = f.readlines()

	clusters = [format_pdb(ent) for clust in lines for ent in clust.split(" ") if ent[:5] != "AF_AF"]

	for cluster in clusters:
		res = defaultdict([])
		for ent in cluster:
			res[ent] = get_res(ent)
		
		repr = min(res, res.get)
	


if __name__ == "__main__":
	main()