from collections import defaultdict
from filt import clusterize
from res import get_res


def main():
	clusters = clusterize()

	reprs = []
	for cluster in clusters:
		res = defaultdict([])
		for ent in cluster:
			res[ent] = get_res(ent)
		
		repr = min(res, res.get)
		reprs.append(repr)
	


if __name__ == "__main__":
	main()