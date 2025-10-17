"""Picks representatives from binding site cluster @ filt/dat/clusters-by-bsite70.json and
stores them in filt/dat/reprs.json
stores binding sites with no resolution in filt/dat/missing_resolution.txt
"""

from collections import defaultdict
import json
from res import get_res
from bio import fwrite


INF = "filt/dat/f-clusters-by-bsite-70.json"
OUTF = "filt/dat/reprs.json"
NORESF = "missing_resolution.txt" # filt/dat/_

def choose_reprs():
	clusters = json.load(open(INF))
	print("number of clusters: ", len(clusters))

	reprs = []
	nores = []
	c= 0
	for cluster in clusters:
		# progress logging
		c+=1
		if not c%500:
			print(c)

		# if there is exactly one binding site in the cluster it must be chosen
		if len(cluster) == 1:
			reprs.append(cluster[0])
			continue

		resl = {}
		for bsite in cluster:
			if bsite in resl.keys():
				continue

			path = f"{bsite[4:6]}/{bsite[:7]}.ent"
			try:
				# get the resolution of the binding site's pdb and store it in the dictionary
				resl[bsite] = get_res(path)

			except ValueError as e:
				nores.append(bsite)
				# print(e)

		if len(resl):
			# pick the lowest resolution binding site and choose it as a representative
			repr = min(resl, key=resl.get)
			reprs.append(repr)
		else:
			# if the cluster has no sites with resolution, pick the first site in the cluster
			print(f"no res-based repr for cluster {c}")
			reprs.append(cluster[0])

	print("---------")
	print(len(reprs))
	print(len(nores))

	json.dump(reprs, open(OUTF, 'w'))
	fwrite(NORESF, nores, "binding sites with no resolution")


if __name__ == "__main__":
	choose_reprs()
