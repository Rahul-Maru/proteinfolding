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

def main():
	clusters = json.load(open(INF))
	print("number of clusters: ", len(clusters))

	reprs = []
	nores = []
	c= 0
	for cluster in clusters:
		c+=1
		if not c%500:
			print(c)

		if len(cluster) == 1:
			reprs.append(cluster[0])
			continue

		resl = defaultdict(list)
		for bsite in cluster:
			path = f"{bsite[4:6]}/{bsite[:7]}.ent"
			try:
				resl[bsite] = get_res(path)

			except ValueError as e:
				nores.append(bsite)
				# print(e)

		if len(resl):	
			repr = min(resl, key=resl.get)
			reprs.append(repr)
		else:
			print(f"no res-based repr for cluster {c}")
			reprs.append(cluster[0])
		
	print("---------")
	print(len(reprs))
	print(len(nores))
	json.dump(reprs, open(OUTF, 'w'))
	fwrite(NORESF, nores, "binding sites with no resolution")



if __name__ == "__main__":
	main()