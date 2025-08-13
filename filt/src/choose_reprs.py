from collections import defaultdict
import json
from res import get_res

# binding site: 
# pdb####_C_****_LIG.pdb
# eg pdb1g29_1_0377_NA.pdb

# pdb file:
# pdb####.ent
# eg pdb1avn.ent

# #### = pdb id, C = chain, **** = residue number, LIG = ligand code

def fwrite(fp, x):
	with open(f"filt/dat/{fp}", 'w') as f:
		f.write(f"{x}")

def main():
	clusters = json.load(open("filt/dat/f-clusters-by-bsite-70.json"))
	print(len(clusters))

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
	json.dump(reprs, open("filt/dat/reprs.json", 'w'))
	fwrite("missing_resolution.txt", nores)



if __name__ == "__main__":
	main()