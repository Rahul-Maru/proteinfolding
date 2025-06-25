from collections import defaultdict
from clusterize import clusterize
from res import get_res

# binding site: 
# pdb####_C_****_LIG.pdb
# eg pdb1g29_1_0377_NA.pdb

# pdb file:
# pdb####.ent
# eg pdb1avn.ent

# #### = pdb id, C = chain, **** = residue number, LIG = ligand code



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