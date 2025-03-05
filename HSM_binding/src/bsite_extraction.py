"""Extracts the binding sites of a given ligand from a PDB file
and stores them in separate files. """

from Protein import Protein
from consts import pdb_ids_all, CHAIN_LET

def main():
	# file path template
	path = "HSM/pdbs/$"

	for p in pdb_ids_all:
		prot = Protein(path.replace('$', p))
		# create a separate file for the binding site of each chain
		for i in range(len(prot.chains)):
			lig = prot.get_ligand('HSM', i)
			if len(lig["atoms"]) == 0:
				continue
			bsite = prot.get_bsite(lig)

			with open(f"HSM/MAPP-3D/MultipleSiteAlignment/HSM/{p[:-4]}_{CHAIN_LET(i)}.pdb", 'w') as f:
				f.writelines(sum([res['atoms'] for res in bsite], []))

if __name__ == "__main__":
	main()
