"""Extracts the binding sites of a given ligand from a PDB file
and stores them in separate files. """

from Protein import Protein
from consts import pdb_ids_all, CHAIN_LET

# whether to consider all ligand molecules together or separately
combined = True

def main():
	# file path template
	path = "hsm/pdbs/"

	for p in pdb_ids_all:
		prot = Protein(path + p)
		if combined:
			lig = prot.get_ligand('HSM')
			bsite = prot.get_bsite(lig)

			with open(f"hsm/bsites_combined/{p[:-4]}.pdb", 'w') as f:
					f.writelines(sum([res['atoms'] for res in bsite], []))

		else:
			# create a separate file for the binding site of each chain
			for i in prot.chains:
				lig = prot.get_ligand('HSM', i)
				if len(lig["atoms"]) == 0:
					continue
				bsite = prot.get_bsite(lig)

				with open(f"hsm/bsites/{p[:-4]}_{i}.pdb", 'w') as f:
					f.writelines(sum([res['atoms'] for res in bsite], []))

if __name__ == "__main__":
	main()
