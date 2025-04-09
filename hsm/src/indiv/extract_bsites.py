"""Extracts the binding sites of a given ligand from a PDB file
and stores them in separate files. """

from Protein import Protein
from consts import pdb_ids_all

# whether to consider all ligand molecules together or separately
combined = False

def extract_bsites():
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
			# TODO ligands in chains with no other atoms not recognized
			x = False
			if p in ['7pbd.pdb', '7pbz.pdb', '7pc0.pdb']:
				print(f"SUP, {p}")
				x = True
			for i in prot.chains:
				lig = prot.get_ligand('HSM', i)

				if len(lig["atoms"]) == 0:
					print("oops", p, i)
					continue
				else:
					print("huh", p, i)
				bsite = prot.get_bsite(lig)

				with open(f"hsm/bsites/{p[:-4]}_{i}.pdb", 'w') as f:
					f.writelines(sum([res['atoms'] for res in bsite], []))

if __name__ == "__main__":
	extract_bsites()
