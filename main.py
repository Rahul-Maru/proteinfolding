import argparse
from consts import DEF_PROT, RAINBOW_OVR, SHOW_HETAMS_OVR
from Protein import Protein
from renderProtein import render


def main():
	# parses the command line to get the protein file and conditional arguments
	parser = argparse.ArgumentParser(description="Renders the atom in a protein")
	parser.add_argument('pdbfile', type=str, nargs='?', default=DEF_PROT,
					 help='the PDB file containing the protein data')

	parser.add_argument('--rainbow', action="store_true", dest="rainbow",
					 help="Using this flag colors the atoms in a rainbow")

	parser.add_argument('--hetatm', action='store_true', dest='show_hetatms',
					 help="Shows Heterogens")

	args = parser.parse_args()
	prot_file = args.pdbfile
	prot_path = f"proteins/{prot_file}.pdb"
	rainbow = args.rainbow or RAINBOW_OVR
	show_hetatms = args.show_hetatms or SHOW_HETAMS_OVR

	# load the protein
	if rainbow: p = Protein(prot_path, "rainbow")
	else: p = Protein(prot_path)

	render(p, [rainbow, show_hetatms])

if __name__ == "__main__":
	main()