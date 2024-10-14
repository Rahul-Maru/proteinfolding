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
					 help="applies a rainbow color scheme to the atoms")

	parser.add_argument('--hetatm', action='store_true', dest='show_hetatms',
					 help="shows heterogens")

	parser.add_argument('--bsite', action="store_true", dest="b_site",
					 help='stores only the binding site of the ligand')

	args = parser.parse_args()
	prot_file = args.pdbfile
	prot_path = f"proteins/{prot_file}.pdb"

	# if the OVR variables are set to true, set the flag to true
	rainbow = args.rainbow or RAINBOW_OVR
	show_hetatms = args.show_hetatms or SHOW_HETAMS_OVR
	b_site = args.b_site
	dp = (rainbow, show_hetatms, b_site)

	# load the protein
	p = Protein(prot_path, dp)

	render(p)

if __name__ == "__main__":
	main()