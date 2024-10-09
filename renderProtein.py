"""Plots the atoms in a protein using VPython"""
import numpy as np
import vpython as vp
import argparse
from Protein import Protein, CHAIN_ID


#——DEFAULTS——
DEF_PROT = "5fil"
# Overrides. If set to true, the respective flag will be assumed to be present
SHOW_HETAMS_OVR = False # Shows Heterogens
RAINBOW_OVR = False # Displays the atoms with a rainbow color scheme

# radius of the atoms
RAD = 0.8
# opacity of HETATMs (WARNING: changing this may cause severe lag)
OPCTY = 1

# colors for depending on atom, chain, and rainbow mode
COLORS = {'C': [vp.vector(0.5, 0.42, 0.42), vp.vector(0.42, 0.5, 0.42), vp.vector(0.42, 0.42, 0.5)],
		'N': [vp.vector(0, 0, 1), vp.vector(0, 0.52, 0.93), vp.vector(0.13, 0, 0.75)],
		'O': [vp.vector(1, 0, 0), vp.vector(0.85, 0.2, 0.1), vp.vector(0.75, 0, 0.2)],
		'S': [vp.vector(1, 0.5, 0), vp.vector(1, 0.58, 0), vp.vector(1, 0.7, 0.1)],
		'HETATM': [vp.vector(0.08, 0.69, 0.08), vp.vector(1, 0.96, 0.85)]} # 0 - Normal, 1 - Rainbow


def main():
	# parses the command line to get the protein file and conditional arguments
	parser = argparse.ArgumentParser(description="Renders the atom in a protein")
	parser.add_argument('pdbfile', type=str, nargs='?', default=DEF_PROT,
					 help='the PDB file containing the protein data')

	parser.add_argument('--rainbow', action="store_true", dest="rainbow",
					 help="Using this flag colors the atoms in a rainbow")

	parser.add_argument('--show_hetatms', action='store_true', dest='show_hetatms',
					 help="Shows Heterogens")

	args = parser.parse_args()
	prot_file = args.pdbfile
	prot_path = f"proteins/{prot_file}.pdb"
	rainbow = args.rainbow or RAINBOW_OVR
	show_hetamts = args.show_hetatms or SHOW_HETAMS_OVR

	# load the protein
	p = Protein(prot_path)

	print(f"Rendering protein at {prot_path}")
	p.info()


	for i, atom in enumerate(p.atoms):
		spr = vp.sphere()
		spr.pos = vp.vector(p.x[i], p.y[i], p.z[i])
		spr.radius = RAD

		elem = atom[77] # get the element of the atom
		chain = CHAIN_ID(atom)

		if rainbow:
			# calculates the relative position of the atom in its chain reduced to the range [0, 2π]
			rb_pos = 2*np.pi * (i - sum(p.chainlens[:chain])) / p.chainlens[chain]
			spr.color = hue_to_RGB(rb_pos)

		else:
			# pick the appropriate color based on element and chain
			if chain <= 2:
				spr.color = COLORS[elem][chain]
			else:
				# if the chain number exceeds the number of available colors,
				# loop through the colors
				spr.color = COLORS[elem][chain % 3]

		# get the appropriate color based on element and side chain

	if show_hetamts:
		for i, j, k in zip(p.hx, p.hy, p.hz):
			spr = vp.sphere()
			spr.radius = RAD
			spr.pos = vp.vector(i, j, k)
			spr.color = COLORS['HETATM'][int(rainbow)]
			spr.opacity = OPCTY

def hue_to_RGB(θ: float) -> vp.vector:
	"""Given a hue value θ ∈ [0, 2π] and converts it to an RGB vector """
	rcos = lambda x : min(1, max(0, np.cos(x) + 0.5))
	gcos = lambda x : min(1, max(0, np.cos(x - 2*np.pi/3) + 0.5))
	bcos = lambda x : min(1, max(0, np.cos(x + 2*np.pi/3) + 0.5))

	return vp.vector(rcos(θ), gcos(θ), bcos(θ))

if __name__ == "__main__":
	main()
