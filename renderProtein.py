"""Plots the atoms in a protein using VPython"""
import numpy as np
import vpython as vp
import argparse
from Protein import Protein, CHAIN_ID


#——DEFAULTS——
DEF_PROT = "1b41"
# Overrides. If set to true, the respective flag will be assumed to be present
SHOW_HETAMS_OVR = False # Shows Heterogens
RAINBOW_OVR = False # Displays the atoms with a rainbow color scheme

# radius of the atoms
RAD = 0.8
# opacity of HETATMs (WARNING: changing this may cause severe lag)
OPCTY = 1

# colors for depending on atom, chain, and rainbow mode
COLORS = {'C': [(0.5, 0.42, 0.42), (0.42, 0.5, 0.42), (0.42, 0.42, 0.5)],
		'N': [(0, 0, 1), (0, 0.52, 0.93), (0.15, 0, 0.78)],
		'O': [(1, 0, 0), (0.85, 0.2, 0.1), (0.75, 0, 0.2)],
		'S': [(1, 0.5, 0), (1, 0.58, 0), (1, 0.7, 0.1)],
		'HETATM': [(0.08, 0.7, 0.08), (1, 0.96, 0.85)],	#   v
		'CENTRD': [(1, 1, 1), (0.2, 0.21, 0.24)],		# 0 - Normal, 1 - Rainbow
		'HETCTD': [(0.2, 1, 0.65), (0.64, 0.6, 0.55)]}	#   ^


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
			# colors the atoms in a rainbow spectrum
			# calculates the relative position of the atom in its chain reduced to the range [0, 2π]
			rel_pos = 2*np.pi * (i - len(sum(p.chains[:chain], []))) / len(p.chains[chain])
			spr.color = vp.vector(*hue_to_RGB(rel_pos))
		else:
			# pick the appropriate color based on element and chain
			# if the chain number exceeds the number of available colors,
			# loop through the colors
			spr.color = vp.vector(*COLORS[elem][chain % 3])

	# display the centroid
	cent = vp.sphere()
	cent.pos = vp.vector(*tuple(p.centroid()))
	cent.radius = RAD*2
	cent.color = vp.vector(*COLORS['CENTRD'][int(rainbow)])
	# vp.scene.center = cent.pos

	if show_hetamts and len(p.hetatms) > 0:
		hx, hy, hz = p.get_xyz(-2)

		for hi, hj, hk in zip(hx, hy, hz):
			spr = vp.sphere()
			spr.pos = vp.vector(hi, hj, hk)
			spr.radius = RAD
			spr.color = vp.vector(*COLORS['HETATM'][int(rainbow)])
			spr.opacity = OPCTY

		# display the hetatm centroid
		hcent = vp.sphere()
		hcent.pos = vp.vector(*tuple(p.centroid()))

		hcent.radius = RAD*1.8
		hcent.color = vp.vector(*COLORS['HETCTD'][int(rainbow)])


def hue_to_RGB(θ: float) -> tuple:
	"""Given a hue value θ ∈ [0, 2π] and converts it to an RGB vector."""

	rcos = lambda x : min(1, max(0, np.cos(x) + 0.5))
	gcos = lambda x : min(1, max(0, np.cos(x - 2*np.pi/3) + 0.5))
	bcos = lambda x : min(1, max(0, np.cos(x + 2*np.pi/3) + 0.5))

	return (rcos(θ), gcos(θ), bcos(θ))

if __name__ == "__main__":
	main()
