"""Plots the atoms in a protein using VPython"""
import numpy as np
import vpython as vp
import argparse

#——UTILS——
# convert the chain letter of an atom to a numerical index
CHAIN_ID = lambda a : ord(a[21]) - 65

#——DEFAULTS——
DEF_PROT = "5fil"
# Overrides. If set to true, the respective flag will be assumed to be present
SHOW_HETAMS_OVR = False # Shows Heterogens
RAINBOW_OVR = False # Displays the atoms with a rainbow color scheme

# radius of the atoms
RAD = 0.8
# opacity of HETATMs (WARNING: changing this may cause severe lag)
OPCTY = 1

# colors for depending on atom and chain
COLORS = {'C': [vp.vector(0.5, 0.42, 0.42), vp.vector(0.42, 0.5, 0.42), vp.vector(0.42, 0.42, 0.5)],
		'N': [vp.vector(0, 0, 1), vp.vector(0, 0.52, 0.93), vp.vector(0.13, 0, 0.75)],
		'O': [vp.vector(1, 0, 0), vp.vector(0.85, 0.2, 0.1), vp.vector(0.75, 0, 0.2)],
		'S': [vp.vector(1, 0.5, 0), vp.vector(1, 0.58, 0), vp.vector(1, 0.7, 0.1)],
		'HETATM': [vp.vector(0.08, 0.69, 0.08), vp.vector(1, 0.96, 0.85)]} # 0 - Normal, 1 - Rainbow


def main():
	# parses the command line to get the protein
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

	with open(f"{prot_path}") as file:
		lines = file.readlines()

	# extract only the ATOM and HETATM lines from the database and their coords
	atoms = list(filter(lambda l : l[:4] == "ATOM", lines))
	hetatms = list(filter(lambda l : l[:6] == "HETATM", lines))
	x, y, z = [[float(atom[i:i+8]) for atom in atoms] for i in range(30, 54, 8)]
	hx, hy, hz = [[float(hetatm[i:i+8]) for hetatm in hetatms] for i in range(30, 54, 8)]
	# sort atoms by element
	c, o, n, s = [list(filter(lambda x: x[77] == e, atoms)) for e in ['C', 'O', 'N', 'S']]

	chainlens = []
	for atom in atoms:
		try:
			chainlens[CHAIN_ID(atom)] += 1
		except:
			chainlens.append(1)

	print(f"Rendering protein at {prot_path}")
	print(f"{len(chainlens)} Chains,")
	print(len(atoms), " Atoms:")
	print(len(c), " Carbon (grey),")
	print(len(n), " Nitrogen (blue),")
	print(len(o), " Oxygen (red),")
	print(len(s), " Sulfur (orange),")
	print(len(hetatms), " Hetero-atoms (magenta)")

	rcos = lambda x : min(1, max(0, np.cos(x) + 0.5))
	gcos = lambda x : min(1, max(0, np.cos(x - 2*np.pi/3) + 0.5))
	bcos = lambda x : min(1, max(0, np.cos(x + 2*np.pi/3) + 0.5))

	for i, atom in enumerate(atoms):
		spr = vp.sphere()
		spr.pos = vp.vector(x[i], y[i], z[i])
		spr.radius = RAD

		elem = atom[77] # get the element of the atom
		chain = CHAIN_ID(atom)

		if rainbow:
			rb_pos = 2*np.pi * (i - sum(chainlens[:chain])) / chainlens[chain]
			spr.color = vp.vector(rcos(rb_pos), gcos(rb_pos), bcos(rb_pos))

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
		for i, j, k in zip(hx, hy, hz):
			spr = vp.sphere()
			spr.radius = RAD
			spr.pos = vp.vector(i, j, k)
			spr.color = COLORS['HETATM'][int(rainbow)]
			spr.opacity = OPCTY

if __name__ == "__main__":
	main()
