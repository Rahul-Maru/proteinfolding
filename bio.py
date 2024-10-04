"""Plots the atoms in a protein using VPython"""
import numpy as np
import vpython as vp

# the protein to be modeled (models file at proteins/{PROT}.pdb)
PROT = "1rwt"

SHOW_HETAMS = True
RAINBOW = False

# convert the chain letter of an atom to a numerical index
CHAIN_ID = lambda a : ord(a[21]) - 65

# radius of the atoms
RAD = 0.8

# opacity of HETATMs
OPCTY = 1

# colors for depending on atom and chain
COLORS = {'C': [vp.vector(0.5, 0.42, 0.42), vp.vector(0.42, 0.5, 0.42), vp.vector(0.42, 0.42, 0.5)],
		  'N': [vp.vector(0, 0, 1), vp.vector(0, 0.52, 0.87), vp.vector(0.13, 0, 0.75)],
		  'O': [vp.vector(1, 0, 0), vp.vector(0.85, 0.2, 0.1), vp.vector(0.75, 0, 0.2)],
		  'S': [vp.vector(1, 0.5, 0), vp.vector(1, 0.58, 0), vp.vector(1, 0.7, 0.1)],
		  'HETATM': [vp.vector(0.08, 0.69, 0.08), vp.vector(1, 0.96, 0.85)]} # 0 - Normal, 1 - Rainbow

with open(f"proteins/{PROT}.pdb") as file:
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

	if RAINBOW:
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

if SHOW_HETAMS:
	for i, j, k in zip(hx, hy, hz):
		spr = vp.sphere()
		spr.radius = RAD
		spr.pos = vp.vector(i, j, k)
		spr.color = COLORS['HETATM'][int(RAINBOW)]
		spr.opacity = OPCTY
