"""Plots the atoms in a protein using VPython"""
from random import randint
import vpython as vp


SHOW_HETAMS = False

# radius of the atoms
RAD = 0.8

# colors for depending on atom and chain
COLORS = {'C': [vp.vector(0.5, 0.42, 0.42), vp.vector(0.42, 0.5, 0.42), vp.vector(0.42, 0.42, 0.5)],
		  'N': [vp.vector(0, 0, 1), vp.vector(0, 0.52, 0.87), vp.vector(0.13, 0, 0.75)],
		  'O': [vp.vector(1, 0, 0), vp.vector(0.85, 0.2, 0.1), vp.vector(0.75, 0, 0.2)],
		  'S': [vp.vector(1, 0.5, 0), vp.vector(1, 0.58, 0), vp.vector(1, 0.7, 0.1)],
		  'HETATM': vp.vector(0.85, 0.14, 0.85)}

with open("proteins/1rwt.pdb") as file:
	lines = file.readlines()

# extract only the ATOM and HETATM lines from the database and their coords
atoms = list(filter(lambda l : l[:4] == "ATOM", lines))
hetatms = list(filter(lambda l : l[:6] == "HETATM", lines))
x, y, z = [[float(atom[i:i+8]) for atom in atoms] for i in range(30, 54, 8)]
hx, hy, hz = [[float(hetatm[i:i+8]) for hetatm in hetatms] for i in range(30, 54, 8)]
# sort atoms by element
c, o, n, s = [list(filter(lambda x: x[77] == e, atoms)) for e in ['C', 'O', 'N', 'S']]

print(len(atoms), " Atoms:")
print(len(c), " Carbon (grey),")
print(len(n), " Nitrogen (blue),")
print(len(o), " Oxygen (red),")
print(len(s), " Sulfur (orange),")
print(len(hetatms), " Hetero-atoms (magenta)")

for i, atom in enumerate(atoms):
	spr = vp.sphere()

	elem = atom[77]
	chain = ord(atom[21]) - 65 # convert the chain letter to index
	# pick the appropriate color based on element and chain
	if chain <= 2:
		spr.color = COLORS[elem][chain]
	else:
		# if the chain number exceeds the number of available colors,
		# loop through the colors
		spr.color = COLORS[elem][chain % 3]

	spr.pos = vp.vector(x[i], y[i], z[i])
	spr.radius = RAD
	# get the appropriate color based on element and side chain

if SHOW_HETAMS:
	for i, j, k in zip(hx, hy, hz):
		spr = vp.sphere()
		spr.radius = RAD
		spr.color = COLORS['HETATM']
		spr.pos = vp.vector(i, j, k)
