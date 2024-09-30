"""Plots the atoms in a protein using VPython"""

import vpython as vp

RAD = 0.8

COLORS = {'C': vp.vector(0.42, 0.42, 0.42),
		  'N': vp.vector(0.04, 0.02, 0.99),
		  'O': vp.vector(0.99, 0.0, 0.01),
		  'S': vp.vector(0.98, 0.5, 0.01)}

with open("2bv6.pdb") as file:
	lines = file.readlines()

# extract only the ATOM lines from the database, and their coords
atoms = list(filter(lambda l : l[:4] == "ATOM", lines))
x, y, z = [[float(atom[i:i+8]) for atom in atoms] for i in range(30, 54, 8)]

# sort atoms by element
c, o, n, s = [list(filter(lambda x: x[77] == e, atoms)) for e in ['C', 'O', 'N', 'S']]

# extract lists of x, y, z coordinates for each element
# cx, cy, cz = [[float(atom[i:i+8]) for atom in c] for i in range(30, 54, 8)]
# nx, ny, nz = [[float(atom[i:i+8]) for atom in n] for i in range(30, 54, 8)]
# ox, oy, oz = [[float(atom[i:i+8]) for atom in o] for i in range(30, 54, 8)]
# sx, sy, sz = [[float(atom[i:i+8]) for atom in s] for i in range(30, 54, 8)]

print(len(atoms), " Atoms:")
print(len(c), " Carbon,")
print(len(n), " Nitrogen,")
print(len(o), " Oxygen,")
print(len(s), " Sulfur")

print(vp.scene.camera.pos, vp.scene.camera.axis)

for i, atom in enumerate(atoms):
	ele = atom[77]
	spr = vp.sphere()
	spr.radius = RAD
	spr.color = COLORS[ele]
	spr.pos = vp.vector(x[i], y[i], z[i])


