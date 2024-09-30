"""Plots the atoms in a protein using VPython"""

import vpython as vp

RAD = 0.8

COLORS = {'C': vp.vector(0.42, 0.42, 0.42),
		  'N': vp.vector(0.04, 0.02, 0.99),
		  'O': vp.vector(0.99, 0.0, 0.01),
		  'S': vp.vector(0.98, 0.5, 0.02),
		  'HETATM': vp.vector(0.85, 0.14, 0.85)}

with open("1rwt.pdb") as file:
	lines = file.readlines()

# extract only the ATOM and HETATM lines from the database and their coords
atoms = list(filter(lambda l : l[:4] == "ATOM", lines))
hetatms = list(filter(lambda l : l[:6] == "HETATM", lines))
x, y, z = [[float(atom[i:i+8]) for atom in atoms] for i in range(30, 54, 8)]
hx, hy, hz = [[float(hetatm[i:i+8]) for hetatm in hetatms] for i in range(30, 54, 8)]
# sort atoms by element
c, o, n, s = [list(filter(lambda x: x[77] == e, atoms)) for e in ['C', 'O', 'N', 'S']]

# extract lists of x, y, z coordinates for each element
# cx, cy, cz = [[float(atom[i:i+8]) for atom in c] for i in range(30, 54, 8)]
# nx, ny, nz = [[float(atom[i:i+8]) for atom in n] for i in range(30, 54, 8)]
# ox, oy, oz = [[float(atom[i:i+8]) for atom in o] for i in range(30, 54, 8)]
# sx, sy, sz = [[float(atom[i:i+8]) for atom in s] for i in range(30, 54, 8)]

print(len(atoms), " Atoms:")
print(len(c), " Carbon (grey),")
print(len(n), " Nitrogen (blue),")
print(len(o), " Oxygen (red),")
print(len(s), " Sulfur (orange),")
print(len(hetatms), " Hetero-atoms (magenta)")

for i, atom in enumerate(atoms):
	ele = atom[77]
	spr = vp.sphere()
	spr.radius = RAD
	spr.color = COLORS[ele]
	spr.pos = vp.vector(x[i], y[i], z[i])

for i, j, k in zip(hx, hy, hz):
	spr = vp.sphere()
	spr.radius = RAD
	spr.color = COLORS['HETATM']
	spr.pos = vp.vector(i, j, k)
