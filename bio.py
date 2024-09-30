from matplotlib import pyplot as plt

with open("1acj.pdb") as file:
	lines = file.readlines()
	# extract only the ATOM lines from the database 
	atoms = list(filter(lambda l : l[:4] == "ATOM", lines))

# sort atoms by element
c, o, n, s = [list(filter(lambda x: x[77] == e, atoms)) for e in ['C', 'O', 'N', 'S']]

# extract lists of x, y, z coordinates for each element
cx, cy, cz = [[float(atom[i:i+8]) for atom in c] for i in range(30, 54, 8)]
ox, oy, oz = [[float(atom[i:i+8]) for atom in o] for i in range(30, 54, 8)]
nx, ny, nz = [[float(atom[i:i+8]) for atom in n] for i in range(30, 54, 8)]
sx, sy, sz = [[float(atom[i:i+8]) for atom in s] for i in range(30, 54, 8)]

fig = plt.figure()
# 3D graph
ax = fig.add_subplot(111, projection='3d')

# black background
fig.patch.set_facecolor('black')
ax.patch.set_facecolor('black')

# plot and color the atoms by element
plt.plot(cx, cy, cz, "o", color="xkcd:dark grey", label="C")
plt.plot(ox, oy, oz, "o", color="xkcd:fire engine red", label="O")
plt.plot(nx, ny, nz, "o", color="xkcd:electric blue", label="N")
plt.plot(sx, sy, sz, "o", color="xkcd:pumpkin orange", label="S")

plt.legend()
plt.show()
