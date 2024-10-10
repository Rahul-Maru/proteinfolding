import numpy as np

#——UTILS——
# convert the chain letter of an atom to a numerical index
CHAIN_ID = lambda a : ord(a[21]) - 65

class Protein:
	"""A class to store data about the spacial information of a protein"""
	def __init__(self, prot_path) -> None:
		# opens the file and splits it into lines
		with open(f"{prot_path}") as file:
			lines = file.readlines()

		# extracts only the ATOM and HETATM lines from the database
		self.atoms = list(filter(lambda l : l[:4] == "ATOM", lines))
		self.hetatms = list(filter(lambda l : l[:6] == "HETATM", lines))

		# extracts coordinates from the atoms lists as separate lists
		self.x, self.y, self.z = [[float(atom[i:i+8]) for atom in self.atoms]
			 for i in range(30, 54, 8)]
		self.hx, self.hy, self.hz = [[float(hetatm[i:i+8]) for hetatm in self.hetatms]
				for i in range(30, 54, 8)]

		# stores the length of each monomeric chain in a list
		self.chainlens = []
		for atom in self.atoms:
			try:
				self.chainlens[CHAIN_ID(atom)] += 1
			except:
				self.chainlens.append(1)

	def elems(self):
		c, n, o, s = [list(filter(lambda x: x[77] == e, self.atoms))
				for e in ['C', 'N', 'O', 'S']]
		return (c, n, o, s)

	def info(self):
		"""Outputs information about the protein to the terminal"""
		c, n, o, s = self.elems()

		print(f"{len(self.chainlens)} Chains,")
		print(len(self.atoms), " Atoms:")
		print(len(c), " Carbon (grey),")
		print(len(n), " Nitrogen (blue),")
		print(len(o), " Oxygen (red),")
		print(len(s), " Sulfur (orange),")
		print(len(self.hetatms), " Hetero-atoms (magenta)")

	def centroid(self, het=False):
		"""Calculates the centroid of the protein or the heterogens"""
		sum = np.zeros(3)
		for x, y, z in (zip(self.hx, self.hy, self.hz) if het else zip(self.x, self.y, self.z)):
			sum += np.array([x, y, z])

		return sum/len(self.hetatms if het else self.atoms)

