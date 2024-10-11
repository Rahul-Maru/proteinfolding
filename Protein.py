from functools import cached_property
import numpy as np

#——UTILS——
# convert the chain letter of an atom to a numerical index
CHAIN_ID = lambda atom : ord(atom[21]) - 65

class Protein:
	"""A class to store data about the spacial information of a protein."""

	def __init__(self, prot_path) -> None:
		# opens the file and splits it into lines
		with open(f"{prot_path}") as file:
			lines = file.readlines()

		# extracts only the ATOM and HETATM lines from the database
		self.atoms = list(filter(lambda l : l[:4] == "ATOM", lines))
		self.hetatms = list(filter(lambda l : l[:6] == "HETATM", lines))

		# store the coordinates to avoid recalculating unnecessarily
		self.x, self.y, self.z = self.get_xyz()

	def get_xyz(self, chain=-1):
		"""Extracts coordinates from the atoms lists as separate lists. -1 is the full list
		-2 denotes the heterogen list."""

		match chain:
			case -1:
				x, y, z = [[float(atom[i:i+8]) for atom in (self.atoms)]
						for i in range(30, 54, 8)]
			case -2:
				x, y, z = [[float(atom[i:i+8]) for atom in (self.hetatms)]
						for i in range(30, 54, 8)]
			case _:
				x, y, z = [[float(atom[i:i+8]) for atom in (self.chains[chain])]
						for i in range(30, 54, 8)]
		return (x, y, z)

	@cached_property
	def chains(self):
		"""Splits self.atoms into a list for each monomeric chain."""
		chains = []
		for atom in self.atoms:
			try:
				chains[CHAIN_ID(atom)].append(atom)
			except IndexError:
				chains.append([atom])

		return chains

	@cached_property
	def elems(self):
		"""Splits self.atoms into different lists based on the element of each atom."""

		c, n, o, s = [list(filter(lambda x: x[77] == e, self.atoms))
				for e in ['C', 'N', 'O', 'S']]
		return (c, n, o, s)

	def info(self):
		"""Outputs information about the protein to the terminal."""

		c, n, o, s = self.elems

		print(f"{len(self.chains)} Chains,")
		print(len(self.atoms), " Atoms:")
		print(len(c), " Carbon (grey),")
		print(len(n), " Nitrogen (blue),")
		print(len(o), " Oxygen (red),")
		print(len(s), " Sulfur (orange),")
		print(len(self.hetatms), " Hetero-atoms (magenta)")

	def centroid(self, chain=-1):
		"""Calculates the centroid of the protein or the heterogens.
		-1 = full protein, -2 = hetatms."""

		coords = self.get_xyz(chain)
		sum = np.zeros(3)
		for x, y, z in zip(*coords):
			sum += np.array([x, y, z])

		return sum/len(coords)
