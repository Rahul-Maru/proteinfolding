from functools import cached_property
from typing import Literal
import numpy as np

#——UTILS——
# convert the chain letter of an atom to a numerical index
CHAIN_ID = lambda atom : ord(atom[21]) - 65
# gets the residue number from an atom
RES_NUM = lambda atom : int(atom[22:26])

# map from 3- to 1-letter codes of the amino acids (+ the blank residue)
AA_MAP = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
		  'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
		  'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
		  'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
		  'MSE': 'SeM', '___': '-'}

class Protein:
	"""A class to store data about the spacial information of a protein."""

	def __init__(self, prot_path, display_mode="normal") -> None:
		"""Initializes a protein."""

		self.display_mode = display_mode

		# opens the file and splits it into lines
		with open(f"{prot_path}") as file:
			self.lines = file.readlines()

		# extracts only the ATOM and HETATM lines from the database
		self.atoms = self.get_record("ATOM")
		self.hetatms = self.get_record("HETATM")

		# store the coordinates to avoid recalculating unnecessarily
		self.x, self.y, self.z = self.get_xyzlist()

	@cached_property
	def chains(self):
		"""self.atoms split by monomeric chain."""

		chains = []
		for atom in self.atoms:
			try:
				# add the new atom to its chain
				chains[CHAIN_ID(atom)].append(atom)
			except IndexError:
				# if this is the first atom, create a new chain entry
				chains.append([atom])

		return chains

	@cached_property
	def residues(self):
		"""self.atoms split by residue."""

		res = []
		for atom in self.atoms:
			n = RES_NUM(atom)
			try:
				# (try to) add the new atom to its residue
				res[n-1]['atoms'].append(atom)
			except IndexError:
				# if this is the first atom, create a new residue entry

				# what should be the number of the next residue, assuming continuity
				next_n = 1 if len(res) == 0 else res[-1]['n'] + 1

				# if there is a gap between this residue and the previous, fill it with
				#   blank residues before adding the next one
				if n > next_n:
					res.extend([{'code': '___', 'n': i, 'atoms': []}
						for i in range(next_n, n)])

				res.append({'code': atom[17:20], 'n': n, 'atoms': [atom]})

		return res

	@cached_property
	def elems(self):
		"""Splits self.atoms into different lists based on the element of each atom."""

		c, n, o, s = [list(filter(lambda x: x[77] == e, self.atoms))
				for e in ['C', 'N', 'O', 'S']]
		return (c, n, o, s)

	def get_xyzlist(self, chain=-1):
		"""Extracts coordinates from a list as 3 separate lists. -1 is the full atoms list,
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


	def get_record(self, record: str):
		"""Returns all the entries in self with the given record."""

		return list(filter(lambda l : l[:4] == record, self.lines))

	def get_ss(self, ss_type: Literal['HELIX', 'SHEET']):
		"""Returns the residues of all the instances of the given secondary structure"""

		if ss_type == 'HELIX':
			helices = [atom for atom in self.atoms if atom[:6] == 'HELIX']
			return sum([self.res_subset(int(ln[21:25]), int(ln[33:37])) for ln in helices], [])
		else:
			sheets = [atom for atom in self.atoms if atom[:6] == 'SHEET']
			return sum([self.res_subset(int(ln[22:26]), int(ln[33:37])) for ln in sheets], [])


	def centroid(self, chain=-1):
		"""Calculates the centroid of the protein or the heterogens.
		-1 = full protein, -2 = hetatms."""

		coords = self.get_xyzlist(chain)
		sum = np.zeros(3)
		for x, y, z in zip(*coords):
			sum += np.array([x, y, z])

		return sum/len(coords[0])

	def res_subset(self, start:int = None, stop:int = None):
		"""Returns a subset of self.residues containing all the residues from
		indices start to stop, inclusive"""

		if start is None and stop is None:
			start = 1
			stop = self.residues[-1]['n']
		elif stop is None:
			start, stop = 1, start

		return [res for res in self.residues if start <= res['n'] <= stop]

	def seq(self, residues=None):
		"""Returns the 1 letter sequence of the residues in the given list"""
		if not residues:
			residues = self.residues

		return ''.join([AA_MAP[res['code']] for res in residues])

	def __str__(self):
		"""Outputs information about the protein to the terminal."""

		c, n, o, s = self.elems
		rb = self.display_mode == 'rainbow'

		s = f'{len(self.chains)} Chains,\n\
{len(self.atoms)} Atoms (centroid: {"grey" if rb else "white"}):\n\
{len(c)} Carbon{"" if rb else" (grey)"},\n\
{len(n)} Nitrogen{"" if rb else" (blue)"},\n\
{len(o)} Oxygen{"" if rb else" (red)"},\n\
{len(s)} Sulfur{"" if rb else" (orange)"},\n\
{len(self.hetatms)} Hetero-atoms ({"white" if rb else "green"}, centroid: light {"grey" if rb else "green"}).\n\
Sequence: {self.seq()}'

		return s
