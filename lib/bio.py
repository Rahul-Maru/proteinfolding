from functools import cached_property
import math
import numpy as np
from typing import Literal

"""CONSTANTS. TO CHANGE PARAMS, MODIFY THIS SECTION."""
#——DEFAULTS——
DEF_PROT = "1a0i"
DEF_LIG = "HSM"
# overrides. if set to true, the respective flag will be assumed to be present
SHOW_HETAMS_OVR = False # Shows Heterogens
RAINBOW_OVR = False # Displays the atoms with a rainbow color scheme

# maximum distance of a residue from a ligand for it to be considered
#	part of the binding site (Å)
BS_DIST = 4.5

#——UTILS——
# converts the chain letter of an atom to a numerical index
CHAIN_ID = lambda atom : ord(atom[21]) - 65
# converts the numerical index of a chain back to a letter
CHAIN_LET = lambda i : chr(i + 65)
# gets the residue number from an atom
RES_NUM = lambda atom : int(atom[22:26])
# gets the distance between 3 triplets of coordinates
DIST = lambda p, q: math.sqrt(sum([(p[i] - q[i])**2 for i in range(3)]))


#——TABULAR DATA——
ELEMS = ['H', 'C', 'N', 'O', 'S']

# map from 3- to 1-letter codes of the amino acids (+ a blank residue and terminator record)
AA_MAP = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
		  'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
		  'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
		  'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
		  'MSE': 'SeM', '___': '-', 'TER' : '|\n'}
# to invert:
# AA_MAP_INV = {v: k for k, v in AA_MAP.items()}

# list of pdb ids of HSM-binding proteins
pdb_ids_all = ['1avn.pdb', '1ike.pdb', '1jqd.pdb', '1kar.pdb', '1np1.pdb', '1qft.pdb',
			   '1qfv.pdb', '1u18.pdb', '2qeb.pdb', '2x45.pdb', '3bu1.pdb', '3g7x.pdb',
			   '3rxh.pdb', '4xmf.pdb', '5rus.pdb', '6dyn.pdb', '6fu4.pdb', '6qfa.pdb',
			   '7a5v.pdb', '7dfl.pdb', '7pbd.pdb', '7pbz.pdb', '7pc0.pdb', '7qn7.pdb',
			   '7qn8.pdb', '7qn9.pdb', '7qnc.pdb', '7qnd.pdb', '7yfc.pdb', '8hn8.pdb',
			   '8jt5.pdb', '8jxt.pdb', '8pok.pdb', '8pvb.pdb', '8tgk.pdb', '8ucn.pdb',
			   '8yn2.pdb', '8yn3.pdb', '8yn4.pdb', '8yn5.pdb', '8yn9.pdb', '8yuu.pdb']

pdb_ids = ['1avn.pdb', '1ike.pdb', '1jqd.pdb', '1kar.pdb', '1np1.pdb', '1qft.pdb', '2qeb.pdb',
		   '2x45.pdb', '3bu1.pdb', '3rxh.pdb', '4xmf.pdb', '5rus.pdb', '6dyn.pdb', '6fu4.pdb',
		   '6qfa.pdb', '7dfl.pdb', '7yfc.pdb', '8jt5.pdb', '8pok.pdb', '8tgk.pdb', '8yuu.pdb']

# bsites = ['1avn_A.pdb', '1ike_A.pdb', '1jqd_A.pdb', '1jqd_B.pdb', '1kar_A.pdb', '1kar_B.pdb', '1np1_A.pdb',
# 		  '1np1_B.pdb', '1qft_A.pdb', '1qft_B.pdb', '1qfv_A.pdb', '1qfv_B.pdb', '1u18_A.pdb', '1u18_B.pdb',
# 		  '2qeb_A.pdb', '2qeb_B.pdb', '2x45_A.pdb', '2x45_B.pdb', '2x45_C.pdb', '3bu1_A.pdb', '3g7x_A.pdb',
# 		  '3g7x_B.pdb', '3rxh_A.pdb', '4xmf_A.pdb', '5rus_A.pdb', '6dyn_A.pdb', '6fu4_A.pdb', '6fu4_B.pdb',
# 		  '6fu4_C.pdb', '6fu4_D.pdb', '6qfa_A.pdb', '6qfa_B.pdb', '6qfa_C.pdb', '6qfa_D.pdb', '6qfa_E.pdb',
# 		  '7a5v_A.pdb', '7dfl_R.pdb', '7qn7_C.pdb', '7qn7_D.pdb', '7qn8_B.pdb', '7qn8_D.pdb', '7qn9_B.pdb',
# 		  '7qn9_C.pdb', '7qn9_D.pdb', '7qnc_C.pdb', '7qnc_D.pdb', '7qnd_B.pdb', '7qnd_C.pdb', '7yfc_R.pdb',
# 		  '8hn8_R.pdb', '8jt5_A.pdb', '8jxt_R.pdb', '8pok_A.pdb', '8pvb_A.pdb', '8tgk_A.pdb', '8ucn_A.pdb',
# 		  '8yn2_R.pdb', '8yn3_R.pdb', '8yn4_R.pdb', '8yn5_R.pdb', '8yn9_R.pdb', '8yuu_R.pdb']


#——DISPLAY PROPERTIES——
# radius of the atoms
RAD = 0.8
# opacity of HETATMs (WARNING: changing this may cause severe lag)
OPCTY = 1

# colors for depending on element, chain, and rainbow mode
COLORS = {'H': [(0.75, 0.7, 0.7), (0.7, 0.75, 0.7), (0.7, 0.7, 0.75)],
		'C': [(0.5, 0.42, 0.42), (0.42, 0.5, 0.42), (0.42, 0.42, 0.5)],
		'N': [(0, 0, 1), (0, 0.52, 0.93), (0.15, 0, 0.78)],
		'O': [(1, 0, 0), (0.85, 0.2, 0.1), (0.75, 0, 0.2)],
		'S': [(1, 0.5, 0), (1, 0.58, 0), (1, 0.7, 0.1)],
		'HETATM': [(0.08, 0.7, 0.08), (1, 0.96, 0.85)],	#  v
		'CENTRD': [(1, 1, 1), (0.2, 0.21, 0.24)],		# 0 - Normal, 1 - Rainbow
		'HETCTD': [(0.2, 1, 0.65), (0.64, 0.6, 0.55)]}	#  ^

def hue_to_RGB(θ: float) -> tuple:
	"""Given a hue value θ ∈ [0, 2π] and converts it to an RGB vector."""

	rcos = lambda x : min(1, max(0, np.cos(x) + 0.5))
	gcos = lambda x : min(1, max(0, np.cos(x - 2*np.pi/3) + 0.5))
	bcos = lambda x : min(1, max(0, np.cos(x + 2*np.pi/3) + 0.5))

	return (rcos(θ), gcos(θ), bcos(θ))


"""PROTEIN CLASS"""
class Protein:
	"""A class to store data about the spacial information of a protein."""

	def __init__(self, prot_path, dp_mode: list[bool] = [False]*3) -> None:
		"""Initializes a protein. `dp_mode` is a list of the following values:
		color_scheme: `False` (normal), `True` (rainbow)
		het_atms: `False` (only shows normal atoms), `True` (shows hetatms as well)."""

		self.path = prot_path
		self.id = prot_path[-8:-4].upper()

		self.display_mode = dp_mode

		# opens the file and splits it into lines
		with open(f"{prot_path}") as file:
			self.lines = file.readlines()

		# extracts only the ATOM and HETATM lines from the database
		self.atoms = self.get_record("ATOM")
		self.hetatms = self.get_record("HETATM")

		# if there are no hetatms, turn off show_hetatms
		self.display_mode[1] = dp_mode[1] and len(self.hetatms)

		# store the coordinates to avoid recalculating unnecessarily
		self.x, self.y, self.z = self.get_xyzlist()


	@cached_property
	def chains(self):
		"""`self.atoms` split by monomeric chain."""

		chains = {}
		for atom in self.atoms:
			try:
				# add the new atom to its chain
				chains[atom[21]].append(atom)
			except KeyError:
				# if this is the first atom, create a new key-list pair for the chain
				chains[atom[21]] = [atom]

		return chains

	@cached_property
	def residues(self):
		"""`self.atoms` split by residue

		Returns:
			dict: keys:
				code (str): the 3-letter code of the residue
				n (int): residue number
				atoms (list): the list of ATOM records comprising the residue
		"""

		res = []

		# yields the atom and residue number of the next TER record
		ter = self.terminations()
		ter_atm_no, ter_res_no = next(ter)
		# start point of the current chain
		st = 0
		for atom in self.atoms:
			# If the loop has crossed the TER record, update st and yield the next TER record,
			#  moving to the next chain
			if int(atom[6:11]) > ter_atm_no:
				res.append({'code': 'TER', 'n': ter_res_no + 1, 'atoms': []})
				st = ter_res_no + 1
				ter_atm_no, ter_res_no = next(ter)
				ter_res_no += st

			# each chain residue number must have st added to it to avoid overlaps between
			#   residues with the same number in different chain
			n = RES_NUM(atom) + st

			try:
				# (try to) add the new atom to its residue
				res[n]['atoms'].append(atom)
			except IndexError:
				# if this is the first atom, create a new residue entry

				# what should be the number of the next residue, assuming continuity
				next_n = 0 if len(res) == 0 else res[-1]['n'] + 1

				# if there is a gap between this residue and the previous, fill it with
				#   blank residues before adding the next one
				if n > next_n:
					res.extend([{'code': '___', 'n': i, 'atoms': []}
						for i in range(next_n, n)])

				res.append({'code': atom[17:20], 'n': n, 'atoms': [atom]})

		return res


	@cached_property
	def elems(self):
		"""Splits `self.atoms` into different lists based on the element of each atom."""

		h, c, n, o, s = [list(filter(lambda x: x[77] == elmt, self.atoms))
				for elmt in ELEMS]
		return (h, c, n, o, s)

	def get_xyzlist(self, query=None, qtype="", triplet = False):
		"""Extracts 3D coordinates from an atom or set of atoms.

		Args:
			query (optional): Residue, record, or chain to get coordinates from.
			qtype (str, optional): Type of `query`. Can be `'all'`, which extracts from all atoms \
				in the protein; `'het'`, which extracts from the hetatms; `'res'`, which extracts \
				from a given residue/molecule; `'lin'`, which extracts from a given line of the PDB file; \
				or a number `n`, which extracts from the n-th chain. \
				Defaults to `"res"` if `query` else `"all"`.
			triplet (bool, optional): If `True`, returns the coordinates in triplets of the form \
				`(x, y, z)`. If `False`, returns the lists of `x`'s, `y`'s, & `z`'s individually. \
				Defaults to `False`.
		"""

		# default value
		# TODO make this work for chains
		if not qtype:
			qtype = "res" if query  else "all"

		match qtype:
			case "all":
				x, y, z = [[float(atom[i:i+8]) for atom in self.atoms]
						for i in range(30, 54, 8)]
			case "het":
				x, y, z = [[float(atom[i:i+8]) for atom in self.hetatms]
						for i in range(30, 54, 8)]
			case "res":
				x, y, z = [[float(atom[i:i+8]) for atom in query['atoms']]
					for i in range(30, 54, 8)]

			case "lin":
				x, y, z = (float(query[i:i+8]) for i in range(30, 54, 8))
			case _:
				x, y, z = [[float(atom[i:i+8]) for atom in self.chains[query]]
						for i in range(30, 54, 8)]

		if triplet:
			return [(i, j, k) for i, j, k, in zip(x, y, z)]
		else:
			return (x, y, z)


	def get_record(self, record: str):
		"""Returns all the entries in the PDB file with the given record."""

		return list(filter(lambda l : l[:6].strip() == record, self.lines))


	def get_ss(self, ss_type: Literal['HELIX', 'SHEET']):
		"""Returns the residues of all the instances of the given secondary structure"""

		if ss_type == 'HELIX':
			helices = self.get_record('HELIX')
			return sum([self.res_subset(int(ln[21:25]), int(ln[33:37])) for ln in helices], [])
		else:
			sheets = self.get_record('SHEET')
			return sum([self.res_subset(int(ln[22:26]), int(ln[33:37])) for ln in sheets], [])


	def centroid(self, query = None, qtype=""):
		"""Calculates the centroid of a given set of atoms."""

		coords = self.get_xyzlist(query, qtype=qtype)
		if n := len(coords[0]):
			sum = np.zeros(3)
			for x, y, z in zip(*coords):
				sum += np.array([x, y, z])

			return sum/n
		else:
			raise ValueError(f"Query (type '{qtype}') has no atoms")


	def res_subset(self, start:int = None, stop:int = None):
		"""Returns a subset of `self.residues` containing all the residues from
		indices start to stop, inclusive"""

		if start is None and stop is None:
			start = 1
			stop = self.residues[-1]['n']
		elif stop is None:
			start, stop = 1, start

		return [res for res in self.residues if start <= res['n'] <= stop]


	def seq(self, residues=None):
		"""Returns the 1-letter sequence of the residues in the given list"""
		if not residues:
			residues = self.residues

		return ''.join([AA_MAP[res['code']] for res in residues])


	def terminations(self):
		"""Successively yield the positions (atom no. and residue no.) of each TER record in the file."""

		ters = self.get_record('TER')
		if ters:
			for ter in ters:
				yield (int(ter[6:11]), int(ter[22:26]))
		else:
			# if there are no TERs, yield an extremely large position
			yield (int(1e15), int(1e15))



	def get_ligand(self, query: str, chain = "all") -> dict:
		"""Extracts a ligand based on a query string

		Args:
			chain (str): chain from which to extract ligands
			query (str): the 3-letter code of the ligand

		Returns:
			dict: The HETATM lines associated with the target molecule
		"""

		lig = {"code": query}

		if chain == "all":
			# each record must be a HETATM record corresponding to the query ligand
			lig["atoms"] = [atm for atm in self.lines
				   if atm[17:20] == query
				   and atm[:6].strip() in ['ATOM', 'HETATM']]
		else:
			lig["atoms"] = [atm for atm in self.lines
				   if atm[17:20] == query
				   and atm[:6].strip() in ['ATOM', 'HETATM']
				   and atm[21] == chain]
	
		lig['n'] = len(lig["atoms"])
	
		return lig

	def get_bsite(self, lig: dict):
		"""Extracts the binding site or sites that bind to the given ligand.

		Args:
			lig (dict): dict containing the atoms of the ligand.

		Returns:
			(list): list of residues that comprise the binding site(s).
		"""

		lig_coords = self.get_xyzlist(lig, triplet=True)

		bsite = []

		for res in self.residues:
			for atom in res['atoms']:
				atom_coords = self.get_xyzlist(atom, "lin")
				# if the atom is close to the ligand add it to the binding
				# 	site and move on to the next residue
				for coords in lig_coords:
					if DIST(coords, atom_coords) <= BS_DIST:
						bsite.append(res)
						break
				else:
					# if the inner loop was broken, break out of this one too
					continue
				break

		return bsite


	def print_res(residues: list[dict]) -> None:
		"""Prints a list of residues in the format: "3-letter Code" - "Number"

		Args:
			residues (list): list of residues to print.
		"""

		for r in residues:
			print(f"{r['code']} - {r['n']}")

	def __str__(self):
		"""Outputs information about the protein to the terminal."""

		h, c, n, o, s = self.elems
		rb, hetatm, _ = self.display_mode

		lc = len(self.chains)
		hetcolor = f"({"white" if rb else "green"}) (centroid: {"light grey" if rb else "teal"})"\
			if hetatm else "(hidden)"

		s = "\n".join([f"———{self.id}———",
			"Sequence (- represents a gap in sequence, | represents the end of a chain):",
			f"{self.seq()}.",
			f"{lc} Chain{"s" if lc != 1 else ""},",
			f"{len(self.atoms)} Atoms (centroid: {"grey" if rb else "white"}):",
			f"{len(h)} Hydrogen{"" if rb else" (light grey)"},",
			f"{len(c)} Carbon{"" if rb else" (grey)"},",
			f"{len(n)} Nitrogen{"" if rb else" (blue)"},",
			f"{len(o)} Oxygen{"" if rb else" (red)"},",
			f"{len(s)} Sulfur{"" if rb else" (orange)"},",
			f"{len(self.hetatms)} Heterogens {hetcolor}."])

		return s

	# Probably can be removed
	def distance_matrix(self, query, qtype="res", query2 = None, qtype2="all") \
						-> list[list[float]]:
		"""OBSOLETE. Returns a bipartite distance matrix of all atoms in query vs all atoms in query2

		Args:
			query: The first set of atoms.
			qtype (str, optional): Type of `query`. Defaults to `"res"`.
			query2 (optional): The second set of atoms. Defaults to `None`.
				If `qtype2` is also left as default, this makes query2 the set of
				all atoms in the protein.
			qtype2 (str, optional): Type of `query2`. Defaults to `"all"`.

		Returns:
			list[list[float]]: Distance Matrix. d[q2 atom #][q1 atom #s]
		"""

		coords1, coords2 = *self.get_xyzlist(query, qtype, True), *self.get_xyzlist(query2, qtype2, True)

		d = []
		for i in coords2:
			d.append([])
			for j in coords1:
				d[-1].append(DIST(i, j))

		return d