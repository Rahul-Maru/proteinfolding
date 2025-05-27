"""PROTEIN CLASS."""

from functools import cached_property
from os.path import basename
import numpy as np
from typing import Literal

from consts import *
from residue import Residue

#TODO add a residue class or at least a typed dict for residues

class Protein:
	"""A class to store data about the spacial information of a protein."""

	def __init__(self, prot_path, dp_mode: list[bool] = [False]*3) -> None:
		"""Initializes a protein. `dp_mode` is a list of the following values:
		color_scheme: `False` (normal), `True` (rainbow)
		het_atms: `False` (only shows normal atoms), `True` (shows hetatms as well)."""

		self.path = prot_path
		self.id = basename(prot_path).upper().split('.')[0]

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
		self.x, self.y, self.z = self.get_coords()


	@cached_property
	def chains(self):
		"""Group the atoms by chain.

		Returns:
			dict: A dictionary mapping chain IDs (str) to lists of ATOM records for that chain.
		"""

		chains = {}
		for atom in self.atoms:
			try:
				# add the new atom to its chain
				chains[atom[CHAIN]].append(atom)
			except KeyError:
				# if this is the first atom, create a new key-list pair for the chain
				chains[atom[CHAIN]] = [atom]

		return chains

	@cached_property
	def residues(self) -> list[Residue]:
		"""Group atoms by residue.

		Returns:
			List of residue dictionaries containing:
				code (str): 3-letter residue code
				id (str): Chain ID + sequence number
				atoms (list): ATOM records for residue
		"""
		residues = []
		last_res_id = ''

		for atom in self.atoms + self.hetatms:
			res_id = f'{atom[CHAIN]}_{atom[RES_SEQ].strip()}'
			if res_id != last_res_id:
				residues.append({
					'code': atom[RESN],
					'id': res_id,
					'type': atom[REC].strip().lower(),
					'atoms': [atom]
				})
				last_res_id = res_id
			else:
				residues[-1]['atoms'].append(atom)

		return residues

	def res_subset(self, start:int, stop:int = None):
		"""Returns a subset of `self.residues` containing all the residues from
		indices start to stop, inclusive"""

		if not stop:
			start, stop = 1, start

		# TODO fix this
		return [res for res in self.residues if start <= res['n'] <= stop]


	@cached_property
	def elems(self):
		"""Splits `self.atoms` into different lists based on the element of each atom."""

		h, c, n, o, s = [[atom[ELEM] for atom in self.atoms if atom[ELEM] == elmt]
				for elmt in ELEMS]
		return (h, c, n, o, s)

	def get_coords(self, query=None, qtype="", triplet = False):
		"""Extracts 3D coordinates from an atom or set of atoms.

		Args:
			query (optional): Residue, record, or chain to get coordinates from.
			qtype (str, optional): Type of `query`. Can be `'atm'`, which extracts from all atoms \
				in the protein; `'het'`, which extracts from the hetatms; `all`, which extracts from both atoms and hetatms;
				`'res'`, which extracts from a given residue/molecule; \
				`'lin'`, which extracts from a given line of the PDB file; \
				or a number `n`, which extracts from the n-th chain. \
				Defaults to `"res"` if `query` else `"all"`.
			triplet (bool, optional): If `True`, returns the coordinates in triplets of the form \
				`(x, y, z)`. If `False`, returns the lists of `x`'s, `y`'s, & `z`'s individually. \
				Defaults to `False`.
		"""

		# default value
		# TODO make this work for chains
		if not qtype:
			qtype = "res" if query  else ("all" if self.display_mode[1] else "atm")

		match qtype:
			case "all":
				x, y, z = [[float(atom[coord]) for atom in (self.atoms + self.hetatms)]
						for coord in [X_COORD, Y_COORD, Z_COORD]]
			case "atm":
				x, y, z = [[float(atom[coord]) for atom in self.atoms]
						for coord in [X_COORD, Y_COORD, Z_COORD]]
			case "het":
				x, y, z = [[float(atom[coord]) for atom in self.hetatms]
						for coord in [X_COORD, Y_COORD, Z_COORD]]
			case "res":
				x, y, z = [[float(atom[coord]) for atom in query['atoms']]
					for coord in [X_COORD, Y_COORD, Z_COORD]]
			case "lin":
				x, y, z = (float(query[coord]) for coord in [X_COORD, Y_COORD, Z_COORD])
			case _:
				x, y, z = [[float(atom[coord]) for atom in self.chains[query]]
						for coord in [X_COORD, Y_COORD, Z_COORD]]

		if triplet:
			return [(i, j, k) for i, j, k, in zip(x, y, z)]
		else:
			return (x, y, z)

	def get_record(self, record: str) -> list[str]:
		"""Get all PDB lines with given record type."""

		return [line for line in self.lines if line[REC].strip() == record]

	def get_secondary_structure(self, ss_type: Literal['HELIX', 'SHEET']):
		"""Returns the residues of all the instances of the given secondary structure"""

		# TODO fix res_subset
		if ss_type == 'HELIX':
			helices = self.get_record('HELIX')
			return sum([self.res_subset(int(lin[21:25]), int(lin[33:37])) for lin in helices], [])
		else:
			sheets = self.get_record('SHEET')
			return sum([self.res_subset(int(lin[22:26]), int(lin[33:37])) for lin in sheets], [])


	def centroid(self, query = None, qtype=""):
		"""Calculates the centroid of a given set of atoms."""

		coords = self.get_coords(query, qtype=qtype)
		if n := len(coords[0]):
			sum = np.zeros(3)
			for x, y, z in zip(*coords):
				sum += np.array([x, y, z])

			return sum/n
		else:
			raise ValueError(f"Query (type '{qtype}') has no atoms")


	def seq(self, residues: list[Residue] = None):
		"""Returns the 1-letter sequence of the residues in the given list,
		or in the whole protein if no list is given."""
		# TODO test this


		if not residues:
			residues = self.residues

		if not residues: # Handle empty list case
			return ""

		seq = []

		# --- Handle leading dashes for the very first chain ---
		first_num_str = residues[0]['id'].split('_')[1]
		first_num = int(first_num_str)
		if first_num > 1:
			seq.extend(['-'] * (first_num - 1))

		for i in range(len(residues)-1):
			# Add current residue
			seq.append(amino_acid_map(residues[i]['code']))

			# Check for gap to next residue
			curr_chain, curr_num_str = residues[i]['id'].split('_')
			next_chain, next_num_str = residues[i+1]['id'].split('_')
			curr_num = int(curr_num_str)
			next_num = int(next_num_str)

			# Only add gaps within same chain
			if curr_chain == next_chain:
				gap_size = next_num - curr_num - 1
				if gap_size > 0:
					seq.extend(['-'] * gap_size)
			else:
				# Chain break
				seq.append('|')
				# --- Handle leading dashes for the new chain ---
				if next_num > 1:
					seq.extend(['-'] * (next_num - 1))

		# Add final residue
		seq.append(amino_acid_map(residues[-1]['code']))

		return ''.join(seq)


	@DeprecationWarning
	def terminations(self):
		"""Successively yield the positions (atom no. and residue no.) of each TER record in the file."""

		# TODO remove this probably
		ters = self.get_record('TER')
		if ters:
			for ter in ters:
				yield (int(ter[6:11]), int(ter[22:26]))
		else:
			# if there are no TERs, yield an extremely large position
			yield (int(1e15), int(1e15))



	def get_ligand(self, query: str = DEF_LIG, chain: str = "all") -> list:
		"""Extracts a ligand based on a query string

		Args:
			chain (str): chain from which to extract ligands
			query (str): the 3-letter code of the ligand

		Returns:
			list: The list of `query` residues in `self.residues`
		"""

		if chain == "all":
			# each record must be a HETATM record corresponding to the query ligand
			lig = [res for res in self.residues if res['code'] == query]

		else:
			lig = [res for res in self.residues if res['code'] == query and res['id'][0] == chain]

		return lig

	def get_bsite(self, lig, include_lig: bool = False, max_dist: float = MAX_BSITE_DIST):
		"""Extracts the binding site(s) that bind to the given ligand.

		Args:
			lig (dict): dict containing the atoms of the ligand.
			include_lig (bool): whether to include the ligand in the binding site.
			max_dist (float): the maximum distance for a residue to be included in the binding site. (default: 4.5Å)
		Returns:
			(list): list of residues that comprise the binding site(s).
		"""
		#TODO try out different cutoff values
		bsite = []

		lig_coords = self.get_coords(lig, triplet=True)

		aminos = [res for res in self.residues if res['type'] == 'atom']

		for res in aminos:
			for atom in res['atoms']:
				atom_coords = self.get_coords(atom, "lin")
				# if the atom is close to the ligand add it to the binding
				# 	site and move on to the next residue
				for coords in lig_coords:
					if dist(coords, atom_coords) <= max_dist:
						bsite.append(res)
						break
				else:
					# if the inner loop was broken, break out of this one too
					continue
				break

		if include_lig:
			bsite.append(lig)

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

	@DeprecationWarning
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

		# TODO remove this probably

		coords1, coords2 = *self.get_coords(query, qtype, True), *self.get_coords(query2, qtype2, True)

		d = []
		for i in coords2:
			d.append([])
			for j in coords1:
				d[-1].append(dist(i, j))

		return d