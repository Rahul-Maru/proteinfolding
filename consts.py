"""CONSTANTS. TO CHANGE THE PARAMS, MODIFY THIS FILE."""
import math
import numpy as np

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
			   '7a5v.pdb', '7dfl.pdb', '7qn7.pdb', '7qn8.pdb', '7qn9.pdb', '7qnc.pdb',
			   '7qnd.pdb', '7yfc.pdb', '8hn8.pdb', '8jt5.pdb', '8jxt.pdb', '8pok.pdb',
			   '8pvb.pdb', '8tgk.pdb', '8ucn.pdb', '8yn2.pdb', '8yn3.pdb', '8yn4.pdb',
			   '8yn5.pdb', '8yn9.pdb', '8yuu.pdb']

pdb_ids = ['1avn.pdb', '1ike.pdb', '1jqd.pdb', '1kar.pdb', '1np1.pdb', '1qft.pdb', '2qeb.pdb',
		   '2x45.pdb', '3bu1.pdb', '3rxh.pdb', '4xmf.pdb', '5rus.pdb', '6dyn.pdb', '6fu4.pdb',
		   '6qfa.pdb', '7dfl.pdb', '7yfc.pdb', '8jt5.pdb', '8pok.pdb', '8tgk.pdb', '8yuu.pdb']


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