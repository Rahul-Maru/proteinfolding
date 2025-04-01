"""CONSTANTS. TO CHANGE PARAMS, MODIFY THIS SECTION."""
import math
import numpy as np

# TODO rearrange

#——DEFAULTS——
DEF_PROT = "sim/proteins/1a0i.pdb"
DEF_LIG = "HSM"
# overrides. if set to true, the respective flag will be assumed to be present
SHOW_HETAMS_OVR = False # Shows Heterogens
RAINBOW_OVR = False # Displays the atoms with a rainbow color scheme

# maximum distance of a residue from a ligand for it to be considered
#	part of the binding site (Å)
BS_DIST = 4.5

#——PDB FIELDS——
# RECORD
REC = slice(0, 6)
# ATOM SERIAL NUMBER
ATNO = slice(6, 11)
# RESIDUE NAME (3-letter)
RESN = slice(17, 20)
# CHAIN ID
CHAIN = 21
# RESIDUE SEQUENCE NUMBER
RES_SEQ = slice(22, 26)
# X COORDINATE
X_COORD = slice(30, 38)
# Y COORDINATE
Y_COORD = slice(38, 46)
# Z COORDINATE
Z_COORD = slice(46, 54)
# ELEMENT
ELEM = 77


#——TABULAR DATA——
ELEMS = ['H', 'C', 'N', 'O', 'S']

# map from 3- to 1-letter codes of the amino acids (+ a blank residue and terminator record)
AA_MAP = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
		  'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
		  'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
		  'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
		  'MSE': 'SeM', '___': '-', 'TER' : '|\n'}
# to invert, uncomment the following line:
# AA_MAP_INV = {v: k for k, v in AA_MAP.items()}

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

def chain_to_color(chain: str, elem: str) -> tuple:
	"""Converts a chain ID to a color index."""
	# TODO why is this red?
	return (ord(chain) - ord('A')) % len(COLORS[elem])

def chainlen(chain: str) -> int:
	"""Returns the length of a chain."""
	# TODO
	pass


def dist(p, q):
    """Gets the distance between two triplets of coordinates.
    
    Args:
        p: First triplet of coordinates
        q: Second triplet of coordinates
        
    Returns:
        float: Euclidean distance between the two points
    """
    return math.sqrt(sum([(p[i] - q[i])**2 for i in range(3)]))

def hue_to_RGB(θ: float) -> tuple:
	"""Given a hue value θ ∈ [0, 2π] and converts it to an RGB vector."""

	rcos = lambda x : min(1, max(0, np.cos(x) + 0.5))
	gcos = lambda x : min(1, max(0, np.cos(x - 2*np.pi/3) + 0.5))
	bcos = lambda x : min(1, max(0, np.cos(x + 2*np.pi/3) + 0.5))

	return (rcos(θ), gcos(θ), bcos(θ))