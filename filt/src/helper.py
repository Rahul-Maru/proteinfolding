"""Constants and the fwrite helper used across the filt pipeline: PDB
fixed-width field offsets, the standard amino-acid codes, and fwrite.
Mirrors the subset of lib/consts.py that filt uses."""

# PDB fixed-width field offsets (see filt/file_format.txt)
RESN = slice(17, 20)      # residue name (3-letter)
CHAIN = 21            				  	 # chain ID
RES_SEQ = slice(22, 26)   # residue sequence number

# the 20 standard amino acids, 3-letter codes
AA_CODES_3L = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE',
			   'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']


def fwrite(fp, dat, text="^", shortcut=True):
	"""Writes to a file.

	Args:
		fp (str): file path.
		dat: data to write
		text (str, optional): text to print. Defaults to "^".
		shortcut (bool, optional): whether to write to dat/fp. Defaults to True.
	"""

	fp = ("dat/" if shortcut else "") + fp
	with open(f"{fp}", 'w') as f:
		print(f"writing to {fp}: {text}\n")
		f.write(f"{dat}")
