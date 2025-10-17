"""BIOINFORMATICS LIBRARY."""

# Constants imported in Protein.py
from Protein import *


def fwrite(fp, dat, text = "^", shortcut=True):
	"""Writes to a file. (Defaults to filt/dat/fp)

	Args:
		fp (str): file path.
		dat: data to write
		text (str, optional): text to print. Defaults to "^".
		shortcut (bool, optional): whether to write to filt/dat/fp. Defaults to True.
	"""
	fp = ("filt/dat/" if shortcut else "") + fp
	with open(f"{fp}", 'w') as f:
		print(f"writing to {fp}: {text}\n")
		f.write(f"{dat}")
