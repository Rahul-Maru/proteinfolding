"""BIO LIBRARY."""

# Constants imported in Protein.py
from Protein import *


def fwrite(fp, x, text = "^", dat=True):
	"""Easier way to write files"""
	fp = ("filt/dat" if dat else "") + fp
	with open(f"filt/dat/{fp}", 'w') as f:
		print(f"writing to {fp}: {text}\n")
		f.write(f"{x}")