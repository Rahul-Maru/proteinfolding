"""Function to extract the resolution of a given"""

import os
import re

PDBDIR = "filt/dat/struct/pdb"

def get_res(pdb, dir=PDBDIR):
	"""Extracts the resolution of a given PDB file"""

	# regex patterns
	RESPATTERN = re.compile(r'REMARK\s{3}2\sRESOLUTION\.\s+(.*)\sANGSTROMS\.')
	REMARKxPATTERN = re.compile(r'REMARK\s{3}[3-5]')

	with open(f'{dir}/{pdb}', 'r') as f:
		for line in f:
			# the remark [x] section is after the resolution section
			# if we get here, then the entry has no resolution
			if REMARKxPATTERN.match(line):
				raise ValueError(f"Resolution not found for {pdb} (0)")

			pattern = RESPATTERN.match(line)
			if pattern:
				try:
					res = float(pattern.group(1))
				except:
					# resolution not defined
					raise ValueError(f"Resolution {pattern.group(1)} for {pdb}")

				break

	return res

def main():
	# gets the resolutions for all proteins in the pdb directory
	for dir in os.listdir(PDBDIR):
		for pdb in os.listdir(f'{PDBDIR}/{dir}'):
			print(f"{pdb}: ", end="")
			print(get_res(f'{dir}/{pdb}'))


if __name__ == "__main__":
	main()