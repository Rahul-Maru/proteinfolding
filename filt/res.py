import os
import re
from sys import float_info

PDBDIR = "filt/dat/struct/pdb"
def get_res(pdb):
	with open(f'{PDBDIR}/{pdb}', 'r') as f:
		lines = ''.join(f.readlines())
		# Find resolution in format "RESOLUTION.    xxx ANGSTROMS."
		resolution_match = re.search(r'RESOLUTION\.\s{4}(\d+\.\d{2})\sANGSTROMS\.', lines)
		if resolution_match:
			res = float(resolution_match.group(1))
		else:
			print("Resolution not found")
			res = float_info.max

	return res

def main():
	for dir in os.listdir(PDBDIR):
		for pdb in os.listdir(f'{PDBDIR}/{dir}'):
			print(f"{pdb}: ", end="")
			print(get_res(f'{dir}/{pdb}'))

if __name__ == "__main__":
	main()