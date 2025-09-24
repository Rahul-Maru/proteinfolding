from multiprocessing import Value
import os
import re
from sys import float_info

PDBDIR = "filt/dat/struct/pdb"

def get_res(pdb, dir=PDBDIR):
	RESPATTERN = re.compile(r'REMARK\s{3}2\sRESOLUTION\.\s+(.*)\sANGSTROMS\.')
	REMARKxPATTERN = re.compile(r'REMARK\s{3}[3-5]')
	with open(f'{dir}/{pdb}', 'r') as f:
		for line in f:
			if REMARKxPATTERN.match(line):
				raise ValueError(f"Resolution not found for {pdb} (0)")

			pattern = RESPATTERN.match(line)
			if pattern:
				try: 
					res = float(pattern.group(1))			
				except:
					raise ValueError(f"Resolution {pattern.group(1)} for {pdb}")

				break
		else:
			raise ValueError(f"Resolution not found for {pdb} (1)")

	return res

def main():
	for dir in os.listdir(PDBDIR):
		for pdb in os.listdir(f'{PDBDIR}/{dir}'):
			print(f"{pdb}: ", end="")
			print(get_res(f'{dir}/{pdb}'))


if __name__ == "__main__":
	main()