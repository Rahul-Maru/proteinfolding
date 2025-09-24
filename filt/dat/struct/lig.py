"""Copies all representative binding sites from filt/dat/struct/representative-binding-sites to
filt/dat/struct/ligand-inclsv-binding-sites, appending to each the HETATM records associated with the
ligand they bind to.
"""

import os
import shutil
from pathlib import Path

def main():
	INDR = "representative-binding-sites"
	OUTDIR = "ligand-inclsv-binding-sites"

	if not os.path.exists(OUTDIR):
		print(OUTDIR)
		os.makedirs(OUTDIR)

	i = 0
	nf = []
	for file_path in Path(INDR).iterdir():
		nam = file_path.name
		ligArr = nam.split('_')
		n = ligArr[2]
		og_nam = f"{ligArr[0]}.ent"
		subdir = og_nam[4:6]

		pdb_file = Path(f"pdb/{subdir}/{og_nam}")

		found = False
		output_file = Path(OUTDIR) / nam

		with open(pdb_file, 'r') as pdb_f:
			for line in pdb_f:
				rec = line[:6].strip()
				if rec == "HETATM":
					num_str = line[22:26].strip()
					num = f"{int(num_str):04d}" if num_str else ""

					# found the match
					if num == n:
						# if this is the first matching line
						if not found:
							found = True
							i += 1
							# Copy the original file if output doesn't exist
							if not output_file.exists():
								shutil.copy2(file_path, output_file)

						# Append the HETATM line to output file
						with open(output_file, 'a') as out_f:
							out_f.write(line)
					elif found:
						# if the number doesn't match and we've already found the ligand, then we are done
						break
		if not found:
			print(f"{i}: not found {nam} ({pdb_file}), {n}, {num}")
			nf.append(nam)
		elif not i % 1000:
			print(i)

	print(len(nf))
	with open("nf_lig.txt", 'w') as f:
		f.write("\n".join(nf))

if __name__ == "__main__":
	main()
