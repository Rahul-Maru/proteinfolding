"""Copies all representative binding sites from filt/dat/struct/representative-binding-sites to
filt/dat/struct/ligand-inclsv-binding-sites, appending to each the HETATM records associated with the
ligand they bind to.
"""

import os
import shutil
from pathlib import Path

def main():
	INDIR = "representative-binding-sites"
	OUTDIR = "ligand-inclsv-binding-sites"

	if not os.path.exists(OUTDIR):
		print(OUTDIR)
		os.makedirs(OUTDIR)

	i = 0
	nf = []
	for f in os.listdir(INDIR):
		st = f[:-4].split("_")
		_, ch, n, lig = st

		found = False
		with open(f"{INDIR}/{f}") as fl:
			for l in fl:
				if l[:6] == "HETATM":
					rch, rn, rlig = (l[21], l[22:26], l[17:20])
					if (ch, int(n), lig) == (rch, int(rn), rlig):
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
			else:
				nf.append(f)
				print("lig not found: ", f)

	print(len(nf))
	with open("nf_lig.txt", 'w') as f:
		f.write("\n".join(nf))

if __name__ == "__main__":
	main()
