"""Copies all filt/dat/struct/representative-binding-sites intto
filt/dat/struct/ligand-inclsv-binding-sites, appending to each site
its ligand data
"""

import os
import shutil

def main():
	INDIR = "representative-binding-sites"
	INDIR = "representative-binding-sites"
	OUTDIR = "ligand-inclsv-binding-sites"

	# make output directory if it doesn't exist
	if not os.path.exists(OUTDIR):
		print(OUTDIR)
		os.makedirs(OUTDIR)

	i = 0
	nf = []
	for f in os.listdir(INDIR):
		i+=1
		if not i%1000:
			print (i)

		inf = f"{INDIR}/{f}"
		outf = f"{OUTDIR}/{f}"

		# parses site name
		st = f[:-4].split("_")
		pdbid, ch, n, lig = st

		id = f"{pdbid}.ent"
		mid = id[4:6]

		found = False
		with open(f"pdb/{mid}/{id}") as fl:
			for l in fl:
				if l[:6] == "HETATM":
					# parses ligand info from record
					rch, rn, rlig = (l[21], l[22:26], l[17:20].strip())

					# checks if the identifying info matches
					if (ch, int(n), lig) == (rch, int(rn), rlig):
						# if this is the first matching line
						if not found:
							found = True
							# Copy the original file if output doesn't exist
							if not os.path.exists(outf):
								shutil.copy2(inf, outf)

						# Append the HETATM line to output file
						with open(outf, 'a') as f2:
							f2.write(l)

					elif found:
						# if the number doesn't match and we've already found the ligand, then we are done
						break
			else:
				if not found:
					nf.append(f)
					print(f"{i}: not found {pdbid} ({f}), {n}, {rn}")

	if nf:
		print(len(nf))
		with open("nf_lig.txt", 'w') as f:
			f.write("\n".join(nf))

if __name__ == "__main__":
	main()

