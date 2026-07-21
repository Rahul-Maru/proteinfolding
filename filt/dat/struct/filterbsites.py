import os
from bio import RESN, AA_CODES_3L, CHAIN, RES_SEQ

def filter_site(lines):
	return [l for l in lines if l.startswith("ATOM") and l[RESN].strip() in AA_CODES_3L]


def count_aa_residues(lines):
	seen = set()
	for line in lines:
		if not line.startswith("ATOM"):
			raise ValueError(f"Line does not start with ATOM:\n{line}")
		seen.add((line[CHAIN], line[RES_SEQ]))

	return len(seen)


def main():
	if os.getcwd().split("/")[-1] != "struct":
		print("Error: must be run from the struct/ directory")
		exit(1)

	bsites_base = 'binding-sites'
	out_base = 'filtered-binding-sites'
	unwanted_ligs = 'unwanted-ligs'

	with open(unwanted_ligs) as f:
		unwanted = set(line.strip() for line in f)

	os.makedirs(out_base, exist_ok=True)

	# make subdirs if missing
	if not any(os.path.isdir(os.path.join(out_base, e)) for e in os.listdir(out_base)):
		for d in sorted(os.listdir(bsites_base)):
			if os.path.isdir(os.path.join(bsites_base, d)):
				os.makedirs(os.path.join(out_base, d), exist_ok=True)

	for d in sorted(os.listdir(out_base)):
		if not os.path.isdir(os.path.join(bsites_base, d)):
			continue

		print(f"d - {d}")
		outd = os.path.join(out_base, d)
		os.makedirs(outd, exist_ok=True)

		for fname in sorted(os.listdir(os.path.join(bsites_base, d))):
			if not fname.endswith(".pdb"):
				continue

			fpath = os.path.join(bsites_base, d, fname)
			outf = os.path.join(outd, fname)

			# get ligand from filename: name_chain_resnum_LIG.pdb
			lig = fname.split("_")[-1].removesuffix(".pdb")

			with open(fpath) as f:
				lines = f.readlines()

			filtered = filter_site(lines)

			if count_aa_residues(filtered) >= 4 and lig not in unwanted:
				with open(outf, "w") as f:
					f.writelines(filtered)
			else:
				if os.path.exists(outf):
					print(f"removing {outf}")
					os.remove(outf)

if __name__ == "__main__":
	main()
