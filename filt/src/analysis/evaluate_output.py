import json
import os
from bio import CHAIN, RES_SEQ

INDIR = "dat/struct/ligand-inclsv-binding-sites/"
OUTF = "dat/repr_error.txt"

def main():
	print(f"---evaluating representative binding sites in {INDIR}---")
	print(f"total number of b-sites: {len(os.listdir(INDIR))}")
	print("--------")

	with open("dat/struct/unwanted-ligs") as f:
		unwanted_ligs = [lig.strip() for lig in f.readlines()]
	if not unwanted_ligs:
		raise ValueError("`unwanted-ligs` is empty.")

	error = {"bad_lig": [], "bad_records": [], "too_few_residues": []}
	for i, f_nam in enumerate(os.listdir(INDIR)):
		if i % 10000 == 0:
			print(f"{i} sites assessed")
	
		f_path = os.path.join(INDIR, f_nam)
		if not os.path.isfile(f_path):
			raise RuntimeError(f"unexpected non-file entry in {INDIR}: {f_nam}")

		lig = f_nam.split('_')[:-4]
		if lig in unwanted_ligs:
			print(f"binding site {f_nam} binds to unwanted ligand {lig}")
			error["bad_lig"].append(f_nam)

		with open(f_path) as f:
			seen = set()
			for line in f:
				if not (line.startswith("ATOM") or line.startswith("HETATM")):
					print(f"binding site {f_nam} contains non-ATOM/HETATM records")
					error["bad_records"].append(f_nam)
					break
				elif line.startswith("ATOM"):
					seen.add((line[CHAIN], line[RES_SEQ]))
				elif line.startswith("HETATM"):
					pass

			if len(seen) < 4:
				print(f"binding site {f_nam} contains too few residues")
				error["too_few_residues"].append(f_nam)

	n = sum([len(error[k]) for k in error])

	print(f"done. {n} problematic sites found. see {OUTF}")

	json.dump(error, open(OUTF, 'w'))



if __name__ == "__main__":
	main()