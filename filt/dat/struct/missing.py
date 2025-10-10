import os

indir = "ligand-inclsv-binding-sites"

nf = []

for f in os.listdir(indir):
	st = f[:-4].split("_")
	_, ch, n, lig = st

	found = False
	with open(f"{indir}/{f}") as fl:
		for l in fl:
			if l[:6] == "HETATM":
				rch, rn, rlig = (l[21], l[22:26], l[17:20])
				if (ch, int(n), lig) == (rch, int(rn), rlig):
					found = True
					break
		else:
			nf.append(f)
			print("lig not found: ", f)

with open("wrong_ligs.txt", 'w') as f2:
	f2.write('\n'.join(nf))
	
