from collections import defaultdict


NRES = 4
MDIST = 0.4

with open("hsm/tools/FLAPP/align_output.txt") as f:
	next(f)

	sites = set()
	sitegrps = defaultdict(list)
	for l in f:
		dat = l.split("\t")
		scores = dat[2].split(" ")
		mdist_min = float(scores[3])
		nres = int(scores[0])
		if mdist_min > MDIST and nres >= NRES:
			s0 = dat[0]
			s1 = dat[1]
			if "HSM" not in s0:
				print("what", l)
				break
			if "HSM" in s1:
				continue

			sites.add(s1)
			sitegrps[s0].append(s1)

sitelist = [f"{k}: {v}" for k, v in sitegrps.items()]
print("\n".join(sitelist))
print(len(sites))