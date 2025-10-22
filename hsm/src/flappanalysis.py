from collections import defaultdict
import json


NRES = 6
MDIST_MIN = 0.6
MDIST_MAX = 0.4
MDIST_SEQ = 1
target = "motif"

inf = "hsm/tools/" + \
		("MAPP-3D/MultipleSiteAlignment" if target == "motif" else "FLAPP") + \
		"/align_output.txt"
suffix = "2" if target == "motif" else ""

with open(inf) as f:

	sites = set()
	sitegrps = defaultdict(list)
	for l in f:
		try:
			dat = l.split("\t")
			scores = dat[2].split(" ")
			
			if target == "motif":
				mdist_min = float(scores[2])
				mdist_max = float(scores[3])
				mdist_seq = float(scores[4])
				nres = int(dat[2].split('/')[0])
			else:
				mdist_min = float(scores[3])
				mdist_max = float(scores[4])
				nres = int(scores[0])
		except:
			continue

		if mdist_min > MDIST_MIN and mdist_max > MDIST_MAX and mdist_seq > MDIST_SEQ and nres >= NRES:
			s0 = dat[0]
			s1 = dat[1]
			if "HSM" not in s0:
				print("what", l)
				break

			sites.add(s1)
			sitegrps[s0].append(s1)
			
			als = dat[3].strip().split("_")
			
			print(any([al.split(" ")[0][:3] != al.split(" ")[1][:3] for al in als]))

# sitelist = "\n".join([f"{k}: {v}" for k, v in sitegrps.items()])
print(len(sites))

json.dump(sitegrps, open(f"hsm/outs/FLAPP/alignedsites{suffix}.txt", 'w'))

with open(f"hsm/outs/FLAPP/alignlist{suffix}.txt", 'w') as f3:
	f3.write('\n'.join(sites))
	f3.write("\n")
