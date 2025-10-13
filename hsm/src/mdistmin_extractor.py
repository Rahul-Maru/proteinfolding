import csv

MD_CUTOFF = 0.4
RES_CUTOFF = 4
REPR = "2x45_C6.pdb"
mode = "directed"
target = "flapp"


def mdistmin_extractor():
	inf = "hsm/tools/" + \
		("MAPP-3D/MultipleSiteAlignment" if target == "sitemotif" else "FLAPP") + \
		"/align_output.txt"
	outdir = "hsm/outs/" + ("SiteMotif" if target == "sitemotif" else "FLAPP")


	final_list = [["Source", "Target", "Score"]]
	with open(inf) as f:
		pairs = f.readlines()
		aligns = []
		for pair in pairs:
			try:
				dat = pair.split("\t")
				scores = dat[2].split(" ")
				if target == "sitemotif":
					mdist_min = float(scores[2])
					nres = int(dat[2].split('/')[0])
				else:
					mdist_min = float(scores[3])
					nres = int(scores[0])


				if mode == "filter":
					if mdist_min > MD_CUTOFF and nres >= RES_CUTOFF:
						if dat[0] < dat[1]:
							final_list.append([dat[0], dat[1], mdist_min])
				elif mode == "directed":
					if mdist_min > MD_CUTOFF and nres >= RES_CUTOFF:
						if dat[0] != dat[1] :
							final_list.append([dat[0], dat[1], mdist_min])
				elif mode == "repr":
					if dat[0] == REPR:
						if nres >= RES_CUTOFF:
							print(pair)
							aligns += ['\n'.join([dat[1], *dat[-1].split('_')])]
							final_list.append([dat[0], dat[1], mdist_min])
							
			except:
				# If the alignment is None, move on
				continue
	
	print(''.join([f"{x}\n" for x in final_list]))

	with open(f"{outdir}/mdist.csv", 'w') as f2:
		writer = csv.writer(f2)
		writer.writerows(final_list)

	if mode == "repr":
		with open(f"{outdir}/aligns.txt", 'w') as f3:
			print(aligns)
			f3.writelines(aligns)


if __name__ == "__main__":
	mdistmin_extractor()
