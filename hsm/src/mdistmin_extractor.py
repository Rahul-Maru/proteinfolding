import csv

MD_CUTOFF = 0
RES_CUTOFF = 4

def mdistmin_extractor():
	mode = "directed"
	final_list = [["Source", "Target", "Score"]]
	with open("hsm/tools/MAPP-3D/MultipleSiteAlignment/align_output.txt") as f:
		pairs = f.readlines()
		for pair in pairs:
			try:
				dat = pair.split("\t")
				mdist_min = float(dat[2].split(" ")[2])
				nres = int(dat[2].split('/')[0])

				if mode == "filter":
					if mdist_min > MD_CUTOFF and nres >= RES_CUTOFF:
						if dat[0] < dat[1]:
							final_list.append([dat[0], dat[1], mdist_min])
				elif mode == "directed":
					if mdist_min > MD_CUTOFF and nres >= RES_CUTOFF:
						if dat[0] != dat[1] :
							final_list.append([dat[0], dat[1], mdist_min])
			except:
				# If the alignment is None, move on
				continue
	
	with open("hsm/outs/SiteMotif/mdist.csv", 'w') as f2:
		writer = csv.writer(f2)
		writer.writerows(final_list)



if __name__ == "__main__":
	mdistmin_extractor()
