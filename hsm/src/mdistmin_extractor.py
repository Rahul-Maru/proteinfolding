import csv

def main():
	final_list = [["Source", "Target", "Score"]]
	with open("hsm/tools/MAPP-3D/MultipleSiteAlignment/align_output.txt") as f:
		pairs = f.readlines()
		for pair in pairs:
			try:
				if (mdist_min := float((dat:=pair.split("\t"))[2].split(" ")[2])) >= 0.5:
					if dat[0] < dat[1]:
						final_list.append([dat[0], dat[1], mdist_min])
			except:
				continue
	
	with open("hsm/outs/SiteMotif/mdist.csv", 'w') as f2:
		writer = csv.writer(f2)
		writer.writerows(final_list)



if __name__ == "__main__":
	main()
