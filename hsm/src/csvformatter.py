import csv

def csvformatter():
	with open("hsm/outs/TMalign/out2.csv") as f:
		reader = csv.reader(f)
		labels = next(reader)[1:]

		rows = [row[1:] for row in reader]
	
	out = [["Source", "Target", "Score"]]
	for i, (row, prot) in enumerate(zip(rows, labels)):
		for score, prot2 in zip(row, labels[:i]):
			if float(score) >= 0.7:
				out.append([prot, prot2, score])
	
	print('\n'.join([' '.join(r) for r in out]))
	print(len(out) -1)

	with open("hsm/outs/TMalign/network2.csv", "w") as f2:
		writer = csv.writer(f2)
		writer.writerows(out)

if __name__ == "__main__":
	csvformatter()