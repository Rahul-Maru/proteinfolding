import csv

def main():
	with open("hsm/outs/Clustal/mat.txt") as f:
		lines = f.readlines()
		labels = [x[:6] for x in lines[1:]]

		rows = [row[9:-1].split(' ') for row in lines[1:]]
		# print (rows)
		# print(rows[0:4])
	
	out = [["Source", "Target", "Score"]]
	for i, (row, prot) in enumerate(zip(rows, labels)):
		for score, prot2 in zip(row, labels[:i]):
			if (s := 1 - float(score)) >= 0.9:
				out.append([prot, prot2, s])
	
	# print('\n'.join([' '.join(r) for r in out]))
	print(len(out) - 1)

	with open("hsm/outs/Clustal/network.csv", "w") as f2:
		writer = csv.writer(f2)
		writer.writerows(out)

if __name__ == "__main__":
	main()