import csv

def main():
	with open("HSM/PDB2Fasta/mat.txt") as f:
		lines = f.readlines()
		labels = [x[:4] + x[5] for x in lines[1:]]

		rows = [row[7:-1].split(' ') for row in lines[1:]]
		# print (rows)
		# print(rows[0:4])
	
	out = [["Source", "Target", "Score"]]
	for i, (row, prot) in enumerate(zip(rows, labels)):
		for score, prot2 in zip(row, labels[:i]):
			if 1 - float(score) >= 1:
				out.append([prot, prot2, score])
	
	# print('\n'.join([' '.join(r) for r in out]))
	print(len(out) - 1)

	with open("HSM/PDB2Fasta/network.csv", "w") as f2:
		writer = csv.writer(f2)
		writer.writerows(out)

if __name__ == "__main__":
	main()