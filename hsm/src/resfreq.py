import csv
import os
import bio
import matplotlib.pyplot as plt

def count_amino_acid_frequencies(pdb_dir):
	"""
	Count the frequency of each amino acid type in all PDB files in the given directory.
	Returns a dictionary mapping 3-letter amino acid codes to their counts.
	"""
	from collections import Counter
	aa_counts = Counter()
	i = 0

	for fname in os.listdir(pdb_dir):

		if not fname.endswith('.pdb'):
			continue

		p = bio.Protein(os.path.join(pdb_dir, fname))
		seq = p.seq().strip('-')
		print(p.id, seq)

		for aa in seq:
			if 65 <= ord(aa) <= 90:
				i += 1
				aa_counts[aa] += 1

	print(i)

	return dict(aa_counts), i

aa_groups = {'polar': 'STNQ',
			'non-polar': 'ILMGPVA',
			'pos-charged': 'KRH',
			'neg-charged': 'DE',
			'aromatic': 'FWY'}

aa_colors = {'polar': 'green',
			'non-polar': 'black',
			'pos-charged': 'blue',
			'neg-charged': 'red',
			'aromatic': 'purple'}

def main():


	pdb_dir = "hsm/bsites_final"
	aa_counts, n = count_amino_acid_frequencies(pdb_dir)
	aa_freq = {k: round(v/n, 3) for k, v in sorted(aa_counts.items(), key=lambda item: item[1], reverse=True)}
	print(aa_freq)

	with open("hsm/outs/resfreq.csv", 'w') as f:
		writer = csv.writer(f)
		writer.writerows(aa_counts.items())

	plt.bar(*zip(*aa_freq.items()), color=[get_color(aa) for aa in aa_freq.keys()])
	plt.ylabel("frequency")
	plt.xlabel("residue")
	plt.title(f"{n=} residues over {len(os.listdir(pdb_dir))} binding sites")
	plt.show()

def get_color(aa):
	for group, chars in aa_groups.items():
		if aa in chars:
			return aa_colors[group]
	return 'gray'

if __name__ == "__main__":
	main()