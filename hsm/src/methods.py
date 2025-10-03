from collections import defaultdict
import csv
from matplotlib import pyplot as plt
import numpy as np
import os
import subprocess
import timeit

from bio import Protein, MAX_BSITE_DIST


def timed(func):
	"""Decorator that times a function."""

	def wrapper(*args, **kwargs):
		t = timeit.Timer(lambda: func(*args, **kwargs))
		print(f"done in {t.timeit(1):.3f}s")
	return wrapper


@timed
def chain_getter():
	print("Extracting chains with binding sites")
	with open('hsm/outs/PDB2Fasta/bchain_identifier.fa') as f:
		chains = [(x[1:-3], x[-2]) for x in f.readlines()[::2]]

	for c in chains:
		with open(f'hsm/bchains/{c[0]}_{c[1]}.pdb', 'w') as f2:
			prot = Protein(f'hsm/pdbs/{c[0]}.pdb')
			chain = prot.chains[c[1]]

			f2.writelines(chain)

@timed
def csv_formatter(mode, cutoff = 0):
	if mode == "seq":
		# sequence similarity network

		with open("hsm/outs/Clustal/sequence_identity_matrix.txt") as f:
			lines = f.readlines()
			labels = [x[:6] for x in lines[1:]]

			rows = [row[9:-1].split(' ') for row in lines[1:]]

		out = [["Source", "Target", "Score"]]
		for i, (row, prot) in enumerate(zip(rows, labels)):
			for score, prot2 in zip(row, labels[:i]):
				if (s := 1 - float(score)) >= cutoff:
					out.append([prot, prot2, s])

		print(len(out) - 1)

		with open("hsm/outs/Clustal/seq_edge_list.csv", "w") as f2:
			writer = csv.writer(f2)
			writer.writerows(out)

	elif mode == "struct":
		# structural similarity network

		with open("hsm/outs/TMalign/structure_identity_matrix.csv") as f:
			reader = csv.reader(f)
			labels = next(reader)[1:]

			rows = [row[1:] for row in reader]

		out = [["Source", "Target", "Score"]]
		for i, (row, prot) in enumerate(zip(rows, labels)):
			for score, prot2 in zip(row, labels[:i]):
				if float(score) >= cutoff:
					out.append([prot, prot2, score])

		print(len(out) -1)

		with open("hsm/outs/TMalign/struct_edge_list.csv", "w") as f2:
			writer = csv.writer(f2)
			writer.writerows(out)

def choose_reprs():
	with open("hsm/outs/Clustal/seq_edge_list.csv") as f:
		reader = csv.reader(f)
		next(reader)

		c_ids = defaultdict(lambda:-1)
		clusters = []
		c = 0

		for row in reader:
			si, sj, _ = row
			if c_ids[si] != -1 and c_ids[sj] != -1:
				ci = c_ids[si]
				cj = c_ids[sj]

				if ci == cj:
					continue

				print(ci, cj, si, sj, clusters[ci], clusters[cj])

				# merge sj's cluster ino si's clusters
				for s in clusters[cj]:
					c_ids[s] = ci
				print(clusters[ci], clusters[cj], '\n')

				if type(clusters[ci]) == list:
					clusters[ci].extend(clusters[cj])
					clusters[cj] = []

			# add the unclustered site to the cluster of the clustered site
			elif c_ids[si] != -1:
				ci = c_ids[si]
				c_ids[sj] = ci
				clusters[ci].append(sj)
			elif c_ids[sj] != -1:
				cj = c_ids[sj]
				c_ids[si] = cj
				clusters[cj].append(si)
			else:
				c_ids[si] = c
				c_ids[sj] = c
				clusters.append([si, sj])
				c += 1
		
		clusters = [c for c in clusters if c]
		print(len(clusters), print([f"{i} - {len(c)}" for i, c in enumerate(clusters)]))


@timed
def extract_bsites(combined=False, dist=MAX_BSITE_DIST):
	"""Extracts the binding sites of a given ligand from a PDB file
	and stores them in separate files.

	Args:
		combined (bool, optional): #whether to consider all ligand molecules together or separately.
			Defaults to False.
	"""

	print(f"Extracting Binding Sites ({"Combined" if combined else "Ligand-wise"})")

	# Get PDB IDs from the directory
	pdb_ids = [f for f in os.listdir("hsm/pdbs")]
	pdb_ids.sort()

	for p in pdb_ids:
		prot = Protein(f"hsm/pdbs/{p}")
		if combined:
			ligs = prot.get_ligand('HSM')
			bsites = [prot.get_bsite(lig, True, dist) for lig in ligs]

			combined_bsite = []
			for site in bsites:
				for res in site:
					if res not in combined_bsite:
						# combines all binding sites into one, removing duplicates
						combined_bsite.append(res)

			combined_bsite.sort(key=lambda x: ((i:=x['id'].split('_'))[0], int(i[1])))

			with open(f"hsm/bsites_combined/{p[:-4]}.pdb", 'w') as f:
				f.writelines(sum([res['atoms'] for res in combined_bsite], []))

		else:
			# create a separate file for the binding site of each chain
			# TODO ligands in chains with no other atoms not recognized

			ligs = prot.get_ligand('HSM')
			bsites = [prot.get_bsite(lig, True, dist) for lig in ligs]

			for i, lig in enumerate(ligs):
				if len(lig["atoms"]) == 0:
					continue

				with open(f"hsm/bsites/{p[:-4]}_{lig['id'][0]}{i}.pdb", 'w') as f:
					f.writelines(sum([res['atoms'] for res in bsites[i]], []))

@timed
def filter_bsites(src="hsm/bsites", ref="hsm/bchains_final", dst="hsm/bsites_final"):
	"""Moves binding sites from the source directory
	to the destination directory if their chain is present in the reference directory."""

	print(f"Filtering binding sites from \033[1m{src}\033[0m to \033[1m{dst}\033[0m if their chain is in \033[1m{ref}\033[0m")

	bsites = os.listdir(src)
	bsites.sort()

	reference = os.listdir(ref)
	reference.sort()

	i = 0
	# for every binding site in src, check if its chain is in reference
	# if so, copy it to dst
	for s in bsites:
		if s[:6] + s[7:] in reference:
			print(s)
			i += 1
			subprocess.run(['cp', f'{src}/{s}', f'{dst}/{s}'])
	print(f"final count: {i}")

@timed
def load_sitemotif_files(path):
    """
    Moves the given directory to the sitemotif directory.
    """

    try:
        subprocess.run(["rm", "-r", "hsm/tools/MAPP-3D/MultipleSiteAlignment/hsm/"],)
    except:
        pass

    subprocess.run(["cp", "-r", path, "hsm/tools/MAPP-3D/MultipleSiteAlignment/hsm/"])

@timed
def mdistmin_extractor():
	final_list = [["Source", "Target", "Score"]]
	with open("hsm/tools/MAPP-3D/MultipleSiteAlignment/align_output.txt") as f:
		pairs = f.readlines()
		for pair in pairs:
			try:
				if (mdist_min := float((dat:=pair.split("\t"))[2].split(" ")[2])) >= 0:
					if dat[0] < dat[1]:
						final_list.append([dat[0], dat[1], mdist_min])
			except:
				continue

	with open("hsm/outs/SiteMotif/mdist.csv", 'w') as f2:
		writer = csv.writer(f2)
		writer.writerows(final_list)

@timed
def pdbs2fasta(dir, out):
	outf = f"hsm/outs/PDB2Fasta/{out}.fa"
	print(f"Converting {dir} → {outf}")

	CMD = './hsm/tools/PDB2Fasta/pdb2fasta.sh'

	files = os.listdir(dir)
	files.sort()

	fasta = ''.join([subprocess.run([CMD, f'{dir}/{f}'], capture_output=True, text=True).stdout for f in files])

	with open(outf, "w") as f:
		f.write(fasta)

def scoreplotter(mode):
	inf = "Clustal/seq_edge_list" if mode == "seq" else "TMalign/struct_edge_list"

	# TODO more efficient way to do this
	with open(f"hsm/outs/{inf}.csv") as f:
		r = csv.reader(f)
		next(r)
		dat = [i[2] for i in r]

	n = len(dat)
	x = np.linspace(0, 1, 101)
	y = [len([j for s in dat if (j:=float(s))<=i])/n for i in x]

	plt.plot(x, y,label=f"cumulative frequency of {'sequence' if mode == 'seq' else 'structural'} identity scores")

	plt.xlabel("sequence identity score")
	plt.ylabel("cumulative frequency")
	plt.title(f"Cumulative frequency of {'sequence' if mode == 'seq' else 'structural'} identity scores")

	plt.show()

@timed
def tmalign():
	dir = "hsm/bchains_seq"
	outf = "structure_identity_matrix"

	print(f"Aligning {dir} using TMalign → {outf}.csv")

	prots = os.listdir(dir)
	prots.sort()

	mat = [["."] + [id[:-4] for id in prots]]
	for p1 in prots:
		mat.append([p1])
		for p2 in prots:
			out = subprocess.run(["hsm/tools/TMalign/TMalign_cpp", "-a", "T", f"{dir}/{p1}", f"{dir}/{p2}"],
						capture_output=True, text=True)
			try:
				x = float(out.stdout.split('\n')[15][9:17].strip())
				mat[-1].append(x)

			except:
				print("ERR", out.stdout)
				mat[-1].append(-1)

	with open(f"hsm/outs/TMalign/{outf}.csv", "w", newline="") as csvfile:
		writer = csv.writer(csvfile)
		writer.writerows(mat)

@timed
def weighted_sum():
	with open("hsm/outs/SiteMotif/mdist.csv") as f:
		reader = csv.reader(f)
		next(reader)

		prots = {}

		for row in reader:
			if row[0] not in prots:
				prots[row[0]] = 0
			if row[1] not in prots:
				prots[row[1]] = 0

			prots[row[0]] += float(row[2])
			prots[row[1]] += float(row[2])

		print(prots)
		print(max(prots, key=prots.get))
