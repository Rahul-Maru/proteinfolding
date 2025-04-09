import argparse
import csv
import os
from subprocess import run

def tmalign():
	parser = argparse.ArgumentParser()
	parser.add_argument('dir', type=str, nargs='?', default="hsm/bchains_final")
	parser.add_argument('out', type=str, nargs='?', default="out2")

	args = parser.parse_args()
	dir = args.dir
	outf = args.out

	prots = os.listdir(dir)
	prots.sort()
	print(prots)

	mat = [["."] + [id[:-4] for id in prots]]
	i=0
	for p1 in prots:
		mat.append([p1])
		for p2 in prots:
			out = run(["hsm/tools/TMalign/TMalign_cpp", "-a", "T", f"{dir}/{p1}", f"{dir}/{p2}"], capture_output=True, text=True)
			print(f"done {i}")
			i += 1
			try:
				x = float(out.stdout.split('\n')[15][9:17].strip())
				mat[-1].append(x)

			except:
				print("ERR", out.stdout)
				mat[-1].append(-1)
	
	with open(f"hsm/outs/TMalign/{outf}.csv", "w", newline="") as csvfile:
		writer = csv.writer(csvfile)
		writer.writerows(mat)


if __name__ == "__main__":
	tmalign()