import csv
from subprocess import run
from consts import pdb_ids

def main():
	mat = [["."] + [id[:-4] for id in pdb_ids]]
	for p1 in pdb_ids:
		mat.append([p1])
		for p2 in pdb_ids:
			out = run(["HSM/TMalign/TMalign_cpp", "-a", "T", f"HSM/pdbs/{p1}", f"HSM/pdbs/{p2}"], capture_output=True, text=True)
			try:
				x = float(out.stdout.split('\n')[15][9:17].strip())
				mat[-1].append(x)

			except:
				print("ERR", out.stdout)
				mat[-1].append(-1)
	
	with open("HSM/TMalign/out.csv", "w", newline="") as csvfile:
		writer = csv.writer(csvfile)
		writer.writerows(mat)


if __name__ == "__main__":
	main()