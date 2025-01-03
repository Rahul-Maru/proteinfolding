import csv
from subprocess import run

pdb_ids = ['1avn.pdb', '1ike.pdb', '1jqd.pdb', '1kar.pdb', '1np1.pdb',
		   '1qft.pdb', '2qeb.pdb', '2x45.pdb', '3bu1.pdb', '3rxh.pdb',
		   '4xmf.pdb', '5rus.pdb', '6dyn.pdb', '6fu4.pdb', '7dfl.pdb',
		   '7yfc.pdb', '8jt5.pdb', '8pok.pdb', '8tgk.pdb', '8yuu.pdb']

def main():
	mat = []
	for p1 in pdb_ids:
		mat.append([])
		for p2 in pdb_ids:
			out = run(["./HSM/TMalign_cpp", "-a", "T", f"HSM/pdbs/{p1}", f"HSM/pdbs/{p2}"], capture_output=True, text=True)
			try:
				x = float(out.stdout.split('\n')[15][9:17].strip())
				mat[-1].append(x)

			except:
				print("ERR", out.stdout)
				mat[-1].append(-1)
	
	with open("out.csv", "w", newline="") as csvfile:
		writer = csv.writer(csvfile)
		writer.writerows(mat)


if __name__ == "__main__":
	main()