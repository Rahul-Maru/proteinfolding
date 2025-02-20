import csv
from subprocess import run

pdb_ids_full = ['1avn.pdb', '1ike.pdb', '1jqd.pdb', '1kar.pdb',
		   		'1np1.pdb', '1qft.pdb', '1qfv.pdb', '1u18.pdb',
		  		'2qeb.pdb', '2x45.pdb', '3bu1.pdb', '3g7x.pdb',
				'3rxh.pdb', '4xmf.pdb', '5rus.pdb', '6dyn.pdb',
				'6fu4.pdb', '6qfa.pdb', '7a5v.pdb', '7dfl.pdb',
				'7qn7.pdb', '7qn8.pdb', '7qn9.pdb', '7qnc.pdb',
				'7qnd.pdb', '7yfc.pdb', '8hn8.pdb', '8jt5.pdb',
				'8jxt.pdb', '8pok.pdb', '8pvb.pdb', '8tgk.pdb',
				'8ucn.pdb', '8yn2.pdb', '8yn3.pdb', '8yn4.pdb',
				'8yn5.pdb', '8yn9.pdb', '8yuu.pdb']

pdb_ids = ['1avn.pdb', '1ike.pdb', '1jqd.pdb', '1kar.pdb', '1np1.pdb', '1qft.pdb', '2qeb.pdb',
		   '2x45.pdb', '3bu1.pdb', '3rxh.pdb', '4xmf.pdb', '5rus.pdb', '6dyn.pdb', '6fu4.pdb',
		   '6qfa.pdb', '7dfl.pdb', '7yfc.pdb', '8jt5.pdb', '8pok.pdb', '8tgk.pdb', '8yuu.pdb']

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