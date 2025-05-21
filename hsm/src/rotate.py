import subprocess
import os
import bio

QUERIES = ['2x45_C6', '2x45_C8', '3bu1_A0', '4xmf_A0', "1u18_A0", "6dyn_A0"]
TEMPLATE = '2x45_C7_HSM.mol2'
IN = 'hsm/ligtest'
OUT = 'hsm/outs/LSalign'

def main():

	for q in QUERIES:
		q_path = f"hsm/bsites/{q}.pdb"
		prot = bio.Protein(q_path, [0, 1])

		if not os.path.isfile(f"{IN}/{q}_HSM.pdb"):
			print(f"extracting HSM to ligtest from {q}")
			with open(f"{IN}/{q}_HSM.pdb", 'w') as hf:
				hf.writelines(prot.hetatms)
		if not os.path.isfile(f"{IN}/{q}_HSM.mol2"):
			# print("hi")
			subprocess.run(["obabel", "-i", "pdb", f"{IN}/{q}_HSM.pdb", "-o", "mol2", "-O" f"{IN}/{q}_HSM.mol2"])

		subprocess.run(["./hsm/tools/LSalign/src/LSalign", f"{IN}/{q}_HSM.mol2",
				  f"{IN}/{TEMPLATE}", "-m", (outpath:=f"{OUT}/mat_{q}"), "-rf", "1"])

		with open(outpath) as f:
			lines = [[c for c in l.split() if c][1:] for l in f.readlines()][1:4]
			t, u = [], []
			for l in lines:
				t.append(float(l[0]))
				u.append([float(n) for n in l[1:]])
		
		coords = prot.get_coords(triplet=1)
		# print(prot.lines)

		coords_new = [[sum([coord[j]*u[i][j] for j in range(3)]) + t[i] for i in range(3)] for coord in coords]
		# print(coords_new[0][0])
		coord_strings = [''.join([f"{c:8.3f}" for c in coord]) for coord in coords_new]
		lines_new = [lin[:30] + cs + lin[54:] for cs, lin in zip(coord_strings, prot.lines)]
		# print(''.join(lines_new))

		with open(f"{OUT}/{q}_HSMr.pdb", 'w') as f2:
			f2.writelines(lines_new)

		

if __name__ == "__main__":
	main()