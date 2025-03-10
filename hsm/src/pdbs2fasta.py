from subprocess import run
from consts import bsites, pdb_ids_all

combined = True
CMD = './hsm/tools/PDB2Fasta/pdb2fasta.sh'

def main():
	if combined:
		fasta = ''.join([run([CMD, f'hsm/bsites_combined/{id}'], capture_output=True, text=True).stdout for id in pdb_ids_all])
	else:
		fasta = ''.join([run([CMD, f'hsm/bsites/{bs}'], capture_output=True, text=True).stdout for bs in bsites])

	print(fasta)
	with open("hsm/outs/PDB2Fasta/fasta.fa", "w") as f:
		f.write(fasta)

if __name__ == "__main__":
	main()