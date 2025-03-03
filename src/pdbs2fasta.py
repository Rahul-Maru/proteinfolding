from subprocess import run
from consts import pdb_ids, pdb_ids_all

def main():
	fasta = ""
	for pdb_id in pdb_ids_all:
		out = run(['./HSM/PDB2Fasta/pdb2fasta.sh', f'HSM/pdbs/{pdb_id}'], capture_output=True, text=True)
		fasta += out.stdout
	with open("HSM/PDB2Fasta/fasta.fa", "w") as f:
		f.write(fasta)

if __name__ == "__main__":
	main()