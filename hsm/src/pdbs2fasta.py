from subprocess import run
from consts import bsites

def main():
	fasta = ""

	for bs in bsites:
		out = run(['./hsm/tools/PDB2Fasta/pdb2fasta.sh', f'hsm/bsites/{bs}'], capture_output=True, text=True)
		fasta += out.stdout
	print(fasta)
	with open("hsm/outs/PDB2Fasta/fasta.fa", "w") as f:
		f.write(fasta)

if __name__ == "__main__":
	main()