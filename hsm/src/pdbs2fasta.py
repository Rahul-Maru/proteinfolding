import argparse
import os
from subprocess import run

combined = True
CMD = './hsm/tools/PDB2Fasta/pdb2fasta.sh'

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('dir', type=str, nargs='?', default="hsm/bsites")
	parser.add_argument('out', type=str, nargs='?', default="bsite")

	args = parser.parse_args()
	dir = args.dir
	out = args.out

	files = os.listdir(dir)
	files.sort()
	print()

	fasta = ''.join([run([CMD, f'{dir}/{f}'], capture_output=True, text=True).stdout for f in files])

	print(fasta)
	with open(f"hsm/outs/PDB2Fasta/{out}.fa", "w") as f:
		f.write(fasta)

if __name__ == "__main__":
	main()