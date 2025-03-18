import argparse
from methods import *

#WIP TODO
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('mode', type=str, nargs='?', default="")

    args = parser.parse_args()
    mode = args.mode

    # sequence alignment
    extract_bsites()
    pdbs2fasta()
    subprocess.run(['clustalo', '-i', 'hsm/outs/PDB2Fasta/fasta.fa', '-o',
                    'hsm/outs/Clustal/out.txt', '--distmat-out=hsm/outs/Clustal/mat.txt', '--full', '--force'])

    # structural alignment


if __name__ == "__main__":
    main()