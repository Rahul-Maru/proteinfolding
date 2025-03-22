import argparse
from methods import *
import timeit

#WIP TODO
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('mode', type=str, nargs='?', default="")

    args = parser.parse_args()
    mode = args.mode

    # sequence alignment
    extract_bsites(False)
    extract_bsites(True)
    pdbs2fasta("hsm/bsites_combined", "bsite")
    chain_getter()
    pdbs2fasta("hsm/bchains", "fasta")
    print("Running clustalo on fasta.fa")
    t = timeit.Timer(lambda: subprocess.run(['clustalo', '-i', 'hsm/outs/PDB2Fasta/fasta.fa', '-o',
                    'hsm/outs/Clustal/out.txt', '--distmat-out=hsm/outs/Clustal/mat.txt', '--full', '--force']))
    print(f"done in {t.timeit(1):.3f}s\n")

    print("Creating plot of similarity scores (close plot to continue)")
    csvformatter("seq", 0)
    scoreplotter("seq")
    cutoff = float(input("cutoff: "))
    print("done\n")

    csvformatter("seq", cutoff)

    # structural alignment


if __name__ == "__main__":
    main()