import argparse
from methods import *
import timeit

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('mode', type=str, nargs='?', default="seq")

    args = parser.parse_args()
    mode = args.mode

    if mode == "seq":
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
        csv_formatter("seq")
        scoreplotter("seq")
        cutoff = float(input("cutoff: "))
        print("done\n")

        csv_formatter("seq", cutoff)

    elif mode == "struct":
        # structural alignment
        tmalign()

        csv_formatter("struct")
        scoreplotter("struct")
        cutoff = float(input("cutoff: "))
        print("done\n")

        csv_formatter("struct", cutoff)

if __name__ == "__main__":
    main()