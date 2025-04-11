import argparse
from methods import *
import timeit

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('mode', type=str, nargs='?', default="seq")

    args = parser.parse_args()
    mode = args.mode

    if mode == "seq":
        # ———sequence alignment———

        extract_bsites(False)
        extract_bsites(True)

        # identifies which chains have binding sites
        pdbs2fasta("hsm/bsites_combined", "bchain_identifier")
        chain_getter()

        # gets the residues of the chains with binding sites
        pdbs2fasta("hsm/bchains", "bchains")

        # sequence alignment of bchains.fa using clustal omega 
        print("Running clustalo on bchains.fa")
        t = timeit.Timer(lambda: subprocess.run(
            ['clustalo', '-i', 'hsm/outs/PDB2Fasta/bchains.fa', '-o',
            'hsm/outs/Clustal/bchains.txt', '--distmat-out=hsm/outs/Clustal/sequence_identity_matrix.txt',
            '--full', '--force']))
        print(f"done in {t.timeit(1):.3f}s\n")

        print("Creating plot of sequence similarity scores (close plot to continue)")
        csv_formatter("seq")
        scoreplotter("seq")
        try:
            cutoff = float(input("cutoff: "))
        except:
            print("\033[91mInvalid cutoff, program terminating\033[0m")
            return

        print("done\n")

        # creates a csv file of the sequence similarity scores
        csv_formatter("seq", cutoff)

        print("\n——INPUT \033[1mseq_edge_list.csv\033[0m INTO CYTOSCAPE——")

    elif mode == "struct":
        # ———structural alignment———

        # all vs all structural alignment of bchain_seq/ using TMalign
        tmalign()

        csv_formatter("struct")
        scoreplotter("struct")

        print("Creating plot of structural similarity scores (close plot to continue)")
        try:
            cutoff = float(input("cutoff: "))
        except:
            print("\033[91mInvalid cutoff, program terminating\033[0m")
            return

        print("done\n")

        # creates a csv file of the structural similarity scores
        csv_formatter("struct", cutoff)

        print("\n——INPUT \033[1mstruct_edge_list.csv\033[0m INTO CYTOSCAPE——")
    elif mode == 'filter':
        filter_bsites()


if __name__ == "__main__":
    main()