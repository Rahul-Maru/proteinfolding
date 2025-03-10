from Protein import Protein
from consts import pdb_ids_all

def main():
    with open('hsm/outs/PDB2Fasta/fasta.fa') as f:
        chains = [(x[1:-3], x[-2]) for x in f.readlines()[::2]]
    print(chains)

    for c in chains:
        print(c)
        with open(f'hsm/bchains/{c[0]}_{c[1]}.pdb', 'w') as f2:
            prot = Protein(f'hsm/pdbs/{c[0]}.pdb')
            chain = prot.chains[c[1]]

            f2.writelines(chain)
if __name__ == "__main__":
    main()
