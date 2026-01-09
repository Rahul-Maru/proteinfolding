from pdbfixer import PDBFixer
from openmm.app import PDBFile
from time import time

inf = "2v30"
pdb_file = f'filt/dat/struct/pdb/{inf[1:3]}/pdb{inf}.ent'

t0 = time()
fixer = PDBFixer(pdb_file)
fixer.findMissingResidues()
fixer.findNonstandardResidues()
fixer.findMissingAtoms()
fixer.replaceNonstandardResidues()
fixer.removeHeterogens(False)
fixer.addMissingHydrogens(7.4)
PDBFile.writeFile(fixer.topology, fixer.positions, open('hsm/src/autodock/output.pdb', 'w'))
print(f"time taken: {time() - t0} s")