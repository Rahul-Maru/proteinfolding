import os


INDIR = "pdb"

for d in os.listdir(INDIR):
    for f in os.listdir(f"{INDIR}/{d}"):
        print(f)