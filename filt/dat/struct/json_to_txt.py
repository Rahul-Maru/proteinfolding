import json
import os

print(os.getcwd())

IN = "../reprs.json"
OUT = "cluster_reprs.txt"

reprs = json.load(open(IN))
with open(OUT, 'w') as f:
	f.writelines('\n'.join(reprs))