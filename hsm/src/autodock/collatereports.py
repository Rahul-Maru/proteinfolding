"""
This script collates the PLIP reports into a single JSON file.
File should be formatted afterwards for readability.
"""
import json
import os

indir = "hsm/outs/autodock/reports"

out = {}

for F in os.listdir(indir):
		with open(f"{indir}/{F}") as f:
			table = False
			pdb = ""
			interaction = ""
			for l in f:
				if 'PDB' in l:
					pdb = l[62:-14]
					print(l)
					out[pdb] = {}

				elif l.startswith('**'):
					table = False
					interaction = l[2:-3]
					print(interaction)
					out[pdb][interaction] = []

				elif l.startswith('+=='):
					table = True

				elif l.startswith('| ') and table:
					out[pdb][interaction].append(l)

json.dump(out, open("hsm/outs/autodock/reports.json", "w"))
