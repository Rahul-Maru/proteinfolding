import json
import os

out = {}
for F in os.listdir("hsm/outs/autodock/reports"):
		with open(F) as f:
			table = False
			pdb = ""
			interaction = ""
			for l in f:
				if 'pdb' in l:
					pdb = l[62:-14]
					print(pdb)
					out[pdb] = {}

				elif l.startswith('**'):
					table = False
					interaction = l[2:-2]
					print(interaction)
					out[pdb][interaction] = []

				elif l.startswith('+=='):
					table = True

				elif l.startswith('| ') and table:
					out[pdb][interaction].append(l)

json.dump(out, open("hsm/outs/autodock/reports.json", "w"))