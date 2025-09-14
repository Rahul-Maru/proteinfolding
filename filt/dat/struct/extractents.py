import os
import re


INDIR = "pdb"
OUTF = "ents2.txt"

SOURCE_PATTERN = re.compile(r"SOURCE")
ENTY_PATTERN = re.compile(r"COMPND\s+\d*\s*MOL_ID:\s*(\d+)")

entys = []
failures = []
i=0

for d in os.listdir(INDIR):
	for nam in os.listdir(f"{INDIR}/{d}"):
		with open(f"{INDIR}/{d}/{nam}", "r") as f:
			c = 0
			for line in f:
				if SOURCE_PATTERN.match(line):
					if not c:
						print(f"FAILED - {INDIR}/{d}/{nam}")
						failures.append(nam[3:7])
					break

				if (ent := ENTY_PATTERN.match(line)):
					enty_n = int(ent.group(1))
					entys.append(f"{nam[3:7]}_{enty_n}")
					# print(f"{nam[:4]}_{enty_n}")
					c+=1
		i+=1
		if not i%2000:
			print(i)

if failures:
	with open("failedents.txt", "w") as f:
		f.write("\n".join(failures))

with open(OUTF, "w") as f:
	f.write("\n".join(entys))
