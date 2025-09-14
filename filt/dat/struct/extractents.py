import os
import re


INDIR = "filt"
OUTF = "ents2.txt"

SOURCE_PATTERN = re.compile(r"SOURCE")
ENTY_PATTERN = re.compile(r"COMPND\s+\d*\s*MOL_ID:\s*(\d+)")
MOL_PATTERN = re.compile(r"COMPND\s+\d*\s*MOLECULE:\s*(.*)")
CHAIN_PATTERN = re.compile(r"COMPND\s+\d*\s*CHAIN:(.*)")

ents=[]

for d in os.listdir(INDIR):
	for nam in os.listdir(f"{INDIR}/{d}"):
		with open(f"{INDIR}/{d}/{nam}") as f:
			for line in f:
				if SOURCE_PATTERN.match(line):
					print(f"FAILED - {INDIR}/{d}/{nam}")
					break

				if (ent := ENTY_PATTERN.match(line)):
					enty_n = int(ent.group(1))
					ents.append(f"{nam}_{enty_n}")

with open(OUTF) as f:
	f.write('\n'.join(ents))