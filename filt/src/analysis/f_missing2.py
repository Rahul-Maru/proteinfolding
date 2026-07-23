""" SECOND VER (supersedes f_mising.py)
Script that investigates why some sites present in the filtered binding sites directory are
not in the final clustered list of binding sites.
"""

import json
import re
from bio import fwrite


CREATE_NONENTS = True

# regex patterns to parse PDB files
SOURCE_PATTERN = re.compile(r"SOURCE")
CHAIN_PATTERN = re.compile(r"COMPND\s+\d*\s*CHAIN:(.*)")
ENTY_PATTERN = re.compile(r"COMPND\s+\d*\s*MOL_ID:\s*(\d+)")
MOL_PATTERN = re.compile(r"COMPND\s+\d*\s*MOLECULE:\s*(.*)")

def main():
	print("––f_bsites.txt––")
	print(" - filtered list of all binding sites")

	with open("dat/struct/f_bsites.txt") as f:
		lins = f.readlines()
		f_bsites = {b.strip() for b in lins}
		print("# of lines w/o whitespace ( = number of binding sites): ", len(f_bsites))

	print("\n––clusters-by-bsite-70.json––")
	print("- clustered list of filtered binding sites")

	clust = json.load(open("dat/f-clusters-by-bsite-70.json"))
	clust_set = {site for c in clust for site in c}
	print("number of clustered sites: ", len(clust_set))

	print() #-----------------

	sites_m = f_bsites - clust_set

	print("filtered binding sites that were not clustered: ", len(sites_m))

	print("\n––nonentities.txt––")
	print(" - list of binding sites whose chains do not belong to a named entity")

	# if CREATE_NONENTS, manually create the set, otherwise read it from the file
	# 	to save time, as creating nonentities.txt is time-consuming
	if CREATE_NONENTS:
		nonentities = create_nonents(sites_m)
	else:
		print("\n––nonentities.txt––")
		print(" - list of binding sites whose chains do not belong to a named entity")
		print("note: this file is created in this script. Set CREATE_NONENTS to True to create it again.")

		with open("dat/nonentities.txt") as f:
			nonentities = set(f.readlines()[0][2:-2].split("', '"))
			print("number of binding sites not in an entity: ", len(nonentities))


	# unclustered binding sites whose chains ARE part of an entity
	# there must be another reason as to why these sites were not clustered
	sites_m2 = sites_m - nonentities
	print("missing bsites not accounted for by nonentities: ", len(sites_m2))

	enty_names = []
	# gets all entities of the second set of missing binding sites
	for bsite in sites_m2:
		tokens = bsite.strip().split("_")
		pdb_id = f"{tokens[0]}.ent"
		enty_n = -1
		enty_nam = ""

		with open(f"dat/struct/pdb/{pdb_id[4:6]}/{pdb_id}", "r") as f:
			for line in f:
				if SOURCE_PATTERN.match(line):
					print("FAILED - ", bsite, chains)
					break

				if (ent := ENTY_PATTERN.match(line)):
					enty_n = int(ent.group(1))
				elif (mol := MOL_PATTERN.match(line)):
					enty_nam = mol.group(1)
				elif (ch := CHAIN_PATTERN.match(line)):
					l = ch.group(1).strip()
					chains = []
					# if the line ends on a character, cut it
					if l[-1] == ';' or l[-1] == ',':
						chains += l[:-1].split(', ')
					else:
						chains += l.split(', ')
					# if the line ends in a comma, there are more chains in the next lines
					if l[-1] == ',':
						for line2 in f:
							l2 = line2[11:].strip()
							if l2[-1] == ';' or l2[-1] == ',':
								chains += l2[:-1].split(', ')
							else:
								chains += l2.split(', ')
							# if there is no comma, then the list has ended
							if l2[-1] != ',':
								break
						else:
							print("NO BREAK -", bsite)

					if tokens[1] in chains:
						enty_names.append((f"{tokens[0][3:].upper()}_{enty_n}", enty_nam))
						break

	# for nam in enty_names:
	# 	print(f"\n---{nam[0]}---")
	# 	print(nam[1])

	print("number of entities represented by ^: ", len(enty_names))
	enties, names = map(list, zip(*enty_names))
	fwrite("entities_missing_names.txt", '\n'.join(f"{e}: {n}" for e, n in enty_names), "list of names of aforementioned entities")


	CLUSTF = "dat/clusters-by-entity-70.txt"
	with open(CLUSTF, "r") as f:
		lines = f.readlines()
		# ignore non-PDB entries
		clusters = [[e for enty in clust.split() if len(e:=enty.strip()) <= 8] for clust in lines]
		clusters = [clust for clust in clusters if len(clust) > 0]
		print("number of clusters: ", len(clusters))

	enties_clust = [enty for c in clusters for enty in c]
	pdbs_clust = {enty[:4] for enty in enties_clust}

	print("number of entities in cluster-file (only PDB): ", len(enties_clust))
	print("sample item: ", enties_clust[0])
	print("number of PDBs in cluster-file: ", len(pdbs_clust))

	print()

	missing_enties = []
	for ent in enties:
		if ent not in enties_clust:
			missing_enties.append(ent)

	print("number of unaccounted-for entities not in the clusterfile: ", len(missing_enties))


def create_nonents(sites_m):
	"""Finds all binding sites whose chains are not part of an entity.
	Takes in set of all filtered binding sites that were not clustered."""

	print("---finding all binding sites whose chains are not part of an entity---")
	print("note: this may take a while. Set CREATE_NONENTS to False to skip this step.")

	removed = []
	nonentities = []

	for c, bsite in enumerate(sites_m):
		# counter
		if c % 10000 == 0:
			print(c)

		# split the site name into its components
		tokens = bsite.strip().split("_")
		pdb_id = f"{tokens[0]}.ent"
		chains = []

		try:
			with open(f"dat/struct/pdb/{pdb_id[4:6]}/{pdb_id}", "r") as f:
				for line in f:
					# if a SOURCE record is found then we have gone past the relevent section
					if SOURCE_PATTERN.match(line):
						break

					if CHAIN_PATTERN.match(line):
						# find all the chains in the entity section
						l = line[18:].strip()
						if l[-1] == ';' or l[-1] == ',':
							chains += l[:-1].split(', ')
						else:
							chains += l.split(', ')
						if l[-1] == ',':
							for line2 in f:
								l2 = line2[11:].strip()
								if l2[-1] == ';' or l2[-1] == ',':
									chains += l2[:-1].split(', ')
								else:
									chains += l2.split(', ')
								if l2[-1] != ',':
									break
							else:
								print("NO BREAK -", bsite)
		except FileNotFoundError as e:
			removed.append(bsite)
			continue

		# check if the chain is not in the list
		if tokens[1] not in chains:
			nonentities.append(bsite)

	fwrite("nonentities.txt", nonentities)
	fwrite("removed.txt", removed)
	print("number of non-entities: ", len(nonentities))
	print("number of sites that should have been removed: ", len(removed))

	return set(nonentities)


if __name__ == "__main__":
	main()