with open("bsites_old.txt") as f:
	lins = f.readlines()
	bsites_old = {b.strip() for b in lins}

with open("bsites.txt") as f:
	lins = f.readlines()
	bsites = {b.strip() for b in lins}

print(len(bsites_old), len(bsites), len(bsites_old-bsites), len(bsites-bsites_old), len(bsites&bsites_old))
