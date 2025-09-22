with open("old/bsites_old.txt") as f:
	lins = f.readlines()
	bsites_old = {b.strip() for b in lins}

with open("bsites.txt") as f:
	lins = f.readlines()
	bsites = {b.strip() for b in lins}

with open("old/f_bsites_old.txt") as f:
    lins = f.readlines()
    fbsites_old = {b.strip() for b in lins}

with open("f_bsites.txt") as f:
    lins = f.readlines()
    fbsites = {b.strip() for b in lins}

print(len(bsites_old), len(bsites), len(old:=bsites_old-bsites), len(new:=bsites-bsites_old), len(bsites&bsites_old))


print(len(fbsites_old), len(fbsites), len(fold:=fbsites_old-fbsites), len(fnew:=fbsites-fbsites_old), len(fbsites&fbsites_old))

print(len(fold-old), len(fnew-new))
print(fold-old)
print("----")
print(fnew-new)
print(len(fbsites-bsites), len(fbsites_old-bsites_old))
