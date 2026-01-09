import bio

inf = "hsm/bsites/2x45_C6.pdb"

prot = bio.Protein(inf)

hsm = [x.strip() for x in prot.hetatms]

print()

tail = [a for a in hsm if a[bio.ATN].strip() == "N"][0]

coords = lambda l : (float(l[bio.X_COORD]), float(l[bio.Y_COORD]), float(l[bio.Z_COORD]))
tailcords = coords(tail)

print(tail)
print(tailcords)

head = [a for a in hsm if a[bio.ATN].strip() in ['CG', 'ND1', 'CD2', 'CE1', 'NE2']]


print('\n'.join(head))