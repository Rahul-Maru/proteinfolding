import bio

def boxsize(fp):
	p = 16
	s = 0.375

	prot = bio.Protein(fp)
	cx, cy, cz = prot.centroid(qtype="atm")
	x, y, z = prot.get_coords()

	dx = round((max(max(x) - cx, cx - min(x)) + p) / (2*s)) * 2
	dy = round((max(max(y) - cy, cy - min(y)) + p) / (2*s)) * 2
	dz = round((max(max(z) - cz, cz - min(z)) + p) / (2*s)) * 2

	return (f"{cx},{cy},{cz}", f"{dx},{dy},{dz}")

if __name__ == "__main__":
	print(boxsize("hsm/bsites/2x45_C7.pdb"))
