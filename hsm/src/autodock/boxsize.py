from math import ceil
import bio

def boxsize(fp, s=0.375, p=0):
	"""Calculates an appropriate size for the autodock grid box given a receptor file path,
	point spacing (defaults to `0.375`), and padding on each side (in Å, defaults to `0`).
	Returns (center, size)."""

	prot = bio.Protein(fp)
	cx, cy, cz = prot.centroid(qtype="atm")
	x, y, z = prot.get_coords()

	def npts(k, ck):
		"""Returns the required number of points in the k-dimension"""
		# finds the max of the farthest +ve dist and -ve dist from the centroid ck,
		#	adds padding dist, divides by spacing, rounds, and doubles to make it even for autogrid
		return ceil((max(max(k) - ck, ck - min(k)) + p) / (s)) * 2

	dx = npts(x, cx)
	dy = npts(y, cy)
	dz = npts(z, cz)

	return (f"{cx},{cy},{cz}", f"{dx},{dy},{dz}")

if __name__ == "__main__":
	print(boxsize("hsm/bsites/2x45_C7.pdb", p=4))
