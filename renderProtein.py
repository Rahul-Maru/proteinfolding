"""Plots the atoms in a protein using VPython."""
import vpython as vp
from consts import *
from Protein import Protein


def render(p: Protein):
	"""Renders a Protein p in vpython using the sphere model."""

	rainbow, show_hetatms, b_site = p.display_mode

	print(f'Rendering protein at "{p.path}"')
	print("———α HELIX RESIDUES———")
	for r in p.get_ss('HELIX'):
		print(f"{r['code']} - {r['n']}")
	print("———β SHEET RESIDUES———")
	for r in p.get_ss('SHEET'):
		print(f"{r['code']} - {r['n']}")
	print(p)

	if b_site:
		ligand = p.get_ligand('ATP')
		x, y, z = p.get_xyzlist(res=ligand)
		for i, atom in enumerate(ligand['atoms']):
			spr = vp.sphere()
			spr.pos = vp.vector(x[i], y[i], z[i])
			spr.radius = RAD

			spr.color = vp.vector(0.42, 0, 0.69)

		cent = vp.sphere()
		cent.pos = vp.vector(*tuple(p.centroid(res=ligand)))
		cent.radius = RAD*1.4
		cent.color = vp.vector(*COLORS['CENTRD'][int(rainbow)])

		print("Position of ligand's centroid (Å):", cent.pos)

		# Center the camera on the centroid
		vp.scene.center = cent.pos

	else:
		for i, atom in enumerate(p.atoms):
			spr = vp.sphere()
			spr.pos = vp.vector(p.x[i], p.y[i], p.z[i])
			spr.radius = RAD

			elem = atom[77] # get the element of the atom
			chain = CHAIN_ID(atom)

			if rainbow:
				# colors the atoms in a rainbow spectrum
				# calculates the relative position of the atom in its chain reduced to the range [0, 2π]
				rel_pos = 2*np.pi * (i - len(sum(p.chains[:chain], []))) / len(p.chains[chain])
				spr.color = vp.vector(*hue_to_RGB(rel_pos))
			else:
				# pick the appropriate color based on element and chain
				# if the chain number exceeds the number of available colors,
				# loop through the colors
				try:
					spr.color = vp.vector(*COLORS[elem][chain % len(COLORS[elem])])
				except KeyError:
					# use the hetatm color if for some reason an element other than the
					#  standard 5 is encountered in an ATOM record
					spr.color = vp.vector(*COLORS['HETATM'][int(rainbow)])

		# display the centroid
		cent = vp.sphere()
		cent.pos = vp.vector(*tuple(p.centroid()))
		cent.radius = RAD*2
		cent.color = vp.vector(*COLORS['CENTRD'][int(rainbow)])
		print("Position of protein's centroid (Å):", cent.pos)

		# Center the camera on the centroid
		vp.scene.center = cent.pos


		if show_hetatms and len(p.hetatms) > 0:
			hx, hy, hz = p.get_xyzlist(-2)

			for hi, hj, hk in zip(hx, hy, hz):
				spr = vp.sphere()
				spr.pos = vp.vector(hi, hj, hk)
				spr.radius = RAD
				spr.color = vp.vector(*COLORS['HETATM'][int(rainbow)])
				spr.opacity = OPCTY

			# display the hetatm centroid
			hcent = vp.sphere()
			hcent.pos = vp.vector(*tuple(p.centroid(-2)))

			hcent.radius = RAD*1.8
			hcent.color = vp.vector(*COLORS['HETCTD'][int(rainbow)])

			print("Position of heterogens' centroid (Å):", hcent.pos)

