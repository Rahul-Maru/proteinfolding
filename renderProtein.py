"""Plots the atoms in a protein using VPython."""
import vpython as vp
from vpython import vector as vect
from consts import *
from Protein import Protein


def render(p: Protein):
	"""Renders a Protein in vpython using the atom-sphere model."""

	rainbow, show_hetatms, show_bsite = p.display_mode

	print(f'Rendering protein at "{p.path}"')

	# lists residues in helices and sheets
	print("———α HELIX RESIDUES———")
	Protein.print_res(p.get_ss('HELIX'))
	print("———β SHEET RESIDUES———")
	Protein.print_res(p.get_ss('SHEET'))

	# prints general info about the protein
	print(p)

	# if b_site:
	# 	ligand = p.get_ligand(DEF_LIG)
	# 	atoms = ligand['atoms']
	# 	x, y, z = p.get_xyzlist(res=ligand)
	# 	colors = [vect(0.42, 0, 0.69)]*len(atoms)

	# else:
	# 	atoms = p.atoms
	# 	x, y, z = p.get_xyzlist()
	# 	if rainbow:
	# 		colors = []
	# 		for i, atom in enumerate(p.atoms):
	# 			# colors the atoms in a rainbow spectrum
	# 			# calculates the relative position of the atom in its chain reduced to the range [0, 2π]
	# 			chain = CHAIN_ID(atom)

	# 			rel_pos = 2*np.pi * (i - len(sum(p.chains[:chain], []))) / len(p.chains[chain])
	# 			colors.append(vect(*hue_to_RGB(rel_pos)))
	# 	else:
	# 		colors = []


	if show_bsite:
		ligand = p.get_ligand('ATP')
		x, y, z = p.get_xyzlist(ligand, qtype="res")
		for i, atom in enumerate(ligand['atoms']):
			spr = vp.sphere()
			spr.pos = vect(x[i], y[i], z[i])
			spr.radius = RAD
			spr.color = vect(0.42, 0, 0.69)

		bsite = p.get_bsite(ligand)
		Protein.print_res(bsite)
		for res in bsite:
			coords = p.get_xyzlist(res, triplet=True)

			for x, y, z in coords:
				spr = vp.sphere()
				spr.pos = vect(x, y, z)
				spr.radius = RAD
				spr.color = vect(1, 0.75, 0)


		try:
			cent = vp.sphere()
			cent.pos = vect(*tuple(p.centroid(ligand)))
			cent.radius = RAD*1.4
			cent.color = vect(*COLORS['CENTRD'][int(rainbow)])

			print("Position of ligand's centroid (Å):", cent.pos)

			# Center the camera on the centroid
			vp.scene.center = cent.pos

		except ValueError as e:
			print("WARNING: ", e)
			s = vp.sphere()
			s.color = vect(1,0,0)


	else:
		for i, atom in enumerate(p.atoms):
			spr = vp.sphere()
			spr.pos = vect(p.x[i], p.y[i], p.z[i])
			spr.radius = RAD

			elem = atom[77] # get the element of the atom
			chain = CHAIN_ID(atom)

			if rainbow:
				# colors the atoms in a rainbow spectrum
				# calculates the relative position of the atom in its chain reduced to the range [0, 2π]
				rel_pos = 2*np.pi * (i - len(sum(p.chains[:chain], []))) / len(p.chains[chain])
				spr.color = vect(*hue_to_RGB(rel_pos))
			else:
				# pick the appropriate color based on element and chain
				# if the chain number exceeds the number of available colors,
				# loop through the colors
				try:
					spr.color = vect(*COLORS[elem][chain % len(COLORS[elem])])
				except KeyError:
					# use the hetatm color if for some reason an element other than the
					#  standard 5 is encountered in an ATOM record
					spr.color = vect(*COLORS['HETATM'][int(rainbow)])

		# display the centroid
		cent = vp.sphere()
		cent.pos = vect(*tuple(p.centroid()))
		cent.radius = RAD*2
		cent.color = vect(*COLORS['CENTRD'][int(rainbow)])
		print("Position of protein's centroid (Å):", cent.pos)

		# Center the camera on the centroid
		vp.scene.center = cent.pos


		if show_hetatms and len(p.hetatms) > 0:
			hx, hy, hz = p.get_xyzlist(qtype="het")

			for hi, hj, hk in zip(hx, hy, hz):
				spr = vp.sphere()
				spr.pos = vect(hi, hj, hk)
				spr.radius = RAD
				spr.color = vect(*COLORS['HETATM'][int(rainbow)])
				spr.opacity = OPCTY

			# display the hetatm centroid
			hcent = vp.sphere()
			hcent.pos = vect(*tuple(p.centroid(qtype='het')))

			hcent.radius = RAD*1.8
			hcent.color = vect(*COLORS['HETCTD'][int(rainbow)])

			print("Position of heterogens' centroid (Å):", hcent.pos)

