"""Plots the atoms in a protein using VPython."""
import vpython as vp
from consts import *
from Protein import Protein


def render(p: Protein, dp_mode: list[bool] = [False, False]):
	"""Renders the Protein p in vpython. dp_mode is a list of the following values:
	color_scheme: 0 (normal), 1 (rainbow)
	het_atms: 0 (only shows normal atoms), 1 (shows hetatms as well)."""

	rainbow, show_hetatms = (*dp_mode,)

	print(f'Rendering protein at "{p.path}"')
	print(p)
	print(p.get_ss('HELIX'))
	print(p.get_ss('SHEET'))

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
				spr.color = vp.vector(*COLORS[elem][chain % 3])
			except KeyError:
				# use the hetatm color if for some reason an element other than the
				#  standard 5 is encountered in an ATOM record
				spr.color = vp.vector(*COLORS['HETATM'][int(rainbow)])

	# display the centroid
	cent = vp.sphere()
	cent.pos = vp.vector(*tuple(p.centroid()))
	cent.radius = RAD*2
	cent.color = vp.vector(*COLORS['CENTRD'][int(rainbow)])
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
