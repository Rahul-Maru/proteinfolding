import subprocess
from boxsize import boxsize

autodir = "hsm/tools/autodock"
mgldir = "hsm/tools/mgltools/MGLToolsPckgs/AutoDockTools/Utilities24"

inf = "hsm/outs/FLAPP/alignlist2.txt"
indir = "hsm/tools/FLAPP/sites"
lig = f"{autodir}/hsm.pdbqt"

def main():
	with open(inf) as sites: 
		for site in sites:
			site = site.strip()
			receptor = f"{indir}/{site}"

			run(f"pythonsh {mgldir}/prepare_receptor4.py -r {receptor} -o {outf} -A checkhydrogens")

			cent, size = boxsize(receptor)
			run(f"pythonsh {mgldir}/prepare_gpf4.py -l {lig} -r {outf} -p npts={size} -p gridcenter={cent} -o {gpf}")

			run(f"./{{autodir}}/autogrid4 -p {gpf} -l {glg}")

			run(f"pythonsh {mgldir}/prepare_dpf4.py -l {lig} -r {outf} -p ga_run={50} -p ga_pop_size={300} -p ga_num_evals={2500000} -o {dpf}")


def run(cmd):
	subprocess.run(cmd.split(" "))


if __name__ == "__main__":
	main()
