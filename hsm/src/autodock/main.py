import getpass
import os
import shlex
import shutil
import subprocess
import time
from boxsize import boxsize


usr = getpass.getuser()

# put rootdir here.
if usr == "rahul":
	ROOT = "~/Desktop/programs/bio/hsm"
elif usr == "aibio":
	ROOT = "~/Desktop/bio/hsm"
else:
	raise Exception("no root dir set. set your root directory and remove this check from the code")
ROOT = os.path.expanduser(ROOT)

autodir = f"{ROOT}/tools/autodock"
mgldir = f"{ROOT}/tools/mgltools/MGLToolsPckgs/AutoDockTools/Utilities24"

inf = f"{ROOT}/outs/FLAPP/alignlist2.txt"
if usr == 'rahul':
	indir = f"{ROOT}/tools/FLAPP/sites"
elif usr == 'aibio':
	indir = f"{ROOT}/../filt/dat/struct/representative-binding-sites"

lig = f"{autodir}/hsm.pdbqt"

# params
pt_space = 0.375
padding = 12

runs = 50
pop = 300
evals = 2500000

def main():
	with open(inf) as sites:
		c = 0
		for site in sites:
			if c == 1: break
			c += 1

			site = site.strip()
			mid = site[4:6]
			receptor = f"{indir}/{site}"

			outdir = f"{ROOT}/outs/autodock/{mid}/{site[:-4]}"
			os.makedirs(outdir, exist_ok=True)

			# output files
			outf = f"{outdir}/site.pdbqt"
			gpf = f"{outdir}/1.gpf"
			glg = f"{outdir}/1.glg"
			dpf = f"{outdir}/1.dpf"
			dlg = f"{outdir}/1.dlg"

			# move on from already-processed site
			if os.path.exists(dlg) and os.path.getsize(dlg) > 0:
				continue

			# navigate to the output directory and temporarily copy the ligand there
			os.chdir(outdir)
			shutil.copy(lig, ".")


			print(f"{c}) processing site: {site}")
			site_start_time = time.time()

			# prepare receptor
			run(f"pythonsh {mgldir}/prepare_receptor4.py -r {receptor} -o {outf} -A checkhydrogens")

			# grid
			cent, size = boxsize(receptor, pt_space, padding)
			print("grid center, size: ", cent, size)
			run(f"pythonsh {mgldir}/prepare_gpf4.py -l {lig} -r {outf} -p npts={size} -p gridcenter={cent} -o {gpf}")

			run(f"{autodir}/autogrid4 -p {gpf} -l {glg}")

			# docking
			run(f"pythonsh {mgldir}/prepare_dpf42.py -l {lig} -r {outf} -p ga_run={runs} -p ga_pop_size={pop} -p ga_num_evals={evals} -o {dpf}")

			run(f"{autodir}/autodock4 -p {dpf} -l {dlg}")

			os.remove("hsm.pdbqt")

			site_end_time = time.time()
			site_elapsed = site_end_time - site_start_time
			print(f"Site {site} completed in {site_elapsed:.2f} seconds\n")

def run(cmd, **kwargs):
	result = subprocess.run(
		shlex.split(cmd), 
		check=True, 
		capture_output=True, 
		text=True,
		**kwargs
	)
	# Print stdout but filter stderr
	if result.stdout and "setting PYTHONHOME" not in result.stdout:
		print(result.stdout)

if __name__ == "__main__":
	main()
