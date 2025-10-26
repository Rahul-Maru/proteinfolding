import getpass
import math
from mpi4py import MPI
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
padding = 12 # Å

runs = 25
pop = 150
evals = 2500000

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

def main():
	if rank == 0:
		print(f"MPI Initiated. Running with {size} cores. \n")

	failed = []
	with open(inf) as f:
		# split up the sites for each core
		sites = f.readlines()
		subsites = sites[rank::size]

		if rank == 0:
			print(f"processing {len(sites)} sites.")

		for i in range(size):
			if rank == i:
				print(f"core {rank} processing {len(subsites)} sites.")
			comm.Barrier()

		if rank == 0:
			print("————————————————\n")

		comm.Barrier()

		for c, site in enumerate(subsites):
			site = site.strip()
			receptor = f"{indir}/{site}"

			outdir = f"{ROOT}/outs/autodock/{site[:-4]}"
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


			print(f"{rank}.{c}) processing site: {site}")
			site_start_time = time.time()

			try:
				x = "preparing receptor"
				# prepare receptor
				run(f"pythonsh {mgldir}/prepare_receptor4.py -r {receptor} -o {outf} -A checkhydrogens")

				x = "preparing grid file"
				# grid
				cent, dims = boxsize(receptor, pt_space, padding)
#				print(f"{site} ({rank}): grid center, size: ", cent, dims)
				run(f"pythonsh {mgldir}/prepare_gpf4.py -l {lig} -r {outf} -p npts={dims} -p gridcenter={cent} -o {gpf}")

				x = "running autogrid"
				run(f"{autodir}/autogrid4 -p {gpf} -l {glg}")

				x = "preparing dock"
				# docking
				run(f"pythonsh {mgldir}/prepare_dpf42.py -l {lig} -r {outf} -p ga_run={runs} -p ga_pop_size={pop} -p ga_num_evals={evals} -o {dpf}")

				x = "docking"
				run(f"{autodir}/autodock4 -p {dpf} -l {dlg}")

				os.remove("hsm.pdbqt")
			except Exception as e:
				os.remove("hsm.pdbqt")
				failed.append(site)
				print(f"\033[093m{type(e).__name__} in {site}, while {x}, skipping\033[0m")
				continue

			site_end_time = time.time()
			site_elapsed = site_end_time - site_start_time
			print(f"{rank}.{c}) Site {site} completed in {site_elapsed:.2f} seconds\n")

	print(f"completed on core {rank}.")

	dat = comm.gather(failed)

	if rank == 0:
		allfailed = sum(dat, [])
		with open(f"{ROOT}/outs/autodock/failed.txt", 'w') as f2:
			f2.write('\n'.join(allfailed))
			print(len(allfailed), "sites failed.")


def run(cmd, **kwargs):
	try:
		result = subprocess.run(
			shlex.split(cmd),
			check=True,
			capture_output=True,
			text=True,
			**kwargs
		)
		if result.stdout and "setting PYTHONHOME" not in result.stdout:
			print(result.stdout)

		if result.returncode != 0:
			print(f"\033[31m[ERROR IN CORE {rank}]\033[0m")
			print(f"command {cmd} failed with error code {result.returncode}")
			print(result.stderr)
#			comm.Abort(42)
			raise Exception

		if "WARNING" in result.stderr:
			print(f"\033[93m[WARNING IN CORE {rank}]\033[0m")
			print(result.stderr)

	except subprocess.CalledProcessError as e:
		print(f"\033[31m[ERROR IN CORE {rank}]\033[0m")
		print(f"command {cmd} failed with error code {e.returncode}")
		if e.stderr:
			# these errors are likely due to malformed geometry and can be skipped
			if "ZeroDivisionError" in e.stderr:
				raise ZeroDivisionError
			elif "ValueError" in e.stderr:
				raise ValueError
			elif "RuntimeError" in e.stderr:
				print(e.stderr)
				raise RuntimeError

			print(f"Error output: {e.stderr}")

		if e.stdout:
			print(f"Standard output: {e.stdout}")
#		comm.Abort(42)
		raise Exception

	return result

if __name__ == "__main__":
	main()
