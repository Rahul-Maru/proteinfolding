import getpass
from concurrent.futures import ProcessPoolExecutor, as_completed
import os
import shlex
import shutil
import subprocess
import time

from openmm.app import PDBFile
from pdbfixer import PDBFixer
from boxsize import boxsize


# params
pt_space = 0.375
padding = 4 # Å

runs = 25
pop = 150
evals = 2500000


usr = getpass.getuser()
# put rootdir here
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
sitedir = f"{ROOT}/../filt/dat/struct/representative-binding-sites"
pdbdir = f"{ROOT}/../filt/dat/struct/pdb"

lig = f"{autodir}/hsm.pdbqt"

currsite = ""
worker_name = "main"

def main():
	global worker_name

	with open(inf) as f:
		sites = [line.strip() for line in f if line.strip()]

	if not sites:
		print("site file empty.")
		return

	max_workers = min(len(sites), os.cpu_count() or 1)
	print(f"running with {max_workers} workers")
	print(f"processing {len(sites)} sites")
	print("————————————————\n")

	abs_st_t = time.time()
	failed = []

	with ProcessPoolExecutor(max_workers=max_workers) as executor:
		futures = {executor.submit(process_site, i, site): site for i, site in enumerate(sites)}
		for future in as_completed(futures):
			site = futures[future]
			try:
				future.result()
			except Exception as e:
				failed.append(site)
				print(f"\033[093mworker failure in {site}: {type(e).__name__}. skipping. \033[0m")
				print("--------------------------------\n")

	abs_end_t = time.time()
	print(f"———all workers completed in {abs_end_t - abs_st_t:.2f}s———\n")

	with open(f"{ROOT}/outs/autodock/failed.txt", 'w') as f2:
		f2.write('\n'.join(failed))
		print(len(failed), "sites failed.")

	worker_name = "main"

def process_site(i, site):
	global currsite, worker_name

	currsite = site
	worker_name = f"pid-{os.getpid()}:{i}"

	sitef = f"{sitedir}/{site}"
	nam = site[:7]
	pdb = f"{pdbdir}/{site[4:6]}/{nam}.ent"

	outdir = f"{ROOT}/outs/autodock/{site[:-4]}"
	os.makedirs(outdir, exist_ok=True)

	recf = f"{outdir}/site.pdbqt"
	gpf = f"{outdir}/1.gpf"
	glg = f"{outdir}/1.glg"
	dpf = f"{outdir}/1.dpf"
	dlg = f"{outdir}/1.dlg"

	root = os.getcwd()
	os.chdir(outdir)

	# if the site has already been processed, the temporary ligand file will not exist
	#   and the dlg file will not be empty, so we can skip the site
	if os.path.exists(dlg) and os.path.getsize(dlg) > 0 and not os.path.exists("lig.py"):
		os.chdir(root)
		return

	# move ligand to output directory temporarily to match with the dpf file
	shutil.copy(lig, ".")
	
	tmp_pdb = f"{nam}.pdb"
	fixpdb(pdb, tmp_pdb)

	print(f"{worker_name}) processing site: {site}")
	site_st_t = time.time()

	curr_action = "initializing"
	try:
		curr_action = "preparing receptor"
		run(f"pythonsh {mgldir}/prepare_receptor4.py -r {nam}.pdb -o {recf} -A checkhydrogens")

		curr_action = "preparing grid file"
		cent, dims = boxsize(sitef, pt_space, padding)
		run(f"pythonsh {mgldir}/prepare_gpf4.py -l {lig} -r {recf} -p npts={dims} -p gridcenter={cent} -o {gpf}")

		curr_action = "running autogrid"
		run(f"{autodir}/autogrid4 -p {gpf} -l {glg}")

		curr_action = "preparing dock"
		run(f"pythonsh {mgldir}/prepare_dpf42.py -l {lig} -r {recf} -p ga_run={runs} -p ga_pop_size={pop} -p ga_num_evals={evals} -o {dpf}")

		curr_action = "docking"
		run(f"{autodir}/autodock4 -p {dpf} -l {dlg}")

		site_end_t = time.time()
		print(f"{worker_name}) Site {site} completed in {site_end_t - site_st_t:.2f}s\n")

	except Exception:
		print(f"\033[093m{worker_name}) Error in {site} while {curr_action}\033[0m")
		raise

	finally:
		for temp_file in ("hsm.pdbqt", tmp_pdb):
			os.remove(temp_file)
		os.chdir(root)


def fixpdb(pdb_file, outf):
	fixer = PDBFixer(pdb_file)
	fixer.findMissingResidues()
	fixer.findNonstandardResidues()
	fixer.findMissingAtoms()
	fixer.addMissingAtoms()
	fixer.replaceNonstandardResidues()
	fixer.removeHeterogens(False)
	fixer.addMissingHydrogens(7.4)
	PDBFile.writeFile(fixer.topology, fixer.positions, open(outf, 'w'))

def run(cmd, **kwargs):
	"""Runs a command and processes error handling.
	"""
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
			print(f"\033[31m[ERROR IN WORKER {worker_name} - {currsite}]\033[0m")
			print(f"{worker_name}) command {cmd} failed with error code {result.returncode}")
			print(result.stderr)
			raise RuntimeError(f"command {cmd} failed with code {result.returncode}")

		if "WARNING" in result.stderr:
			print(f"\033[93m[WARNING IN WORKER {worker_name} - {currsite}]\033[0m")
			print(result.stderr)

	except subprocess.CalledProcessError as e:
		print(f"\033[31m[ERROR IN WORKER {worker_name} - {currsite}]\033[0m")
		print(f"command {cmd} failed with error code {e.returncode}")
		if e.stderr:
			# these errors are likely due to malformed geometry and can be skipped
			if "ZeroDivisionError" in e.stderr:
				print(f"{worker_name}) {e.stderr}")
				raise ZeroDivisionError
			elif "ValueError" in e.stderr:
				print(f"{worker_name}) {e.stderr}")
				raise ValueError
			elif "RuntimeError" in e.stderr:
				print(f"{worker_name}) {e.stderr}")
				raise RuntimeError

			print(f"{worker_name}) Error output: {e.stderr}")

		if e.stdout:
			print(f"{worker_name}) Standard output: {e.stdout}")

		raise

	return result

if __name__ == "__main__":
	main()
