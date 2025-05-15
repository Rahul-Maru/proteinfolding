import os
import subprocess
import glob

REPR = "2x45_C6.pdb"

def run_mapp3d_comparisons(query_site, sites_dir, output_dir):
    """
    Run MAPP-3D pairwise comparisons between query site and all sites in directory

    Args:
        query_site: Path to query PDB file
        sites_dir: Directory containing all binding site PDB files
        output_dir: Directory to store output files
    """
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print("Made directory")
    # Delete all non-hidden files in output directory
    for file in os.listdir(output_dir):
        if not file.startswith('.'):  # Skip hidden files
            file_path = os.path.join(output_dir, file)
            if os.path.isfile(file_path):
                os.remove(file_path)
                print(f"Deleted {file}")

    print("----------------------")

    # Get absolute path to MAPP-3D script
    mapp3d_script = os.path.abspath(os.path.join("hsm", "tools", "MAPP-3D", "PairWiseComparison", "pocket_matrix7.py"))

    # Get all PDB files in sites directory
    site_pdbs = glob.glob(os.path.join(sites_dir, "*.pdb"))
    print(f"Found {len(site_pdbs)} PDB files in {sites_dir}")
    print(f"Query site: {os.path.basename(query_site)}")

    # Run comparison against each site
    for site_pdb in site_pdbs:
        print(f"\nProcessing: {os.path.basename(site_pdb)}")
        # Skip comparing against self
        if os.path.basename(site_pdb) == os.path.basename(query_site):
            print(f"Skipping self-comparison with {os.path.basename(site_pdb)}")
            continue

        # Run MAPP-3D comparison

        print(f"Running comparison...")
        result = subprocess.run([
            "python2.7",
            mapp3d_script,
            site_pdb,
            query_site
        ], cwd=output_dir, capture_output=True, text=True)

        if result.returncode != 0:
            print(f"Error comparing {os.path.basename(query_site)} with {os.path.basename(site_pdb)}:")
            print(result.stderr)
        else:
            print(f"Successfully completed comparison")

        try:
            repr_path = os.path.join(output_dir, os.path.basename(query_site))
            if not os.path.isfile(repr_path):
                os.rename(os.path.join(output_dir, 'fixed.pdb'), repr_path)
                print(f"Renamed fixed.pdb of {os.path.basename(site_pdb)}")
            else:
                os.remove(os.path.join(output_dir, 'fixed.pdb'))

            for file in ['align.txt', 'frag.pdb', 'site1.pdb', 'site2.pdb']:
                os.rename(os.path.join(output_dir, file), os.path.join(output_dir,
                                                        f"{os.path.basename(site_pdb)[:-4]}_{file}"))
            
            
            
        except Exception:
            print(f"Comparison failed for {os.path.basename(site_pdb)}")


def main():
	# Run comparisons
	query_site = os.path.abspath(os.path.join("hsm", "bsites_final", REPR))
	sites_dir = os.path.abspath(os.path.join("hsm", "bsites_final"))
	output_dir = os.path.abspath(os.path.join("hsm", "outs", "mapp3d_results"))

	run_mapp3d_comparisons(query_site, sites_dir, output_dir)

if __name__ == "__main__":
	main()
