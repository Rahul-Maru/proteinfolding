import os
import shutil
from pathlib import Path

def main():
    # Directory paths
    INDR = "representative-binding-sites"
    OUTDIR = "ligand-inclsv-binding-sites"
    
    # Create output directory if it doesn't exist
    if not os.path.exists(OUTDIR):
        print(OUTDIR)
        os.makedirs(OUTDIR)
    
    # Process each file in the input directory
    for file_path in Path(INDR).iterdir():
        if file_path.is_file():
            nam = file_path.name
            ligArr = nam.split('_')
            n = ligArr[2]
            og_nam = f"{ligArr[0]}.ent"
            subdir = og_nam[4:6]

            print(n, og_nam, subdir)

            # Check if the corresponding PDB file exists
            pdb_file = Path(f"pdb/{subdir}/{og_nam}")
            if not pdb_file.exists():
                print(f"PDB file not found: {pdb_file}")
                continue
            
            found = False
            output_file = Path(OUTDIR) / nam
            
            # Read the PDB file and process HETATM records
            with open(pdb_file, 'r') as pdb_f:
                for line in pdb_f:
                    rec = line[:6].strip()
                    if rec == "HETATM":
                        num_str = line[22:26].strip()
                        num = f"{int(num_str):04d}" if num_str else ""
                        
                        if num == n:
                            if not found:
                                found = True
                                # Copy the original file if output doesn't exist
                                if not output_file.exists():
                                    shutil.copy2(file_path, output_file)
                            
                            # Append the HETATM line to output file
                            with open(output_file, 'a') as out_f:
                                out_f.write(line)
            
            if not found:
                print(f"not found {nam}")

if __name__ == "__main__":
    main()
