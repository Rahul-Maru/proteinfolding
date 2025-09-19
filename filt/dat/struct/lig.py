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
    i=0
    nf=[]
    for file_path in Path(INDR).iterdir():
        nam = file_path.name
        ligArr = nam.split('_')
        n = ligArr[2]
        og_nam = f"{ligArr[0]}.ent"
        subdir = og_nam[4:6]

        # print(n, og_nam, subdir)

        pdb_file = Path(f"pdb/{subdir}/{og_nam}")

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
                            i+=1
                            # Copy the original file if output doesn't exist
                            if not output_file.exists():
                                shutil.copy2(file_path, output_file)

                        # Append the HETATM line to output file
                        with open(output_file, 'a') as out_f:
                            out_f.write(line)
                    elif found:
                        break
        if not found:
            print(f"{i}: not found {nam} ({pdb_file}), {n}, {num}")
            nf.append(nam)
        elif not i%1000:
            print(i)

    print(len(nf))
    with open("nf_lig.txt", 'w') as f:
        f.write("\n".join(nf))

if __name__ == "__main__":
    main()
