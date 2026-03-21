import json
from subprocess import run
from bio import ATNO


incl_receptor = False

clusts = json.load(open("hsm/outs/autodock/docking_clusts.json"))

poses = {}
for prot, clust in clusts.items():
    ranks, runs = zip(*clust)
    with open(f"hsm/outs/autodock/results/{prot}/1.dlg") as f:
        for l in f:
            if l.startswith("Run"):
                r = int(l.split()[1])
                if r in ranks:
                    # get site data
                    with open(f"hsm/outs/autodock/results/{prot}/site.pdbqt", "r") as sitef:
                        site = sitef.readlines()
                        maxn = int(site[-1][ATNO])

                        # check if extraneous hetatms in receptor
                        if "HETATM" in ''.join(site):
                            print("HETATM in site: ", prot)
                    
                    poses[f"{prot}_{r}"] = []
                    for l2 in f:
                        if l2.startswith("DOCKED: ATOM"):
                            # replace DOCKED: ATOM with HETATM, fix line no, and set chain identifier
                            l2 = l2[8:]
                            n = int(l2[ATNO]) + maxn
                            poses[f"{prot}_{r}"].append(f"HETATM{n:5d}" + l2[11:21] + "A" + l2[22:])
                        elif l2.startswith("DOCKED: ENDMDL"):
                            break

                    with open(f"hsm/outs/autodock/poses/{prot}_{r}.pdbqt", "w") as outf:
                        if incl_receptor:
                            # write receptor lines to pose file
                            outf.writelines(site)

                        outf.writelines(poses[f"{prot}_{r}"])

# print('\n'.join([f"{k}:\n{''.join(v)}" for k, v in poses.items()]))
print(len(poses))
print(sum([len(v) for _, v  in poses.items()]))
json.dump(poses, open("hsm/outs/autodock/docking_poses.json", "w"))