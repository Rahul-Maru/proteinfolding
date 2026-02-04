import json
import os
import sh

clusts = json.load(open("hsm/outs/autodock/docking_clusts.json"))

poses = {}
for p, c in clusts.items():
    with open(f"hsm/outs/autodock/{p}/1.dlg") as f:
        for l in f:
                if l.startswith("USER    Cluster Rank = "):
                    r = int(l.split()[-1])
                    if r in c:
                        poses[f"{p}_{r}"] = []
                        for l2 in f:
                            if l2.startswith("ATOM"):
                                poses[f"{p}_{r}"].append("HETATM" + l2[6:])
                            elif l2.startswith("ENDMDL"):
                                break

                        with open(f"hsm/outs/autodock/poses/{p}_{r}.pdb", "w") as outf, \
                              open(f"hsm/outs/autodock/{p}/site.pdbqt", "r") as sitef:
                            site = sitef.readlines()

                            outf.writelines(site)
                            outf.writelines(poses[f"{p}_{r}"])

print('\n'.join([f"{k}:\n{''.join(v)}" for k, v in poses.items()]))
print(len(poses))
print(sum([len(v) for _, v  in poses.items()]))
json.dump(poses, open("hsm/outs/autodock/docking_poses.json", "w"))