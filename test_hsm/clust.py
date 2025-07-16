import os

mode = 'ver'
hsm = [f[:-4].upper() for f in os.listdir("hsm/pdbs")]

if mode == 'wrt':
    with open("filt/dat/clusters-by-entity-70.txt") as f:
        clusters = f.readlines()

    # # ignore non-PDB entries
    clusters = [[enty.strip() for enty in clust.split(" ") if len(enty) <= 8] for clust in clusters]
    clusters = [clust for clust in clusters if len(clust) > 0]

    cl2 = [[]]

    for l in clusters:
        added = False
        for e in l:
            if e[:4] in hsm:
                cl2[-1].append(e)
                added = True
        if added:
            cl2.append([])

    print ('\n'.join([' '.join(l) for l in cl2]))
else:
    with open("test_hsm/clusters-by-entity-70.txt") as f:
        lines = f.read()
    lines = lines.replace('\n', '')
    for h in hsm:
        if h not in lines:
            print(f"{h} NOT OKAY")
        else:
            print(f"{h} OKAY")