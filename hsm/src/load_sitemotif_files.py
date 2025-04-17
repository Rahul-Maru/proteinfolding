import argparse
from subprocess import run

def load_sitemotif_files(path):
    """
    Moves the given directory to the sitemotif directory.
    """

    try:
        run(["rm", "-r", "hsm/tools/MAPP-3D/MultipleSiteAlignment/hsm/"])
    except:
        pass

    run(["cp", "-r", path, "hsm/tools/MAPP-3D/MultipleSiteAlignment/hsm/"])
    

def main():
    DEF = "hsm/bsites_split"

    parser = argparse.ArgumentParser()
    parser.add_argument("--path", type=str, default=DEF, dest="path")
    args = parser.parse_args()

    path = args.path

    load_sitemotif_files(path)

if __name__ == "__main__":
    main()
