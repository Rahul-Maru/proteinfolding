import sys
from bio import fixpdb

# usage: python fixreceptor.py input.pdb output.pdb
if __name__ == "__main__":
	fixpdb(sys.argv[1], sys.argv[2])

