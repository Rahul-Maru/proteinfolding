Author: Rahul Maru
Date: October 2024
Language: Python 3.10.12

A collection of projects related to structural biology
1) A program that models protein folding using VPython
2) The programs used to derive a binding site motif for Histamine

Done as part of an internship with Dr. Nagasuma Chandra at IISc's department of Biochemistry 



To run the code, run main.py. Constants can be changed in consts.py
The flags --rainbow and --hetatm can be appended to 'python3 main.py' to display the atoms in a rainbow color scheme, and to show the HETATMs/heterogens, respectively. You can also append the name of the file to the command to make it open that file (ONLY the 4 letter file name, not the .pdb or the path). Otherwise it will open the file specified by consts.DEF_PROT.
e.g. >> python3 1acj --rainbow --hetatm

To be able to import lib/bio.py and use the programs herein, add the following line to ~/.bashrc:
export PYTHONPATH=${PYTHONPATH}:$[PATH_TO_REPOSITORY]/lib
then open a new command window or run `source ~/.bashrc`
