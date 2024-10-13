Author: Rahul Maru
Date: October 2024
Language: Python 3.10.12

A program that models protein folding using VPython

Done as part of an internship with Dr. Nagasuma Chandra at IISc's department of Biochemistry 

To run the code, run main.py. Constants can be changed in consts.py
The flags --rainbow and --hetatm can be appended to 'python3 main.py' to display the atoms in a rainbow color scheme, and to show the HETATMs/heterogens, respectively. You can also append the name of the file to the command to make it open that file (ONLY the 4 letter file name, not the .pdb or the path). Otherwise it will open the file specified by consts.DEF_PROT.
e.g. >> python3 1acj --rainbow --hetatm