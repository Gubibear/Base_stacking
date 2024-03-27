# Base_stacking
Repository for files and scripts related to the stack scoring algorithm.

The planes.py script has and option to read pdb or xyz files and calculate stacking score between adjacent bases.
For the pdb file it is enough to just call the script with -p option followed by the name of the file containing the structure in pdb format.

$ python planes.py -p _filename_.pdb

Current version of the script requires more involvement from the user. 
Firstly, it is necessary to provide file containing the atom numbers of atoms from the xyz file, that are forming the aromatic rings of nucleotides. The file should be named bases.txt and it needs to be formated in the same way as the file bases.txt located in the example directory.
It is also necessary to provide number of bases that re to be read from the xyz file while calling the script.

$ python planes.py -x xyz _filename_.xyz -n _number_of_bases_

The calculated stacking scores with additional information will be gathered in the raport.txt file. 

The planes.py script can be tested with db_0 files which are avaiable in the example directory. The example usage of planes.py script for these files would be as follows:

$ python planes.py -p db_0.pdb 

or 

$ python planes.py -x db_0.xyz -n 3

