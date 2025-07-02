# Base_stacking
Repository for files and scripts related to the stack scoring algorithm.

#########IMPORTANT#########
##########UPDATE###########
New version of stacking score algorithm is available since 2025-06-28

The main novelty is an updated definition of perfect overlap of two nucleobase rings. The algorithm is now tuned to give a higher overlap value to structures in which at least two rings (or subunits of purine rings) completely overlap.

The new algorithm is available in the stack_2.0.py script.

This new version does not work with pdb files (for now - we are working on restoring this feature) and requires a different file format that specifies the atoms in aromatic rings.
For purine rings, the script now requires splitting the ring into two subunits. The order of atoms in the rings is very important for the algorithm's workflow. The last two atoms of each subunit must be the same carbon atoms shared by each purine ring subunit. The atom numbers must also be specified in the same order (i.e. clockwise/counterclockwise) for the two subunits.

An example of this new format is available in the base_spec file in the example_2.0 directory.

To use the new version of the script, use the following command:

$ python3 stack_2.0.py -b bases-spec -x filename.xyz

The user no longer needs to specify the number of bases in the xyz file, but the exact number of lines in the rings specification file must equal the total number of bases.
##########################
##########################

The stack.py script has an option to read pdb or xyz files and calculate stacking score between adjacent bases.
For the pdb file it is enough to just call the script with -p option followed by the name of the file containing the structure in pdb format.

$ python stack.py -p _filename_.pdb

Current version of the script requires more involvement from the user. 
Firstly, it is necessary to provide file containing the atom numbers of atoms from the xyz file, that are forming the aromatic rings of nucleotides. The file should be named bases.txt and it needs to be formated in the same way as the file bases.txt located in the example directory. The order of atom numbers is very crucial, as the shape of polygon representing the base will be constructed according to this order. The numbers should be written in the order that corresponds to the atom numbering when going around the aromatic ring.
It is also necessary to provide number of bases that re to be read from the xyz file while calling the script.

$ python stack.py -x xyz _filename_.xyz -n _number_of_bases_

The calculated stacking scores with additional information will be gathered in the raport.txt file. 

The stack.py script can be tested with db_0 files which are avaiable in the example directory. The example usage of planes.py script for these files would be as follows:

$ python stack.py -p db_0.pdb 

or 

$ python stack.py -x db_0.xyz -n 3

The stacking score for pair of bases is calculated based on algorithm developed by Mikolaj Gurba within research group led by Rafal Szabla. 
The work of Taghavi et al. (J. Chem. Theory Comput. 2022, 18, 3637−3653) was inspiration for some elements of the algorithm.
The stacking score is calculated as a function of three parameters - the distance betwenn the bases, the angle between the planes that describe the bases and the overlap area of the aromatic rings of the bases. 
These parameters are all transformed into values ranging from 0 to 1. 
The stacking score is simple multiplication of those parameters. 
Thus, the stacking score of 1 is corresponding to the perfectly stacked bases, i.e. pair of bases with following features: the distance between the bases is less than 3.5 A, the bases are almost parallel and are compleatly overlaping with each other.
As the value of stacking score approaches the 1, the pair of bases more resembles the perfectly scaked one.
