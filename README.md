# Base_stacking
Repository for files and scripts related to the stack scoring algorithm.

The planes.py script has and option to read pdb or xyz files and calculate stacking score between adjacent bases.
For the pdb file it is enough to just call the script with -p option followed by the name of the file containing the structure in pdb format.

$ python planes.py -p _filename_.pdb

Current version of the script requires more involvement from the user. 
Firstly, it is necessary to provide file containing the atom numbers of atoms from the xyz file, that are forming the aromatic rings of nucleotides. The file should be named bases.txt and it needs to be formated in the same way as the file bases.txt located in the example directory. The order of atom numbers is very crucial, as the shape of polygon representing the base will be constructed according to this order. The numbers should be written in the order that corresponds to the atom numbering when going around the aromatic ring.
It is also necessary to provide number of bases that re to be read from the xyz file while calling the script.

$ python planes.py -x xyz _filename_.xyz -n _number_of_bases_

The calculated stacking scores with additional information will be gathered in the raport.txt file. 

The planes.py script can be tested with db_0 files which are avaiable in the example directory. The example usage of planes.py script for these files would be as follows:

$ python planes.py -p db_0.pdb 

or 

$ python planes.py -x db_0.xyz -n 3

The stacking score for pair of bases is calculated based on algorithm developed by Mikolaj Gurba within research group led by Rafal Szabla. 
The work of Taghavi et al. (J. Chem. Theory Comput. 2022, 18, 3637−3653) was inspiration for some elements of the algorithm.
The stacking score is calculated as a function of three parameters - the distance betwenn the bases, the angle between the planes that describe the bases and the overlap area of the aromatic rings of the bases. 
These parameters are all transformed into values ranging from 0 to 1. 
The stacking score is simple multiplication of those parameters. 
Thus, the stacking score of 1 is corresponding to the perfectly stacked bases, i.e. pair of bases with following features: the distance between the bases is less than 3.5 A, the bases are almost parallel and are compleatly overlaping with each other.
As the value of stacking score approaches the 1, the pair of bases more resembles the perfectly scaked one.
