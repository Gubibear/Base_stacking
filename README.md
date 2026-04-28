Stack-3.1.py is the current version of code developed by Mikołaj Gurba in Rafał Szabla's research group. This code allows for the calculation of the geometric parameters of a system containing stacked aromatic rings and the estimation of the stacking score. The value of the stacking coefficient can be related to the strength of the interaction between these rings, as well as to the properties of the excited states present in such systems.

#############  Description of the algorithm  #############

For a given pair of aromatic rings, three main parameters are calculated:
- the distance between the rings,
- the angle at which the rings are oriented relative to each other,
- the area of ​​geometric overlap between these rings.

A description of these parameters is provided below:

1. Distance
The distance between rings can be defined in various ways. The most common definition in the literature is the distance between the centers of mass or the centroids of the rings. Because this distance does not precisely describe the distance by which the different parts of the stacked rings are separated, a different definition is considered in this code. Each ring is described by a plane fitted to the positions of the heavy atoms present in the aromatic rings using the least-squares method. From these two planes, the mean plane of the two rings is constructed midway between their centroids. The distance between the rings is calculated as the sum of the distances from the center of mass of each ring to the mean plane.

2. Angle
The angle parameter is straightforward, it is calculated as the angle between the normal vectors of the planes describing the rings. For each ring in the system the equation of the plane describing them is obtained from the leas-squares approach. Only the heavy atoms of the aromatic rings are considered in the fitting.

4. Overlap Area
The code currently uses two methods to calculate the overlap area. Both project the rings onto a mean plane and calculate their intersection area. The first method calculates the intersection area of ​​the full aromatic ring shapes. To obtain the overlap coefficient, this area is then divided by the area of ​​the smaller of the stacked rings. The second method first divides polycyclic aromatic rings into their smallest subunits (e.g., a purine ring is split into pyrimidine and imidazole) and then calculates the overlap area of ​​all the subunits. The overlap coefficient of each pair of subunits is calculated as the ratio of their intersection area to the area of ​​the smaller subunit. The second method allows to calculate the degree of aromatic ring overlap (DARO), which is defined as the sum of the squared overlap ratios of all subunits divided by the minimum number of subunits in one of the stacked rings.

For each of these three parameters, a score value ranging from 0 to 1 is calculated. Then, as the product of these parameter values, the stacking score for the ring pair is calculated. 
The code now provides two score values: Old and New. 
Old Score value is calculated based on the overlap parameter obtained by calculating the overlap of the entire aromatic rings. 
New Score value uses the DARO definition as the overlap parameter.
Both stacking scores take values ​​from 0 to 1. A value of 1 is obtained for two aromatic rings that are no more than 3.4 angstroms apart, are oriented almost parallel to each other, and completely overlap (or two subunits from different rings completely overlap - this only describes the New Score).

#############  How to use the Stack-3.1.py  #############

Current version of the code can calculate the stacking score (together with the distance, angle and overlap area) between two aromatic rings present in the xyz file. Since the driving incentive for the developement of the algorithm and code was studies of stacking between the nucleobases in nucleic acid structures, the code can also read the pdb file which contains single strand of DNA or RNA, but only with canonical bases. The description of the usage of the code for the xyz and pdb files will be provided bellow.


XYZ files 
The code *should* work for all structures with stacked aromatic rings. In case of problems please reach out to the contributor through github issue or directly at: mikolaj.gurba@pwr.edu.pl.

To use the Stack-3.1.py code for calulation of stacking score between rings found in the xyz file 1 additional file must be prepared by the user beforehand. In this txt file the user must specify the atom numbers (corresponding to those in the xyz file) for each atom present in the aromatic rings that are to be inculed in the calculation. The structure of this file and the order of the atom numbers must be provided as requested or the code will crash. Each line of this file corresponds to different aromatic ring in the system and must start with the "ring_n = " phrase, where n is the number of the line (e. g. 1, 2, etc.). Each ring is defined with atom numbers provided in currly brackets. The order of these numbers is very important, as it will define the shape of the ring. The atoms must be provided in order which allows to draw the shape of the ring without any line crossing. For polycyclic rings brackets separated by coma must be provided for each monocyclic subunit. Additionally, for polycyclic rings the full ring must be specified in the last currly bracket. To better explain the structure of this additional file two examples will be detailed.

Aromatic ring specification for adenine:
<img width="960" height="965" alt="image" src="https://github.com/user-attachments/assets/c33dffc9-8ab2-4274-be75-473a5542f186" />
For the adenine molecule the line in the txt file specyfing the aromatic ring will be as follows:
  ring_1 = {1, 2, 3, 4, 5, 6}, {5, 4, 9, 8, 7}, {1, 2, 3, 4, 9, 8, 7, 5, 6}


pdb files: this works for all pdb files containing canonical nucleobases 
xyz files: all molecules will work, but additional file is necessary to specify the order of atoms in aromatic rings. further details in example description.



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
