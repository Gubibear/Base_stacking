This directory contains the input and output files for the analysis of the 5'-AC-3' DNA dinucleotide using the Stack-3.1.py code.
The dna.xyz file contains the xyz coordinates for the AC dinucleotide, while the aromatic rings are specified in the dna_ring_atom_numbers file.
The following command was used to execute the code:

python Stack-3.1.py -x dna.xyz -b dna_ring_atom_numbers -l full

The analysis results are saved in the Stack.log file. The program also provides PDF files containing approximate contours of the overlapping rings, visualized using the matplotlib package, and a tcl file that can be used to visualize planes and vectors using the VMD package.

Additionally, temporary files are saved – xyz files with the coordinates of the atoms in each aromatic ring.
