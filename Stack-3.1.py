#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np 
import math
import argparse
import matplotlib.pyplot as plt
import shapely 
import re
import itertools
from itertools import tee, islice, chain
import Bio
from Bio import PDB
from openbabel import pybel as pb 

# for iteration on list of files 
def previous_and_next(some_iterable):
    prevs, items, nexts = tee(some_iterable, 3)
    prevs = chain([None], prevs)
    nexts = chain(islice(nexts, 1, None), [None])
    return zip(prevs, items, nexts)

# PDB file parser 
# Tool for parsing PDB files
def read_pdb(pdb_filename):
    parser = PDB.PDBParser()
    io = PDB.PDBIO()
    structure_id = f"{pdb_filename}-structure"
    filename = pdb_filename
    structure = parser.get_structure(structure_id, filename)
    return structure
    
# Tool for gathering restypes present in PDB structure
def get_restype(structure):
    restype=[]
    for model in structure:
        for chain in model:
            for residue in chain:
                l = residue.__repr__()
                l = ''.join(l.split())
                l = l.translate({ord(i): None for i in '<Residuht=rsqoc>'})
                restype.append(l)
    return restype

def get_bases_and_rings_xyz(pdb_filename):
    base = read_pdb(pdb_filename)
    extract(base)

# Tool for extracting atom coordinates of bases within the oligonucleotide 
def extract(structure_name):
    substring_1 = "'"
    substring_2 = "P"
    substring_3 = "H"
    full_coor = []
    full_E_names = []
    ring_1_coor = []
    ring_1_E_names = []    
    substring_A = "A"
    substring_T = "T"
    substring_G = "G"
    substring_C = "C"
    substring_U = "U"
    for model in structure_name:
        for chain in model:
            for residue in chain:
                l = residue.__repr__()
                l = ''.join(l.split())
                l = l.translate({ord(i): None for i in '<Residuht=rsqoc>'})
                if substring_A in residue.get_resname():
                    substring_4 = "N6"
                    purine_atom_1 = 'N1'
                    purine_atom_2 = 'C2'
                    purine_atom_3 = 'N3'
                    purine_atom_4 = 'C4'
                    purine_atom_5 = 'C5'
                    purine_atom_6 = 'C6'
                    imidazole_atom_1 = 'C4'
                    imidazole_atom_2 = 'C5'
                    imidazole_atom_3 = 'N7'
                    imidazole_atom_4 = 'C8'
                    imidazole_atom_5 = 'N9'
                    ring_2_coor = []
                    ring_2_E_names = []
                    for atom in residue:
                        if (purine_atom_1 in atom.get_fullname() or purine_atom_2 in atom.get_fullname() or purine_atom_3 in atom.get_fullname() or purine_atom_4 in atom.get_fullname() or purine_atom_5 in atom.get_fullname() or purine_atom_6 in atom.get_fullname()):
                            ring_1_E_names.append(atom.element)
                            x,y,z = atom.get_coord()
                            ring_1_coor.append(x)
                            ring_1_coor.append(y)
                            ring_1_coor.append(z)
                            lst = np.array(ring_1_coor)
                            n = int((len(ring_1_coor)/3))
                            ring_1_xyz = np.reshape(lst, (n,3))
                        if (imidazole_atom_1 in atom.get_fullname() or imidazole_atom_2 in atom.get_fullname() or imidazole_atom_3 in atom.get_fullname() or imidazole_atom_4 in atom.get_fullname() or imidazole_atom_5 in atom.get_fullname()):
                            ring_2_E_names.append(atom.element)
                            x,y,z = atom.get_coord()
                            ring_2_coor.append(x)
                            ring_2_coor.append(y)
                            ring_2_coor.append(z)
                            lst = np.array(ring_2_coor)
                            n = int((len(ring_2_coor)/3))
                            ring_2_xyz = np.reshape(lst, (n,3))
                        if (substring_1 not in atom.get_fullname() and substring_2 not in atom.get_fullname() 
                            and substring_3 not in atom.get_fullname() and substring_4 not in atom.get_fullname()):
                            full_E_names.append(atom.element)
                            x,y,z = atom.get_coord()
                            full_coor.append(x)
                            full_coor.append(y)
                            full_coor.append(z)
                            lst = np.array(full_coor)
                            n = int((len(full_coor)/3))
                            full_xyz = np.reshape(lst, (n,3))
                    write_xyz(f"{l}_ring_1", ring_1_E_names, ring_1_xyz)
                    write_xyz(f"{l}_ring_2", ring_2_E_names, ring_2_xyz)
                    write_xyz(f"{l}_full_ring", full_E_names, full_xyz)
                     
                elif substring_T in residue.get_resname():
                    substring_4 = "C7"
                    substring_5 = "O"
                    purine_atom_1 = 'N1'
                    purine_atom_2 = 'C2'
                    purine_atom_3 = 'N3'
                    purine_atom_4 = 'C4'
                    purine_atom_5 = 'C5'
                    purine_atom_6 = 'C6'
                    for atom in residue:
                        if (purine_atom_1 in atom.get_fullname() or purine_atom_2 in atom.get_fullname() or purine_atom_3 in atom.get_fullname() or purine_atom_4 in atom.get_fullname() or purine_atom_5 in atom.get_fullname() or purine_atom_6 in atom.get_fullname()):
                            ring_1_E_names.append(atom.element)
                            x,y,z = atom.get_coord()
                            ring_1_coor.append(x)
                            ring_1_coor.append(y)
                            ring_1_coor.append(z)
                            lst = np.array(ring_1_coor)
                            n = int((len(ring_1_coor)/3))
                            ring_1_xyz = np.reshape(lst, (n,3))
                        if (substring_1 not in atom.get_fullname() and substring_2 not in atom.get_fullname() 
                            and substring_3 not in atom.get_fullname() and substring_4 not in atom.get_fullname() 
                            and substring_5 not in atom.get_fullname()):
                            full_E_names.append(atom.element)
                            x,y,z = atom.get_coord()
                            full_coor.append(x)
                            full_coor.append(y)
                            full_coor.append(z)
                            lst = np.array(full_coor)
                            n = int((len(full_coor)/3))
                            full_xyz = np.reshape(lst, (n,3))
                    write_xyz(f"{l}_ring_1", ring_1_E_names, ring_1_xyz)
                    write_xyz(f"{l}_full_ring", full_E_names, full_xyz)
                elif substring_G in residue.get_resname():
                    substring_4 = "N2"
                    substring_5 = "O"
                    purine_atom_1 = 'N1'
                    purine_atom_2 = 'C2'
                    purine_atom_3 = 'N3'
                    purine_atom_4 = 'C4'
                    purine_atom_5 = 'C5'
                    purine_atom_6 = 'C6'
                    imidazole_atom_1 = 'C4'
                    imidazole_atom_2 = 'C5'
                    imidazole_atom_3 = 'N7'
                    imidazole_atom_4 = 'C8'
                    imidazole_atom_5 = 'N9'
                    ring_2_coor = []
                    ring_2_E_names = []
                    for atom in residue:
                        if (purine_atom_1 in atom.get_fullname() or purine_atom_2 in atom.get_fullname() or purine_atom_3 in atom.get_fullname() or purine_atom_4 in atom.get_fullname() or purine_atom_5 in atom.get_fullname() or purine_atom_6 in atom.get_fullname()):
                            ring_1_E_names.append(atom.element)
                            x,y,z = atom.get_coord()
                            ring_1_coor.append(x)
                            ring_1_coor.append(y)
                            ring_1_coor.append(z)
                            lst = np.array(ring_1_coor)
                            n = int((len(ring_1_coor)/3))
                            ring_1_xyz = np.reshape(lst, (n,3))
                        if (imidazole_atom_1 in atom.get_fullname() or imidazole_atom_2 in atom.get_fullname() or imidazole_atom_3 in atom.get_fullname() or imidazole_atom_4 in atom.get_fullname() or imidazole_atom_5 in atom.get_fullname()):
                            ring_2_E_names.append(atom.element)
                            x,y,z = atom.get_coord()
                            ring_2_coor.append(x)
                            ring_2_coor.append(y)
                            ring_2_coor.append(z)
                            lst = np.array(ring_2_coor)
                            n = int((len(ring_2_coor)/3))
                            ring_2_xyz = np.reshape(lst, (n,3))
                        if (substring_1 not in atom.get_fullname() and substring_2 not in atom.get_fullname() 
                            and substring_3 not in atom.get_fullname() and substring_4 not in atom.get_fullname() and substring_5 not in atom.get_fullname()):
                            full_E_names.append(atom.element)
                            x,y,z = atom.get_coord()
                            full_coor.append(x)
                            full_coor.append(y)
                            full_coor.append(z)
                            lst = np.array(full_coor)
                            n = int((len(full_coor)/3))
                            full_xyz = np.reshape(lst, (n,3))
                    write_xyz(f"{l}_ring_1", ring_1_E_names, ring_1_xyz)
                    write_xyz(f"{l}_ring_2", ring_2_E_names, ring_2_xyz)
                    write_xyz(f"{l}_full_ring", full_E_names, full_xyz)

                elif substring_C in residue.get_resname():
                    substring_4 = "N4"
                    substring_5 = "O"
                    purine_atom_1 = 'N1'
                    purine_atom_2 = 'C2'
                    purine_atom_3 = 'N3'
                    purine_atom_4 = 'C4'
                    purine_atom_5 = 'C5'
                    purine_atom_6 = 'C6'
                    for atom in residue:
                        if (purine_atom_1 in atom.get_fullname() or purine_atom_2 in atom.get_fullname() or purine_atom_3 in atom.get_fullname() or purine_atom_4 in atom.get_fullname() or purine_atom_5 in atom.get_fullname() or purine_atom_6 in atom.get_fullname()):
                            ring_1_E_names.append(atom.element)
                            x,y,z = atom.get_coord()
                            ring_1_coor.append(x)
                            ring_1_coor.append(y)
                            ring_1_coor.append(z)
                            lst = np.array(ring_1_coor)
                            n = int((len(ring_1_coor)/3))
                            ring_1_xyz = np.reshape(lst, (n,3))
                        if (substring_1 not in atom.get_fullname() and substring_2 not in atom.get_fullname() 
                            and substring_3 not in atom.get_fullname() and substring_4 not in atom.get_fullname() 
                            and substring_5 not in atom.get_fullname()):
                            full_E_names.append(atom.element)
                            x,y,z = atom.get_coord()
                            full_coor.append(x)
                            full_coor.append(y)
                            full_coor.append(z)
                            lst = np.array(full_coor)
                            n = int((len(full_coor)/3))
                            full_xyz = np.reshape(lst, (n,3))
                    write_xyz(f"{l}_ring_1", ring_1_E_names, ring_1_xyz)
                    write_xyz(f"{l}_full_ring", full_E_names, full_xyz)
                elif substring_U in residue.get_resname():
                    substring_4 = "O"
                    purine_atom_1 = 'N1'
                    purine_atom_2 = 'C2'
                    purine_atom_3 = 'N3'
                    purine_atom_4 = 'C4'
                    purine_atom_5 = 'C5'
                    purine_atom_6 = 'C6'
                    for atom in residue:
                        if (purine_atom_1 in atom.get_fullname() or purine_atom_2 in atom.get_fullname() or purine_atom_3 in atom.get_fullname() or purine_atom_4 in atom.get_fullname() or purine_atom_5 in atom.get_fullname() or purine_atom_6 in atom.get_fullname()):
                            ring_1_E_names.append(atom.element)
                            x,y,z = atom.get_coord()
                            ring_1_coor.append(x)
                            ring_1_coor.append(y)
                            ring_1_coor.append(z)
                            lst = np.array(ring_1_coor)
                            n = int((len(ring_1_coor)/3))
                            ring_1_xyz = np.reshape(lst, (n,3))
                        if (substring_1 not in atom.get_fullname() and substring_2 not in atom.get_fullname() 
                            and substring_3 not in atom.get_fullname() and substring_4 not in atom.get_fullname()):
                            full_E_names.append(atom.element)
                            x,y,z = atom.get_coord()
                            full_coor.append(x)
                            full_coor.append(y)
                            full_coor.append(z)
                            lst = np.array(full_coor)
                            n = int((len(full_coor)/3))
                            full_xyz = np.reshape(lst, (n,3))
                    write_xyz(f"{l}_ring_1", ring_1_E_names, ring_1_xyz)
                    write_xyz(f"{l}_full_ring", full_E_names, full_xyz)                           

# Tool for writting xyz files from pdb files
def pdb2xyz(file_name, ele_name, x_coor, y_coor, z_coor):
    n = len(ele_name)
    g = []
    for i in range(0,n):
        g.append(float(x_coor[i]))
        g.append(float(y_coor[i]))
        g.append(float(z_coor[i]))
    lst = np.array(g)
    xyz = np.reshape(lst, (n,3))
    sym = np.array(ele_name)
    e = np.reshape(sym, (n,1))              
    exyz = np.concatenate((e,xyz), axis=1)
    with open(f"{file_name}.xyz","w") as f:
        f.write(f"{n}" + '\n')
        f.write("opt" + '\n')
        for k in range(0,len(ele_name)):
            f.write('\t'.join(map(str, exyz[k])) + '\n')


# XYZ file parser
def get_xyz(file_name):
    with open(file_name, "r") as xyz_file:
        lines = xyz_file.readlines()[2:] # Skipping the first two lines
        atomic_symbols = []
        for line in lines:
            atomic_symbols.append(line.split()[0]) #Extracting atomic symbols
        atomic_coordinates = np.array([line.split()[1:4] for line in lines], dtype=float) #Extracting atomic coordinates as an anrray
    return atomic_symbols, atomic_coordinates 

# Tool to write xyz files from atomic symbol list and atomic coordinate array
def write_xyz(file_name, atom_symbols, atom_coors):
    n = len(atom_symbols)
    sym = np.array(atom_symbols)
    e = np.reshape(sym, (n,1))              
    exyz = np.concatenate((e,atom_coors), axis=1)
    with open(f"{file_name}.xyz","w") as f:
        f.write(f"{n}" + '\n')
        f.write("opt" + '\n')
        for k in range(0, n):
            f.write('\t'.join(map(str, exyz[k])) + '\n')  

# Tool to extract aromatic rings as separate xyz files
def parse_rings_into_xyz(structure, name, structure_len, out_name, xyz_corr, ele_list):
    el = []
    xyz = []
    for j in structure:
        el.append(ele_list[j - 1])
        xyz.append(xyz_corr[j - 1][0])
        xyz.append(xyz_corr[j - 1][1])
        xyz.append(xyz_corr[j - 1][2])
    xyz = np.reshape(xyz, (structure_len,3))
    write_xyz(f"{name}_{out_name}", el, xyz)
    
# Tool to read aromatic ring specifications for xyz file:
def get_rings_from_list(file_name):
    data_dict = {}
    
    with open(file_name, 'r') as file:
        lines = file.readlines()
    
    n  = len(lines)

    for line in lines:
        # Extract ring name (e.g., ring_1)
        ring_match = re.match(r'(\w+)\s*=', line)
        if ring_match:
            ring_name = ring_match.group(1)
        else:
            print("Wrong format of the ring specification file")
            continue  # Skip if no ring name found

        # Find all {...} groups
        groups = re.findall(r'\{([^}]+)\}', line)
        parsed_groups = []
        full_ring = []

        for group in groups:
            # Convert space-separated string of numbers to list of integers
            numbers = list(map(int, group.strip().split()))
            parsed_groups.append(numbers)
            #full = numbers
            if len(groups) == 1:
                full = numbers
                parsed_groups.append(full)
            
        # Add to dictionary
        data_dict[ring_name] = parsed_groups
        
    return data_dict, n

# An algorithm for obtaining the best-fit plane equation for a given set of atom coordinates
def find_best_fitting_plane_and_centroid(points):
    
    # Calculate the coordinates of centroid of given points  
    n = len(points)
    sum_points = np.sum(points, 0)
    cent = sum_points/n
    
    # Formulate the design matrix X and target vector y
    X = np.column_stack((points[:, 0], points[:, 1], np.ones(len(points))))
    y = points[:, 2]

    # Use least squares to find the coefficients vector w
    w, _, _, _ = np.linalg.lstsq(X, y, rcond=None)

    # Extract coefficients from w
    A, B, D = w
    C = -1

    # Return the coefficients of the plane equation
    return A, B, C, D, cent

# TCL filling the aromatic rings with triangles starting form the ring centroid
def tcl_ring_filler(points, cent):
    n = len(points)
    p0 = cent
    for i in range(n -1):
        p1 = points[i]
        p2 = points[i +1]
        line = f"draw_triangle {p0} {p1} {p2} blue 0.5"
        line = line.replace('[', '{')
        line = line.replace(']', '}')
        write_tcl(line)
    p1 = points[n -1]
    p2 = points[0]
    line = f"draw_triangle {p0} {p1} {p2} blue 0.5"
    line = line.replace('[', '{')
    line = line.replace(']', '}')
    write_tcl(line)

# TCL plane
def tcl_plane(plane, cent):
    line = f"plane_from_normal {cent} [{plane[0]} {plane[1]} {plane[2]}] 2.5 green 0.75"
    line = line.replace('[', '{')
    line = line.replace(']', '}')
    write_tcl(line)

# TCL vector
def tcl_vector(point_1, point_2):
    line = f"draw_arrow {point_1} {point_2} 0.03 red"
    line = line.replace('[', '{')
    line = line.replace(']', '}')
    write_tcl(line)

# TCL vectors for angle
def tcl_angle(point, plane_1, plane_2):
    A = plane_1[0]
    B = plane_1[1]
    C = plane_1[2]
    E = plane_2[0]
    F = plane_2[1]
    G = plane_2[2]

    end_1 = []
    for i in range(3):
        end_1.append(point[i] - 2 * plane_1[i])
    
    end_2 = []
    for i in range(3):
        end_2.append(point[i] - 2 * plane_2[i])

    vec_end_1 = np.array(end_1)
    vec_end_2 = np.array(end_2)

    tri_1 = []
    for i in range(3):
        tri_1.append(point[i] - plane_1[i])
    
    tri_2 = []
    for i in range(3):
        tri_2.append(point[i] - plane_2[i])

    tri_end_1 = np.array(tri_1)
    tri_end_2 = np.array(tri_2)

    line_1 = f"draw_arrow {point} {vec_end_1} 0.03 red"
    line_1 = line_1.replace('[', '{')
    line_1 = line_1.replace(']', '}')
    write_tcl(line_1)
    line_2 = f"draw_arrow {point} {vec_end_2} 0.03 blue"
    line_2 = line_2.replace('[', '{')
    line_2 = line_2.replace(']', '}')
    write_tcl(line_2)
    line_3 = f"draw_triangle {point} {tri_end_1} {tri_end_2} purple 0.75"
    line_3 = line_3.replace('[', '{')
    line_3 = line_3.replace(']', '}')
    write_tcl(line_3)


# An algorithm for obtaining the "average" plane equation for two given planes
def middle_plane(A, B, C, D, A2, B2, C2, D2, point_1, point_2):
    S_a = A + A2
    S_b = B + B2
    S_c = C + C2
    a = S_a/(np.sqrt(S_a**2 + S_b**2 + S_c**2))
    b = S_b/(np.sqrt(S_a**2 + S_b**2 + S_c**2))
    c = S_c/(np.sqrt(S_a**2 + S_b**2 + S_c**2))
    point = (point_1 + point_2)/2
    d = -(a*point[0] + b*point[1] + c*point[2])
    return(a, b, c, d, point)

#Tool for calculating the angle between two planes
def find_angle(plane_1_parameter_1, plane_1_parameter_2, plane_1_parameter_3, plane_2_parameter_1, plane_2_parameter_2, plane_2_parameter_3):

    A = plane_1_parameter_1
    B = plane_1_parameter_2
    C = plane_1_parameter_3
    A2 = plane_2_parameter_1
    B2 = plane_2_parameter_2
    C2 = plane_2_parameter_3

    a = A * A2 + B * B2 + C * C2 
    b = np.sqrt(A**2 + B**2 + C**2) 
    c = np.sqrt(A2**2 + B2**2 + C2**2) 
    angle = np.rad2deg(np.arccos(a/(b*c)))
    # Calculation of an acute angle if the above calculation gives an obtuse angle:
    if angle > 90.000:
        acute_angle = 180.000 - angle
    else:
        acute_angle = angle
        
    return(acute_angle) 

#Tool for projecting point on a plane
def project_point(point, plane_cent, plane):
    A = plane[0]
    B = plane[1]
    C = plane[2]
    D = plane[3]
    sqrs = (A**2 + B**2 + C**2) 
    t = (A * (point[0]-plane_cent[0]) + B * (point[1]-plane_cent[1]) + C * (point[2]-plane_cent[2]))/(sqrs)
    translated_point = []
    translated_point.append(point[0] - t * A)
    translated_point.append(point[1] - t * B)
    translated_point.append(point[2] - t * C)
    projected_point = np.array(translated_point)
    return(projected_point)

#Tool for calculating the distance between two rings
def calculate_cent_plane_distance(A, B, C, D, point):
    dis = (np.absolute(A * point[0] + B * point[1] + C * point[2] + D))/(np.sqrt(A**2 + B**2 + C**2)) 
    return(dis)

def find_distance(plane_1_parameter_1, plane_1_parameter_2, plane_1_parameter_3, plane_1_parameter_4, center_of_plane_2, center_of_plane_3):
        
    A = plane_1_parameter_1
    B = plane_1_parameter_2
    C = plane_1_parameter_3
    D = plane_1_parameter_4
    E = center_of_plane_2
    F = center_of_plane_3
        
    cent_plane_cent_dis = (calculate_cent_plane_distance(A, B, C, D, E) + calculate_cent_plane_distance(A, B, C, D, F))
    cent_cent_dis = np.sqrt((E[0]-F[0])**2 + (E[1]-F[1])**2 + (E[2]-F[2])**2)
    return(cent_plane_cent_dis, cent_cent_dis)

def find_distance_parallel(plane_1_parameter_1, plane_1_parameter_2, plane_1_parameter_3, plane_1_parameter_4, plane_2_parameter_4, center_of_plane_2, center_of_plane_3):
        
    A = plane_1_parameter_1
    B = plane_1_parameter_2
    C = plane_1_parameter_3
    D = plane_1_parameter_4
    H = plane_2_parameter_4
    E = center_of_plane_2
    F = center_of_plane_3
        
    cent_plane_cent_dis = (np.absolute(D-H))/(np.sqrt(A**2 + B**2 + C**2))
    cent_cent_dis = np.sqrt((E[0]-F[0])**2 + (E[1]-F[1])**2 + (E[2]-F[2])**2)
    return(cent_plane_cent_dis, cent_cent_dis)

#Tool for calculating the overlap area between two rings
def find_overlap(atoms_in_ring_1, atoms_in_ring_2, plane_1_parameter_1, plane_1_parameter_2, plane_1_parameter_3, plane_1_parameter_4, center_of_plane_1, output_file_name, xyz_file_name):
    ring_xyz = atoms_in_ring_1
    nxt_ring_xyz = atoms_in_ring_2
    A = plane_1_parameter_1
    B = plane_1_parameter_2
    C = plane_1_parameter_3
    D = plane_1_parameter_4
    cent_1 = center_of_plane_1
    name = output_file_name
    
    overlap = calculate_overlap(ring_xyz, nxt_ring_xyz, A, B, C, D, cent_1, name, xyz_file_name)
    plane_3 = [A, B, C, D]
    return overlap
    
#Tool for calculating the overlap area between two rings
# Aligning the coordinates of two rings on the same plane
def calculate_overlap(ring_xyz, nxt_ring_xyz, A, B, C, D, point, name, xyz_file_name):
    c = (A**2 + B**2 + C**2) 
    ring_1 = np.empty((len(ring_xyz),3))
    ring_2 = np.empty((len(nxt_ring_xyz),3))
    for i in range(0,len(ring_xyz)):
        t = (A * (point[0]-ring_xyz[i][0]) + B * (point[1]-ring_xyz[i][1]) + C * (point[2]-ring_xyz[i][2]))/(c)
        ring_1[i][0] = ring_xyz[i][0] + t * A
        ring_1[i][1] = ring_xyz[i][1] + t * B
        ring_1[i][2] = ring_xyz[i][2] + t * C
    # print("######")
    # print("ring_1")
    # print(ring_1)
    # print("######")
    polygon_1 = shapely.Polygon(ring_1)
    for i in range(0,len(nxt_ring_xyz)):
        t = (A * (point[0]-nxt_ring_xyz[i][0]) + B * (point[1]-nxt_ring_xyz[i][1]) + C * (point[2]-nxt_ring_xyz[i][2]))/(c)
        ring_2[i][0] = nxt_ring_xyz[i][0] + t * A
        ring_2[i][1] = nxt_ring_xyz[i][1] + t * B
        ring_2[i][2] = nxt_ring_xyz[i][2] + t * C
    polygon_2 = shapely.Polygon(ring_2)
    # print("######")
    # print("ring_2")
    # print(ring_2)
    # print("######")
    # Calculating intersection area between two planes 
    if polygon_2.area >= polygon_1.area:       
        intersect = shapely.intersection(polygon_1, polygon_2)
        overlap = (polygon_1.intersection(polygon_2).area/polygon_1.area)
        overlap_p = ((polygon_1.intersection(polygon_2).area/polygon_1.area)*100)
    else: 
        intersect = shapely.intersection(polygon_2, polygon_1)
        overlap = (polygon_2.intersection(polygon_1).area/polygon_2.area)
        overlap_p = ((polygon_2.intersection(polygon_1).area/polygon_2.area)*100)
    try:
        x,y = polygon_1.exterior.xy
        u,v = polygon_2.exterior.xy
        r,t = intersect.exterior.xy    
    
        plt.figure()
        plt.plot(x,y, 'r')
        plt.plot(u,v, 'b')
        plt.plot(r,t, 'g')
        plt.fill_between(r,t, alpha=0.3)
        #plt.fill_between(r,t)
        plt.xlabel("Contour of the overlap area (green) between two rings (blue and red)")
        plt.show
        plt.savefig(f"{xyz_file_name}_{name}_overlap_shape.pdf")

    except:
        print("An error occured while drawing overlap scheme")

    return(overlap)

#Tool for calculating the value of the D parameter based on the ring-ring distance
def dis2D(distance):
    if distance <= 3.371959:
        D = 1
    elif distance >= 5.00:
        D = 0
    else:
        D = -3*(0.39109629*(np.exp(-2*1.51146861*(distance-3.37195972))-2*np.exp(-1.51146861*(distance-3.37195972))) + 0.06449501)
    return D

#Tool for calculating the value of the A parameter based on the ring-ring angle
def ang2A(angle):
    r_angle = np.deg2rad(angle)
    A = np.exp(-1 * r_angle**4) + np.exp(-(r_angle - np.pi)**4) + 0.1 * np.exp(-(r_angle - 0.5*np.pi)**4)
    return A
    
#Tool for calculating the stacking score based on the full ring overlap (old score)
def old_stack_score(distance, angle, overlap):
    D = dis2D(distance)
    A = ang2A(angle)
    O = overlap

    stack_score = D*A*O
    return stack_score

#Tool for calculating the stacking score based on the sum of squared subring overlap values (new score)
def new_stack_score(distance, angle, divided_sum):
    D = dis2D(distance)
    A = ang2A(angle)
    O = divided_sum

    stack_score = D*A*O
    return stack_score

#Tool for writing the final program report
def write_report(prtlvl, base_1, base_2, name, full_O, DARO, reversed_O, angle, plane_1, plane_2, plane_3, cent_plane_cent_dis, cent_cent_dis, old_score, new_score, rev_score):
    with open("Stack.log","a") as f:
        f.write(f"Results for structure {name}" + '\n')       
        f.write(f"Calculating the stacking parameters between {base_1} and {base_2}" + '\n')
        if prtlvl > 0:
            f.write(f"The parameters of the plane equation for the first ring are: {plane_1}" + '\n')
            f.write(f"The parameters of the plane equation for the second ring are: {plane_2}" + '\n')
            f.write(f"The parameters of the mean plane equation are as follows: {plane_3}" + '\n' )
        f.write(f"The distance between the two aromatic rings is: {cent_plane_cent_dis}" + '\n')
        f.write(f"The distance between the two aromatic ring centroids is: {cent_cent_dis}" + '\n')
        f.write(f"The angle between the first and the second plane is: {angle}" + '\n')
        f.write(f"The degree of full ring overlap between {base_1} and {base_2} is: {full_O}" + '\n')
        f.write(f"The degree of aromatic ring overlap (DARO) between {base_1} and {base_2} is: {DARO}" + '\n')
        # f.write(f"The reversed and refined value of the aromatic ring overlap between {base_1} and {base_2} is: {reversed_O}" + '\n')
        if prtlvl > 0:
            f.write(f"The old stacking score for this base pair is: {old_score}" + '\n')
        f.write(f"The new stacking score for this base pair is: {new_score}" + '\n')
        # f.write(f"The reversed O new stacking score for this base pair is: {rev_score}" + '\n')
        f.write('\n')

#Tool for generating tcl script
def write_tcl(line):
    with open("tcl_figure_parm.tcl", "a") as f:
        f.write(line + '\n')
        f.write('\n')

# Main program
def find_stacking_xyz(prtlvl, xyz_file_name, ring_file_name):

    rings, number_of_rings = get_rings_from_list(ring_file_name)
    elements, xyz_corr = get_xyz(xyz_file_name)
    
    nn = []
    for name, values in rings.items():
        number_of_subrings_in_line = len(values)
        nn.append(number_of_subrings_in_line)
        for i in range(number_of_subrings_in_line - 1):
            lenn = len(rings[name][i])
            parse_rings_into_xyz(rings[name][i], name, lenn, f"subring_{i}", xyz_corr, elements)
        for z in range(number_of_subrings_in_line - 1, number_of_subrings_in_line):
            lenz = len(rings[name][z])
            parse_rings_into_xyz(rings[name][z], name, lenz, "full_ring", xyz_corr, elements)

    for i in range(1, number_of_rings):
        # Obtaining the plane equation for the first ring
        e, full_ring_xyz = get_xyz(f"ring_{i}_full_ring.xyz")
        A, B, C, D , cent_1 = find_best_fitting_plane_and_centroid(full_ring_xyz)
        plane_1 = [A, B, C, D]
        # btaining the plane equation for the second ring
        e, nxt_full_ring_xyz = get_xyz(f"ring_{i+1}_full_ring.xyz")
        E, F, G, H, cent_2 = find_best_fitting_plane_and_centroid(nxt_full_ring_xyz)        
        plane_2 = [E, F, G, H]
        # Calculating the angle between planes
        angle = find_angle(A, B, C, E, F, G)
        tcl_angle(cent_1, plane_1, plane_2)
        # Calculating overlap between full rings 
        if math.isclose(angle, 0.0, rel_tol=1e-5) == True or math.isclose(angle, 180.0, rel_tol=1e-5) == True:
            cent_plane_cent_dis, cent_cent_dis = find_distance_parallel(A, B, C, D, H, cent_1, cent_2)
            full_O = find_overlap(full_ring_xyz, nxt_full_ring_xyz, A, B, C, D, cent_1, f"full_ring_overlap_{i}_{i+1}", xyz_file_name)
            plane_3 = [E, F, G, H]
            tcl_plane(plane_3, cent_2)
            pp = project_point(cent_1, cent_2, plane_3)
            tcl_vector(cent_1, pp)
        else:
            I, J, K, L, point = middle_plane(A, B, C, D, E, F, G, H, cent_1, cent_2)
            plane_3 = [I, J, K, L]
            tcl_plane(plane_3, point)
            pp_1 = project_point(cent_1, point, plane_3)
            tcl_vector(cent_1, pp_1)
            pp_2 = project_point(cent_2, point, plane_3)
            tcl_vector(cent_2, pp_2)
            cent_plane_cent_dis, cent_cent_dis = find_distance(I, J, K, L, cent_1, cent_2)
            full_O = find_overlap(full_ring_xyz, nxt_full_ring_xyz, I, J, K, L, point, f"full_ring_overlap_{i}_{i+1}", xyz_file_name)

        # Calculating stack score with old algorithm
        old_score = old_stack_score(cent_plane_cent_dis, angle, full_O)
        
        LK = 0
        Sum_of_DARO_squared = 0
        for j in range(nn[i-1] - 1):
            for k in range(nn[i] - 1):
                # Obtaining the XYZ coordinates of first subring 
                e, ring_xyz = get_xyz(f"ring_{i}_subring_{j}.xyz")
                # Obtaining plane equation for first subring
                A, B, C, D , cent_1 = find_best_fitting_plane_and_centroid(ring_xyz)
                tcl_ring_filler(ring_xyz, cent_1)
                plane_1 = [A, B, C, D]
                # Obtaining XYZ coordinates of second subring 
                e, nxt_ring_xyz = get_xyz(f"ring_{i+1}_subring_{k}.xyz")
                # Obtaining plane equation that describes the centroid of second subring 
                E, F, G, H, cent_2 = find_best_fitting_plane_and_centroid(nxt_ring_xyz)
                tcl_ring_filler(nxt_ring_xyz, cent_2)
                plane_2 = [E, F, G, H]
                # Calculating overlap between individual subrings 
                if math.isclose(angle, 0.0, rel_tol=1e-5) == True or math.isclose(angle, 180.0, rel_tol=1e-5) == True:
                    DARO = find_overlap(ring_xyz, nxt_ring_xyz, A, B, C, D, cent_1, f"ring_{i}_subring_{j}_{k}", xyz_file_name)
                else:
                    DARO = find_overlap(ring_xyz, nxt_ring_xyz, I, J, K, L, point, f"ring_{i}_subring_{j}_{k}", xyz_file_name)
                    
                Sum_of_DARO_squared = Sum_of_DARO_squared + DARO**2
                #print(Sum_of_DARO_squared)
                if DARO > 0.0001:
                    LK = LK + 1

        # Calculating stack score with new algorithm
        div_sum = Sum_of_DARO_squared/(min(nn) - 1)
        #print(nn, min(nn), div_sum)
        new_score = new_stack_score(cent_plane_cent_dis, angle, div_sum)

        # Calculating reversed and refined overlap
        # if LK > 1:
        #     reversed_O = 1 - div_sum
        # else:
        #     reversed_O = full_O
        
        reversed_O = 1 - div_sum
        rev_score = new_stack_score(cent_plane_cent_dis, angle, reversed_O)

        write_report(prtlvl, f"ring_{i}", f"ring_{i+1}", xyz_file_name, full_O, div_sum, reversed_O, angle, plane_1, plane_2, plane_3, cent_plane_cent_dis, cent_cent_dis, old_score, new_score, rev_score)

def find_stacking_pdb(prtlvl, pdb_file_name):

#    Read PDB file and extract residue names 
    Select = PDB.Select
    io = PDB.PDBIO()
    sequence = read_pdb(pdb_file_name)
    res_list = get_restype(sequence)    

    # Extract bases and save them in individual PDB files 
    for i in res_list:
        class name(Select):
            def accept_residue(self,residue):
                name = residue.__repr__()
                name = ''.join(name.split())
                name = name.translate({ord(i): None for i in '<Residuht=rsqoc>'})
                if name == f"{i}":
                    return True
                else:
                    return False
            def accept_atom(self,atom):
                substring_1 = "'"
                substring_2 = "P"
                if (substring_1 not in atom.get_fullname()) and (substring_2 not in atom.get_fullname()):
                    return True
                else:
                    return False
        inp = f"{i}.pdb"
        io.set_structure(sequence)
        io.save(inp, name())

    for i in res_list:
        inp = f"{i}.pdb"
        get_bases_and_rings_xyz(inp)

    for previous, item, nxt in previous_and_next(res_list):
        if nxt != None:
            # Obtaining plane equation for first ring
            e, full_ring_xyz = get_xyz(f"{item}_full_ring.xyz")
            A, B, C, D , cent_1 = find_best_fitting_plane_and_centroid(full_ring_xyz)
            plane_1 = [A, B, C, D]
            # Obtaining XYZ coordinates of second ring 
            e, nxt_full_ring_xyz = get_xyz(f"{nxt}_full_ring.xyz")
            # Obtaining plane equation that describes the centroid of second ring 
            E, F, G, H, cent_2 = find_best_fitting_plane_and_centroid(nxt_full_ring_xyz)
            plane_2 = [E, F, G, H]
            # Calculating angle between planes 
            angle = find_angle(A, B, C, E, F, G)
            tcl_angle(cent_1, plane_1, plane_2)
            # Calculating overlap between full rings 
            if math.isclose(angle, 0.0, rel_tol=1e-5) == True or math.isclose(angle, 180.0, rel_tol=1e-5) == True:
                cent_plane_cent_dis, cent_cent_dis = find_distance_parallel(A, B, C, D, H, cent_1, cent_2)
                full_O = find_overlap(full_ring_xyz, nxt_full_ring_xyz, A, B, C, D, cent_1, f"full_ring_overlap_{item}_{nxt}", pdb_file_name)
                plane_3 = [E, F, G, H]
                tcl_plane(plane_3, cent_2)
                pp = project_point(cent_1, cent_2, plane_3)
                tcl_vector(cent_1, pp)
            else:
                I, J, K, L, point = middle_plane(A, B, C, D, E, F, G, H, cent_1, cent_2)
                plane_3 = [I, J, K, L]
                tcl_plane(plane_3, point)
                pp_1 = project_point(cent_1, point, plane_3)
                tcl_vector(cent_1, pp_1)
                pp_2 = project_point(cent_2, point, plane_3)
                tcl_vector(cent_2, pp_2)
                cent_plane_cent_dis, cent_cent_dis = find_distance(I, J, K, L, cent_1, cent_2)
                full_O = find_overlap(full_ring_xyz, nxt_full_ring_xyz, I, J, K, L, point, f"full_ring_overlap_{item}_{nxt}", pdb_file_name)
            #Clauclate old stack score 
            old_score = old_stack_score(cent_plane_cent_dis, angle, full_O)
        
            #Claculate new stack score        
            Sum_of_DARO_squared = 0
            LK = 0
            nn = []
            if ("A" in item or "G" in item):
                number_of_rings_in_first_base = 2
            elif ("C" in item or "T" in item or "U" in item):
                number_of_rings_in_first_base = 1
            nn.append(number_of_rings_in_first_base)

            if ("A" in nxt or "G" in nxt):
                number_of_rings_in_second_base = 2
            elif ("C" in nxt or "T" in nxt or "U" in nxt):
                number_of_rings_in_second_base = 1
            nn.append(number_of_rings_in_second_base)

            for j in range(number_of_rings_in_first_base):
                for k in range(number_of_rings_in_second_base):
                    # Obtaining the XYZ coordinates of first ring 
                    e, ring_xyz = get_xyz(f"{item}_ring_{j+1}.xyz")
                    # Obtaining plane equation for first ring
                    A, B, C, D , cent_1 = find_best_fitting_plane_and_centroid(ring_xyz)
                    tcl_ring_filler(ring_xyz, cent_1)
                    plane_1 = [A, B, C, D]
                    # Obtaining XYZ coordinates of second ring 
                    e, nxt_ring_xyz = get_xyz(f"{nxt}_ring_{k+1}.xyz")
                    # Obtaining plane equation that describes the centroid of second ring 
                    E, F, G, H, cent_2 = find_best_fitting_plane_and_centroid(nxt_ring_xyz)
                    tcl_ring_filler(nxt_ring_xyz, cent_2)
                    plane_2 = [E, F, G, H]
                    # Calculating overlap between individual rings 
                    if math.isclose(angle, 0.0, rel_tol=1e-5) == True or math.isclose(angle, 180.0, rel_tol=1e-5) == True:
                        DARO = find_overlap(ring_xyz, nxt_ring_xyz, A, B, C, D, cent_1, f"{item}_ring_{j+1}_{k+1}", pdb_file_name)
                    else:                 
                        DARO = find_overlap(ring_xyz, nxt_ring_xyz, I, J, K, L, point, f"{item}_ring_{j+1}_{k+1}", pdb_file_name)

                    Sum_of_DARO_squared = Sum_of_DARO_squared + DARO**2
                    #print(Sum_of_DARO_squared)
                    if DARO > 0.0001:
                        LK = LK + 1
    
            # Calculating stack score with new algorithm
            div_sum = Sum_of_DARO_squared/(min(nn))
            #print(nn, min(nn), div_sum)
            new_score = new_stack_score(cent_plane_cent_dis, angle, div_sum)
            reversed_O = 1 - div_sum
            rev_score = new_stack_score(cent_plane_cent_dis, angle, reversed_O)

            write_report(prtlvl, f"{item}", f"{nxt}", pdb_file_name, full_O, div_sum, reversed_O, angle, plane_1, plane_2, plane_3, cent_plane_cent_dis, cent_cent_dis, old_score, new_score, rev_score)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Claculate the stacking score between adjacent aromatic rings from an XYZ or PDB file.")
    # parser.add_argument( "-n", "--name", type=str, metavar = "", required = True, help = "name of uppdated PDB file")
    parser.add_argument( "-p", "--pdb_filename", type=str, metavar = "" ,required = False, help = "name of the PDB file")

    parser.add_argument( "-l", "--print_level", type=str, metavar = "", required = False, help = "set to full to print additional information in the output file")
    parser.add_argument( "-b", "--ring_specification_filename", type=str, metavar = "", required = False, help = "name of the file containing the atom numbers defining the aromatic rings in the system")
    parser.add_argument( "-x", "--xyz_filename", type=str, metavar = "", required = False, help = "name of the XYZ file")
    args = parser.parse_args()

    if args.print_level == "full":
        prtlvl = 1
    else:
        prtlvl = 0

    #find_stacking_xyz(args.xyz_filename, args.base_specification_filename)

    if args.pdb_filename != None:
        find_stacking_pdb(prtlvl, args.pdb_filename)
    elif args.xyz_filename != None:
        find_stacking_xyz(prtlvl, args.xyz_filename, args.ring_specification_filename)
    else:
        print("The file name and type were not specified. Please, use -h option to get help")


        

