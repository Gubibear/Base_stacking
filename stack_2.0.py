#!/usr/bin/env python
# coding: utf-8

# In[7]:


import numpy as np 
import math
import urllib.request
import Bio
from Bio import PDB
from openbabel import pybel as pb 
import argparse
import matplotlib.pyplot as plt
import shapely 
import itertools
from itertools import tee, islice, chain
import re

# Tool for reading XYZ files
def get_xyz(file_name):
    d = []
    e = []
    x_coor = []
    y_coor = []
    z_coor = []
    g = []
    with open(file_name, "r") as z:
        for line in z:
            d = line.split()
            if len(d) > 0:
                e.append(d[0])
            if len(d) > 1 and len(d[0]) == 1:
                g.append(float(d[1]))
                g.append(float(d[2]))
                g.append(float(d[3]))
    num = int(e.pop(0))
    if len(e) > num:
        e.pop(0)
    lst = np.array(g)
    xyz = np.reshape(lst, (num, 3))
    return e, xyz 

# Tool for writting xyz files from xyz files
def write_xyz(file_name, ele_name, coors):
    n = len(ele_name)
    sym = np.array(ele_name)
    e = np.reshape(sym, (n,1))              
    exyz = np.concatenate((e,coors), axis=1)
    with open(f"{file_name}.xyz","w") as f:
        f.write(f"{n}" + '\n')
        f.write("opt" + '\n')
        for k in range(0,len(ele_name)):
            f.write('\t'.join(map(str, exyz[k])) + '\n')   

def parse_rings_into_xyz(structure, name, structure_len, out_name, xyz_corr, ele_list):
    # print(name)
    el = []
    xyz = []
    for j in structure:
        el.append(ele_list[j - 1])
        xyz.append(xyz_corr[j - 1][0])
        xyz.append(xyz_corr[j - 1][1])
        xyz.append(xyz_corr[j - 1][2])
    xyz = np.reshape(xyz, (structure_len,3))
    write_xyz(f"{name}_{out_name}", el, xyz)

#Tool for reading aromatic rings specifications for xyz file:
def get_rings_from_list(file_name):
    data_dict = {}
    
    with open(file_name, 'r') as file:
        lines = file.readlines()
    
    n  = len(lines)

    for line in lines:
        # Extract base name (e.g., base_1)
        base_match = re.match(r'(\w+)\s*=', line)
        if base_match:
            base_name = base_match.group(1)
        else:
            print("Wrong format of the base specification file")
            continue  # Skip if no base name found

        # Find all {...} groups
        groups = re.findall(r'\{([^}]+)\}', line)
        parsed_groups = []
        full_ring = []

        for group in groups:
            # Convert space-separated string of numbers to list of integers
            numbers = list(map(int, group.strip().split()))
            parsed_groups.append(numbers)
            full = numbers
            if len(groups) > 1:
                for i in range(len(full) - 1):
                    full_ring.append(full[i])
            else:
                full_ring = full
        parsed_groups.append(full_ring)
            
        # Add to dictionary
        data_dict[base_name] = parsed_groups
        
    return data_dict, n

# Algorithm for obtaining plane equation that fits best for given base 
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

def calculate_angle(A, B, C, A2, B2, C2):
    a = A * A2 + B * B2 + C * C2 
    b = np.sqrt(A**2 + B**2 + C**2) 
    c = np.sqrt(A2**2 + B2**2 + C2**2) 
    teta = np.rad2deg(np.arccos(a/(b*c)))
    return(teta)

def find_angle(plane_1_parameter_1, plane_1_parameter_2, plane_1_parameter_3, plane_2_parameter_1, plane_2_parameter_2, plane_2_parameter_3):

    A = plane_1_parameter_1
    B = plane_1_parameter_2
    C = plane_1_parameter_3
    E = plane_2_parameter_1
    F = plane_2_parameter_2
    G = plane_2_parameter_3

    angle = calculate_angle(A, B, C, E, F, G)
    r_angle = np.deg2rad(angle)
    rot = np.exp(-1 * r_angle**4) + np.exp(-(r_angle - np.pi)**4) + 0.1 * np.exp(-(r_angle - 0.5*np.pi)**4)
    return(angle, rot)

def calculate_distance(A, B, C, D, point):
    dis = (np.absolute(A * point[0] + B * point[1] + C * point[2] + D))/(np.sqrt(A**2 + B**2 + C**2)) 
    return(dis)

def find_distance(plane_1_parameter_1, plane_1_parameter_2, plane_1_parameter_3, plane_1_parameter_4, center_of_plane_2, center_of_plane_3):
        
    A = plane_1_parameter_1
    B = plane_1_parameter_2
    C = plane_1_parameter_3
    D = plane_1_parameter_4
    E = center_of_plane_2
    F = center_of_plane_3
        
    dis = (calculate_distance(A, B, C, D, E) + calculate_distance(A, B, C, D, F))
    D = dis2D(dis)
    return(dis, D)

def find_distance_parallel(plane_1_parameter_1, plane_1_parameter_2, plane_1_parameter_3, plane_1_parameter_4, plane_2_parameter_4):
        
    A = plane_1_parameter_1
    B = plane_1_parameter_2
    C = plane_1_parameter_3
    D = plane_1_parameter_4
    H = plane_2_parameter_4
        
    dis = (np.absolute(D-H))/(np.sqrt(A**2 + B**2 + C**2))
    D = dis2D(dis)
    return(dis, D)

def dis2D(d):
    if d <= 3.350:
        D = 1
    elif d >= 5.00:
        D = 0
    elif d > 3.35 and d < 3.37195972:
        D = -0.9196911436029287*d + 4.080965331069811
    else:
        D = -3*(0.39109629*(np.exp(-2*1.51146861*(d-3.37195972))-2*np.exp(-1.51146861*(d-3.37195972))) + 0.06449501)
    return D

def find_overlap(atoms_in_ring_1, atoms_in_ring_2, plane_1_parameter_1, plane_1_parameter_2, plane_1_parameter_3, plane_1_parameter_4, center_of_plane_1, output_file_name):
    # Calculating overlap between bases
    ring_xyz = atoms_in_ring_1
    nxt_ring_xyz = atoms_in_ring_2
    A = plane_1_parameter_1
    B = plane_1_parameter_2
    C = plane_1_parameter_3
    D = plane_1_parameter_4
    cent_1 = center_of_plane_1
    name = output_file_name
    
    overlap = calculate_overlap(ring_xyz, nxt_ring_xyz, A, B, C, D, cent_1, name)
    plane_3 = [A, B, C, D]
    return overlap

# Alligning coordinates of the two bases in the same plane and gettin their shapes in shapely format
def calculate_overlap(base_xyz, nxt_base_xyz, A, B, C, D, point, name):
    c = (A**2 + B**2 + C**2) 
    base_1 = np.empty((len(base_xyz),3))
    base_2 = np.empty((len(nxt_base_xyz),3))
    for i in range(0,len(base_xyz)):
        t = (A * (point[0]-base_xyz[i][0]) + B * (point[1]-base_xyz[i][1]) + C * (point[2]-base_xyz[i][2]))/(c)
        base_1[i][0] = base_xyz[i][0] + t * A
        base_1[i][1] = base_xyz[i][1] + t * B
        base_1[i][2] = base_xyz[i][2] + t * C
    polygon_1 = shapely.Polygon(base_1)
    for i in range(0,len(nxt_base_xyz)):
        t = (A * (point[0]-nxt_base_xyz[i][0]) + B * (point[1]-nxt_base_xyz[i][1]) + C * (point[2]-nxt_base_xyz[i][2]))/(c)
        base_2[i][0] = nxt_base_xyz[i][0] + t * A
        base_2[i][1] = nxt_base_xyz[i][1] + t * B
        base_2[i][2] = nxt_base_xyz[i][2] + t * C
    polygon_2 = shapely.Polygon(base_2)
    
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
        plt.plot(x,y)
        plt.plot(u,v)
        plt.plot(r,t)
        plt.fill_between(r,t, alpha=0.3)
        #plt.fill_between(r,t)
        plt.xlabel("Contour of the overlap area (green) between two bases (blue and orange)")
        plt.show
        plt.savefig(f"{name}_overlap_shape.pdf")

    except:
        print("An error occured while drawing overlap scheme")

    return(overlap)

#Tool for writting the final raport of the program
def write_raport(base_1, base_2, new_stack, old_stack, name):
    with open("stack.out","a") as f:
        f.write(f"Results for structure {name}" + '\n')       
        f.write(f"Stacking score calculation between {base_1} and {base_2}" + '\n')
        # f.write(f"The parameters of the plane equation for the first base are: {plane_1}" + '\n')
        # f.write(f"The parameters of the plane equation for the second base are: {plane_2}" + '\n')
        # f.write(f"The parameters of the average plane equation are: {plane_3}" + '\n' )
        # f.write(f"The distance between the two aromatic rings is: {dis}" + '\n')
        # f.write(f"The angle between the first and the second plane is: {angle}" + '\n')
        # f.write(f"The percentage of the shared area between two two bases is: {overlap}" + '\n')
        f.write(f"The old stacking score between {base_1} and {base_2} is: {old_stack}" + '\n')
        f.write(f"The new stacking score between {base_1} and {base_2} is: {new_stack}" + '\n')
        f.write('\n')

def find_stacking_xyz(xyz_file_name, base_file_name):

    bases, number_of_bases = get_rings_from_list(base_file_name)
    elements, xyz_corr = get_xyz(xyz_file_name)
    
    for name, values in bases.items():
        number_of_rings_in_line = len(values)
        for i in range(number_of_rings_in_line - 1):
            lenn = len(bases[name][i])
            parse_rings_into_xyz(bases[name][i], name, lenn, f"ring_{i}", xyz_corr, elements)
        for z in range(number_of_rings_in_line - 1, number_of_rings_in_line):
            lenz = len(bases[name][z])
            parse_rings_into_xyz(bases[name][z], name, lenz, "full_ring", xyz_corr, elements)

    for i in range(1, number_of_bases):
        # Obtaining plane equation for first ring
        e, full_ring_xyz = get_xyz(f"base_{i}_full_ring.xyz")
        A, B, C, D , cent_1 = find_best_fitting_plane_and_centroid(full_ring_xyz)
        plane_1 = [A, B, C, D]
        # Obtaining XYZ coordinates of second ring 
        e, nxt_full_ring_xyz = get_xyz(f"base_{i+1}_full_ring.xyz")
        # Obtaining plane equation that describes the centroid of second ring 
        E, F, G, H, cent_2 = find_best_fitting_plane_and_centroid(nxt_full_ring_xyz)
        plane_2 = [E, F, G, H]
        # Calculating angle between planes 
        angle, parameter_A = find_angle(A, B, C, E, F, G)
        # Calculating overlap between full rings 
        if math.isclose(angle, 0.0, rel_tol=1e-5) == True or math.isclose(angle, 180.0, rel_tol=1e-5) == True:
            dis, parameter_D  = find_distance_parallel(A, B, C, D, H)
            full_overlap = find_overlap(full_ring_xyz, nxt_full_ring_xyz, A, B, C, D, cent_1, f"full_ring_overlap_{i}_{i+1}")

        else:
            I, J, K, L, point = middle_plane(A, B, C, D, E, F, G, H, cent_1, cent_2)
            plane_3 = [I, J, K, L]
            dis, parameter_D = find_distance(I, J, K, L, cent_1, cent_2)
            full_overlap = find_overlap(full_ring_xyz, nxt_full_ring_xyz, I, J, K, L, point, f"full_ring_overlap_{i}_{i+1}")
        #Clauclate old stack score 
        old_stack_score = parameter_D*parameter_A*full_overlap
        #Claculate new stack score        
        new_stack_score = 0
        for j in range(number_of_rings_in_line - 1):
            for k in range(number_of_rings_in_line - 1):
                # Obtaining the XYZ coordinates of first ring 
                e, ring_xyz = get_xyz(f"base_{i}_ring_{j}.xyz")
                # Obtaining plane equation for first ring
                A, B, C, D , cent_1 = find_best_fitting_plane_and_centroid(ring_xyz)
                plane_1 = [A, B, C, D]
                # Obtaining XYZ coordinates of second ring 
                e, nxt_ring_xyz = get_xyz(f"base_{i+1}_ring_{k}.xyz")
                # Obtaining plane equation that describes the centroid of second ring 
                E, F, G, H, cent_2 = find_best_fitting_plane_and_centroid(nxt_ring_xyz)
                plane_2 = [E, F, G, H]
                # Calculating overlap between individual rings 
                if math.isclose(angle, 0.0, rel_tol=1e-5) == True or math.isclose(angle, 180.0, rel_tol=1e-5) == True:
                    ring_overlap = find_overlap(ring_xyz, nxt_ring_xyz, A, B, C, D, cent_1, f"base_{i}_ring_{j}_{k}")
                else:
                    ring_overlap = find_overlap(ring_xyz, nxt_ring_xyz, I, J, K, L, point, f"base_{i}_ring_{j}_{k}")

                # Calculating stack score
                stack_score = parameter_D*parameter_A*ring_overlap
                new_stack_score = new_stack_score + stack_score**2

        write_raport(f"base_{i}", f"base_{i+1}", new_stack_score, old_stack_score, xyz_file_name)

parser = argparse.ArgumentParser(description = "Claculate stacking score between adjacent bases in XYZ file.")
# parser.add_argument( "-n", "--name", type=str, metavar = "", required = True, help = "name of uppdated PDB file")
#parser.add_argument( "-p", "--pdb_filename", type=str, metavar = "" ,required = False, help = "name of PDB file")

#parser.add_argument( "-n", "--number_of_bases", type=int, metavar = "", required = False, help = "number of bases to be read from XYZ file")
parser.add_argument( "-b", "--base_specification_filename", type=str, metavar = "", required = False, help = "name of file containing the number of atoms for each ring in each base")
parser.add_argument( "-x", "--xyz_filename", type=str, metavar = "", required = False, help = "name of XYZ file")
args = parser.parse_args()

find_stacking_xyz(args.xyz_filename, args.base_specification_filename)


# In[ ]:




