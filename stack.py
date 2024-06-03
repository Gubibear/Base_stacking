#!/usr/bin/env python
# coding: utf-8

# In[1]:


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

#Tool for reading aromatic rings specifications for xyz file:
def get_rings_from_list(file_name, n):
    l = []
    with open(file_name, "r") as z:
        for i, line in enumerate(z):
            if i == n:
                d = line.split()
                for i in range(2,len(d)):
                    l.append(int(d[i]))
            elif i > n:
                break
    return(l)

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

# Tool that allows to take previous and next elements in loop
def previous_and_next(some_iterable):
    prevs, items, nexts = tee(some_iterable, 3)
    prevs = chain([None], prevs)
    nexts = chain(islice(nexts, 1, None), [None])
    return zip(prevs, items, nexts)

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

# Tool for extracting atom coordinates of bases within the oligonucleotide 
def extract(structure_name):
    substring_1 = "'"
    substring_2 = "P"
    substring_3 = "H"
    coor = []
    E_names = []
    substring_A = "A"
    substring_T = "T"
    substring_G = "G"
    substring_C = "C"
    substring_U = "U"
    for model in structure_name:
        for chain in model:
            for residue in chain:
                if substring_A in residue.get_resname():
                    substring_4 = "N6"
                    for atom in residue:
                        if (substring_1 not in atom.get_fullname() and substring_2 not in atom.get_fullname() 
                            and substring_3 not in atom.get_fullname() and substring_4 not in atom.get_fullname()):
                            E_names.append(atom.element)
                            x,y,z = atom.get_coord()
                            coor.append(x)
                            coor.append(y)
                            coor.append(z)
                elif substring_T in residue.get_resname():
                    substring_4 = "C7"
                    substring_5 = "O"
                    for atom in residue:
                        if (substring_1 not in atom.get_fullname() and substring_2 not in atom.get_fullname() 
                            and substring_3 not in atom.get_fullname() and substring_4 not in atom.get_fullname() 
                            and substring_5 not in atom.get_fullname()):
                            E_names.append(atom.element)
                            x,y,z = atom.get_coord()
                            coor.append(x)
                            coor.append(y)
                            coor.append(z)
                elif substring_G in residue.get_resname():
                    substring_4 = "N2"
                    substring_5 = "O"
                    for atom in residue:
                        if (substring_1 not in atom.get_fullname() and substring_2 not in atom.get_fullname() 
                            and substring_3 not in atom.get_fullname() and substring_4 not in atom.get_fullname() 
                            and substring_5 not in atom.get_fullname()):
                            E_names.append(atom.element)
                            x,y,z = atom.get_coord()
                            coor.append(x)
                            coor.append(y)
                            coor.append(z)
                elif substring_C in residue.get_resname():
                    substring_4 = "N4"
                    substring_5 = "O"
                    for atom in residue:
                        if (substring_1 not in atom.get_fullname() and substring_2 not in atom.get_fullname() 
                            and substring_3 not in atom.get_fullname() and substring_4 not in atom.get_fullname() 
                            and substring_5 not in atom.get_fullname()):
                            E_names.append(atom.element)
                            x,y,z = atom.get_coord()
                            coor.append(x)
                            coor.append(y)
                            coor.append(z)
                elif substring_U in residue.get_resname():
                    substring_4 = "O"
                    for atom in residue:
                        if (substring_1 not in atom.get_fullname() and substring_2 not in atom.get_fullname() 
                            and substring_3 not in atom.get_fullname() and substring_4 not in atom.get_fullname()):
                            E_names.append(atom.element)
                            x,y,z = atom.get_coord()
                            coor.append(x)
                            coor.append(y)
                            coor.append(z)                            
    lst = np.array(coor)
    n = int((len(coor)/3))
    xyz = np.reshape(lst, (n,3))
    return E_names, xyz

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

# # Algorithm for obtaining equation of a plane that is parallel to the first one and is describing the centroid of second base 
# def find_second_plane(points, A, B, C, D):
    
#     # Calculate the coordinates of centroid of points from the second base 
#     n = len(points)
#     sum_points = np.sum(points, 0)
#     cent = sum_points/n
    
#     # Calculate distance between two planes   
#     dis = (np.absolute(A * cent[0] + B * cent[1] + C * cent[2] + D))/(np.sqrt(A**2 + B**2 + C**2))
#     new_D = -(A * cent[0] + B * cent[1] + C * cent[2])
    
#     # Return the coefficients of the plane equation
#     return A, B, C, new_D, cent, dis

# Tool for getting xyz coordinates of bases from PDB file
def get_base_xyz(pdb_filename):
    base = read_pdb(pdb_filename)
    E, XYZ = extract(base)
    return XYZ

def calculate_angle(A, B, C, A2, B2, C2):
    a = A * A2 + B * B2 + C * C2 
    b = np.sqrt(A**2 + B**2 + C**2) 
    c = np.sqrt(A2**2 + B2**2 + C2**2) 
    teta = np.rad2deg(np.arccos(a/(b*c)))
    return(teta)

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

def calculate_distance(A, B, C, D, point):
    dis = (np.absolute(A * point[0] + B * point[1] + C * point[2] + D))/(np.sqrt(A**2 + B**2 + C**2))
    #dis = np.absolute(A * point[0] + B * point[1] + C * point[2] + D)  
    return(dis)

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

# Tool for writting the final raport of the program
def write_raport(base_1, base_2, plane_1, plane_2, plane_3, dis, angle, overlap, stack, name):
    with open("raport.txt","a") as f:
        f.write(f"Results for structure {name}" + '\n')       
        f.write(f"Stacking score calculation between {base_1} and {base_2}" + '\n')
        f.write(f"The parameters of the plane equation for the first base are: {plane_1}" + '\n')
        f.write(f"The parameters of the plane equation for the second base are: {plane_2}" + '\n')
        f.write(f"The parameters of the average plane equation are: {plane_3}" + '\n' )
        f.write(f"The distance between the two aromatic rings is: {dis}" + '\n')
        f.write(f"The angle between the first and the second plane is: {angle}" + '\n')
        f.write(f"The percentage of the shared area between two two bases is: {overlap}" + '\n')
        f.write(f"The final stacking score between {base_1} and {base_2} is: {stack}" + '\n')
        f.write('\n')

# Algortithm for calculating stacking score for two adjacent bases in sequence 
def find_stacking_pdb(pdb_filename):
    
    # Read PDB file and extract residue names 
    Select = PDB.Select
    io = PDB.PDBIO()
    sequence = read_pdb(pdb_filename)
    res_list = get_restype(sequence)
    print(res_list)
    
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
        
    # Stacking score calculator 
    # Parsing two adjacent bases in sequence 
    for previous, item, nxt in previous_and_next(res_list):
        #print("Item is now", item, "next is", nxt, "previous is", previous)
        
        # Obtaining the XYZ coordinates of first base 
        inp = f"{item}.pdb"
        base_xyz = get_base_xyz(inp)
        #print(base_xyz)
        
        # Obtaining plane equation of the first base 
        A, B, C, D, cent_1 = find_best_fitting_plane_and_centroid(base_xyz)
        plane_1 = [A, B, C, D]
        
        # Parsing second base 
        if nxt != None:
            
            # Obtaining XYZ coordinates of second base 
            nxt_inp = f"{nxt}.pdb"
            nxt_base_xyz = get_base_xyz(nxt_inp)
            
            # Obtaining plane equation of the second base 
            E, F, G, H, cent_2 = find_best_fitting_plane_and_centroid(nxt_base_xyz)
            plane_2 = [E, F, G, H]
            
            # Calculating angle between planes 
            angle = calculate_angle(A, B, C, E, F, G)
            
            # Calculating overlap between bases 
            if math.isclose(angle, 0.0, rel_tol=1e-5) == True or math.isclose(angle, 180.0, rel_tol=1e-5) == True:
                dis = (np.sqrt((D-H)**2))/(np.sqrt(A**2 + B**2 + C**2))
                overlap = calculate_overlap(base_xyz, nxt_base_xyz, A, B, C, D, cent_1, item)
                plane_3 = [A, B, C, D]
            else:
                I, J, K, L, point = middle_plane(A, B, C, D, E, F, G, H, cent_1, cent_2)
                plane_3 = [I, J, K, L]
                # dis = (calculate_distance(E, F, G, H, cent_1) + calculate_distance(A, B, C, D, cent_2))/2
                dis = (calculate_distance(I, J, K, L, cent_1) + calculate_distance(I, J, K, L, cent_2))
                overlap = calculate_overlap(base_xyz, nxt_base_xyz, I, J, K, L, point, item)

            # Calculating stack score
            r_angle = np.deg2rad(angle)
            rot = np.exp(-1 * r_angle**4) + np.exp(-(r_angle - np.pi)**4) + 0.1 * np.exp(-(r_angle - 0.5*np.pi)**4)
            D = dis2D(dis)
            Stack_score = D*rot*overlap
            write_raport(item, nxt, plane_1, plane_2, plane_3, dis, angle, overlap, Stack_score, pdb_filename)
                
def find_stacking_xyz(file_name):
    elements, xyz_corr = get_xyz(file_name)
    number_of_bases = args.number_of_bases 
    res_list = []
    for i in range(1,number_of_bases + 1):
        res_list.append(f"base_{i}")
        n = i - 1
        bases_file_name = "bases.txt"
        atoms_in_base = get_rings_from_list(bases_file_name, n)
        el = []
        xyz = []
        for j in atoms_in_base:
            el.append(elements[j - 1])
            xyz.append(xyz_corr[j -1][0])
            xyz.append(xyz_corr[j -1][1])
            xyz.append(xyz_corr[j -1][2])
        xyz = np.reshape(xyz, (len(atoms_in_base),3))
        write_xyz(f"base_{i}", el, xyz)
    
    for previous, item, nxt in previous_and_next(res_list):
        #print("Item is now", item, "next is", nxt, "previous is", previous)
        
        # Obtaining the XYZ coordinates of first base 
        inp = f"{item}.xyz"
        e, base_xyz = get_xyz(inp)
        #print(base_xyz)
        
        # Obtaining plane equation for first base 
        A, B, C, D , cent_1 = find_best_fitting_plane_and_centroid(base_xyz)
        plane_1 = [A, B, C, D]
        
        # Parsing second base 
        if nxt != None:
            
            # Obtaining XYZ coordinates of second base 
            nxt_inp = f"{nxt}.xyz"
            e, nxt_base_xyz = get_xyz(nxt_inp)
            
            # Obtaining plane equation that describes the centroid of second base 
            E, F, G, H, cent_2 = find_best_fitting_plane_and_centroid(nxt_base_xyz)
            plane_2 = [E, F, G, H]
            
            # Calculating angle between planes 
            angle = calculate_angle(A, B, C, E, F, G)

            # Calculating overlap between bases 
            if math.isclose(angle, 0.0, rel_tol=1e-5) == True or math.isclose(angle, 180.0, rel_tol=1e-5) == True:
                dis = (np.absolute(D-H))/(np.sqrt(A**2 + B**2 + C**2))
                overlap = calculate_overlap(base_xyz, nxt_base_xyz, A, B, C, D, cent_1, item)
                plane_3 = [A, B, C, D]
            else:
                I, J, K, L, point = middle_plane(A, B, C, D, E, F, G, H, cent_1, cent_2)
                plane_3 = [I, J, K, L]
                dis = (calculate_distance(I, J, K, L, cent_1) + calculate_distance(I, J, K, L, cent_2))
                overlap = calculate_overlap(base_xyz, nxt_base_xyz, I, J, K, L, point, item)

            # Calculating stack score
            r_angle = np.deg2rad(angle)
            rot = np.exp(-1 * r_angle**4) + np.exp(-(r_angle - np.pi)**4) + 0.1 * np.exp(-(r_angle - 0.5*np.pi)**4)
            D = dis2D(dis)
            Stack_score = D*rot*overlap
            write_raport(item, nxt, plane_1, plane_2, plane_3, dis, angle, overlap, Stack_score, file_name)
    
parser = argparse.ArgumentParser(description = "Claculate stacking score between adjacent bases in XYZ or PDB file.")
# parser.add_argument( "-n", "--name", type=str, metavar = "", required = True, help = "name of uppdated PDB file")
parser.add_argument( "-p", "--pdb_filename", type=str, metavar = "" ,required = False, help = "name of PDB file")
parser.add_argument( "-x", "--xyz_filename", type=str, metavar = "", required = False, help = "name of XYZ file")
parser.add_argument( "-n", "--number_of_bases", type=int, metavar = "", required = False, help = "number of bases to be read from XYZ file")
args = parser.parse_args()

if args.pdb_filename != None:
    find_stacking_pdb(args.pdb_filename)
elif args.xyz_filename != None:
    find_stacking_xyz(args.xyz_filename)
else:
    print("The file name and type were not specified. Please, use -h option to get help")


# In[ ]:




