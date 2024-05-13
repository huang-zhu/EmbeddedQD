# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 09:06:20 2023

@author: huangzhu
"""

import argparse
import numpy as np
import mdtraj as md

def getDistancesWithPBC_XY(positions_i, positions_j, size_box):
    """ ====================================================\n
        === Provide coordinates of atom i (positions_i). ===\n
        === Provide coordinates of atom j (positions_j). ===\n
        === Provide size of box (size_box).              ===\n
        === Return distances array between i and j.      ===\n
        ====================================================\n """
    # Take distance between particle i to every other particle j.
    distances = positions_i[:,np.newaxis,:] - positions_j[np.newaxis,:,:]
    # Apple PBC in X and Y dimensions.
    distances[:,:,0] -= ( distances[:,:,0]  >=  size_box[0]/2 )*size_box[0]
    distances[:,:,0] += ( distances[:,:,0]  < -size_box[0]/2 )*size_box[0]
    distances[:,:,1] -= ( distances[:,:,1]  >=  size_box[1]/2 )*size_box[1]
    distances[:,:,1] += ( distances[:,:,1]  < -size_box[1]/2 )*size_box[1]
    return(distances)

class parserNP:
    pass
core = parserNP()
parser = argparse.ArgumentParser()
parser.add_argument('--input_GRO')
parser.add_argument('--output_GRO')
args = parser.parse_args()

input_GRO  = args.input_GRO
output_GRO = args.output_GRO

print('=== WILL NOW CLEAN THE CORE ===')

traj        = md.load(input_GRO)
top         = traj.topology
md_size_box = traj.unitcell_lengths

resnames    = np.unique([res.name for res in top._residues if res.name != 'QD'])
top_QD_CSZS = top.select('resname QD and name CSZS')
pos_QD_CSZS = traj.xyz[:, top_QD_CSZS, :]

com_QD_CSZS   = np.mean(pos_QD_CSZS, axis=1)[:,np.newaxis,:]
distances     = getDistancesWithPBC_XY(com_QD_CSZS[0], pos_QD_CSZS[0], md_size_box[0])
norms         = np.linalg.norm(distances, axis=2).T
norms_rounded = np.array([np.round(i, decimals=2) for i in norms]) 
NP_radius     = np.round( np.mean(norms_rounded), decimals=2)
solvents      = ['PW']

atom_data = np.array( [ [str(atom.residue.resSeq)+atom.residue.name, atom.name, atom.index, atom.index+1 ] for atom in top._atoms] )
gro_data  = np.array( [ [str(atom.residue.resSeq), atom.residue.name, atom.name, atom.index+1, traj.xyz[:,atom.index,:][0][0], traj.xyz[:,atom.index,:][0][1], traj.xyz[:,atom.index,:][0][2] ] for atom in top._atoms] )

resNumNames_to_remove = []
resNumNames_to_keep   = np.full(np.shape(atom_data[:,0]), True)
total_atoms_to_remove = 0
for res in resnames: 

    print('Screening residue: ' + res)
    res_top      = top.select('resname ' + res)
    res_pos      = traj.xyz[:, res_top, :]
    radius       = [ NP_radius*(res not in solvents) + NP_radius*(res in solvents)*1.25 ][0]
    dist         = getDistancesWithPBC_XY(com_QD_CSZS[0], res_pos[0], md_size_box[0])
    norms        = np.array([ np.round(norm, decimals=2) for norm in np.linalg.norm(dist, axis=2).T])
    boolean      = norms <= radius
    indices      = np.argwhere(norms <= radius)[:, 0]

    marked_atoms_raw = np.isin( atom_data[:,2].astype(int), res_top[indices] )
    marked_residues  =  np.unique( atom_data[:,0][marked_atoms_raw])
    marked_atoms     =  np.isin( atom_data[:,0], marked_residues )
    resNumNames_to_keep[(np.argwhere(marked_atoms))] = res not in solvents
    resNumNames      = [ str(r) for r in marked_residues ]
    total_atoms_to_remove += sum(marked_atoms)*(res in solvents)

    if resNumNames != []:
        print('  --> Found ' + str(len(resNumNames)) + ' ' + res + ' within ' + str(radius) + ' nm the core COM. Processing...\n')
        resNumNames_to_remove += resNumNames

input_GRO_OPENED = open(input_GRO, 'r')
lines            = input_GRO_OPENED.readlines() 
lines            = [x[:-1] for x in lines] 
title            = lines[0]
lines            = [' '.join(line.split()).split(" ") for line in lines] 
total_old_atoms  = int(lines[1][0])
total_new_atoms  = total_old_atoms-total_atoms_to_remove
box_size         = lines[-1]
input_GRO_OPENED.close()

output_GRO_OPENED = open(output_GRO, 'w')
output_GRO_OPENED.write(title)
output_GRO_OPENED.write( '\n{:}\n'.format( str(total_new_atoms) ))
for lineIndex, line in enumerate( gro_data[resNumNames_to_keep] ): 

    if (int(line[3])/100000 >= 1): 
        line[3] = str(int(line[3]) - 100000*int(np.floor((int(line[3])/100000))))
    if (line[0]+line[1] in resNumNames_to_remove) and (line[1] not in solvents): 
        line[4] = float(line[4]) + [md_size_box[0][0]/2+1.5][0]*(float(line[4])>=md_size_box[0][0]/2) - [md_size_box[0][0]/2+1.5][0]*(float(line[4])<md_size_box[0][0]/2)

    output_GRO_OPENED.write( '{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}\n'.format( int(line[0]), str(line[1]), str(line[2]), int(line[3]), float(line[4]), float(line[5]), float(line[6]) ) )
output_GRO_OPENED.write('{:10.5f}{:10.5f}{:10.5f}\n'.format( float(box_size[0]),float(box_size[1]),float(box_size[2])))
output_GRO_OPENED.close()
print('=== GRO UPDATED ===')
print('=== THE CORE OF YOUR PARTICLE SHOULD NOW BE EMPTY ===')
print('=== USE GMX GENCONF WITH -renumber FLAG ===')
print('=== WATER WAS REMOVED, ALL ELSE WAS DISPLACED TO THE OUTSIDE OF THE BOX ===\n')




