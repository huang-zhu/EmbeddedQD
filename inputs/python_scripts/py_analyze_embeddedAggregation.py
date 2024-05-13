# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 15:08:16 2024

@author: huangzhu
"""

import os
import numpy as np
import mdtraj as md
import pickle
import scipy.linalg
from scipy.optimize import curve_fit
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from matplotlib import rc
import matplotlib as mpl
from matplotlib import mathtext
plt.rcParams.update({
    'text.usetex':False,
    'font.family':'Arial',
    'font.sans-serif':['Arial'],
    'mathtext.default':'regular',
    })

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

#######################################################################################################################################
### COMPUTE & PLOT 2D HISTOGRAM DATA AND AREA DENSITIES
#######################################################################################################################################
inputPath  = os.getcwd() + '/md'
outputPath = os.getcwd() + '/md'

# inputPath = '//wsl.localhost/Ubuntu-20.04/home/huangzhu/github/EmbeddedQD/QD_CSZS_5nm/Native/rep_0/md'
# outputPath = '//wsl.localhost/Ubuntu-20.04/home/huangzhu/github/EmbeddedQD/QD_CSZS_5nm/Native/rep_0/md'

pickling_histos          = True
saving_fig_histos        = True
pickling_areaDensities   = True
saving_fig_areaDensities = True

NP_name  = [i for i in inputPath.split('/') if 'QD_CSZS' in i][0]
membrane = [i for i in inputPath.split('/') if i in ['Native','Thick']][0]
rep      = [i for i in inputPath.split('/') if 'rep_' in i][0].split('_')[-1]

NP_prod_dict = { 
                'QD_CSZS_5nm'  : np.array([np.round(i, decimals=1) for i in np.arange(100.2,200.2,0.2)], dtype='float32'),
                'QD_CSZS_7nm'  : np.array([np.round(i, decimals=1) for i in np.arange(100.2,200.2,0.2)], dtype='float32'),
                'QD_CSZS_11nm' : np.array([np.round(i, decimals=1) for i in np.arange(300.2,400.2,0.2)], dtype='float32'),
                }

input_GRO     = inputPath + '/md_dry_centered.gro'
input_XTC     = inputPath + '/md_dry_centered.xtc'
traj          = md.load(input_XTC, top=input_GRO)
top           = traj.topology
md_sim_time   = traj.time
md_num_frames = traj.n_frames
md_size_box   = traj.unitcell_lengths

topologies_headgroups = {
                         'POPC'  : top.select('resname POPC and name PO4'),
                         'DOTAP' : top.select('resname DOTAP and name TAP'),
                         'CHOL'  : top.select('resname CHOL and name ROH'),
                         'PC24'  : top.select('resname DNPC and name PO4')
                         }

timesteps_wanted         = NP_prod_dict[NP_name]
timesteps_wanted_indices = np.argwhere( np.isin(md_sim_time/1000, timesteps_wanted) )[:,0]

com_QD = np.mean(traj.xyz[:,top.select('resname QD and name CSZS'),:], axis=1)[:,np.newaxis,:]

fig_areaDensities, ax_areaDensities = plt.subplots(figsize=(3,3))
fig_histos, ax_histos = plt.subplots(nrows=4, ncols=1, figsize = (3.25,8))
ax_histos_index = 0
for headgroupName,headgroupTop in topologies_headgroups.items(): 
    
    ### EXTRACT HEADGROUP POSITIONS
    if membrane == 'Native' and headgroupName == 'PC24':        
        positions_headgroups = traj.xyz[:,topologies_headgroups['POPC'],:]
    else:
        positions_headgroups = traj.xyz[:,headgroupTop,:]
            
    ### PREPARE 2D HISTOS 
    box_length  = [(NP_name in ['QD_CSZS_5nm'])*(30) + (NP_name in ['QD_CSZS_7nm'])*(31) + (NP_name in ['QD_CSZS_11nm'])*(32)][0] 
    box_height  = [(NP_name in ['QD_CSZS_5nm'])*(20) + (NP_name in ['QD_CSZS_7nm'])*(21) + (NP_name in ['QD_CSZS_11nm'])*(29)][0] 
    
    bin_width_x = 0.375 # nm
    bin_width_y = 0.375 # nm
    
    bin_count_x = int(np.ceil(box_length / bin_width_x)) # add an additional bin to account for min box size truncation
    bin_count_y = int(np.ceil(box_height / bin_width_y)) # add an additional bin to account for min box size truncation
    
    hist_xedges = np.arange(0, box_length + bin_width_x, bin_width_x) # add one bin to account for the edge
    hist_yedges = np.arange(0, box_height + bin_width_y, bin_width_y) # add one bin to account for the edge
    
    hist_xedges = np.tile( hist_xedges, (len(positions_headgroups[1]), 1) ) # matrix of the xedges to identify the bin to which each headgroup belongs
    hist_yedges = np.tile( hist_yedges, (len(positions_headgroups[1]), 1) ) # matrix of the yedges to identify the bin to which each headgroup belongs
    
    ### LOOP THROUGH FRAMES TO GENERATE HIST COUNTS
    hist_counts = np.zeros([bin_count_y, bin_count_x]) # generate zero matrix to add counts of each bin 
    for index_frame, frame in enumerate(timesteps_wanted_indices): 
        
        print("\n=== Analyzing %s %s rep_%s %s frame %d - timestep %.1f ns ===" % (NP_name, membrane, rep, headgroupName, frame, timesteps_wanted[index_frame]))

        ### 2D Histogram Operations
        distances_com_headgroups = getDistancesWithPBC_XY(com_QD[frame], positions_headgroups[frame], md_size_box[frame])[:,:,:2] # XY coordinates
        norms_com_headgroups     = np.linalg.norm(distances_com_headgroups, axis=2).T
        
        hist_rValues = norms_com_headgroups.copy() # radial values to bin
        hist_hValues = positions_headgroups[frame][:,2][:,np.newaxis] # height values to bin
        
        hist_binned_r = ( (hist_rValues < hist_xedges+bin_width_x).cumsum(axis=1).cumsum(axis=1) == 1) # mark F if bin radius < radial value, then compare the boolean array with itself and mark T the first True == True for each frame (row) 
        hist_binned_h = ( (hist_hValues < hist_yedges+bin_width_y).cumsum(axis=1).cumsum(axis=1) == 1 )# mark F if bin height < height value, then compare the boolean array with itself and mark T the first True == True for each frame (row) 
        hist_binned_r[:,-1]+=(np.sum(hist_binned_r,axis=1)==False) # every frame that has radii greater than the min box radius has all F, so convert the last radial value to True (this radial value is discarded, this is only done to match array sizes)
        hist_binned_h[:,-1]+=(np.sum(hist_binned_h,axis=1)==False) # every frame that has heights greater than the min box height has all F, so convert the last height value to True (this height value is discarded, this is only done to match array sizes)
    
        hist_binned_rh  = np.vstack([ np.where(hist_binned_r == True)[1] , np.where(hist_binned_h == True)[1] ]).T # stack arrays as [bin_radius, bin_height]
        
        hist_shellAreas = np.zeros([bin_count_x]) # generate zero matrix to populate with shell areas
        hist_shellRadii = hist_xedges[0].copy()  # radial values to calculate shell areas
        
        for shellRadius in range(len(hist_shellRadii)-1):
            hist_shellAreas[shellRadius] = np.pi*(hist_shellRadii[shellRadius + 1]**2 - hist_shellRadii[shellRadius]**2) 
            
        for rh in hist_binned_rh:
            if (rh[0]<(np.shape(hist_counts)[1])) and (rh[1]<(np.shape(hist_counts)[0])): 
                hist_counts[rh[1], rh[0]] += 1 / hist_shellAreas[rh[0]] / bin_width_y
    hist_counts /= len(timesteps_wanted_indices) 
    
    ### CONVERT HISTOS FROM COUNTS TO ENRICHMENT
    if membrane == 'Native' and headgroupName == 'PC24':  
        hist_counts[:,:] = 0
    else: 
        pickle_file         = os.getcwd() + '/../../../membrane_%s/rep_0/md/data_volumetricDensities.p' % (membrane)
        # pickle_file         = '//wsl.localhost/Ubuntu-20.04/home/huangzhu/github/EmbeddedQD/membrane_%s/rep_0/md/data_volumetricDensities.p' % (membrane)
        pickle_file_OPENED  = open(pickle_file, 'rb')
        bulk_densities_dict = pickle.load(pickle_file_OPENED)
        pickle_file_OPENED.close()
        hist_counts /= bulk_densities_dict[headgroupName] 
                
    ### PLOT HISTOS
    vmax, vmin = 2, 0
    cutoff = (NP_name=='QD_CSZS_5nm')*[27] + (NP_name=='QD_CSZS_7nm')*[28] + (NP_name=='QD_CSZS_11nm')*[40]
    im_straight = ax_histos[ax_histos_index].pcolormesh(hist_xedges[0], hist_yedges[0]-hist_yedges[0][cutoff[0]], hist_counts, cmap='binary', vmin=vmin, vmax=vmax)
    im_mirror   = ax_histos[ax_histos_index].pcolormesh(-hist_xedges[0], hist_yedges[0]-hist_yedges[0][cutoff[0]], hist_counts,cmap='binary', vmin=vmin, vmax=vmax)
    ax_histos[ax_histos_index].set_title(headgroupName, x=0.2, y=0.8)
    ax_histos[ax_histos_index].set(
                                   xlabel='$\mathit{D{_r}}$ (nm)', ylabel='$Z$ ($nm$)', 
                                   xlim=(-13,13), ylim=(-13,13),
                                   xticks=[-12,-6,0,6,12], yticks=[-12,-6,0,6,12]
                                   )
    fig_histos.colorbar(im_straight,ax=ax_histos[ax_histos_index], label='Fractional Enrichment' )
    fig_histos.tight_layout()
    ax_histos_index += 1
                    
    ### PLOT AREA DENSITIES
    if membrane == 'Native' and headgroupName == 'PC24':
        continue
    else:
        ### COMPUTE AREA DENSITY 
        areaDensity = np.array([ 
                                hist_shellRadii[1:],
                                np.sum(hist_counts * bulk_densities_dict[headgroupName] * bin_width_y, axis=0)
                                ]).T
        ax_areaDensities.plot(
                              areaDensity[:,0],
                              areaDensity[:,1],
                              label=headgroupName
                              )
        ax_areaDensities.set(
                             xlabel='$\mathit{D{_r}}$ (nm)', ylabel='Area Density ($lipids/nm^{2}$)', 
                             xlim=(0,13), ylim=(0,6),
                             xticks=np.arange(0,14), yticks=np.arange(0,7)
                             )
        ax_areaDensities.legend(frameon=False)
        fig_areaDensities.tight_layout()
       
        ### PICKLE HIST COUNTS FOR HEADGROUP -- REQUIRED FOR MEMBRANE FITTING AND P2 CALCULATIONS
        if pickling_histos:
            print('** Pickling data_histCounts for %s **' % (headgroupName))
            pickle_file        = outputPath + '/data_histCounts_%s.p' % (headgroupName)
            pickle_file_OPENED = open(pickle_file, 'wb')
            pickle.dump(hist_counts, pickle_file_OPENED)
            pickle_file_OPENED.close()

        ### PICKLE AREA DENSITY FOR HEADGROUP -- USEFUL FOR MEAN/ERROR CALCULATIONS/PLOTTING
        if pickling_areaDensities:
            print('** Pickling data_areaDensities for %s **' % (headgroupName))
            pickle_file        = outputPath + '/data_areaDensities_%s.p' % (headgroupName)
            pickle_file_OPENED = open(pickle_file, 'wb')
            pickle.dump(areaDensity, pickle_file_OPENED)
            pickle_file_OPENED.close()

### SAVE FIGURES
if saving_fig_histos: fig_histos.savefig(outputPath + '/plot_histos.png', dpi=2400)
print('\n=======================================')
print(' GENERATED HISTOGRAM ENRICHMENT PLOT   ')
print('=======================================')
if saving_fig_areaDensities: fig_areaDensities.savefig(outputPath + '/plot_areaDensities.png', dpi=2400)
print('\n=======================================')
print(' GENERATED AREA DENSITY PLOT           ')
print('=======================================')
  


    




