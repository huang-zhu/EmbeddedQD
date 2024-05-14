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

#######################################################################################################################################
### PURE MEMBRANE CALCULATIONS
#######################################################################################################################################
inputPath  = os.getcwd() + '/md'
outputPath = os.getcwd() + '/md'

# inputPath = '//wsl.localhost/Ubuntu-20.04/home/huangzhu/github/EmbeddedQD/membrane_Thick/rep_0/md'
# outputPath = '//wsl.localhost/Ubuntu-20.04/home/huangzhu/github/EmbeddedQD/membrane_Thick/rep_0/md'

pickling_volumetricDensities = True
saving_fig_convergence       = True

input_GRO     = inputPath + '/md.gro'
input_XTC     = inputPath + '/md.xtc'
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

membrane = 'Native'*(len(topologies_headgroups['PC24']) == 0) + 'Thick'*(len(topologies_headgroups['PC24']) != 0)

bin_height = 0.375
timesteps_wanted         = np.arange(51,101,1) # frame in ns
timesteps_wanted_indices = np.argwhere( np.isin(md_sim_time/1000, timesteps_wanted) )[:,0]

data_volumetricDensities = {}
for headgroupName in topologies_headgroups.keys():

    bulk_density_avg = np.zeros(len(timesteps_wanted_indices)) # average across number of sampled configs
    for index_frame, frame in enumerate(timesteps_wanted_indices): 
        
        print("Analyzing %s membrane - frame %d - timestep %d ns" % (membrane, frame, timesteps_wanted[index_frame]))

        bulk_numLipidsLeaflet = len(topologies_headgroups[headgroupName]) / 2 # divide by 2 assuming symmetry
        bulk_boxArea          = md_size_box[frame][0] * md_size_box[frame][1]
        bulk_density          = bulk_numLipidsLeaflet / bulk_boxArea / bin_height # per unit volume
        bulk_density_avg[index_frame] = bulk_density
    data_volumetricDensities[headgroupName] = np.mean(bulk_density_avg)
    print("*** Mean volumetric density in pure bilayer: %s = %.2f lipids/nm^3 " % (headgroupName, data_volumetricDensities[headgroupName]) )
    
    ### PICKLE VOL DENSITY FOR HEADGROUP -- REQUIRED FOR HISTOGRAMS (ENRICHMENT)
    if membrane == 'Native' and headgroupName == 'PC24':
        continue
    else:
        if pickling_volumetricDensities:
            print('** Pickling data_volumetricDensities **')
            pickle_file        = outputPath + '/data_volumetricDensities.p'
            pickle_file_OPENED = open(pickle_file, 'wb')
            pickle.dump(data_volumetricDensities, pickle_file_OPENED)
            pickle_file_OPENED.close()
#######################################################################################################################################

#######################################################################################################################################
### COMPUTE & PLOT PURE MEMBRANE CONVERGENCE
#######################################################################################################################################
nrows,ncols = 1,1
fig_width,fig_height = 4,4 # inches

fig_convergence, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(fig_width,fig_height))
ax.plot(
        md_sim_time/1000, 
        md_size_box[:,0] * md_size_box[:,1] / (len([res for res in top.residues if res.name not in ['PW','ION', 'CHOL']]) / 2)
        )
ax.set(xlabel='Time (ns)', ylabel='Area per lipid ($nm{^2}$)', ylim=(0.66,0.71))
fig_convergence.tight_layout()
if saving_fig_convergence: fig_convergence.savefig(outputPath + '/plot_convergence.png', dpi=1200)
print('\n=======================================')
print(' GENERATED BILAYER CONVERGENCE PLOT    ')
print('=======================================')

#######################################################################################################################################


