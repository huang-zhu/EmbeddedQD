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
### COMPUTE & PLOT EMBEDDED QD CONVERGENCE
#######################################################################################################################################
inputPath  = os.getcwd() + '/md'
outputPath = os.getcwd() + '/md'

# inputPath = '//wsl.localhost/Ubuntu-20.04/home/huangzhu/github/EmbeddedQD/QD_CSZS_5nm/Thick/rep_0/md'
# outputPath = '//wsl.localhost/Ubuntu-20.04/home/huangzhu/github/EmbeddedQD/QD_CSZS_5nm/Thick/rep_0/md'

saving_fig_convergence = True

NP_name  = [i for i in inputPath.split('/') if 'QD_CSZS' in i][0]
membrane = [i for i in inputPath.split('/') if i in ['Native','Thick']][0]
rep      = [i for i in inputPath.split('/') if 'rep_' in i][0].split('_')[-1]

input_GRO     = inputPath + '/md_dry_centered.gro'
input_XTC     = inputPath + '/md_dry_centered.xtc'
traj          = md.load(input_XTC, top=input_GRO)
top           = traj.topology
md_sim_time   = traj.time / 1000 # ns
md_num_frames = traj.n_frames
md_size_box   = traj.unitcell_lengths

fig_convergence, ax_convergence = plt.subplots(figsize=(3,3))

ax_convergence.plot( 
                    md_sim_time ,
                    md_size_box[:,0] * md_size_box[:,1] / (len([ res for res in top.residues if res.name not in ['QD','CHOL']]) / 2),
                    )

x_min, x_max, x_step = 0, md_sim_time[-1], 50
x_minor, x_major     = x_step/2, x_step
plot_xticks          = np.arange(x_min, x_max + x_step, x_step)
if NP_name == 'QD_CSZS_5nm':
    y_min, y_max, y_step = 0.64, 0.69, 0.01
    y_minor, y_major     = 0.01, y_step
    plot_yticks          = np.array([ y_min + y_step*step for step in range(int(round((y_max-y_min)/y_step + 1))) ])
elif NP_name == 'QD_CSZS_7nm':
    y_min, y_max, y_step = 0.63, 0.70, 0.02
    y_minor, y_major     = 0.01, y_step
    plot_yticks          = np.array([ y_min + y_step*step for step in range(int(round((y_max-y_min)/y_step + 1))) ]) 
elif NP_name == 'QD_CSZS_11nm':
    y_min, y_max, y_step = 0.58, 0.72, 0.04
    y_minor, y_major     = 0.02, y_step
    plot_yticks          = np.array([ y_min + y_step*step for step in range(int(round((y_max-y_min)/y_step + 1))) ]) 
ax_convergence.set(
                   xlabel='Time (ns)', ylabel='Area per lipid ($nm^{2}$)', 
                   xlim=(x_min,x_max), ylim=(y_min,y_max),
                   xticks=plot_xticks, yticks=plot_yticks
                   )
fig_convergence.tight_layout()

### SAVE FIGURES
if saving_fig_convergence: fig_convergence.savefig(outputPath + '/plot_convergence.png', dpi=2400)
print('\n=======================================')
print(' GENERATED EMBEDDED CONVERGENCE PLOT   ')
print('=======================================')



    




