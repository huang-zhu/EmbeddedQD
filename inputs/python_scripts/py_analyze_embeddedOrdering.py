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

def getTopologiesFromNDX( input_NDX, groupString_ndx ):
    """ ====================================================================\n
        === Provide path to the NDX (input_NDX).                         ===\n
        === Provide string/name of the NDX group (groupString_ndx).      ===\n
        === Return an array of the topologies of the atoms in the group. ===\n
        ====================================================================\n """
    path_ndx_OPENED   = open(input_NDX, 'r')
    path_ndx_lines    = [ line[:-1] for line in path_ndx_OPENED.readlines() ]
    path_ndx_lines.append('[')
    for row, string in enumerate(path_ndx_lines): 
        if '[ '+groupString_ndx+' ]' in string:
            lines_start = row 
            lines_end   = row + 1
            break 
    while ('[' not in path_ndx_lines[lines_end]) or ('' == path_ndx_lines[lines_end]): 
        lines_end += 1
    
    topologies = []
    for string in path_ndx_lines[lines_start + 1 :lines_end]:
        for atom in string.split():
            topologies.append(int(atom)-1) #topology in mdtraj starts at 0
    print("\n** Topologies recovered for " + groupString_ndx + " **")
    print("--> Remember this list starts at 0 instead of 1 for mdtraj. ")
    return(topologies)
def computeDistancesPBC(ref, sel, box):
    """\n
    ==============
    === INPUTS ===
    ==============
    Reference positions (ref): 
        array(numFrames, numPositions, numCoordinates)
    Selection positions (sel): 
        array(numFrames, numPositions, numCoordinates)
    Box size (box): 
        array(numFrames, numCoordinates)
    Number of axes (axes, optional, default='XYZ'): 
        string(axes)
    ===============
    === OUTPUTS ===
    ===============
    Distance vectors (deltas):
        array(numFrames, numRefPositions, numSelPositions, numCoordinates)
    Box size array (boxArray):
        array(numFrames, numRefPositions, numSelPositions, numCoordinates) 
    \n"""
    # Take distance between particle i to every other particle j.
    deltas = ref[:,:,np.newaxis,:] - sel[:,np.newaxis,:,:]
    boxArray = np.repeat(box[:,np.newaxis,:], 
                         np.shape(deltas)[1], 
                         axis=1)
    boxArray = np.repeat(boxArray[:,:,np.newaxis,:], 
                         np.shape(deltas)[2], 
                         axis=2)
    
    deltas[:,:,:,0] -= (deltas[:,:,:,0] >= boxArray[:,:,:,0]/2)*boxArray[:,:,:,0]
    deltas[:,:,:,0] += (deltas[:,:,:,0] < -boxArray[:,:,:,0]/2)*boxArray[:,:,:,0]
    deltas[:,:,:,1] -= (deltas[:,:,:,1] >= boxArray[:,:,:,1]/2)*boxArray[:,:,:,1]
    deltas[:,:,:,1] += (deltas[:,:,:,1] < -boxArray[:,:,:,1]/2)*boxArray[:,:,:,1]
    deltas[:,:,:,2] -= (deltas[:,:,:,2] >= boxArray[:,:,:,2]/2)*boxArray[:,:,:,2]
    deltas[:,:,:,2] += (deltas[:,:,:,2] < -boxArray[:,:,:,2]/2)*boxArray[:,:,:,2]
    return(deltas,boxArray)
def gaussian(x, a, x0, sigma, c):
    return a * np.exp(-(x-x0)**2/(2*sigma)**2) + c
def LL4(x,b,c,d,e):
    return(c+(d-c)/(1+np.exp(b*(np.log(x)-np.log(e)))))
def stretchedExp(x,a,b,c,d):
    return a*np.exp(-(x/b)**c)+d
def gradient_stretchedExp(x,a,b,c,d):
    dfdx = a*c/x * (x/b)**c * np.exp(-(x/b)**c) 
    dfdy = np.ones(len(x))
    unit = np.linalg.norm(np.array([dfdx,dfdy]).T,axis=1)
    dfdx /= unit
    dfdy /= unit
    gradient = np.array([a*c/x * (x/b)**c * np.exp(-(x/b)**c) , np.ones(len(x))])
    gradient /= np.linalg.norm(gradient, axis=0)  
    return dfdx,dfdy,gradient
def set_axes_equal(ax: plt.Axes):
    """Set 3D plot axes to equal scale.

    Make axes of 3D plot have equal scale so that spheres appear as
    spheres and cubes as cubes.  Required since `ax.axis('equal')`
    and `ax.set_aspect('equal')` don't work on 3D.
    """
    limits = np.array([
        ax.get_xlim3d(),
        ax.get_ylim3d(),
        ax.get_zlim3d(),
    ])
    origin = np.mean(limits, axis=1)
    radius = 0.5 * np.max(np.abs(limits[:, 1] - limits[:, 0]))
    _set_axes_radius(ax, origin, radius)
def _set_axes_radius(ax, origin, radius):
    x, y, z = origin
    ax.set_xlim3d([x - radius, x + radius])
    ax.set_ylim3d([y - radius, y + radius])
    ax.set_zlim3d([z - radius, z + radius])

#######################################################################################################################################
### EXTRACT INITIAL DATA FOR ANALYSES
#######################################################################################################################################
inputPath  = os.getcwd() + '/md'
outputPath = os.getcwd() + '/md'

# inputPath = '//wsl.localhost/Ubuntu-20.04/home/huangzhu/github/EmbeddedQD/QD_CSZS_5nm/Native/rep_0/md'
# outputPath = '//wsl.localhost/Ubuntu-20.04/home/huangzhu/github/EmbeddedQD/QD_CSZS_5nm/Native/rep_0/md'

saving_fig_membraneFits          = True
saving_fig_vectorRepresentations = True
saving_fig_radialP2              = True
saving_fig_segmentalP2           = True

NP_name  = [i for i in inputPath.split('/') if 'QD_CSZS' in i][0]
membrane = [i for i in inputPath.split('/') if i in ['Native','Thick']][0]
rep      = [i for i in inputPath.split('/') if 'rep_' in i][0].split('_')[-1]

NP_prod_dict = { 
                'QD_CSZS_5nm'  : np.array([np.round(i, decimals=1) for i in np.arange(100.2,200.2,0.2)], dtype='float32'),
                'QD_CSZS_7nm'  : np.array([np.round(i, decimals=1) for i in np.arange(100.2,200.2,0.2)], dtype='float32'),
                'QD_CSZS_11nm' : np.array([np.round(i, decimals=1) for i in np.arange(300.2,400.2,0.2)], dtype='float32'),
                }

### EXTRACT BULK DENSITIES FROM PURE BILAYER
pickle_file         = os.getcwd() + '/../../../membrane_%s/rep_0/md/data_volumetricDensities.p' % (membrane)
# pickle_file         = '//wsl.localhost/Ubuntu-20.04/home/huangzhu/github/EmbeddedQD/membrane_%s/rep_0/md/data_volumetricDensities.p' % (membrane)
pickle_file_OPENED  = open(pickle_file, 'rb')
bulk_densities_dict = pickle.load(pickle_file_OPENED)
pickle_file_OPENED.close()

input_GRO     = inputPath + '/md_dry_centered.gro'
input_XTC     = inputPath + '/md_dry_centered.xtc'
input_NDX     = inputPath + '/../input_files/bilayer.ndx'

traj          = md.load(input_XTC, top=input_GRO)
top           = traj.topology
md_sim_time   = traj.time
md_num_frames = traj.n_frames
md_size_box   = traj.unitcell_lengths

timesteps_wanted         = NP_prod_dict[NP_name]
timesteps_wanted_indices = np.argwhere( np.isin(md_sim_time/1000, timesteps_wanted) )[:,0]

topologies_headgroups = {
                         'POPC'  : top.select('resname POPC and name PO4'),
                         'DOTAP' : top.select('resname DOTAP and name TAP'),
                         'CHOL'  : top.select('resname CHOL and name ROH'),
                         'PC24'  : top.select('resname DNPC and name PO4')
                         }

membraneFits_dict = {}
fittingData = []
for headgroupName,headgroupTop in topologies_headgroups.items(): 
    
    if (membrane == 'Native' and headgroupName in ['POPC']) or (membrane == 'Thick' and headgroupName in ['POPC','PC24']):
        ### EXTRACT HIST COUNTS FROM EMBEDDED QD
        pickle_file         = inputPath + '/data_histCounts_%s.p' % (headgroupName)
        pickle_file_OPENED  = open(pickle_file, 'rb')
        hist_counts         = pickle.load(pickle_file_OPENED) * bulk_densities_dict[headgroupName]
        pickle_file_OPENED.close()
        
        if len(fittingData) == 0:
            fittingData = hist_counts.copy()
        else: 
            fittingData += hist_counts

bin_width_x = 0.375
bin_width_y = 0.375

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
### MEMBRANE FITTING: UPPER LEAFLET CALCULATIONS
#######################################################################################################################################
fig_membraneFits, ax_membraneFits = plt.subplots(nrows=2, ncols=1, figsize=(3,5))
ax_membraneFits_index = 0
data = fittingData[:, np.sum(fittingData, axis=0)>0][int(len(fittingData)/2):,:]

x_bins = np.arange(0, bin_width_x*np.shape(data)[1], bin_width_x) + bin_width_x/2
y_bins = np.arange(0, bin_width_y*np.shape(data)[0], bin_width_y) + bin_width_y/2

### SELECT ONLY MAX DENSITIES FOR FIT
y0     = np.argwhere(data.T == np.repeat( np.max(data, axis=0)[np.newaxis,:], len(data), axis=0).T ) 
y_dump = []
y      = []
for arr in y0:
    if not np.isin(arr[0],y_dump):
        y_dump.append(arr[0])
        y.append(arr)
y     = np.array(y)
xdata = x_bins.copy()
ydata = y_bins[y[:,1]]

### SCATTER DATA
ax_membraneFits[ax_membraneFits_index].scatter(xdata, ydata, color='none', edgecolor='black', label='Data')

### FITTING TO GAUSSIAN
popt, pcov = curve_fit(gaussian, xdata, ydata)
residuals  = ydata- gaussian(xdata, *popt)
ss_res     = np.sum(residuals**2)
ss_tot     = np.sum((ydata-np.mean(ydata))**2)
r_squared  = 1 - (ss_res / ss_tot)
ax_membraneFits[ax_membraneFits_index].plot(xdata, gaussian(xdata, *popt), linewidth=1.5, color='green', label='(1) R$^2$ = %.4f' %(r_squared))

### FITTING TO LL4
popt, pcov = curve_fit(LL4, xdata, ydata)
residuals  = ydata- LL4(xdata, *popt)
ss_res     = np.sum(residuals**2)
ss_tot     = np.sum((ydata-np.mean(ydata))**2)
r_squared  = 1 - (ss_res / ss_tot)
ax_membraneFits[ax_membraneFits_index].plot(xdata, LL4(xdata, *popt), linewidth=1.5, color='blue', label='(1) R$^2$ = %.4f' %(r_squared))

### FITTING TO STRETCHED EXPONENTIAL
popt, pcov = curve_fit(stretchedExp, xdata, ydata, p0=[1,2,0,1], maxfev=10000) 
residuals  = ydata- stretchedExp(xdata, *popt)
ss_res     = np.sum(residuals**2)
ss_tot     = np.sum((ydata-np.mean(ydata))**2)
r_squared  = 1 - (ss_res / ss_tot)
ax_membraneFits[ax_membraneFits_index].plot(xdata, stretchedExp(xdata, *popt), linewidth=1.5, color='red', label='(1) R$^2$ = %.4f' %(r_squared))

### COMPUTE GRADIENT (VECTOR NORMAL TO SURFACE) FOR BEST FIT MODEL (STRETCHED EXPONENTIAL)
dfdx,dfdy,gradient = gradient_stretchedExp(xdata,*popt)
ax_membraneFits[ax_membraneFits_index].quiver(xdata, stretchedExp(xdata, *popt), dfdx, dfdy, width=0.005) 

### SAVE MEMBRANE FIT FOR UPPER LEAFLET - USED FOR P2
membraneFits_dict['UL'] = np.array([dfdx,dfdy]).T

### ADJUST FIGURE FOR VISUAL PURPOSES
ax_membraneFits[ax_membraneFits_index].set_xlabel('$\mathit{D}{_r}$ [$\mathit{nm}$]')
ax_membraneFits[ax_membraneFits_index].set_ylabel('$\mathit{Z}$ [$\mathit{nm}$]')
ax_membraneFits[ax_membraneFits_index].set_xlim(-0.5, 13)
ax_membraneFits[ax_membraneFits_index].set_ylim(-1, 15)
ax_membraneFits[ax_membraneFits_index].set_xticks(np.arange(0,14,2))
ax_membraneFits[ax_membraneFits_index].set_yticks(np.arange(0,16,2))
leg = ax_membraneFits[ax_membraneFits_index].legend(loc='upper right', frameon=False, labelspacing=0.3, handlelength=1.5)
for line in leg.get_lines():
    line.set_linewidth(2)            

#######################################################################################################################################
### MEMBRANE FITTING: LOWER LEAFLET CALCULATIONS
#######################################################################################################################################
ax_membraneFits_index += 1
data = fittingData[:, np.sum(fittingData, axis=0)>0][:int(len(fittingData)/2),:]

x_bins = np.arange(0, bin_width_x*np.shape(data)[1], bin_width_x) + bin_width_x/2
y_bins = np.arange(-bin_width_y*np.shape(data)[0], 0, bin_width_y) + bin_width_y/2

### SELECT ONLY MAX DENSITIES FOR FIT
y0     = np.argwhere(data.T == np.repeat( np.max(data, axis=0)[np.newaxis,:], len(data), axis=0).T ) 
y_dump = []
y      = []
for arr in y0:
    if not np.isin(arr[0],y_dump):
        y_dump.append(arr[0])
        y.append(arr)
y     = np.array(y)
xdata = x_bins.copy()
ydata = y_bins[y[:,1]]

### SCATTER DATA
ax_membraneFits[ax_membraneFits_index].scatter(xdata, ydata, facecolor='none', edgecolor='black', label='Data')

### FITTING TO GAUSSIAN
popt, pcov = curve_fit(gaussian, xdata, ydata)
residuals  = ydata- gaussian(xdata, *popt)
ss_res     = np.sum(residuals**2)
ss_tot     = np.sum((ydata-np.mean(ydata))**2)
r_squared  = 1 - (ss_res / ss_tot)
ax_membraneFits[ax_membraneFits_index].plot(xdata, gaussian(xdata, *popt), linewidth=1.5, color='green', label='(1) R$^2$ = %.4f' %(r_squared))

### FITTING TO LL4
popt, pcov = curve_fit(LL4, xdata, ydata)
residuals  = ydata- LL4(xdata, *popt)
ss_res     = np.sum(residuals**2)
ss_tot     = np.sum((ydata-np.mean(ydata))**2)
r_squared  = 1 - (ss_res / ss_tot)
ax_membraneFits[ax_membraneFits_index].plot(xdata, LL4(xdata, *popt), linewidth=1.5, color='blue', label='(1) R$^2$ = %.4f' %(r_squared))

### FITTING TO STRETCHED EXPONENTIONAL
popt, pcov = curve_fit(stretchedExp, xdata, ydata, p0=[-1,2,0,1], maxfev=10000) 
residuals  = ydata- stretchedExp(xdata, *popt)
ss_res     = np.sum(residuals**2)
ss_tot     = np.sum((ydata-np.mean(ydata))**2)
r_squared  = 1 - (ss_res / ss_tot)
plt.plot(xdata, stretchedExp(xdata, *popt), linewidth=1.5, color='red', label='(1) R$^2$ = %.4f' %(r_squared))

### COMPUTE GRADIENT (VECTOR NORMAL TO SURFACE) FOR BEST FIT MODEL (STRETCHED EXPONENTIAL)
dfdx,dfdy,gradient = gradient_stretchedExp(xdata,*popt)
ax_membraneFits[ax_membraneFits_index].quiver(xdata, stretchedExp(xdata, *popt), -dfdx, -dfdy, width=0.005)

### SAVE MEMBRANE FIT FOR UPPER LEAFLET - USED FOR P2
membraneFits_dict['LL'] = np.array([-dfdx,-dfdy]).T

### ADJUST FIGURE FOR VISUAL PURPOSES
ax_membraneFits[ax_membraneFits_index].set_xlabel('$\mathit{D}{_r}$ [$\mathit{nm}$]')
ax_membraneFits[ax_membraneFits_index].set_ylabel('$\mathit{Z}$ [$\mathit{nm}$]')
ax_membraneFits[ax_membraneFits_index].set_xlim(-0.5, 13)
ax_membraneFits[ax_membraneFits_index].set_ylim(-15, 1)
ax_membraneFits[ax_membraneFits_index].set_xticks(np.arange(0,14,2))
ax_membraneFits[ax_membraneFits_index].set_yticks(np.arange(-14,2,2))
leg = ax_membraneFits[ax_membraneFits_index].legend(loc='lower right', frameon=False, labelspacing=0.3, handlelength=1.5) 
for line in leg.get_lines():
    line.set_linewidth(2)   
         
fig_membraneFits.tight_layout()
fig_membraneFits.align_ylabels() 
if saving_fig_membraneFits: fig_membraneFits.savefig(outputPath + '/plot_membraneFits.png', dpi=2400)
plt.close(saving_fig_membraneFits)
print('\n=======================================')
print(' GENERATED MEMBRANE FITS PLOT          ')
print('=======================================')

#######################################################################################################################################
#######################################################################################################################################

#######################################################################################################################################
### COMPUTE ANGLES & GENERATE EXAMPLE VECTOR REPRESENTATIONS
### COMPUTE SEGMENTAL P2 & PLOT SEGMENTAL P2
#######################################################################################################################################
com_QD = np.mean(traj.xyz[:,top.select('resname QD and name CSZS'),:], axis=1)[:,np.newaxis,:]

dict_headgroups = { 
                    'Native': ['POPC_PO4'], 
                    'Thick' : ['POPC_PO4','DNPC_PO4']
                    }

for headgroup in dict_headgroups[membrane]:

    box_length   = 14 
    bin_width_x  = 0.375 
    bin_count_x  = int(np.ceil(box_length / bin_width_x))           

    ### PERFORM CALCULATIONS FOR UPPER LEAFLET
    leaflet                 = 'UL'
    indexGroup              = headgroup + '_' + leaflet
    resName                 = headgroup.split('_')[0]
    headgroupLabel          = 'POPC'*(resName == 'POPC') + 'PC24'*(resName == 'DNPC')
    topologies_indexGroup   = getTopologiesFromNDX(input_NDX, indexGroup) 
    coreTop_residues        = [ residue for residue in top.residues if ((residue.name == resName) and any(np.isin([atom.index for atom in residue.atoms],topologies_indexGroup))) ]
    coreTop_atoms           = [ [atom for atom in residue.atoms] for residue in coreTop_residues ]
    topologies_atomsIndices = [ [atom.index for atom in residue.atoms] for residue in coreTop_residues ]
    resStructure            = [atom.name for atom in coreTop_residues[0].atoms ]
    resTailStructure        = resStructure[np.where(headgroup.split('_')[-1] == np.array(resStructure))[0][0]+1:]
    resTailLength           = int(len(resTailStructure)/2)
    if len(resTailStructure) % 2 == 0: 
        resTailStructure.insert(resTailLength,resTailStructure.pop(1)) 
        resTailStructure=np.array(resTailStructure).reshape(2,int(len(resTailStructure)/2)).T
    else:
        assert not len(resTailStructure) % 2 != 0, "====== ERROR: TAILS ARE ASYMMETRIC | CODE UNTESTED ======"
            
    ### BIN LIPIDS RADIALLY
    hist_xedges = np.arange(0, box_length + bin_width_x, bin_width_x) 
    
    distances,_ = computeDistancesPBC(
                                       com_QD[timesteps_wanted_indices],
                                       traj.xyz[:,topologies_indexGroup,:][timesteps_wanted_indices],
                                       md_size_box[timesteps_wanted_indices]
                                       )
    norms       = np.linalg.norm(distances[:,:,:,:2], axis=3).reshape(np.shape(distances[:,0,:,0]))
    norms       = np.repeat(norms[:,:,np.newaxis],len(hist_xedges),axis=2)

    hist_xedges = np.tile( hist_xedges, (np.shape(norms)[1], 1) )
    hist_xedges = np.repeat(hist_xedges[np.newaxis,:,:],np.shape(norms)[0],axis=0)

    hist_binned_r          = ((norms < hist_xedges+bin_width_x).cumsum(axis=2).cumsum(axis=2) == 1)
    hist_binned_r[:,:,-1] += (np.sum(hist_binned_r,axis=2)==False) 
    binned_headgroup       = np.argwhere(hist_binned_r)[:,2].reshape(np.shape(hist_binned_r)[:2])

    shift = 1
    dummy_angles = [] 
    for binned_index, binned in enumerate(binned_headgroup):
        
        frame  = timesteps_wanted_indices[binned_index]
        binned = [ np.array(topologies_indexGroup)[np.where(binned_headgroup[binned_index] == b)] for b,r in enumerate(hist_xedges[0,0]) ][shift:]
        
        binned_angles = []
        for bin_index, bin_atoms in enumerate(binned[:-1]):
        
            bin_atoms = binned[bin_index]
            
            ### EXTRACT MEMBRANE FIT AND COMPUTE 3D NORMAL VECTOR IN CARTESIAN COORDINATES
            bin_headgroup_distances,_ = computeDistancesPBC(
                                                             com_QD[timesteps_wanted_indices][bin_index][np.newaxis,:,:], 
                                                             traj.xyz[:,bin_atoms,:][frame][np.newaxis,:,:], 
                                                             md_size_box[frame][np.newaxis,:]
                                                             )
            
            bin_headgroup_theta        = np.arctan2(bin_headgroup_distances[:,0,:,1],bin_headgroup_distances[:,0,:,0]) 
            bin_rzVector               = membraneFits_dict[leaflet][bin_index][np.newaxis,:][0]
            bin_headgroupNormalVectors = np.array([ [-1*bin_rzVector[0]*np.cos(theta),-1*bin_rzVector[0]*np.sin(theta),bin_rzVector[1]] for theta in bin_headgroup_theta[0] ])
            
            ### COMPUTE ANGLES IF LIPIDS PRESENT IN BIN
            if len(bin_atoms) == 0:
                binned_angles.append(np.array([]))
            else:
                
                ### SELECT TOPOLOGIES OF BINNED HEADGROUPS
                bin_sel        = np.argwhere(np.isin(np.array(topologies_atomsIndices),bin_atoms))
                bin_topologies = []
                for i,l in enumerate(topologies_atomsIndices):
                    if i in bin_sel[:,0]:
                        a = l.copy()
                        a.insert(2+resTailLength,a.pop(2+1))
                        bin_topologies.append(a)
                
                ### EXTRACT TOP & POS TO COMPUTE MOLECULAR VECTORS
                bin_topologies_tails    = np.array(bin_topologies)[:,bin_sel[0,1]+1:]
                bin_topologies_tails    = np.swapaxes(bin_topologies_tails.reshape(len(bin_topologies_tails),2,resTailLength),1,2)
                bin_positions_tails     = traj.xyz[:,bin_topologies_tails.flatten(order='K'),:][frame] 
                bin_moleculeVectors     = np.swapaxes(bin_positions_tails.reshape(len(bin_topologies_tails),2,len(bin_topologies_tails[0]),3),1,2) #[numLipids, numVectorsPerTail, numTails, numCoords]
                bin_moleculeVectors     = bin_moleculeVectors[:,:-1,:,:] - bin_moleculeVectors[:,1:,:,:]
                
                ### COMPUTE MOLECULAR UNIT VECTORS 
                bin_moleculeUnitVectors = np.zeros(np.shape(bin_moleculeVectors))
                for lipid in range(len(bin_moleculeUnitVectors)):
                    for tail in range(2):
                        for vector in range(len(bin_moleculeUnitVectors[lipid])):
                            bin_moleculeUnitVectors[lipid,vector,tail] = bin_moleculeVectors[lipid,vector,tail] / np.linalg.norm(bin_moleculeVectors[lipid,vector,tail])
                
                ### COMPUTE HEADGROUP NORMAL UNIT VECTORS
                bin_headgroupNormalUnitVectors = bin_headgroupNormalVectors / np.linalg.norm(bin_headgroupNormalVectors,axis=1)[:,np.newaxis]
                
                ### COMPUTE INNER PRODUCT BETWEEN BINNED LIPIDS AND LOCAL NORMAL VECTOR
                bin_innerProduct = np.zeros(np.shape(bin_moleculeUnitVectors)[:-1])
                for lipid in range(len(bin_innerProduct)):
                    for tail in range(2):
                        for vector in range(len(bin_innerProduct[lipid])):
                            bin_innerProduct[lipid,vector,tail] = np.inner( bin_moleculeUnitVectors[lipid,vector,tail],bin_headgroupNormalUnitVectors[lipid] )
                
                ### COMPUTE ANGLES FROM INNER PRODUCTS
                bin_angles     = np.arccos(bin_innerProduct)
                bin_angles_deg = np.rad2deg(bin_angles)   
                
                ### PLOT 3D VECTOR REPRESENTATIONS OF THIS BIN ALONG THE RADIAL AXIS FOR LAST FRAME
                if (frame == timesteps_wanted_indices[-1]) and (hist_xedges[0,0][shift:][bin_index] in [0.375, 3.75, 7.875, 11.25]):
                    
                    ### INITIALIZE 3D FIG AND ADJUST FOR VISUAL PURPOSES
                    fig = plt.figure() 
                    ax  = fig.add_subplot(111, projection='3d')
                    ax.set(xlabel=('X'), ylabel=('Y'), zlabel=('Z'))
                    ax.xaxis.pane.set_edgecolor('black')
                    ax.yaxis.pane.set_edgecolor('black')
                    ax.zaxis.pane.set_edgecolor('black')
                    ax.view_init(elev=15, azim=135)
                    ax.tick_params(axis='x', which='major', pad=-3)
                    ax.tick_params(axis='y', which='major', pad=-3)
                    ax.tick_params(axis='z', which='major', pad=-3)
                    ax.zaxis.labelpad=-8
                    
                    ### SCATTER SHELL BEADS
                    NP_pos = traj.xyz[:, top.select('resname QD and name CSZS'), :]
                    ax.scatter(NP_pos[frame,:,0], NP_pos[frame,:,1], NP_pos[frame,:,2],
                               linewidth=0, s=10, c='black')
                    
                    ### PLOT LOCAL NORMAL VECCTORS
                    origins_normalVectors = traj.xyz[:,bin_atoms,:][frame]
                    ax.scatter(origins_normalVectors[:,0], origins_normalVectors[:,1], origins_normalVectors[:,2], 
                               s=10, c='black', linewidth=1, marker='x') 
                    ax.quiver(origins_normalVectors[:,0], origins_normalVectors[:,1], origins_normalVectors[:,2], 
                              bin_headgroupNormalVectors[:,0],bin_headgroupNormalVectors[:,1],bin_headgroupNormalVectors[:,2],
                              linewidth=1, length=3, arrow_length_ratio=0.4, color='black')
                    set_axes_equal(ax)
                    
                    ### SAVE FIGURE
                    filename = outputPath + '/normalVectors_'+headgroupLabel+'_'+leaflet+'_'+str(hist_xedges[0,0][shift:][bin_index])+'nm.png'
                    if saving_fig_vectorRepresentations: fig.savefig(filename, dpi=2400)
                    plt.close(fig)
                    
                    ### INITIALIZE 3D FIG AND ADJUST FOR VISUAL PURPOSES
                    fig = plt.figure(figsize=(5, 5)) 
                    ax  = fig.add_subplot(111, projection='3d')
                    ax.set(xlabel=('X'), ylabel=('Y'), zlabel=('Z'))
                    ax.xaxis.pane.set_edgecolor('black')
                    ax.yaxis.pane.set_edgecolor('black')
                    ax.zaxis.pane.set_edgecolor('black')
                    ax.view_init(elev=20, azim=135)
                    ax.tick_params(axis='x', which='major', pad=-3)
                    ax.tick_params(axis='y', which='major', pad=-3)
                    ax.tick_params(axis='z', which='major', pad=-3)
                    ax.zaxis.labelpad=-8
                    
                    ### PLOT LIPID TAIL VECTORS AND HEADGROUP NORMAL VECTOR
                    origins_headgroups = traj.xyz[:,bin_atoms,:][frame]
                    origins_tails      = traj.xyz[:,bin_topologies_tails,:][frame]
                    output_list        = []
                    for lipid in range(len(bin_atoms[:3])):
                        ax.scatter(origins_headgroups[lipid,0], origins_headgroups[lipid,1], origins_headgroups[lipid,2],
                                   c='black', linewidth=1, marker='x', label='test')
                        ax.quiver(origins_headgroups[lipid,0], origins_headgroups[lipid,1], origins_headgroups[lipid,2],
                                  bin_headgroupNormalVectors[lipid,0],bin_headgroupNormalVectors[lipid,1],bin_headgroupNormalVectors[lipid,2],
                                  linewidth=0.75, length=0.8, arrow_length_ratio=0.4, color='black')
                        tail = 0
                        ax.scatter(origins_tails[lipid,:,tail,0],origins_tails[lipid,:,tail,1], origins_tails[lipid,:,tail,2],
                                   c='darkred', linewidth=2, marker='o')
                        ax.quiver(origins_tails[lipid,1:,tail,0],origins_tails[lipid,1:,tail,1], origins_tails[lipid,1:,tail,2], 
                                  bin_moleculeVectors[lipid,:,tail,0], bin_moleculeVectors[lipid,:,tail,1], bin_moleculeVectors[lipid,:,tail,2],
                                  linewidth=0.75, length=0.8, arrow_length_ratio=0.4, color='red')
                        tail = 1
                        ax.scatter(origins_tails[lipid,:,tail,0],origins_tails[lipid,:,tail,1], origins_tails[lipid,:,tail,2],
                                   c='darkblue', linewidth=2, marker='o')
                        ax.quiver(origins_tails[lipid,1:,tail,0],origins_tails[lipid,1:,tail,1], origins_tails[lipid,1:,tail,2], 
                                  bin_moleculeVectors[lipid,:,tail,0], bin_moleculeVectors[lipid,:,tail,1], bin_moleculeVectors[lipid,:,tail,2],
                                  linewidth=0.75, length=0.8, arrow_length_ratio=0.4, color='blue')
                        ### SAVE LIPID TOPOLOGY FOR USE IN VMD
                        output_list.append( str(np.sort(bin_topologies_tails[lipid].flatten())) )

                    ### SAVE FIGURE
                    filename = outputPath + '/molecularVectors_'+headgroupLabel+'_'+leaflet+'_'+str(hist_xedges[0,0][shift:][bin_index])+'nm.png'
                    if saving_fig_vectorRepresentations: 
                        
                        fig.savefig(filename, dpi=2400)
                        
                        ### WRITE output_list WITH INDEX NUMBER FOR PLOTTED LIPID TAILS
                        filename = outputPath + '/molecularVectors_'+headgroupLabel+'_'+leaflet+'_'+str(hist_xedges[0,0][shift:][bin_index])+'nm.txt'
                        with open(filename, "w") as output_file:
                            output_file.write('indices for molecularVectors_'+headgroupLabel+'_'+leaflet+'_'+str(hist_xedges[0,0][shift:][bin_index])+'nm\n')
                            for indices in output_list:
                                output_file.writelines(indices.replace('[', '').replace(']', '')+'\n')
                    plt.close(fig)            
                binned_angles.append(bin_angles)
        print("=== Analyzing %s frame %d: %.2f ps ===" % (headgroupLabel+'_'+leaflet,frame,timesteps_wanted[binned_index]))
        dummy_angles.append(binned_angles)

    ### PERFORM CALCULATIONS FOR LOWER LEAFLET
    leaflet                 = 'LL'    
    indexGroup              = headgroup + '_' + leaflet
    resName                 = headgroup.split('_')[0]
    headgroupLabel          = 'POPC'*(resName == 'POPC') + 'PC24'*(resName == 'DNPC')
    topologies_indexGroup   = getTopologiesFromNDX(input_NDX, indexGroup) 
    coreTop_residues        = [ residue for residue in top.residues if ((residue.name == resName) and any(np.isin([atom.index for atom in residue.atoms],topologies_indexGroup))) ]
    coreTop_atoms           = [ [atom for atom in residue.atoms] for residue in coreTop_residues ]
    topologies_atomsIndices = [ [atom.index for atom in residue.atoms] for residue in coreTop_residues ]
    resStructure            = [atom.name for atom in coreTop_residues[0]._atoms]
    resTailStructure        = resStructure[np.where(headgroup.split('_')[-1]==np.array(resStructure))[0][0]+1:]
    resTailLength           = int(len(resTailStructure)/2)
    if len(resTailStructure) % 2 == 0: 
        resTailStructure.insert(resTailLength,resTailStructure.pop(1)) 
        resTailStructure=np.array(resTailStructure).reshape(2,int(len(resTailStructure)/2)).T
    else:
        assert not len(resTailStructure) % 2 != 0, "====== ERROR: TAILS ARE ASYMMETRIC | CODE UNTESTED ======"

    ### BIN LIPIDS RADIALLY
    hist_xedges = np.arange(0, box_length + bin_width_x, bin_width_x) 
    
    distances,_ = computeDistancesPBC(
                                       com_QD[timesteps_wanted_indices],
                                       traj.xyz[:,topologies_indexGroup,:][timesteps_wanted_indices], 
                                       md_size_box[timesteps_wanted_indices]
                                       )
    norms       = np.linalg.norm(distances[:,:,:,:2], axis=3).reshape(np.shape(distances[:,0,:,0]))
    norms       = np.repeat(norms[:,:,np.newaxis],len(hist_xedges),axis=2)

    hist_xedges = np.tile(hist_xedges, (np.shape(norms)[1], 1)) 
    hist_xedges = np.repeat(hist_xedges[np.newaxis,:,:],np.shape(norms)[0],axis=0)

    hist_binned_r          = ((norms < hist_xedges+bin_width_x).cumsum(axis=2).cumsum(axis=2) == 1)
    hist_binned_r[:,:,-1] +=(np.sum(hist_binned_r,axis=2)==False) 
    binned_headgroup       = np.argwhere(hist_binned_r)[:,2].reshape(np.shape(hist_binned_r)[:2])

    for binned_index, binned in enumerate(binned_headgroup): 
        
        frame  = timesteps_wanted_indices[binned_index]
        binned = [ np.array(topologies_indexGroup)[np.where(binned_headgroup[binned_index] == b)] for b,r in enumerate(hist_xedges[0,0]) ][shift:]
        
        binned_angles = []
        for bin_index, bin_atoms in enumerate(binned[:-1]):
        
            bin_atoms = binned[bin_index]
            
            ### EXTRACT MEMBRANE FIT AND COMPUTE 3D NORMAL VECTOR IN CARTESIAN COORDINATES
            bin_headgroup_distances,_ = computeDistancesPBC(
                                                             com_QD[timesteps_wanted_indices][bin_index][np.newaxis,:,:], 
                                                             traj.xyz[:,bin_atoms,:][frame][np.newaxis,:,:], 
                                                             md_size_box[frame][np.newaxis,:]
                                                             )
            
            bin_headgroup_theta        = np.arctan2(bin_headgroup_distances[:,0,:,1],bin_headgroup_distances[:,0,:,0]) 
            bin_rzVector               = membraneFits_dict[leaflet][bin_index][np.newaxis,:][0]
            bin_headgroupNormalVectors = np.array([ [-1*bin_rzVector[0]*np.cos(theta),-1*bin_rzVector[0]*np.sin(theta),bin_rzVector[1]] for theta in bin_headgroup_theta[0] ])
            
            ### COMPUTE ANGLES IF LIPIDS PRESENT IN BIN
            if len(bin_atoms) == 0:
                binned_angles.append(np.array([]))
            else:
                
                ### SELECT TOPOLOGIES OF BINNED HEADGROUPS
                bin_sel        = np.argwhere(np.isin(np.array(topologies_atomsIndices),bin_atoms))
                bin_topologies = []
                for i,l in enumerate(topologies_atomsIndices):
                    if i in bin_sel[:,0]:
                        a = l.copy()
                        a.insert(2+resTailLength,a.pop(2+1))
                        bin_topologies.append(a)
                
                ### EXTRACT TOP & POS TO COMPUTE MOLECULAR VECTORS
                bin_topologies_tails    = np.array(bin_topologies)[:,bin_sel[0,1]+1:]
                bin_topologies_tails    = np.swapaxes(bin_topologies_tails.reshape(len(bin_topologies_tails),2,resTailLength),1,2)
                bin_positions_tails     = traj.xyz[:,bin_topologies_tails.flatten(order='K'),:][frame] 
                bin_moleculeVectors     = np.swapaxes(bin_positions_tails.reshape(len(bin_topologies_tails),2,len(bin_topologies_tails[0]),3),1,2) #[numLipids, numVectorsPerTail, numTails, numCoords]
                bin_moleculeVectors     = bin_moleculeVectors[:,:-1,:,:] - bin_moleculeVectors[:,1:,:,:]
                
                ### COMPUTE MOLECULAR UNIT VECTORS 
                bin_moleculeUnitVectors = np.zeros(np.shape(bin_moleculeVectors))
                for lipid in range(len(bin_moleculeUnitVectors)):
                    for tail in range(2):
                        for vector in range(len(bin_moleculeUnitVectors[lipid])):
                            bin_moleculeUnitVectors[lipid,vector,tail] = bin_moleculeVectors[lipid,vector,tail] / np.linalg.norm(bin_moleculeVectors[lipid,vector,tail])
                
                ### COMPUTE HEADGROUP NORMAL UNIT VECTORS
                bin_headgroupNormalUnitVectors = bin_headgroupNormalVectors / np.linalg.norm(bin_headgroupNormalVectors,axis=1)[:,np.newaxis]
                
                ### COMPUTE INNER PRODUCT BETWEEN BINNED LIPIDS AND LOCAL NORMAL VECTOR
                bin_innerProduct = np.zeros(np.shape(bin_moleculeUnitVectors)[:-1])
                for lipid in range(len(bin_innerProduct)):
                    for tail in range(2):
                        for vector in range(len(bin_innerProduct[lipid])):
                            bin_innerProduct[lipid,vector,tail] = np.inner( bin_moleculeUnitVectors[lipid,vector,tail],bin_headgroupNormalUnitVectors[lipid] )
                
                ### COMPUTE ANGLES FROM INNER PRODUCTS
                bin_angles     = np.arccos(bin_innerProduct)
                bin_angles_deg = np.rad2deg(bin_angles)   
                
                ### PLOT 3D VECTOR REPRESENTATIONS OF THIS BIN ALONG THE RADIAL AXIS FOR LAST FRAME
                if (frame == timesteps_wanted_indices[-1]) and (hist_xedges[0,0][shift:][bin_index] in [0.375, 3.75, 7.875, 11.25]):
                    
                    ### INITIALIZE 3D FIG AND ADJUST FOR VISUAL PURPOSES
                    fig = plt.figure() 
                    ax  = fig.add_subplot(111, projection='3d')
                    ax.set(xlabel=('X'), ylabel=('Y'), zlabel=('Z'))
                    ax.xaxis.pane.set_edgecolor('black')
                    ax.yaxis.pane.set_edgecolor('black')
                    ax.zaxis.pane.set_edgecolor('black')
                    ax.view_init(elev=15, azim=135)
                    ax.tick_params(axis='x', which='major', pad=-3)
                    ax.tick_params(axis='y', which='major', pad=-3)
                    ax.tick_params(axis='z', which='major', pad=-3)
                    ax.zaxis.labelpad=-8
                    
                    ### SCATTER SHELL BEADS
                    NP_pos = traj.xyz[:, top.select('resname QD and name CSZS'), :]
                    ax.scatter(NP_pos[frame,:,0], NP_pos[frame,:,1], NP_pos[frame,:,2],
                               linewidth=0, s=10, c='black')
                    
                    ### PLOT LOCAL NORMAL VECCTORS
                    origins_normalVectors = traj.xyz[:,bin_atoms,:][frame]
                    ax.scatter(origins_normalVectors[:,0], origins_normalVectors[:,1], origins_normalVectors[:,2], 
                               s=10, c='black', linewidth=1, marker='x') 
                    ax.quiver(origins_normalVectors[:,0], origins_normalVectors[:,1], origins_normalVectors[:,2], 
                              bin_headgroupNormalVectors[:,0],bin_headgroupNormalVectors[:,1],bin_headgroupNormalVectors[:,2],
                              linewidth=1, length=3, arrow_length_ratio=0.4, color='black')
                    set_axes_equal(ax)
                    
                    ### SAVE FIGURE
                    filename = outputPath + '/normalVectors_'+headgroupLabel+'_'+leaflet+'_'+str(hist_xedges[0,0][shift:][bin_index])+'nm.png'
                    if saving_fig_vectorRepresentations: fig.savefig(filename, dpi=2400)
                    plt.close(fig)
                    
                    ### INITIALIZE 3D FIG AND ADJUST FOR VISUAL PURPOSES
                    fig = plt.figure(figsize=(5, 5)) 
                    ax  = fig.add_subplot(111, projection='3d')
                    ax.set(xlabel=('X'), ylabel=('Y'), zlabel=('Z'))
                    ax.xaxis.pane.set_edgecolor('black')
                    ax.yaxis.pane.set_edgecolor('black')
                    ax.zaxis.pane.set_edgecolor('black')
                    ax.view_init(elev=20, azim=135)
                    ax.tick_params(axis='x', which='major', pad=-3)
                    ax.tick_params(axis='y', which='major', pad=-3)
                    ax.tick_params(axis='z', which='major', pad=-3)
                    ax.zaxis.labelpad=-8
                    
                    ### PLOT LIPID TAIL VECTORS AND HEADGROUP NORMAL VECTOR
                    origins_headgroups = traj.xyz[:,bin_atoms,:][frame]
                    origins_tails      = traj.xyz[:,bin_topologies_tails,:][frame]
                    output_list        = []
                    for lipid in range(len(bin_atoms[:3])):
                        ax.scatter(origins_headgroups[lipid,0], origins_headgroups[lipid,1], origins_headgroups[lipid,2],
                                   c='black', linewidth=1, marker='x', label='test')
                        ax.quiver(origins_headgroups[lipid,0], origins_headgroups[lipid,1], origins_headgroups[lipid,2],
                                  bin_headgroupNormalVectors[lipid,0],bin_headgroupNormalVectors[lipid,1],bin_headgroupNormalVectors[lipid,2],
                                  linewidth=0.75, length=0.8, arrow_length_ratio=0.4, color='black')
                        tail = 0
                        ax.scatter(origins_tails[lipid,:,tail,0],origins_tails[lipid,:,tail,1], origins_tails[lipid,:,tail,2],
                                   c='darkred', linewidth=2, marker='o')
                        ax.quiver(origins_tails[lipid,1:,tail,0],origins_tails[lipid,1:,tail,1], origins_tails[lipid,1:,tail,2], 
                                  bin_moleculeVectors[lipid,:,tail,0], bin_moleculeVectors[lipid,:,tail,1], bin_moleculeVectors[lipid,:,tail,2],
                                  linewidth=0.75, length=0.8, arrow_length_ratio=0.4, color='red')
                        tail = 1
                        ax.scatter(origins_tails[lipid,:,tail,0],origins_tails[lipid,:,tail,1], origins_tails[lipid,:,tail,2],
                                   c='darkblue', linewidth=2, marker='o')
                        ax.quiver(origins_tails[lipid,1:,tail,0],origins_tails[lipid,1:,tail,1], origins_tails[lipid,1:,tail,2], 
                                  bin_moleculeVectors[lipid,:,tail,0], bin_moleculeVectors[lipid,:,tail,1], bin_moleculeVectors[lipid,:,tail,2],
                                  linewidth=0.75, length=0.8, arrow_length_ratio=0.4, color='blue')
                        ### SAVE LIPID TOPOLOGY FOR USE IN VMD
                        output_list.append( str(np.sort(bin_topologies_tails[lipid].flatten())) )

                    ### SAVE FIGURE
                    filename = outputPath + '/molecularVectors_'+headgroupLabel+'_'+leaflet+'_'+str(hist_xedges[0,0][shift:][bin_index])+'nm.png'
                    if saving_fig_vectorRepresentations: 
                        
                        fig.savefig(filename, dpi=2400)
                        
                        ### WRITE output_list WITH INDEX NUMBER FOR PLOTTED LIPID TAILS
                        filename = outputPath + '/molecularVectors_'+headgroupLabel+'_'+leaflet+'_'+str(hist_xedges[0,0][shift:][bin_index])+'nm.txt'
                        with open(filename, "w") as output_file:
                            output_file.write('indices for molecularVectors_'+headgroupLabel+'_'+leaflet+'_'+str(hist_xedges[0,0][shift:][bin_index])+'nm\n')
                            for indices in output_list:
                                output_file.writelines(indices.replace('[', '').replace(']', '')+'\n')
                    plt.close(fig)            
                binned_angles.append(bin_angles)
        print("=== Analyzing %s frame %d: %.2f ps ===" % (headgroupLabel+'_'+leaflet,frame,timesteps_wanted[binned_index]))
        dummy_angles.append(binned_angles)

    dummy_angles = np.array(dummy_angles, dtype='object')
    radii  = hist_xedges[0,0][1:-1]
    angles = []

    for i in range(len(dummy_angles[0])):
        angles.append(np.concatenate([j.flatten(order='F') for j in dummy_angles[:,i] ]))

    #########################################
    ### PLOT RADIAL P2
    #########################################    
    fig_radialP2, ax_radialP2 = plt.subplots(figsize=(3,3))
    ax_radialP2.plot( 
                      radii,
                      np.array([np.mean(1.5*np.cos(ang)**2-0.5) for ang in angles]),
                      label="%s %s" %(membrane,headgroupLabel)
                      )
    ax_radialP2.set(
                     xlabel='$\mathit{D}{_r}$ [$\mathit{nm}$]', ylabel='$\mathit{P}{_2}$',
                     xlim=(0,13), ylim=(0.2,0.5),
                     xticks=(np.arange(0,13+1,1)), yticks=(np.arange(0.2,0.5+0.1,0.1))
                     )
    ax_radialP2.legend(frameon=False, loc='lower right')
    fig_radialP2.tight_layout()
    if saving_fig_radialP2: fig_radialP2.savefig(outputPath + '/plot_radialP2_%s_%s.png' % (membrane,headgroupLabel), dpi=2400)
    plt.close(fig_radialP2)
    print('\n=======================================')
    print(' GENERATED RADIAL P2 PLOT              ')
    print('=======================================')
    #########################################    
    
    #########################################
    ### ANGLE SELECTION FOR BULK LIPIDS
    #########################################
    regionCutoff_start = 11
    regionCutoff_end   = 13
    dummy_angles_bulk  = dummy_angles[:,(hist_xedges[0,0][shift:-1] >= regionCutoff_start) & (hist_xedges[0,0][shift:-1] <= regionCutoff_end)]
    angles_bulk        = []
    for i in range(len(dummy_angles_bulk[0])):
        angles_bulk.append(np.concatenate([j for j in dummy_angles_bulk[:,i] if len(j)>0]))
    angles_bulk = np.concatenate(angles_bulk) 
    #########################################

    #########################################
    ### ANGLE SELECTION FOR INTERFACE LIPIDS
    #########################################
    rdfData  = np.loadtxt(inputPath + '/rdf.xvg',comments=['@','#'])[:,:2][::-1]
    maxIndex = 0
    prevMax_y, newMax_y = rdfData[maxIndex,1], 0
    while maxIndex <= len(rdfData):
        prevMax_y, newMax_y = newMax_y, rdfData[maxIndex,1]
        if newMax_y < prevMax_y:
            rdfMax_x,rdfMax_y = rdfData[maxIndex-1,0], rdfData[maxIndex-1,1]
            break
        maxIndex += 1
        
    pickle_file            = os.getcwd() + '/../../../membrane_Thick/rep_0/md/data_volumetricDensities.p'
    # pickle_file            = '//wsl.localhost/Ubuntu-20.04/home/huangzhu/github/EmbeddedQD/membrane_Thick/rep_0/md/data_volumetricDensities.p'
    pickle_file_OPENED     = open(pickle_file, 'rb')
    bulk_density_headgroup = pickle.load(pickle_file_OPENED)['PC24']
    pickle_file_OPENED.close()
    
    pickle_file         = os.getcwd() + '/../../../%s/Thick/rep_%s/md/data_histCounts_PC24.p' % (NP_name,rep)
    # pickle_file         = '//wsl.localhost/Ubuntu-20.04/home/huangzhu/github/EmbeddedQD/%s/Thick/rep_%s/md/data_histCounts_PC24.p' % (NP_name,rep)
    pickle_file_OPENED  = open(pickle_file, 'rb')
    verticalLine_data   = pickle.load(pickle_file_OPENED) * bulk_density_headgroup * 0.375
    pickle_file_OPENED.close()

    verticalLine_data      = np.sum(verticalLine_data, axis=0) 
    regionCutoff_start     = radii[np.where(abs(radii-rdfMax_x)==np.min(abs(radii-rdfMax_x)))][0] #extract from C5 rdf
    regionCutoff_end       = radii[np.where(verticalLine_data==np.max(verticalLine_data))[0][0]]
    dummy_angles_interface = dummy_angles[:,(hist_xedges[0,0][shift:-1] >= regionCutoff_start) & (hist_xedges[0,0][shift:-1] <= regionCutoff_end)]
    
    angles_interface = []
    for i in range(len(dummy_angles_interface[0])):
        angles_interface.append(np.concatenate([j for j in dummy_angles_interface[:,i] if len(j)>0]))
    angles_interface = np.concatenate(angles_interface) 
    #########################################

    #########################################
    ### PLOT SEGMENTAL P2
    #########################################
    fig_segmentalP2, ax_segmentalP2= plt.subplots(figsize=(3,3))
    
    p2 = 1.5 * np.cos(angles_bulk)**2 - 0.5
    ax_segmentalP2.plot(np.arange(0,len(p2[0]))+1, np.mean(p2,axis=0)[:,0], linestyle='solid', color='black', marker='o', markerfacecolor='black', label='Bulk - sn1')
    ax_segmentalP2.plot(np.arange(0,len(p2[0]))+1, np.mean(p2,axis=0)[:,1], linestyle='solid', color='black', marker='o', markerfacecolor='white', label='Bulk - sn2')
    
    p2 = 1.5 * np.cos(angles_interface)**2 - 0.5
    ax_segmentalP2.plot(np.arange(0,len(p2[0]))+1, np.mean(p2,axis=0)[:,0], linestyle='dashed', color='black', marker='o', markerfacecolor='black', label='Interface - sn1')
    ax_segmentalP2.plot(np.arange(0,len(p2[0]))+1, np.mean(p2,axis=0)[:,1], linestyle='dashed', color='black', marker='o', markerfacecolor='white', label='Interface - sn2')
    
    x_min, x_max, x_step = 1, int(len(np.mean(p2,axis=0)[:,0])), 1
    x_minor, x_major     = 0, x_step        
    plot_xticks          = np.arange(x_min, x_max + x_step, x_step)
    plot_xlim            = (x_min-0.2, x_max+0.2)

    y_min, y_max, y_step = 0.0,0.6,0.1
    y_minor, y_major     = 0.05, y_step
    plot_yticks          = np.arange(y_min, y_max + y_step, y_step)
    plot_ylim            = (y_min-0.05, y_max+0.05)
    if headgroup == 'POPC_PO4': plot_ylim = (0.15, 0.55)
    ax_segmentalP2.set( 
                        xlabel='Tail Bead Number', ylabel='$\mathit{P}{_2}$', 
                        xlim=(plot_xlim), ylim=(plot_ylim),
                        xticks=(plot_xticks), yticks=(plot_yticks),
                        )
    ax_segmentalP2.legend(loc='lower left', frameon=False)
    fig_segmentalP2.tight_layout()
    if saving_fig_segmentalP2: fig_segmentalP2.savefig(outputPath + '/plot_segmentalP2_%s_%s.png' % (membrane,headgroupLabel), dpi=2400)
    plt.close(fig_segmentalP2)
    print('\n=======================================')
    print(' GENERATED SEGMENTAL P2 PLOT           ')
    print('=======================================')
    #########################################

    



































