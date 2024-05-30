
<p align="center">
  <img src=https://github.com/huang-zhu/EmbeddedQD/assets/98200265/aadde84b-1242-4c6e-ae1c-249fa62749b9>
</p>

# Long-Chain Lipids Facilitate Insertion of Large Nanoparticles into Membranes of Small Unilamellar Vesicles

[<img src="https://img.shields.io/badge/DOI-10.1021%2Facs.langmuir.3c03471-blue">](https://doi.org/10.1021/acs.langmuir.3c03471)

This Github Repository was created to replicate the methods and results of the classical molecular dynamics simulations from this article. It includes all necessary files:
-	Coarse-grained models for nanoparticles (GRO and ITP)
-	Force field topology files
-	MDP files
-	Structure and converged trajectory files for one replica of the Native and Thick pure bilayers
-	Python scripts for system preparation
-	Python scripts for analyses and plotting
-	Bash scripts for running system preparation and simulations
-	Link to a Docker Container that contains the necessary software to run the simulations and perform the analyses
These are very large systems (400,000 – 700,000 coarse-grained beads) so the use of high-performance computing (HPC) resources is strongly suggested. 
The example shown below will follow the smallest system (5 nm QD) for the purpose of this tutorial. In this case, ``Singularity``/``Apptainer`` will be used due to software availability on the HPC resource used. {CREDIT TO TAMU FASTER/ACCESS}. The tutorial assumes you will be running this on an HPC resource with an NVIDIA GPU (for CUDA acceleration) using the ``slurm`` and that you have access to a ``/scratch`` directory.

## SET-UP
Pull the Docker container.

```
CONTAINER_IMAGE='huangzhu/github:EmbeddedQD_1.0'
cd ${SCRATCH}
singularity pull docker://${CONTAINER_IMAGE}
```

Clone the GitHub repository.
```
mkdir github
cd github
git clone https://github.com/huang-zhu/EmbeddedQD.git
```
```
### NOTE: IF YOUR HPC RESOURCE DOES NOT HAVE GIT, YOU CAN USE THE CONTAINER
CONTAINER="singularity exec ${SCRATCH}/github_EmbeddedQD_1.0.sif"
GIT="${CONTAINER} git"
${GIT} clone https://github.com/huang-zhu/EmbeddedQD.git
```

Inside the ``EmbeddedQD/bash_scripts`` directory you will find SH scripts that will (1) prepare all files, and run energy minimization and NPT equilibration, and (2) run the production run. The production run will also use ``gmx rdf`` starting from 100,000 (ps), keep in mind that this will change according to (a) your converged trajectory and (b) the QD diameter (i.e., with default settings, the production run for the 11 nm QD is twice as long as the smaller QDs). Plan accordingly or modify the scripts to account for this. You will also find the ``input_files`` directory, which contains all necessary files to run the simulations, analyze them, and plot results. Finally, you will find the ``membrane_Native`` and ``membrane_Thick`` directories that contain production trajectories for the pure bilayers of the lipid composition as described in the article and as followed in this tutorial. For the purpose of this tutorial, for this tutorial you will only run the 5 nm QD with default parameters.

Start by entering the account you will use to run the simulations and the email where you want email notifications to arrive. Replace ``ACCOUNT NUMBER`` and ``EMAIL`` below with your information (don't delete the single quotes) and run the lines. For this tutorial, I am setting a time limit of 5 days (120 hours), but you should change this according to your HPC resource's limits. 
```
DEF_ACCOUNT="Enter your account here"
DEF_EMAIL="Enter your email here"
cd EmbeddedQD/bash_scripts
sed -i s/DEF_ACCOUNT/"$DEF_ACCOUNT"/g bash_equil.sh 
sed -i s/DEF_ACCOUNT/"$DEF_ACCOUNT"/g bash_prod.sh 
sed -i s/DEF_EMAIL/"$DEF_EMAIL"/g bash_equil.sh 
sed -i s/DEF_EMAIL/"$DEF_EMAIL"/g bash_prod.sh 
sed -i s/DEF_TIME/'120:00:00'/g bash_equil.sh 
sed -i s/DEF_TIME/'120:00:00'/g bash_prod.sh 
```

## EQUILIBRATION
You are now ready to submit the equilibration script. This script will prepare the system requested, energy minimize it, and equilibrate. This script will require you to input the nanoparticle, the number of the replica, and the membrane type. Let's run both Native and Thick bilayers.
```
sbatch bash_equil.sh QD_CSZS_5nm 0 Native
sbatch bash_equil.sh QD_CSZS_5nm 0 Thick
```

## PRODUCTION
After the equilibration finishes, you will have to submit the production run. Below, replace ``BILAYER`` with the system you want to submit the production for.
```
REP_DIR=${SCRATCH}/github/EmbeddedQD/QD_CSZS_5nm/"BILAYER"/rep_0
cd ${REP_DIR}
sbatch bash_prod.sh 
```

## DATA ANALYSIS
You can now perform data analysis on the pure bilayers and the embedded quantum dots to reproduce the main results of the paper. You will first need to calculate the bulk density of each lipid in the pure bilayers. Navigate into the replica directory and run the python analysis script. This will generate pickle files containing the bulk densities which will be used for the following set of data analyses. In addition, convergence will be assessed by plotting the area per lipid (APL) as a function of simulation time. You need to do this for both Native and Thick pure bilayers. 

```
CONTAINER="singularity exec ${SCRATCH}/github_EmbeddedQD_1.0.sif"
PYTHON3="${CONTAINER} python3.10" 
cd ${SCRATCH}/github/EmbeddedQD/membrane_Native/rep_0
${PYTHON3} ../../inputs/python_scripts/py_analyze_pureBilayer.py
cd ../../membrane_Thick/rep_0
${PYTHON3} ../../inputs/python_scripts/py_analyze_pureBilayer.py
```
<p align="center">
  <img width="600" src=https://github.com/huang-zhu/EmbeddedQD/assets/98200265/ae942348-ffc2-4d12-8388-0d8fe07e99f1>
</p>

You are now ready to analyze the embedded quantum dots simulations. You will start by doing a convergence check by plotting the area per lipid (APL) as a function of time for both bilayers.
```
cd ${SCRATCH}/github/EmbeddedQD/QD_CSZS_5nm/Native/rep_0
${PYTHON3} ../../../inputs/python_scripts/py_analyze_embeddedConvergence.py
cd ${SCRATCH}/github/EmbeddedQD/QD_CSZS_5nm/Thick/rep_0
${PYTHON3} ../../../inputs/python_scripts/py_analyze_embeddedConvergence.py
```
<p align="center">
  <img width="600" src=https://github.com/huang-zhu/EmbeddedQD/assets/98200265/074d856f-a534-4181-819a-bac1e5d82abd>
</p>

You can see that both simulations converge at around 100 ns. You can then proceed with analyzing the lipid aggregation around the embedded quantum dots. 

These next scripts will assume that trajectory analysis will begin at 100 ns. If you simulate for longer, and/or change the analysis timeframe, you need to edit these scripts (see ``NP_prod_dict`` in each script). For this analysis, lipid aggregation will be quantified by plotting (A) histograms of lipid head group enrichment and (B) area density as a function of the lipid radial distance from the center of the embedded quantum dot. This data will also be stored as pickle files for subsequent analyses.
```
cd ${SCRATCH}/github/EmbeddedQD/QD_CSZS_5nm/Native/rep_0
${PYTHON3} ../../../inputs/python_scripts/py_analyze_embeddedAggregation.py
cd ../../Thick/rep_0
${PYTHON3} ../../../inputs/python_scripts/py_analyze_embeddedAggregation.py
```

<p align="center">
  <img width="600" src=https://github.com/huang-zhu/EmbeddedQD/assets/98200265/6b7ea5cc-c84e-4008-940c-4abbdde61c1a>
</p>

<p align="center">
  <img width="600" src=https://github.com/huang-zhu/EmbeddedQD/assets/98200265/ac9386dc-8911-46e1-8700-98942036be7c>
</p>

You can see that POPC enriched at the nanoparticle-membrane interface in the Native membrane. In the Thick membrane you can see that PC24 becomes enriched at this interface, displacing the POPC. Now that you have analyzed the aggregation in both cases, you can analyze lipid tail ordering along the membrane.

For this analysis, each leaflet will be fitted with non-linear regressions as described in the supplementary information of the article. A plot will be generated showing these regressions and the 2D vectors normal to the membrane fits. Afterwards, the lipid tail order $P{_2}$ will be quantified as a function of the lipid radial distance from the center of the embedded quantum dot using the local normal vectors computed from the membrane fits. Additionally, the segmental $P{_2}$  will also be plotted for both tails of lipids found at the nanoparticle-membrane interface and in the bulk (refer the main article for details). 

This data will not be stored as a pickle file, but you can easily implement pickle saves (take a look at the previous analysis scripts). Furthermore, a series of plots will be generated to show (A) 3D vectors normal to the membrane with the nanoparticle shell and (B) 3D molecular vectors of lipids with their lipid tails. These plots are generated only for some lipids at some radial distances, and only if lipids are present in the bin of the last frame (i.e., radial distance of 0.375 nm may or may not yield plots). The lipid indices for the plotted lipid tails are also saved on a text file for visualization. 
```
cd ${SCRATCH}/github/EmbeddedQD/QD_CSZS_5nm/Native/rep_0
${PYTHON3} ../../../inputs/python_scripts/py_analyze_embeddedOrdering.py
cd ../../Thick/rep_0
${PYTHON3} ../../../inputs/python_scripts/py_analyze_embeddedOrdering.py
```
<p align="center">
  <img width="600" src=https://github.com/huang-zhu/EmbeddedQD/assets/98200265/99562d47-67a4-47c3-85a5-f0533c8fba9e>
</p>

You can see that the membranes are well fit.

<p align="center">
  <img width="600" src=https://github.com/huang-zhu/EmbeddedQD/assets/98200265/7a060822-35e7-43a8-8d00-2c8d4ccd4dd0>
</p>

The $P{_2}$ order analyses, show that POPC possesses a higher order at $D{_r}$ corresponding to the nanoparticle-membrane interface. Furthermore, the segmental $P{_2}$ shows that the lipids at the interface exhibit a lower order in comparison to the lipids in the bulk membrane.

In contrast, the data for the Thick membrane shows that PC24 is more ordered at $D{_r}$ corresponding to the nanoparticle-membrane interface. While POPC shows a decreased order at the interface when compared to lipids in the bulk, PC24 lipids at the interface and the bulk show minimal difference between interfacial and bulk lipids. 

<p align="center">
  <img width="800" src=https://github.com/huang-zhu/EmbeddedQD/assets/98200265/965479e4-1127-4a92-971e-ce4015beabcf>
</p>

Plotting vectors normal to the membrane fits and molucular vectors of certain lipids at certain distances from the nanoparticle shows that fitting the membrane using 2D data allows for appropriate lipid tail order analyses.


