#!/bin/bash
#SBATCH --account=DEF_ACCOUNT
#SBATCH --time DEF_TIME
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --gres=gpu:1
#SBATCH --mem=16G
#SBATCH -J EmbQD_PROD
#SBATCH --mail-user=DEF_EMAIL
#SBATCH --mail-type=all

CONTAINER="singularity exec --nv ${SCRATCH}/github_EmbeddedQD_1.0.sif"
GMX="${CONTAINER} gmx"

SLURM_PATH=$SLURM_SUBMIT_DIR/slurm-$SLURM_JOB_ID.out

start_analysis=0 # ps

OMP_PER_MPI=$SLURM_NTASKS
MPI=1 #$(($SLURM_NTASKS/$OMP_PER_MPI))
export OMP_NUM_THREADS=$OMP_PER_MPI
MDRUN_FLAGS="-ntmpi $MPI -ntomp $OMP_NUM_THREADS -nb gpu -bonded gpu -update cpu -pme gpu -pmefft gpu -pin on -pinstride 1"
##############################################################
##### PROD
##############################################################
mkdir md

${GMX} grompp -f input_files/md.mdp \
              -c npt/npt.gro \
              -p input_files/topology.top \
              -n input_files/index.ndx \
              -po md/mdout_md.mdp \
              -o md/md.tpr \
              -maxwarn 1
${GMX} mdrun -v -deffnm md/md $MDRUN_FLAGS

if [[ ! -f md/md.gro ]]; then
    echo "****************************"
    echo "*** MD PRODUCTION FAILED ***"
    echo "****************************"
else
    ##### PROCESS THE TRAJ
    if [[ "$(pwd)" == *"Native"* ]]; then
        ${GMX} make_ndx -f md/md.gro \
                        -o md/md_dry.ndx << INPUTS
r QD POPC CHOL DOTAP
r QD & a CSZS
q
INPUTS
        ${GMX} trjconv -f md/md.gro \
                       -s md/md.tpr \
                       -n md/md_dry.ndx \
                       -o md/md_dry_centered.gro \
                       -pbc mol \
                       -conect \
                       -center << INPUTS
QD
QD_POPC_CHOL_DOTAP
INPUTS
        ${GMX} trjconv -f md/md.xtc \
                       -s md/md.tpr \
                       -n md/md_dry.ndx \
                       -o md/md_dry_centered.xtc \
                       -pbc mol \
                       -conect \
                       -center << INPUTS
QD
QD_POPC_CHOL_DOTAP
INPUTS
        ${GMX} rdf -f md/md_dry_centered.xtc \
                   -s md/md.tpr \
                   -b $start_analysis \
                   -bin 0.02 \
                   -rmax 14 \
                   -norm number_density \
                   -o md/rdf.xvg << INPUTS
com of name CSZS
resname QD and name S1 C2 C3 C4 C5
resname POPC
resname DOTAP
resname CHOL
INPUTS
    elif [[ "$(pwd)" == *"Thick"* ]]; then
        ${GMX} make_ndx -f md/md.gro \
                        -o md/md_dry.ndx << INPUTS
r QD POPC CHOL DOTAP DNPC
r QD & a CSZS
q
INPUTS
        ${GMX} trjconv -f md/md.gro \
                       -s md/md.tpr \
                       -n md/md_dry.ndx \
                       -o md/md_dry_centered.gro \
                       -pbc mol \
                       -conect \
                       -center << INPUTS
QD
QD_POPC_CHOL_DOTAP_DNPC
INPUTS
        ${GMX} trjconv -f md/md.xtc \
                       -s md/md.tpr \
                       -n md/md_dry.ndx \
                       -o md/md_dry_centered.xtc \
                       -pbc mol \
                       -conect \
                       -center << INPUTS
QD
QD_POPC_CHOL_DOTAP_DNPC
INPUTS
        ${GMX} rdf -f md/md_dry_centered.xtc \
                   -s md/md.tpr \
                   -b $start_analysis \
                   -bin 0.02 \
                   -rmax 14 \
                   -norm number_density \
                   -o md/rdf.xvg << INPUTS
com of name CSZS
resname QD and name S1 C2 C3 C4 C5
resname POPC
resname DOTAP
resname CHOL
resname DNPC
INPUTS
    fi
fi
##############################################################

