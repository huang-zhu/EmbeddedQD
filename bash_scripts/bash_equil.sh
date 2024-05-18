#!/bin/bash
#SBATCH --account=DEF_ACCOUNT
#SBATCH --time DEF_TIME
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --gres=gpu:1
#SBATCH --mem=16G
#SBATCH -J EmbQD_EQ
#SBATCH --mail-user=DEF_EMAIL
#SBATCH --mail-type=all

CONTAINER="singularity exec --nv ${SCRATCH}/github_EmbeddedQD_1.0.sif"
PYTHON2="${CONTAINER} python2.7"
PYTHON3="${CONTAINER} python3.10"
GMX="${CONTAINER} gmx"

SLURM_PATH=$SLURM_SUBMIT_DIR/slurm-$SLURM_JOB_ID.out

NP=$1 #QD_CSZS_5nm
repetition=$2 #0
membrane_type=$3 #Native

NP_resname='QD'
if [[ "$membrane_type" == *"Native"* ]]; then
    composition_DNPC=0 #"$1"
    composition_POPC=70 #"$2"
    composition_CHOL=20 #"$3"
    composition_DOTAP=10 #"$4"
    echo "**GENERATING NATIVE MEMBRANE WITH COMPOSITIONS:"
    echo "DNPC  = $composition_DNPC %"
    echo "POPC  = $composition_POPC %"
    echo "CHOL  = $composition_CHOL %"
    echo "DOTAP = $composition_DOTAP %"

elif [[ "$membrane_type" == *"Thick"* ]]; then
    composition_DNPC=25 #"$1"
    composition_POPC=45 #"$2"
    composition_CHOL=20 #"$3"
    composition_DOTAP=10 #"$4"
    echo "**GENERATING THICK MEMBRANE WITH COMPOSITIONS:"
    echo "DNPC  = $composition_DNPC %"
    echo "POPC  = $composition_POPC %"
    echo "CHOL  = $composition_CHOL %"
    echo "DOTAP = $composition_DOTAP %"
fi

lipid_resname="LIPIDS"
solvent='PW'
if [[ "$NP" == *"5nm"* ]]; then
    NP_radius=2.5 #nm
    box_length=36 #nm
    prod_time=200 #ns
    box_z=16
elif [[ "$NP" == *"7nm"* ]]; then
    NP_radius=3.5 #nm
    box_length=38 #nm
    prod_time=200 #ns
    box_z=16
elif [[ "$NP" == *"11nm"* ]]; then
    NP_radius=5.5 #nm
    box_length=42 #nm
    prod_time=400 #ns
    box_z=20
fi
echo "-- NP_radius = $NP_radius nm --"

BASE_PATH=${SCRATCH}/github/EmbeddedQD
mkdir -p ${BASE_PATH}/${NP}/$membrane_type

rep=${BASE_PATH}/${NP}/$membrane_type/rep_$repetition
mkdir ${rep}

inputs=${rep}/input_files
mkdir ${inputs}

cp -rv ${BASE_PATH}/inputs/nanoparticles/"$NP".* ${inputs}/ # Copy NP GRO and ITP
cp -rv ${BASE_PATH}/inputs/mdp/* ${inputs}/ # Copy MDP
cp -rv ${BASE_PATH}/inputs/python_scripts/* ${inputs}/ # Copy PY
cp -rv ${BASE_PATH}/inputs/ff/* ${inputs}/ # Copy necessary force field ITP
cp -rv ${BASE_PATH}/bash_scripts/*.sh ${rep}/ # Copy bash script

bilayer_x=$(awk "BEGIN { print 1*$box_length }")
bilayer_y=$(awk "BEGIN { print 1*$box_length }")
echo "*** BILAYER SIZE: $bilayer_x nm x $bilayer_y nm ***"
box_x=$bilayer_x
box_y=$bilayer_x
echo "*** BOX SIZE: $box_x nm x $box_y nm x $box_z nm ***"

${PYTHON2} ${inputs}/insane_NP.py -f ${inputs}/$NP.gro \
                                  -l DNPC:"$composition_DNPC" \
                                  -l POPC:"$composition_POPC" \
                                  -l CHOL:"$composition_CHOL" \
                                  -l DOTAP:"$composition_DOTAP" \
                                  -sol $solvent \
                                  -a 0.7 \
                                  -x $bilayer_x \
                                  -y $bilayer_y \
                                  -z $box_z \
                                  -exclNP $NP_radius \
                                  -o ${inputs}/init.gro

${GMX} editconf -f ${inputs}/init.gro \
                -o ${inputs}/init.gro \
                -translate 0 0 -6 \
                -box $box_x $box_y $box_z

${PYTHON3} ${inputs}/tool_cleanShell.py --input_GRO ${inputs}/init.gro \
                                        --output_GRO ${inputs}/init.gro

NUM_DNPC=$(grep -c 'DNPC   PO4' ${inputs}/init.gro)
NUM_POPC=$(grep -c 'POPC   PO4' ${inputs}/init.gro)
NUM_CHOL=$(grep -c 'CHOL   ROH' ${inputs}/init.gro)
NUM_DOTAP=$(grep -c 'DOTAP  TAP' ${inputs}/init.gro)
sed -i -e '/'"DNPC "'/{H;d}' -e '${x;s/^\n//p;x}' ${inputs}/init.gro
sed -i -e '/'"POPC "'/{H;d}' -e '${x;s/^\n//p;x}' ${inputs}/init.gro
sed -i -e '/'"CHOL "'/{H;d}' -e '${x;s/^\n//p;x}' ${inputs}/init.gro
sed -i -e '/'"DOTAP "'/{H;d}' -e '${x;s/^\n//p;x}' ${inputs}/init.gro
sed -i -e '/'"PW "'/{H;d}' -e '${x;s/^\n//p;x}' ${inputs}/init.gro

${GMX} genconf -f ${inputs}/init.gro -o ${inputs}/init.gro -renumber

NUM_SOLVENT=$(grep 'WP' ${inputs}/init.gro -c)
NUM_CATIONS=$(grep 'DOTAP  TAP' ${inputs}/init.gro -c)

if [[ "$membrane_type" == *"Native"* ]]; then
    > ${inputs}/topology.top
    tee -a ${inputs}/topology.top << EOF
#include "martini_v2.3refP_PEO.itp"
#include "lipids_martiniv2.0.itp"
#include "${NP}.itp"
#include "ions.itp"
[ system ]
; name
INSANE!

[ molecules ]
; name  number
$NP_resname 1
POPC $NUM_POPC
CHOL $NUM_CHOL
DOTAP $NUM_DOTAP
$solvent $NUM_SOLVENT
EOF
elif [[ "$membrane_type" == *"Thick"* ]]; then
    > ${inputs}/topology.top
    tee -a ${inputs}/topology.top << EOF
#include "martini_v2.3refP_PEO.itp"
#include "lipids_martiniv2.0.itp"
#include "${NP}.itp"
#include "ions.itp"
[ system ]
; name
INSANE!

[ molecules ]
; name  number
$NP_resname 1
DNPC $NUM_DNPC
POPC $NUM_POPC
CHOL $NUM_CHOL
DOTAP $NUM_DOTAP
$solvent $NUM_SOLVENT
EOF
fi

${GMX} make_ndx -f ${inputs}/init.gro -o ${inputs}/index.ndx << INPUT
q
INPUT
NUM_GROUPS=$(grep -oP '\[.*?\]' ${inputs}/index.ndx -c)
if [[ "$membrane_type" == *"Native"* ]]; then
    ${CONTAINER} \
    gmx make_ndx -f ${inputs}/init.gro -n ${inputs}/index.ndx -o ${inputs}/index.ndx << INPUT
r POPC CHOL DOTAP
name $NUM_GROUPS LIPIDS
q
INPUT
elif [[ "$membrane_type" == *"Thick"* ]]; then
    ${CONTAINER} \
    gmx make_ndx -f ${inputs}/init.gro -n ${inputs}/index.ndx -o ${inputs}/index.ndx << INPUT
r DNPC POPC CHOL DOTAP
name $NUM_GROUPS LIPIDS
q
INPUT
fi
NUM_GROUPS=$(grep -oP '\[.*?\]' ${inputs}/index.ndx -c)
echo "*** TOTAL NUMBER OF GROUPS BEFORE ADDING IONS: $NUM_GROUPS ***"

cd ${rep}

OMP_PER_MPI=$(($SLURM_NNODES*$SLURM_CPUS_ON_NODE)) # TOTAL NUMBER OF CPUs AVAILABLE
MPI=1
export OMP_NUM_THREADS=$OMP_PER_MPI
MDRUN_FLAGS="-ntmpi $MPI -ntomp $OMP_NUM_THREADS"
##############################################################
##### ENERGY MINIMIZATION
##############################################################
mkdir ${rep}/em

${GMX} grompp -f input_files/em.mdp \
              -c input_files/init.gro \
              -p input_files/topology.top \
              -n input_files/index.ndx \
              -po em/mdout_em.mdp \
              -o em/em.tpr \
              -maxwarn 1
${GMX} genion -s em/em.tpr \
              -p input_files/topology.top \
              -o em/em_init.gro \
              -nn $NUM_CATIONS << INPUTS
PW
INPUTS

${GMX} make_ndx -f em/em_init.gro \
                -o input_files/index.ndx << INPUTS
del $NUM_GROUPS-$(($NUM_GROUPS * 2))
q
INPUTS

if [[ "$membrane_type" == *"Native"* ]]; then
    ${GMX} make_ndx -f em/em_init.gro \
                    -n input_files/index.ndx \
                    -o input_files/index.ndx << INPUTS
r POPC CHOL DOTAP
name $NUM_GROUPS LIPIDS
r PW CL NA
name $(($NUM_GROUPS+1)) PW_ION
q
INPUTS
elif [[ "$membrane_type" == *"Thick"* ]]; then
    ${GMX} make_ndx -f em/em_init.gro \
                    -n input_files/index.ndx \
                    -o input_files/index.ndx << INPUTS
r DNPC POPC CHOL DOTAP
name $NUM_GROUPS LIPIDS
r PW CL NA
name $(($NUM_GROUPS+1)) PW_ION
q
INPUTS
fi

${GMX} grompp -f input_files/em.mdp \
              -c em/em_init.gro \
              -p input_files/topology.top \
              -n input_files/index.ndx \
              -po em/mdout_em.mdp \
              -o em/em.tpr -maxwarn 1
${GMX} mdrun -v -deffnm em/em $MDRUN_FLAGS
##############################################################

##############################################################
##### PREPARE bilayer.ndx FUTURE LEAFLET-WISE ANALYSES
##############################################################
cutoff=$(($box_z / 2))
if [[ "$membrane_type" == *"Native"* ]]; then
    ${GMX} select -f em/em.gro \
                  -s em/em.tpr \
                  -on input_files/bilayer.ndx << INPUTS
resname POPC and name PO4 and z > $cutoff
resname POPC and name PO4 and z < $cutoff
resname CHOL and name ROH and z > $cutoff
resname CHOL and name ROH and z < $cutoff
resname DOTAP and name TAP and z > $cutoff
resname DOTAP and name TAP and z < $cutoff
INPUTS
    ${GMX} make_ndx -n input_files/bilayer.ndx \
                    -o input_files/bilayer.ndx << INPUTS
name 0 POPC_PO4_UL
name 1 POPC_PO4_LL
name 2 CHOL_ROH_UL
name 3 CHOL_ROH_LL
name 4 DOTAP_TAP_UL
name 5 DOTAP_TAP_LL
q
INPUTS
elif [[ "$membrane_type" == *"Thick"* ]]; then
    ${GMX} select -f em/em.gro \
                  -s em/em.tpr \
                  -on input_files/bilayer.ndx << INPUTS
resname POPC and name PO4 and z > $cutoff
resname POPC and name PO4 and z < $cutoff
resname CHOL and name ROH and z > $cutoff
resname CHOL and name ROH and z < $cutoff
resname DOTAP and name TAP and z > $cutoff
resname DOTAP and name TAP and z < $cutoff
resname DNPC and name PO4 and z > $cutoff
resname DNPC and name PO4 and z < $cutoff
INPUTS
    ${GMX} make_ndx -n input_files/bilayer.ndx \
                    -o input_files/bilayer.ndx << INPUTS
name 0 POPC_PO4_UL
name 1 POPC_PO4_LL
name 2 CHOL_ROH_UL
name 3 CHOL_ROH_LL
name 4 DOTAP_TAP_UL
name 5 DOTAP_TAP_LL
name 6 DNPC_PO4_UL
name 7 DNPC_PO4_LL
q
INPUTS
fi
##############################################################

OMP_PER_MPI=$SLURM_NTASKS
MPI=1 #$(($SLURM_NTASKS/$OMP_PER_MPI))
export OMP_NUM_THREADS=$OMP_PER_MPI
MDRUN_FLAGS="-ntmpi $MPI -ntomp $OMP_NUM_THREADS -nb gpu -bonded gpu -update cpu -pme gpu -pmefft gpu -pin on -pinstride 1"
##############################################################
##### NPT
##############################################################
mkdir ${rep}/npt

### NPT WARMUP #1
${GMX} grompp -f input_files/npt_warmup_1.mdp \
              -c em/em.gro \
              -p input_files/topology.top \
              -n input_files/index.ndx \
              -po npt/mdout_npt_warmup_1.mdp \
              -o npt/npt_warmup_1.tpr \
              -maxwarn 1
${GMX} mdrun -v -deffnm npt/npt_warmup_1 $MDRUN_FLAGS
if [[ ! -f npt/npt_warmup_1.gro ]]; then
    echo "****************************"
    echo "*** NPT WARMUP #1 FAILED ***"
    echo "****************************"
fi

### NPT WARMUP #2
${GMX} grompp -f input_files/npt_warmup_2.mdp \
              -c npt/npt_warmup_1.gro \
              -p input_files/topology.top \
              -n input_files/index.ndx \
              -po npt/mdout_npt_warmup_2.mdp \
              -o npt/npt_warmup_2.tpr \
              -maxwarn 1
${GMX} mdrun -v -deffnm npt/npt_warmup_2 $MDRUN_FLAGS
if [[ ! -f npt/npt_warmup_2.gro ]]; then
    echo "****************************"
    echo "*** NPT WARMUP #2 FAILED ***"
    echo "****************************"
fi

### NPT
${GMX} grompp -f input_files/npt.mdp \
              -c npt/npt_warmup_2.gro \
              -p input_files/topology.top \
              -n input_files/index.ndx \
              -po npt/mdout_npt.mdp \
              -o npt/npt.tpr \
              -maxwarn 1
${GMX} mdrun -v -deffnm npt/npt $MDRUN_FLAGS
if [[ ! -f npt/npt.gro ]]; then
    echo "********************************"
    echo "*** NPT EQUILIBRATION FAILED ***"
    echo "********************************"
fi
##############################################################

mv ${SLURM_PATH} ${rep}
