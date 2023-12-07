#!/bin/bash

#SBATCH --job-name=MRChem
#SBATCH --output=mrchem_mpi_vecmult.log
#SBATCH --error=mrchem_mpi_vecmult.err

#SBATCH --account=nn4654k
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=8
#SBATCH --time=00-00:01:00
###SBATCH --qos=devel

## README 
## run as :
## sbatch submitsript.sh <input> <path/to/mrchem/build> <tag>
## input is input file in keyword format without extenstion, the path 
## is self explanatory and the tag is whatever tag you want to append at 
## the end of the output files e.g. "mrchem_omp.out", where the tag was "omp".
##
## Additionaly, line 71 and 74 need to be uncommented wrt. if you are running
## mpi or omp. more info at the lines.
##

echo "Starting job"

echo "defining variables"
## make basic input names
input=$1
ext="inp"

echo input ${input}.${ext}
machine_name="betzy"

mrchem_path=$2
echo path ${mrchem_path}

type=$3
## gotta define this so it always starts from sharpVsmooth
##run_dir="/cluster/projects/nn4654k/gabrielgerez/calculations/PCM_MW_2020/"
## jq binary executable
##jq="/cluster/projects/nn4654k/gabrielgerez/jq"

echo "loading modules"
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

module --quiet purge  # Reset the modules to the system default
module load foss/2021b
module load CMake/3.21.1-GCCcore-11.2.0
# module load Eigen/3.4.0-GCCcore-11.2.0  # maybe not needed
module load Python/3.9.6-GCCcore-11.2.0


export UCX_LOG_LEVEL=ERROR
export OMP_NUM_THREADS=16
echo 'SLURM_CPUS_PER_TASK'  ${OMP_NUM_THREADS}

# Creating aliases:
workdir=${USERWORK}/${SLURM_JOB_ID}
echo $workdir
echo creating copy of file in ${workdir}
mkdir -p ${workdir}

echo ${SLURM_SUBMIT_DIR}/${input}.${ext} ${workdir}/mrchem_${type}.inp
cp ${SLURM_SUBMIT_DIR}/${input}.${ext} ${workdir}/mrchem_${type}.inp
cd ${workdir}


echo Running
## use this for mpi 
time ${mrchem_path}/bin/mrchem --launcher="mpirun --rank-by node --map-by numa --bind-to numa --oversubscribe" --executable=${mrchem_path}/bin/mrchem.x  mrchem_${type} #this should already print out the necesary stuff

## use this for omp
#time ${mrchem_path}/bin/mrchem  --executable=${mrchem_path}/bin/mrchem.x  mrchem_${type} #this should already print out the necesary stuff

echo calculation finished, copying back
## move stuff back to submitdir

echo ${workdir}/. ${SLURM_SUBMIT_DIR}
cp -r ${workdir}/. ${SLURM_SUBMIT_DIR}

## cleanup
echo "cleaning up"
if [ -d ${workdir} ]; then
    rm -rf ${workdir}
fi

echo "job done"
exit 0
