#!/bin/bash
#SBATCH --job-name="sparsematmult"
#SBATCH --output="sparsematmult.%j.%N.out"
#SBATCH --partition=cmp
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=56
#SBATCH --mem=32G
#SBATCH --export=ALL
#SBATCH -t 4:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<youremail>@scu.edu

module load GCC

export OMP_PLACES=cores
export OMP_PROC_BIND=close

DIR=${HOME}/csen-145-fa24/labs/lab04
LOGDIR="${DIR}/logs"

for f in 0.20 0.15 0.10 0.05; do
  for t in 28 24 20 16 14 12 8 4 2 1; do
    LF="sparsematmult-1000-${f}-${t}.log"
    if [ ! -f "${LOGDIR}/${LF}" ]; then
      echo ${LF%4}
      ${DIR}/sparsematmult 1000 1000 1000 ${f} -t ${t} > ${LOGDIR}/${LF}
    fi
  done
  for t in 28 24 20 16 14 12 8 4 2 1; do
    LF="sparsematmult-2000-${f}-${t}.log"
    if [ ! -f "${LOGDIR}/${LF}" ]; then
      echo ${LF%4}
      ${DIR}/sparsematmult 2000 530 5000 ${f} -t ${t} > ${LOGDIR}/${LF}
    fi
  done
  for t in 28 24 20 16 14 12 8 4 2 1; do
    LF="sparsematmult-900-${f}-${t}.log"
    if [ ! -f "${LOGDIR}/${LF}" ]; then
      echo ${LF%4}
      ${DIR}/sparsematmult 900 3500 575 ${f} -t ${t} > ${LOGDIR}/${LF}
    fi
  done
done


