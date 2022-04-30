#!/bin/bash

#SBATCH --job-name="emolecules_tautomers_wbo_calcs"
#SBATCH --partition=ilg2.3
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20gb
# one hour
#SBATCH --time=48:00:00
#SBATCH --distribution=block:cyclic
#--------------

# Informational output
echo "=================================== SLURM JOB ==================================="
echo
echo "The job will be started on the following node(s):"
echo $SLURM_JOB_NODELIST
echo
echo "Slurm User:         $SLURM_JOB_USER"
echo "Run Directory:      $(pwd)"
echo "Job ID:             $SLURM_JOB_ID"
echo "Job Name:           $SLURM_JOB_NAME"
echo "Partition:          $SLURM_JOB_PARTITION"
echo "Number of nodes:    $SLURM_JOB_NUM_NODES"
echo "Number of tasks:    $SLURM_NTASKS"
echo "Submitted From:     $SLURM_SUBMIT_HOST"
echo "Submit directory:   $SLURM_SUBMIT_DIR"
echo "=================================== SLURM JOB ==================================="
echo


cd $SLURM_SUBMIT_DIR
echo 'Working Directory:'
pwd

date

. ~/.bashrc
conda activate danielenv
python emolecules_tautomer_wbo_calcs.py

date

