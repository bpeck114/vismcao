#!/bin/sh
# Job name
#SBATCH --account=ulens_g
#SBATCH --qos=regular
#SBATCH --constraint=gpu
#SBATCH --nodes=1
#SBATCH --time=00:45:00
#SBATCH --job-name=grid_act_study
#SBATCH --output=line_act_study.out
echo "---------------------------"
echo "Job id = $SLURM_JOBID"
echo "Proc id = $SLURM_PROCID"
hostname
date
echo "---------------------------"

cd /global/homes/b/bpeck/work/vismcao


srun -N 1 -n 1 -G 1 python /global/homes/b/bpeck/work/vismcao/scripts/grid.py
 
date
echo "All done!"
exit $exitcode
