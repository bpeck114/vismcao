#!/bin/sh
# Job name
#SBATCH --account=ulens_g
#SBATCH --qos=regular
#SBATCH --constraint=gpu
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --job-name=dm2_study
#SBATCH --output=dm2_study.out
echo "---------------------------"
echo "Job id = $SLURM_JOBID"
echo "Proc id = $SLURM_PROCID"
hostname
date
echo "---------------------------"

cd /global/homes/b/bpeck/work/vismcao


srun -N 1 -n 1 -G 1 python /global/homes/b/bpeck/work/vismcao/vismcao_DM_script.py
 
date
echo "All done!"
exit $exitcode
