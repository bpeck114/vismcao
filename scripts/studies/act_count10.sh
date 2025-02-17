#!/bin/sh
# Job name
#SBATCH --account=ulens_g
#SBATCH --qos=regular
#SBATCH --constraint=gpu
#SBATCH --nodes=1
#SBATCH --time=02:00:00
#SBATCH --job-name=dm_order_study
#SBATCH --output=dm_order.out
echo "---------------------------"
echo "Job id = $SLURM_JOBID"
echo "Proc id = $SLURM_PROCID"
hostname
date
echo "---------------------------"

cd /global/homes/b/bpeck/work/vismcao


srun -N 1 -n 1 -G 1 python /global/homes/b/bpeck/ulens/vismcao/scripts/studies/act_count10.py
 
date
echo "All done!"
exit $exitcode
