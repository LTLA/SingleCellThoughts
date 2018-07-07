#################################################################
# This script submits all simulation jobs to the SLURM scheduler.
# The use of separate files makes it trivial to set up multiple jobs on different cores.
# This allows timely execution of all of the scripts, yet with easily altered parameters.
#################################################################

mkdir -p logs
R=/Users/lun01/software/R/devel/bin/R

# Submitting the gaussian cluster simulations with varying parameters.

for npops in 5 10 20
do
    for sd in 0.2 0.5 1
    do
sbatch << EOT
#!/bin/bash
#SBATCH -o logs/out-gauss-${npops}-${sd}
#SBATCH -e logs/err-gauss-${npops}-${sd}
#SBATCH -n 1
#SBATCH --mem 8000

echo 'NPOPS <- ${npops}; SD <- ${sd}; source("sim_gaussclust.R")' | ${R} --slave --vanilla
EOT
    done
done

# Submitting the trajectory simulations with varying parameters.

for npops in 5 10 20
do
    for sd in 0.2 0.5 1
    do
sbatch << EOT
#!/bin/bash
#SBATCH -o logs/out-traj-${npops}-${sd}
#SBATCH -e logs/err-traj-${npops}-${sd}
#SBATCH -n 1
#SBATCH --mem 8000

echo 'NPOPS <- ${npops}; SD <- ${sd}; source("sim_trajectory.R")' | ${R} --slave --vanilla
EOT
    done
done

# Submitting the gene-wise variance simulations with varying parameters.

for npops in 5 10 20
do
    for sd in 0.2 0.5 1
    do
sbatch << EOT
#!/bin/bash
#SBATCH -o logs/out-gene-${npops}-${sd}
#SBATCH -e logs/err-gene-${npops}-${sd}
#SBATCH -n 1
#SBATCH --mem 8000

echo 'NPOPS <- ${npops}; SD <- ${sd}; source("sim_genevar.R")' | ${R} --slave --vanilla
EOT
    done
done

