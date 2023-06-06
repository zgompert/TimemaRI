#!/bin/sh 
#SBATCH --time=144:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --account=gompert-np
#SBATCH --partition=gompert-np
#SBATCH --job-name=bcfCall
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu


module load samtools
module load bcftools
## bcftools 1.16

cd /scratch/general/nfs1/u6000989

perl VarCallFork.pl

