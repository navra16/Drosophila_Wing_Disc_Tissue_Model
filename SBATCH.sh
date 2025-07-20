#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=3G
#SBATCH --output=YB_recalibed_L0d75A0d75B0d135_GrowthFreq25_restiff4d0Rmin_budNeck0d135# This affects the print out of the "std::cout" in the script, make sure this is changed for different jobs.
#SBATCH --mail-user=useremail@address.com
#SBATCH --mail-type=ALL
#SBATCH --job-name="budNeck"
#SBATCH -p gpu # This is the default partition, you can use any of the following; intel, batch, highmem, gpu

module load extra
module load GCC
module load cuda/9.1


srun -p gpu --gres=gpu:1 ./virus-model -dt=0.001 Data_structure.xml

