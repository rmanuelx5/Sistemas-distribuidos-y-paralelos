#!/bin/bash
#SBATCH -N 1 
#SBATCH --exclusive
#SBATCH --partition=Blade
#SBATCH -o directorioSalida/output.txt 
#SBATCH -e directorioSalida/errors.txt 
#SBATCH --time=00:05:00  #Tiempo límite (HH:MM:SS)
./nreinas_pthreads $1 $2