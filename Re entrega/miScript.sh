#!/bin/bash
#SBATCH --exclusive
#SBATCH --partition=Blade
#SBATCH -o directorioSalida/output_%j.txt
#SBATCH -e directorioSalida/errors_%j.txt

mpirun --bind-to none ./miAplicacion $1 $2
