#!/bin/bash
#SBATCH --exclusive
#SBATCH --partition=Blade
#SBATCH -o directorioSalida/output_%j.txt
#SBATCH -e directorioSalida/errors_%j.txt

mpiexec --bind-to none ./miAplicacion $1 $2

# Se ejecuta en el cluster con la linea de comando:
# sbatch -N 2 --tasks-per-node=1 ./miScript.sh 18 8