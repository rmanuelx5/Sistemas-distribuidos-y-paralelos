#!/bin/bash
#SBATCH --exclusive
#SBATCH --partition=Blade
#SBATCH -o directorioSalida/output_%j.txt
#SBATCH -e directorioSalida/errors_%j.txt
#SBATCH --time=00:20:00

./nreinasSec $1

# Se ejecuta en consola con el comando:
# sbatch ./scriptSec.sh 15

# El script para ejecutar el programa nreinas secuencial provisto por la catedra