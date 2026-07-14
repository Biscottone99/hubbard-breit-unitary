#!/bin/bash
set -e

#============ PARAMETRI =====================================
dimension=70
rho_file="rho.bin"       
eigen_file="eigen.bin"  
n_prop=1

# Array con i valori su cui iterare

u_values=10

# Array Bash per gestire la lista dei file delle proprietà
props_files=("spin-density.bin") # Aggiunto ../

deltat=700
points=3000

#==============================================================



ifx basis.f90 -o basis.e
./basis.e
ifx geometria.f90 -o prova.e
./prova.e
#python3 rotate.py
ulimit -s unlimited
ifx newmodule.f90 ppp-breit-pauli.f90 -o vb.e -qmkl -qopenmp -traceback 
./vb.e


pushd unitary > /dev/null
ifx Unitaria.f90 -o unitaria.e -qmkl -qopenmp
popd > /dev/null

echo "--- Compilazioni completate ---"

# Entro in unitary per l'esecuzione
pushd unitary > /dev/null

{
    echo "$dimension"
    echo "$rho_file"
    echo "$eigen_file"
    echo "$n_prop"

    for file in "${props_files[@]}"; do
        echo "$file"
    done

    echo "$deltat"
    echo "$points"
} | ./unitaria.e > /dev/null

# Rinomina dei file di output della dinamica (Spin Polarization e S^2)
# Aggiunto un controllo di esistenza per evitare crash con "set -e"
if [[ -f prop1.dat ]]; then
    mv -f prop1.dat spin-pol.dat
else
    echo "Attenzione: prop1.dat non generato!" >&2
fi

if [[ -f prop2.dat ]]; then
    mv -f prop2.dat s2_evolution.dat
else
    echo "Attenzione: prop2.dat non generato!" >&2
fi

popd > /dev/null

echo "--- Esecuzione completata con successo ---"
