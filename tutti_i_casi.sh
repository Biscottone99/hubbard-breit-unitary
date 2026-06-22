#!/bin/bash
set -e # Ferma l'esecuzione se un qualsiasi comando (es. compilazione) fallisce

#============ PARAMETRI =====================================
dimension=70
rho_file="rho.bin"
eigen_file="eigen.bin"
n_prop=2

# Array con i valori su cui iterare
hoppings=("ncn" "nce")
u_values=(10)

# Array Bash per gestire la lista dei file delle proprietà
props_files=("spin-density.bin" "s2rot.bin")

deltat=700
points=3000

# Path base per evitare percorsi hardcoded ripetuti e facilitare la manutenzione
BASE_DIR="/home/biscottone/ricerca/ciss/unitaria"
#==============================================================

echo "--- Compilazione dei codici Fortran ---"
ifx basis.f90 -o basis.e
./basis.e

ifx geometria.f90 -o prova.e
./prova.e

ulimit -s unlimited
ifx newmodule.f90 ppp-breit-pauli.f90 -o breit.e -qmkl -qopenmp -traceback

# Compila unitaria.e una volta sola nella sua cartella
pushd unitary > /dev/null
ifx Unitaria.f90 -o unitaria.e -qmkl -qopenmp
popd > /dev/null

echo "--- Compilazioni completate ---"

# ==============================================================================
# FUNZIONE DI PROPAGAZIONE UNITARIA
# ==============================================================================
run_model() {
    local input_template=$1
    local out_name=$2
    local cur_hopping=$3
    local cur_u=$4

    echo -e "\n  [>] Esecuzione: $out_name"
    cp "input/${input_template}_${cur_hopping}.inp" input.inp
    
    # Esecuzione breit.e passando il valore di u aggiornato
    ./breit.e <<< "$cur_u"
    
    cp output.out "model_${out_name}_u${cur_u}_${cur_hopping}.out"
    cp psi0.dat "model_${out_name}_u${cur_u}_${cur_hopping}.psi"
    
    # pushd cambia directory ricordando dove eravamo
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
    mv prop1.dat "model_${out_name}_u${cur_u}_${cur_hopping}.dat"
    mv prop2.dat "s2_evolution_model_${out_name}_u${cur_u}_${cur_hopping}.dat"

    # popd torna esattamente alla cartella di partenza
    popd > /dev/null
}

# ==============================================================================
# CICLO PRINCIPALE SUI PARAMETRI
# ==============================================================================

for hopping in "${hoppings[@]}"; do
    for u in "${u_values[@]}"; do
        
        # Aggiorna la variabile 'fine' per il ciclo corrente
        fine="u${u}_${hopping}"
        
        echo -e "\n================================================================="
        echo " INIZIO CALCOLI: hopping = $hopping | U = $u"
        echo "================================================================="

        # Esecuzione dei 4 modelli passando hopping e u in modo sicuro
        run_model "input-mono-model_a" "a_mono" "$hopping" "$u"
        run_model "input-tot-model_a" "a_tot" "$hopping" "$u"
        run_model "input-mono-model_b" "b_mono" "$hopping" "$u"
        run_model "input-tot-model_b" "b_tot" "$hopping" "$u"

        # Generazione grafico per il set corrente
        echo "  [>] Generazione dei grafici per $fine..."
        pushd unitary > /dev/null
        python3 figure.py <<< "$fine"
        popd > /dev/null

    done
done

# ==============================================================================
# SPOSTAMENTO FINALE DEI FILE
# ==============================================================================
echo -e "\n[>] Organizzazione dei file nelle directory di output..."

# Assicuriamoci che le cartelle esistano prima di spostarci i file
mkdir -p "$BASE_DIR/unitary/spin-polarization/0_perc_gs"
mkdir -p "$BASE_DIR/output_file/0_perc_gs/nearest-neighbor"
mkdir -p "$BASE_DIR/output_file/0_perc_gs/non-coherent-nearest-neighbor"
mkdir -p "$BASE_DIR/output_file/0_perc_gs/non-coherent-exponential_decay"
mkdir -p "$BASE_DIR/output_file/0_perc_gs/exponential_decay"
mkdir -p "$BASE_DIR/unitary/0_perc_gs"

# Sposta i plot PNG generati in unitary
mv "$BASE_DIR"/unitary/*.png "$BASE_DIR/unitary/spin-polarization/0_perc_gs/" 2>/dev/null || true

# Sposta i log di output generati nella root (filtrati strettamente per non toccare gli eseguibili)
mv "$BASE_DIR"/model_*_nn.* "$BASE_DIR/output_file/0_perc_gs/nearest-neighbor/" 2>/dev/null || true
mv "$BASE_DIR"/model_*_ncn.* "$BASE_DIR/output_file/0_perc_gs/non-coherent-nearest-neighbor/" 2>/dev/null || true
mv "$BASE_DIR"/model_*_nce.* "$BASE_DIR/output_file/0_perc_gs/non-coherent-exponential_decay/" 2>/dev/null || true
mv "$BASE_DIR"/model_*_e.* "$BASE_DIR/output_file/0_perc_gs/exponential_decay/" 2>/dev/null || true

# Sposta i file di dinamica temporale (.dat) che si trovano dentro unitary
mv "$BASE_DIR"/unitary/*.dat "$BASE_DIR/unitary/0_perc_gs/" 2>/dev/null || true

echo -e "\n================================================================="
echo " Script completato con successo per tutte le combinazioni! "
echo "================================================================="
