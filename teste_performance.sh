#!/bin/bash

# Compilando executáveis em Rust e Fortran
cd ./Rust/Duration/Sequencial || exit
cargo build --release

cd ../Parallel || exit
cargo build --release


cd ../../General/Sequencial || exit
cargo build --release

cd ../Parallel || exit
cargo build --release

cd ../../General_without_Print/Sequencial || exit
cargo build --release

cd ../Parallel || exit
cargo build --release

cd ../../../Fortran/Duration/Sequencial
gfortran -g -O3 -march=native ModuleGlobalData.F90 TCCModule.F90 TCCMain.F90 -o tcc

cd ../../General/Sequencial
gfortran -g -O3 -march=native ModuleGlobalData.F90 TCCModule.F90 TCCMain.F90 -o tcc

cd ../../General_without_Print
gfortran -g -O3 -march=native ModuleGlobalData.F90 TCCModule.F90 TCCMain.F90 -o tcc

cd ../Prepare_Data
gfortran -g -O3 -march=native ModuleGlobalData.F90 TCCModule.F90 TCCMain.F90 -o tcc

cd ../../../


# Preparação dos dados intermediários
FORTRAN_PREPARE_DATA="./Fortran/Prepare_Data/tcc"


# Compilado em Rust e Fortran
RUST_SEQUENCIAL_DURATION_EXEC="./Rust/Duration/Sequencial/target/release/Sequencial"
RUST_PARALLEL_DURATION_EXEC="./Rust/Duration/Parallel/target/release/Parallel"
RUST_SEQUENCIAL_GENERAL_EXEC="./Rust/General/Sequencial/target/release/Sequencial"
RUST_PARALLEL_GENERAL_EXEC="./Rust/General/Parallel/target/release/Parallel"
RUST_SEQUENCIAL_GENERAL_WITHOUT_PRINT_EXEC="./Rust/General_without_Print/Sequencial/target/release/Sequencial"
RUST_PARALLEL_GENERAL_WITHOUT_PRINT_EXEC="./Rust/General_without_Print/Parallel/target/release/Parallel"

FORTRAN_DURATION_EXEC="./Fortran/Duration/Sequencial/tcc"
FORTRAN_GENERAL_EXEC="./Fortran/General/Sequencial/tcc"
FORTRAN_GENERAL_WITHOUT_PRINT_EXEC="./Fortran/General_without_Print/tcc"

# Parâmetros de teste
INPUT_PARAMS=(
    "1_000_000.txt"
    "500_000.txt"
    "250_000.txt"
    "100_000.txt"
    "75_000.txt"
    "50_000.txt"
    "25_000.txt"
    "10_000.txt"
    "5_000.txt"
    1_000.txt
)

CASE_NUMBER=(1 2 3 4 5 6 7 8 9 10)

NUMBER_OF_THREADS=(2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20)


# Criando arquivos intermediários
for INPUT in "${INPUT_PARAMS[@]}"; do
    $FORTRAN_PREPARE_DATA $INPUT
done



# Executando benchmark
for CASE in "${CASE_NUMBER[@]}"; do
    for INPUT in "${INPUT_PARAMS[@]}"; do
        NUMERO="${INPUT%.*}"
        ID="${NUMERO}_${CASE}"

        $RUST_SEQUENCIAL_DURATION_EXEC $CASE $INPUT >> "Data/Output/Rust/Duration/Sequencial/Execution_Result_File/${ID}.txt" 2>&1

    done
done

for CASE in "${CASE_NUMBER[@]}"; do
    for INPUT in "${INPUT_PARAMS[@]}"; do
        NUMERO="${INPUT%.*}"
        ID="${NUMERO}_${CASE}"

        $FORTRAN_DURATION_EXEC $CASE $INPUT >> "Data/Output/Fortran/Duration/Execution_Result_File/${ID}.txt" 2>&1

    done
done


for CASE in "${CASE_NUMBER[@]}"; do
    for THREAD in "${NUMBER_OF_THREADS[@]}"; do
        for INPUT in "${INPUT_PARAMS[@]}"; do
            NUMERO="${INPUT%.*}"
            ID="${THREAD}_${NUMERO}_${CASE}"

            $RUST_PARALLEL_DURATION_EXEC $CASE $THREAD $INPUT >> "Data/Output/Rust/Duration/Parallel/Execution_Result_File/${ID}.txt" 2>&1
        done
    done
done


for CASE in "${CASE_NUMBER[@]}"; do
    for INPUT in "${INPUT_PARAMS[@]}"; do
        NUMERO="${INPUT%.*}"
        ID="${NUMERO}_${CASE}"

        perf stat -e cpu-clock,cycles,cache-misses,cache-references,branch-misses,branches,power/energy-pkg/,power/energy-ram/,page-faults \
            $RUST_SEQUENCIAL_GENERAL_EXEC $CASE $INPUT >> "Data/Output/Rust/General/Sequencial/Execution_Result_File/${ID}.txt" 2>&1

    done
done


for CASE in "${CASE_NUMBER[@]}"; do
    for INPUT in "${INPUT_PARAMS[@]}"; do
        NUMERO="${INPUT%.*}"
        ID="${NUMERO}_${CASE}"

        perf stat -e cpu-clock,cycles,cache-misses,cache-references,branch-misses,branches,power/energy-pkg/,power/energy-ram/,page-faults \
            $FORTRAN_GENERAL_EXEC $CASE $INPUT >> "Data/Output/Fortran/General/Execution_Result_File/${ID}.txt" 2>&1

    done
done


for CASE in "${CASE_NUMBER[@]}"; do
    for THREAD in "${NUMBER_OF_THREADS[@]}"; do
        for INPUT in "${INPUT_PARAMS[@]}"; do
            NUMERO="${INPUT%.*}"
            ID="${THREAD}_${NUMERO}_${CASE}"

            perf stat -e cpu-clock,cycles,cache-misses,cache-references,branch-misses,branches,power/energy-pkg/,power/energy-ram/,page-faults \
                $RUST_PARALLEL_GENERAL_EXEC $CASE $THREAD $INPUT >> "Data/Output/Rust/General/Parallel/Execution_Result_File/${ID}.txt" 2>&1
        done
    done
done


for CASE in "${CASE_NUMBER[@]}"; do
    for INPUT in "${INPUT_PARAMS[@]}"; do
        NUMERO="${INPUT%.*}"
        ID="${NUMERO}_${CASE}"

        perf stat -e cpu-clock,cycles,cache-misses,cache-references,branch-misses,branches,power/energy-pkg/,power/energy-ram/,page-faults \
            $RUST_SEQUENCIAL_GENERAL_WITHOUT_PRINT_EXEC $CASE $INPUT >> "Data/Output/Rust/General_without_Print/Sequencial/Execution_Result_File/${ID}.txt" 2>&1

    done
done


for CASE in "${CASE_NUMBER[@]}"; do
    for INPUT in "${INPUT_PARAMS[@]}"; do
        NUMERO="${INPUT%.*}"
        ID="${NUMERO}_${CASE}"

        perf stat -e cpu-clock,cycles,cache-misses,cache-references,branch-misses,branches,power/energy-pkg/,power/energy-ram/,page-faults \
            $FORTRAN_GENERAL_WITHOUT_PRINT_EXEC $CASE $INPUT >> "Data/Output/Fortran/General_without_Print/Execution_Result_File/${ID}.txt" 2>&1

    done
done


for CASE in "${CASE_NUMBER[@]}"; do
    for THREAD in "${NUMBER_OF_THREADS[@]}"; do
        for INPUT in "${INPUT_PARAMS[@]}"; do
            NUMERO="${INPUT%.*}"
            ID="${THREAD}_${NUMERO}_${CASE}"

            perf stat -e cpu-clock,cycles,cache-misses,cache-references,branch-misses,branches,power/energy-pkg/,power/energy-ram/,page-faults \
                $RUST_PARALLEL_GENERAL_WITHOUT_PRINT_EXEC $CASE $THREAD $INPUT >> "Data/Output/Rust/General_without_Print/Parallel/Execution_Result_File/${ID}.txt" 2>&1
        done
    done
done