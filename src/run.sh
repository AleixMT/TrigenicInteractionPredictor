#!/bin/bash

python_interpreter='python3'
trigenicInteractionPredictor='/home/aleixmt/Escritorio/TrigenicInteractionPredictor/src/TrigenicInteractionPredictor.py'
traintest_datasets_folder='/home/aleixmt/Escritorio/TrigenicInteractionPredictor/data/folds/'
output_base_path=''
arg_k=('2' '3' '4')
trainfile_basename='train'
testfile_basename='test'

# create folder structure to run in parallel
rm -rf results
mkdir results
cd results
for k in ${arg_k[@]}; do
    mkdir K$k
    cd K$k
    for ((fold=0; fold<5; fold++)); do
        mkdir fold$fold
    done
    cd ..
done

# At this point working directory is .../Results/
output_base_path=$(pwd)
cd ..
rm batch_commands.sh
# Generate commands
for k in ${arg_k[@]}; do
    for ((fold=0; fold<5; fold++)); do
         echo -e $python_interpreter $trigenicInteractionPredictor --train=$traintest_datasets_folder$trainfile_basename$fold.dat --test=$traintest_datasets_folder$testfile_basename$fold.dat --k=$k --out=$output_base_path/K$k/fold$fold >> batch_commands.sh
    done
done

parallel --eta --bar --jobs 12 :::: batch_commands.sh

