#!/bin/bash

# Generate independent calls to main script in order to run in parallel

python_interpreter='/home/aleixmt/Desktop/pypy3.5-v7.0.0-linux64/bin/pypy3'
trigenicInteractionPredictor='/home/aleixmt/Desktop/GIT/TrigenicInteractionPredictor/src/TrigenicInteractionPredictor.py'
traintest_datasets_folder='/home/aleixmt/Desktop/GIT/TrigenicInteractionPredictor/data/folds/'
output_base_path=''
arg_k=('2')
trainfile_basename='train'
testfile_basename='test'

# create folder structure to run in parallel
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
         echo -e $python_interpreter $trigenicInteractionPredictor --train=$traintest_datasets_folder$trainfile_basename$fold.dat --test=$traintest_datasets_folder$testfile_basename$fold.dat --k=$k --out=$output_base_path/K$k/fold$fold/ >> batch_commands.sh
    done
done

parallel --eta --bar --jobs 0 :::: batch_commands.sh

