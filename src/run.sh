#!/bin/bash

pypy_bin='/export/home/shared/Projects/TrigenicInteractionPredictor/TrigenicInteractionPredictor/src/pypy3.5-v7.0.0-linux64/bin/pypy3'
trigenicInteractionPredictor='/export/home/shared/Projects/TrigenicInteractionPredictor/TrigenicInteractionPredictor/src/TrigenicInteractionPredictor.py'
traintest_datasets_folder='/export/home/shared/Projects/TrigenicInteractionPredictor/TrigenicInteractionPredictor/src/'
arg_k=('2' '3' '4' '5')
trainfile_basename='train'
testfile_basename='test'

rm -rf results
mkdir results
cd results
for k in ${arg_k[@]}; do
    mkdir K$k
    cd K$k
    for ((fold=0; fold<5; fold++)); do
        mkdir fold$fold
        cd fold$fold
        echo '$pypy_bin $trigenicInteractionPredictor --train=$traintest_datasets_folder$trainfile_basename$fold.dat --test=$traintest_datasets_folder$testfile_basename$fold.dat --k=$k'
        $pypy_bin $trigenicInteractionPredictor --train=$traintest_datasets_folder$trainfile_basename$fold.dat --test=$traintest_datasets_folder$testfile_basename$fold.dat --k=$k
        cd ..
    done
    cd ..
done

