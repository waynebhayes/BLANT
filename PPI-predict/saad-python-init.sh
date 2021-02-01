#!/bin/bash

module load anaconda3/2020.11

eval "$(/pkg/anaconda3/2020.11/bin/conda shell.bash hook)";

conda init;

# conda create -n py38 python=3.8.5 -y;

conda activate py38;

conda install -c anaconda numpy -y;

conda install -c conda-forge matplotlib -y;

conda install -c conda-forge scikit-learn -y;

conda deactivate py38;

echo "Done";
