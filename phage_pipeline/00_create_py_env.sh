#!/bin/bash

## VARIABLES

## Directory containing the virtual environment
py_env_dir=~/envs/coev

## url for the desired python3 executable
py_version=$(readlink $(command -v python3))

echo 'Environment to be installed in: ' $py_env_dir
echo 'Version of python to be used: ' $py_version


## SCRIPT

mkdir -p  $py_env_dir
virtualenv  -p $py_version  $py_env_dir


## INSTALL  PACKAGES

# Activate environment
source  ${py_env_dir}/bin/activate

# Update pip if necessary
pip install --upgrade pip

## 
# pip install  -r scripts/requirements_py-env.txt
pip install biopython multiqc numpy pandas ipython pyvcf
# pip  install  biopython  pandas   matplotlib  multiqc  pyvcf


deactivate
