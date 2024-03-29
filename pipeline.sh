#!/bin/bash

echo Please enter the name of the project

read varname
env_loc=$(conda info --base)/etc/profile.d/conda.sh

source $env_loc
conda activate natrix_repair

varname+=".yaml"
cores=$(grep "cores : " $varname | awk '{print $3}')
pip install pybloomfiltermmap3
python setup.py install
python create_dataframe.py "$varname"

screen -S $varname bash -c "source $env_loc;conda activate natrix_repair;snakemake --use-conda --cores $cores --configfile $varname; exec sh"
