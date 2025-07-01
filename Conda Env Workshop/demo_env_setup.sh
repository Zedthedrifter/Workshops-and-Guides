#!/bin/bash
#SBATCH --job-name="your project"
#SBATCH --export=ALL
#SBATCH --partition=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G 

#
#--------------------------------------------------------------------------------------------------------------------
#use srun after navigating to scripts directory on /scratch
#run this after activating the env 'salix'


function setup_env {

conda config --add channels dnachun
conda config --add channels micromacro
#Set channel priority:
conda config --set channel_priority flexible
#Create environment:
conda create --name $env_name r-base=4.2.2 python=3.10.12 -y
#Install PanSyn:
mamba install -c bioconda -c conda-forge \
    bioconductor-genomeinfodbdata \
    bioconductor-annotationdbi \
    r-base=4.4.2 --name $env_name
mamba install -c micromacro pansyn –y --name $env_name

}

function main {

env_name=env_demo
#USAGES
#=============================================
setup_env #the function that actualy set up env
}

main


