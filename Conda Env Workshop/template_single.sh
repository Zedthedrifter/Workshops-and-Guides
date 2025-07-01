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


#Install samtools, NanoPlot, porechop, chopper
#these are the old QC packages from amplicon-seq project, some might be useful, such as samtools, biopython,etc
function install_basic {    

#install to specific env
conda install -c bioconda samtools -y --name $env_name
conda install pip -y --name $env_name #install pip
}

#hacking with pip
#sometimes when the path is not recognized: can always give a full path
#can also add this to .bashrc permanently... but keeping it within the environment is a safer choice
#use pip to install a lot of pure python packages (when you just need python)
function pip_install {

pip=$HOME/projects/rbge/$USER/env/$env_name/bin/pip #install with pip
$pip install psutil 
$pip install biopython
$pip install matplotlib
$pip install requests 
$pip install beautifulsoup4
$pip install amas
$pip install numpy pandas
$pip install seaborn 
}

#--------------------------------------------------------------------------------------------------------------------
#these 
#occasionally when you have to install directly from the sourcecode on github
function install_trimal { 

mkdir $HOME/projects/rbge/$USER/env/$env_name/bin/trimal
git clone https://github.com/inab/trimal.git $HOME/projects/rbge/$USER/env/$env_name/bin/trimal
cd $HOME/projects/rbge/$USER/env/$env_name/bin/trimal/source
make

}

function install_tree {

conda install -c bioconda iqtree -y --name $env_name
conda install -c bioconda figtree -y --name $env_name
conda install -c conda-forge mafft -y --name $env_name
conda install -c bioconda orthofinder -y --name $env_name
conda install -c anaconda scipy -y --name $env_name
}


function setup_env {

conda create -n $env_name -y #make an env 
install_basic #conda install
install_tree #can organize different groups of packages into different funtions (later if you want to set up a different env, you can copy it over easily)
pip_install #for special installing method
#remember to add path
install_trimal #trim alignment: remove regions of poor alignment: https://vicfero.github.io/trimal/index.html#installation_sec
}

function main {

env_name=env_test
USER=zedchen 
####manual configuration 
easy353path=$HOME/projects/rbge/$USER/env/$env_name/bin/Easy353 #this way you don't need to config in .bashrc and thus contain all the packages/path within the env
trimalpath=$HOME/projects/rbge/$USER/env/$env_name/bin/trimal/source
#USAGES
#=============================================
setup_env #the function that actualy set up env
}

main


