#!/bin/bash
### Created by Achim Klittich
################################################################################
### SLURM CONFIG
################################################################################
#SBATCH --job-name="windowtrees"
### runtime
#SBATCH --time=2-00:00:00
### changedir / output
#SBATCH --chdir=/hpc-cloud/user/windowTrees
#SBATCH --output=/hpc-cloud/user/windowTrees/WindowTrees_%j.log
#SBATCH --error=/hpc-cloud/user/windowtTees/WindowTrees_%j.log
### Account / Partition
#SBATCH --partition=All ### add partition here
#SBATCH --account=all ### add account here if necessary
### NODE / CPU settings
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
### MEM per node
#SBATCH --mem=32G
### Notifications
#SBATCH --mail-type=END
#SBATCH --mail-user=email@email.de ### add your email address here to get notified
################################################################################
### Print Date and get time
################################################################################
date
start_time="$(date -u +%s)"
################################################################################
### Print Slurmscript to log
################################################################################
echo "################################################################################"
echo -e "Used SLURM script:\n"
cat $0
echo -e "\n################################################################################"
echo "### Output:"
echo "################################################################################"

################################################################################
### Define variables
################################################################################
windowtree_py=/hpc-cloud/user/WindowTrees/windowtrees.py
inputfile=/hpc-cloud/user/WindowTrees/inputfile_test.tab
outpath=/hpc-cloud/user/WindowTrees/run1
outgroup="BELUGA"
windowsize=1000000
nthreshold=0.2
ncpu=10


################################################################################
### Running code
################################################################################
# print input file to log
echo "Utilized sample input file:"
cat $inputfile
echo "Starting WindowTrees..."
# create output dir
mkdir $outpath
# start WindowTree py with arguments defined above
python3 -u $windowtree_py -o $outpath --outgroup $outgroup -w $windowsize --cpu $ncpu $inputfile

################################################################################
### Print elapsed time and date
################################################################################
end_time="$(date -u +%s)"
elapsed="$(($end_time-$start_time))"
echo "Job finished after $elapsed sec"
date
