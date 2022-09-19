#!/bin/bash
#SBATCH --job-name=fimo_scan                          # Job name
#SBATCH --mail-type=ALL                          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=tajo5912@colorado.edu                # Where to send mail
#SBATCH --nodes=1                                        # Number of cores job will run on
#SBATCH --ntasks=1                                      # Number of CPU (processers, tasks)
#SBATCH --time=10:45:00                                  # Time limit hrs:min:sec
#SBATCH --partition short                                # Job queue
#SBATCH --mem=12gb                                        # Memory limit
#SBATCH --output=/scratch/Users/tajo5912/eofiles/%x_%j.out
#SBATCH --error=/scratch/Users/tajo5912/eofiles/%x_%j.err
                           
module purge
module load python/3.6.3
module load bedtools/2.25.0
module load samtools/1.8
module load meme/5.0.3

fimo_verbosity=$1
threshold_fimo=$2
background_file=$3
motif=$4
full_outdir=$5
motifs=$6
genome=$7

fimo --verbosity $fimo_verbosity --max-stored-scores 10000000 --thresh $threshold_fimo --bfile $background_file --motif $motif --oc $full_outdir $motifs $genome
