#!/bin/bash
#SBATCH --job-name=fimo_scanner                           # Job name
#SBATCH --mail-type=ALL                          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=tajo5912@colorado.edu                # Where to send mail
#SBATCH --nodes=1                                        # Number of cores job will run on
#SBATCH --ntasks=1                                       # Number of CPU (processers, tasks)
#SBATCH --time=01:50:00                                  # Time limit hrs:min:sec
#SBATCH --partition short                                # Job queue
#SBATCH --mem=8gb                                        # Memory limit
#SBATCH --output=/scratch/Shares/dowell/tajo/hoco_flat/eofiles/%x_%j.out
#SBATCH --error=/scratch/Shares/dowell/tajo/hoco_flat/eofiles/%x_%j.err
                           
module purge
module load python/3.6.3
module load meme/5.0.3

BIN=/Users/tajo5912/fimo_scanner/fimo_scanner/whole_genome_fimo_scan.py

python3 ${BIN} \
-g /scratch/Shares/dowell/genomes/hg38/hg38.fa \
-o /scratch/Shares/dowell/tajo/hoco_flat \
-m /scratch/Shares/dowell/tajo/hoco_flat/scripts/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme \
-v \
-t 1e-5 \
-b /scratch/Shares/dowell/tajo/hoco_flat/scripts/enhancer_background_flat
