#!/bin/bash
#SBATCH --job-name=distances_hoco_mouse                          # Job name
#SBATCH --mail-type=ALL                          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=tajo5912@colorado.edu                # Where to send mail
#SBATCH --nodes=1                                        # Number of cores job will run on
#SBATCH --ntasks=12                                       # Number of CPU (processers, tasks)
#SBATCH --time=06:00:00                                  # Time limit hrs:min:sec
#SBATCH --partition short                                # Job queue
#SBATCH --mem=150gb                                        # Memory limit
#SBATCH --output=/scratch/Users/tajo5912/eofiles/%x_%j.out
#SBATCH --error=/scratch/Users/tajo5912/eofiles/%x_%j.err
                           
module purge
module load python/3.6.3
module load bedtools/2.25.0
module load samtools/1.8
module load meme/5.0.3

BIN=/Users/tajo5912/rbg/distance/
PRESCAN_MOTIFS=/scratch/Shares/dowell/tajo/mouse_hoco_flat/motifs

NUM=$1
mkdir -p /scratch/Shares/dowell/tajo/bidirectionals_dreg_tfit_rbg/mouse_hocomoco_full/split_${NUM}

python3 ${BIN} \
-a /scratch/Shares/dowell/tajo/bidirectionals_dreg_tfit_rbg/annotations/mm10_qc4_gc50perc_tss50perc_MASTER_tfit_dreg_${NUM}.sorted.bed \
-o /scratch/Shares/dowell/tajo/bidirectionals_dreg_tfit_rbg/mouse_hocomoco_full/split_${NUM} \
-s mm10_hoco_${NUM} \
-v \
-c 12 \
-x False \
-q ${PRESCAN_MOTIFS} \
-w 1500
