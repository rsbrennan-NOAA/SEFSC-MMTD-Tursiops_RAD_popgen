#!/bin/bash
#SBATCH -D /home/rbrennan/Tursiops-RAD-popgen/logout
#SBATCH --mail-type=END
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --partition=standard
#SBATCH --cpus-per-task=16
#SBATCH --mem=70G
#SBATCH --time=40:00:00
#SBATCH --job-name=markdups
#SBATCH --output=%x.%A.%a.out
#SBATCH --error=%x.%A.%a.err
#SBATCH --array=[1-3]

indir=~/Tursiops-RAD-popgen/analysis/pop_structure/treemix/aduncus/aligned
outdir=~/Tursiops-RAD-popgen/analysis/pop_structure/treemix/aduncus/aligned/mkdup

cd $indir
inbam=$(ls *.bam | cut -f 1 -d "." |  sed -n $(echo $SLURM_ARRAY_TASK_ID)p)

echo $inbam

~/bin/sambamba-1.0.1 markdup --nthreads=16  ${indir}/${inbam}.bam ${outdir}/${inbam}_mkdup.bam
