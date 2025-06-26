#!/bin/bash
#SBATCH -A GREEN-SL2-CPU
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:59:59
#SBATCH --mail-type=FAIL
#SBATCH --partition=cclake

module load singularity
SLX=slx26226
/usr/bin/singularity exec \
   --bind bam_$SLX:/input \
  --bind ${SLX}_evidence:/output \
  --bind ${SLX}_insertions.xlsx:/input_excel.xlsx \
  --bind collect_evidence.py:/collect_evidence.py \
  /home/jd862/rds/rds-sjl83-green/jd862/py13_pandas_pysam.sif \
  python /collect_evidence.py $1