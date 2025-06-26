#!/bin/bash

for i in {1..96}; do
  echo "submitting $i"
  bc=xGenUDI$i-UMI
  sbatch -J $bc -o slx26173_evidence/$bc.%j.out -e slx26173_evidence/$bc.%j.err collect_evidence.slurm $bc
done