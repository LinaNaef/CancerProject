##!/bin/bash

source /scratch/IAS/AnisGroup/naefpau1/env_tral/bin/activate
python /home/naefpau1/TRAL_Analytics/TR_per_Chromosome.py 23

for i in {1..23};
do echo /home/naefpau1/TRAL_Analytics/TR_per_Chromosome.py $i;
done
