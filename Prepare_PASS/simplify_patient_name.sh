#!/bin/bash

# Simplify file names


cd "/home/lina/Desktop/Unassembled_full_copy"

# patient_ids="/home/lina/Desktop/TRAL_Masterthesis/TRAL_Pipeline_Analytics/Prepare_PASS/patient_ids.txt"
patient_ids_key="/home/lina/Desktop/TRAL_Masterthesis/TRAL_Pipeline_Analytics/Prepare_PASS/patient_ids_key.txt"
# ls > "$patient_ids"

# https://askubuntu.com/questions/318611/using-text-list-to-batch-rename-files

# awk '{print $0" Pat",NR}' "$patient_ids" > $patient_ids_key
#manually deleted the spaces between pat and number

# rename files in a directory according to a list
xargs -a $patient_ids_key -n 2 mv