#!/bin/bash

# Simplify file names


cd "/home/lina/SynologyDrive/TRAL_Masterthesis/IBM_files/Unassembled_genes_fasta_full"
patient_ids="/home/lina/Desktop/TRAL_Masterthesis/TRAL_Pipeline_Analytics/Prepare_PASS/patient_ids.txt"
patient_ids_key="/home/lina/Desktop/TRAL_Masterthesis/TRAL_Pipeline_Analytics/Prepare_PASS/patient_ids_key.txt"
ls > "$patient_ids"


awk '{print $0" Pat",NR}' "$patient_ids" > $patient_ids_key
#manually deleted the spaces between pat and number

xargs -a $patient_ids_key -n 2 mv