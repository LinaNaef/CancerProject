#!/bin/bash

# directory with unassembled genes in fasta (sorted for each patient and each tissue)
ROOT="/home/lina/SynologyDrive/TRAL_Masterthesis/IBM_files/Unassembled_full_copy"
CASE_TYPE=(
	"blood_derived_normal"
	"primary_tumor"
	"solid_tissue_normal"
)
# output directory
OUT_ROOT="/home/lina/SynologyDrive/TRAL_Masterthesis/IBM_files/Assembled_genes"
#GENES_FILE="/ibm/gpfs-dataT/uye/tcga/tcga/colorectal/tmp/colorectal_msi_genes.txt"

for patient in "${ROOT}/"*; do
	case_id=$(basename "$patient")
	echo "$case_id"

    # create a new directory per patient
	out_case="${OUT_ROOT}/${case_id}"
	mkdir -p "$out_case"

	for tissue in "${CASE_TYPE[@]}"; do
		
        # create a new directory per tissue
		out_case_sample="${out_case}/${tissue}"
		mkdir -p "${out_case_sample}"
        
        # echo "${CASE_TYPE[@]}"

        for unassembled_gene in "${ROOT}/${case_id}/$tissue/"*.fasta; do

            # only if fasta file contains sequences pass it to PASS
            if [ -s "$unassembled_gene" ]; then
                gene=$(basename -s .fasta "$unassembled_gene")
                echo "$gene"
                
                assembled_gene="${out_case_sample}/$gene" # filename
                mkdir -p "$assembled_gene"
                echo "$assembled_gene"
                /home/lina/SynologyDrive/TRAL_Masterthesis/programs/pass_v0.3/PASS -f "$unassembled_gene" -b "$assembled_gene/$gene" -m 4 -w 1 -o 1 -r 0.51 # these are just default variables, should be adapted
                echo $(cat "$assembled_gene/$gene.contigs" "$assembled_gene/$gene.singlets" > "$assembled_gene.fasta") # merge singlets and contigs
            fi

            # -f  File containing all the peptide reads (required)
            # -b  Basename for your output files (optional)
            # -m  Minimum number of overlapping amino acids with the seed/contig during overhang consensus build up (default -m 10)
            # -w  Minimum depth of coverage allowed for contigs (e.g. -w 1 = process all reads, required)
                # *The assembly will stop when 50+ contigs with coverage < -w have been seen.*
            # -o  Minimum number of peptide reads needed to call a amino acid during an extension (default -o 2)
            # -r  Minimum ratio used to accept a overhang consensus amino acid (default -r 0.7)


        done
	done
done
