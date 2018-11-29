#!/bin/bash

# Extract all genes listed in a txt file from SAM-file and convert it to FASTA
# Routine used by An-phi (20181126)

ROOT="/ibm/gpfs-dataT/uye/tcga/tcga/colorectal/protein_alignment"
CASE_TYPE=(
	"blood_derived_normal"
	"primary_tumor"
	"solid_tissue_normal"
)
OUT_ROOT="/ibm/gpfs-dataT/uye/tcga/tcga/colorectal/extracted_genes"
GENES_FILE="/ibm/gpfs-dataT/uye/tcga/tcga/colorectal/tmp/colorectal_msi_genes.txt" # plain text file with genes

for i in ${ROOT}/*; do
	case_id=$(basename $i)
	echo $case_id

	out_case=${OUT_ROOT}/${case_id}
	mkdir -p $out_case

	for st in ${CASE_TYPE[@]}; do
		
		out_case_sample=${out_case}/${st}
		mkdir -p ${out_case_sample}

		while read gene; do
			
			#out_file=${out_case_sample}/${gene}.sam			
			out_fasta=${out_case_sample}/${gene}.fasta

			echo $gene

			grep -hw "^\([^|]*|\)\{2\}${gene}_HUMAN" ${i}/${st}/f*.sam | sed '/^@/ d' | awk 'NR%2==0 {print ">"$1"\n"$10}' > $out_fasta


		done < ${GENES_FILE}
	done 
done
