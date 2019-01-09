##########################################################################
### Importing required modules
##########################################################################

import os

# modules in this folder
import get_sequences
import find_TRs_in_genes

##########################################################################
### Defining paths and variables
##########################################################################

genes = ["APC","MLH1"]
# genes = ["APC"]
# patients = [ patient for patient in os.listdir(sequences_path)]
patients = ["Pat1","Pat10","Pat34"]

working_directory = "/home/lina/SynologyDrive/TRAL_Masterthesis/TRAL_Pipeline_Analytics/TRAL_Analytics/working_files"
sequences_path = "/home/lina/SynologyDrive/TRAL_Masterthesis/IBM_files/Assembled_genes"

##########################################################################
### Get the sequences from files
##########################################################################

# sequences_per_gene = sequences_per_gene(genes, sequences_path, show_time=True, as_pickle=True, patient=True)
# sequences_per_gene["MLH1"]["primary_tumor"]["Pat29"]
# sequences_per_gene["APC"]["primary_tumor"]["Pat29"]

pickle_sequences = get_sequences.retrieve_sequences_pickle(working_directory, genes)
# pickle_sequences["APC"]["primary_tumor"]["Pat34"]
# pickle_sequences["MLH1"]["solid_tissue_normal"]["Pat28"]

##########################################################################       
### use detector for sequences and get tandem repeats
##########################################################################


test_TRs = find_TRs_in_genes.get_TRs("XSTREAM", pickle_sequences["APC"]["primary_tumor"]["Pat34"])
find_TRs_in_genes.loop_list_TR(test_TRs)