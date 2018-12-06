##########################################################################


### Importing required modules
##########################################################################

import os
import pickle
import time
import importlib

# to use this modules, TRAL has to be installed properly
from tral.sequence import sequence
from tral.repeat import repeat
from tral.repeat_list import repeat_list

##########################################################################
### Defining paths and variables
##########################################################################

sequences_path = "/home/lina/SynologyDrive/TRAL_Masterthesis/IBM_files/Assembled_genes"

## tissues to analyze
tissues = ["blood_derived_normal","primary_tumor","solid_tissue_normal"]

## genes to analyze
# genes_list = "/home/lina/SynologyDrive/TRAL_Masterthesis/TRAL_Pipeline_Analytics/Prepare_PASS/colorectal_msi_genes.txt"
# with open(genes_list) as g:
#     genes = g.readlines()
# genes = [x.strip() for x in genes]
# or define list with only few genes:
genes = ["APC","MLH1"]
patients = [ patient for patient in os.listdir(sequences_path)]

##########################################################################
### Iteration through files
##########################################################################

# sequences["primary_tumor"]["Pat1"]
# sequences["blood_derived_normal"]["Pat14"]
# sequences["solid_tissue_normal"]["Pat43"]

# create one dictionary per gene
for gene in genes:
    # gene name gives dictionary of all patients for all tissues
    sequences = iterate_patients(sequences_path, gene)
    # vars()[gene] = iterate_patients(sequences_path, gene)
    # print(gene)
gene = "APC"
sequences = iterate_patients(sequences_path, gene)

def iterate_patients(sequences_path, gene, show_time=False):
    start = time.time()

    # initialize dictionary for tissues
    sequences = {tissue: {patient: None for patient in patients} for tissue in tissues}

    for patient in patients:
        for tissue in tissues:
            tissue_path = os.path.join(sequences_path, tissue, patient)
            sequences[tissue][patient] = get_sequences(tissue_path, gene)

            print("Got sequences from tissue {} and patient {}".format(tissue, patient))

    end = time.time()  
    if show_time:
        print("Function get_sequences() took", end - start,"seconds")
    
    return sequences

def get_sequences(tissue_path,gene):
    sequence_file = os.path.join(tissue_path, gene + ".fasta")
    try:
        sequences_gene = "Test {} in {}.".format(gene,sequence_file)
        # sequences_gene = sequence.Sequence.create(file = sequence_file, input_format = 'fasta')
    except FileNotFoundError:
        print("Did not found {} in {}.".format(gene,tissue_path))
        sequences_gene = "Did not found {} in {}.".format(gene,tissue_path)
    except:
        print("Unexpected Error while trying to get the sequence from {}.".format(sequence_file))
        sequences_gene = "Unexpected Error while trying to get the sequence from {}.".format(sequence_file)
    print("sequences_gene", sequences_gene)
    return sequences_gene