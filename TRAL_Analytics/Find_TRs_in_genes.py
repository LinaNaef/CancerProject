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
genes_list = "/home/lina/SynologyDrive/TRAL_Masterthesis/TRAL_Pipeline_Analytics/Prepare_PASS/colorectal_msi_genes.txt"
with open(genes_list) as g:
    genes = g.readlines()
# genes = [x.strip() for x in genes]
# or define list with only few genes:
genes = ["APC","MLH1"]

##########################################################################
### Iteration through files
##########################################################################
# iterate_patients(sequences_path)
def iterate_patients(sequences_path,show_time=False):
    start = time.time()
    sequences = dict()

    for patient in os.listdir(sequences_path):
        patient_path = os.path.join(sequences_path, patient)
        print("patient_path",patient_path)
        for tissue in os.listdir(patient_path):
            tissue_path = os.path.join(patient_path, tissue)
            for gene in genes:
                sequences[patient][tissue][gene] = get_sequences(tissue_path,gene)

            print("Got sequences from tissue {} and patient {}".format(patient, tissue))

        print("Got sequences from {}".format(patient))
    end = time.time()  
    if show_time:
        print("Function get_sequences() took",end - start,"seconds")  

def get_sequences(tissue_path,gene):
        print(gene)
        sequence_file = os.path.join(tissue_path, gene + ".fasta")
        print(sequence_file)
        try:
            sequences_gene = "Test {} in {}.".format(gene,tissue_path)
            # sequences_gene = sequence.Sequence.create(file = sequence_file, input_format = 'fasta')  
        except FileNotFoundError:
            print("Did not found {} in {}.".format(gene,tissue_path))
            sequences_gene = "Did not found {} in {}.".format(gene,tissue_path)
        except:
            print("Unexpected Error while trying to get the sequence from {}.".format(sequence_file))
            sequences_gene = "Unexpected Error while trying to get the sequence from {}.".format(sequence_file)
        return sequences_gene
    

        
    # return sequences_per_TR_type







