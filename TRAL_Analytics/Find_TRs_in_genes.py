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
genes = [x.strip() for x in genes]
# or define list with only few genes:
# genes = ["APC","MLH1"]

##########################################################################
### Iteration through files
##########################################################################
# get_sequences(sequences_path)
def get_sequences(sequences_path,show_time=False):
    
    """ returns a dictionary with detected TRs
            key in form: [l_n]
    
    
    Args:
        sequences_path (str) 
        """
     
    start = time.time()
    sequences_per_TR_type = dict()

    for patient in os.listdir(sequences_path):
        patient_path = os.path.join(sequences_path, patient)
        for tissue in os.listdir(patient_path):
            tissue_path = os.path.join(patient_path, tissue)
            for gene in genes:
                print(gene)
                sequence_file = os.path.join(tissue_path, gene + ".fasta")
                print(sequence_file)
                    


                # if sequence_file.endswith(".fasta"): 
                #     print(tissue_path)
                # else:
                #     continue


    
    # for patient in list_l:
        
    #     for n in list_n:
            
    #         TR_type = "_".join([str(l),str(n)])
    #         filename = TR_type + ".faa"
    #         path = os.path.join(directory_in_str, filename)
    #         sequences = sequence.Sequence.create(file = path, input_format = 'fasta')
    #         sequences_per_TR_type[TR_type] = sequences
    #         print("Got sequences from TR_type",TR_type)
            
    # end = time.time()   
    
    # if show_time:
    #     print("Function get_sequences_in_directory_per_TR_type() took",end - start,"seconds")  
        
    # return sequences_per_TR_type






