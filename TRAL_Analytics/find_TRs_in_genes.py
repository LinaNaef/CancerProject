##########################################################################
### Importing required modules
##########################################################################

import os
import pickle
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

working_directory = "/home/lina/SynologyDrive/TRAL_Masterthesis/TRAL_Pipeline_Analytics/TRAL_Analytics/working_files"
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
# genes = ["APC"]
patients = [ patient for patient in os.listdir(sequences_path)]
# patients = ["Pat1","Pat10","Pat34"]

##########################################################################
### Get the sequences from files
##########################################################################

# sequences_per_gene = sequences_per_gene(genes, sequences_path,show_time=True,as_pickle=True)
# sequences_per_gene["MLH1"]["primary_tumor"]["Pat29"]
# sequences_per_gene["APC"]["primary_tumor"]["Pat29"]

# pickle_sequences = retrieve_sequences_pickle(genes)
# pickle_sequences["APC"]["primary_tumor"]["Pat34"]
# pickle_sequences["MLH1"]["solid_tissue_normal"]["Pat28"]

def retrieve_sequences_pickle(genes):

    """ returns a dictionary of sequences for each gene in the input list
            sequences_per_gene[gene][tissue][patient] if saved as pickle in working_directory

    Args:
        genes (list): defines for which genes the sequences should be taken """

    sequences_per_gene = {}
    for gene in genes:
        with open(os.path.join(working_directory, "sequences", gene + ".pkl"), 'rb') as handle:
            sequences = pickle.load(handle)        
        sequences_per_gene.update({gene:sequences})
    return sequences_per_gene

def sequences_per_gene(genes, sequences_path, show_time=False, as_pickle=False):
    
    """ returns a dictionary of sequences for each gene in the input list
            sequences_per_gene[gene][tissue][patient]

    Args:
        genes (list): defines for which genes the sequences should be taken
        sequences_path (string) 
        show_time: prints time if true
        as_pickle: save sequences as pickle if true """

    start = time.time()
# create one dictionary per gene
    sequences_per_gene = {}
    for gene in genes:
        # gene name gives dictionary of all patients for all tissues
        sequences = iterate_patients(sequences_path, gene)
        if as_pickle:
            save_path = os.path.join(working_directory, "sequences", gene + ".pkl")
            try:
                with open(save_path, "wb") as handle:
                    pickle.dump(sequences, handle)
            except (IOError, OSError) as e:
                print("Was not able to create a pickle from the sequences.\n Error:{}".format(e))
        sequences_per_gene.update({gene:sequences})
    
    end = time.time()  
    if show_time:
        print("It took {} seconds to get the sequences.".format(end - start))
    
    return(sequences_per_gene)


def iterate_patients(sequences_path, gene):

    """ returns a dictionary of sequences the gene in the input list
            sequences[tissue][patient]

    Args:
        gene (string)
        sequences_path (string) """

    

    # initialize dictionary for tissues
    sequences = {tissue: {patient: None for patient in patients} for tissue in tissues}

    for patient in patients:
        for tissue in tissues:
            tissue_path = os.path.join(sequences_path, patient, tissue)
            try:
                sequences[tissue][patient] = get_sequences(tissue_path, gene)
                print("Got sequences from tissue {} and patient {}".format(tissue, patient))
            except:
                print("Unexpected Error while trying to call get_sequences({}, {}) .".format(tissue_path, gene))
                exit
    
    return sequences


def get_sequences(tissue_path, gene):

    """ returns TRAL sequences from a tissue

    Args:
        gene (string)
        tissue_path (string) """
      
    sequence_file = os.path.join(tissue_path, gene + ".fasta")
    try:
        sequences_gene = sequence.Sequence.create(file = sequence_file, input_format = 'fasta')
    except FileNotFoundError:
        print("Did not found {} in {}.".format(gene,tissue_path))
        sequences_gene = "Did not found {} in {}.".format(gene,tissue_path)
    except:
        print("Unexpected Error while trying to get the sequence from {}.".format(sequence_file))
        sequences_gene = "Unexpected Error while trying to get the sequence from {}.".format(sequence_file)
    # print("sequences_gene", sequences_gene)
    return sequences_gene