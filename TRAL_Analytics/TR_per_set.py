import sys

##########################################################################
### Importing required modules
##########################################################################

import os
import pickle

# # modules in this folder
import get_sequences
import find_TRs_in_genes
from pyfaidx import Fasta

import logging
import logging.config
import os

from tral.paths import config_file, PACKAGE_DIRECTORY
from tral import configuration
from tral.sequence import sequence
from tral.hmm import hmm

logging.config.fileConfig(config_file("logging.ini"))
log = logging.getLogger('root')

# define this as in the config
CONFIG_GENERAL = configuration.Configuration.instance().config
CONFIG = CONFIG_GENERAL["repeat_list"]
score = CONFIG["model"]

##########################################################################
######### Defining Paths and Parameters

## AA reference
working_dir = "/home/lina/Desktop/TRAL_Masterthesis/references/Uniprot_data/proteom/pickles"
sequences_path = "/home/lina/Desktop/TRAL_Masterthesis/references/Uniprot_data/proteom"
output_path = "/home/lina/SynologyDrive/TRAL_Masterthesis/TRAL_Pipeline_Analytics/test_output/proteom"

# Thresholds for filtering
pvalue_threshold = 0.05
divergence_threshold = 0.8
n_threshold = 2.5 # minimun repeat unit count
l_threshold = 3 # maximum repeat unit length

# set_file = "AUP000005640_Mitochondrion.fasta"
set_file = "AUP000005640_Chromosome21.fasta"
set_name = set_file.split(".fasta")[0]
sequences_file = os.path.join(sequences_path, set_file)
result_dir = os.path.join(output_path, set_name)
try:
    if not os.path.isdir(result_dir):
        os.makedirs(result_dir)
except:
    raise Exception(
        "Could not create path to result directory: {}".format(
            os.path.dirname(result_dir)))


##########################################################################
######### Getting sequences

## obtaining sequences from fasta format
## Pyfaix Documentation (https://pythonhosted.org/pyfaidx/#installation)
proteins = Fasta(sequences_file)
# print(proteins.keys())
# proteins['sp|P03886|NU1M_HUMAN'][:].seq

for pyfaidx in proteins:
    seq_name = pyfaidx.name.split("|")[1]
    sequence_pkl = os.path.join(working_dir, seq_name + "_sequence.pkl")

    if os.path.exists(sequence_pkl):
        # Getting back the sequences:
        with open(sequence_pkl,'rb') as f: # files saved before  
            seq = pickle.load(f)
    else:
        seq = sequence.Sequence(seq=str(pyfaidx), name=seq_name) # name is protein identifier
        # Saving this sequences as binary files:
        with open(sequence_pkl, 'wb') as f: 
            pickle.dump(seq, f)

    log.debug("Work on sequence {}".format(seq_name))
    ##########################################################################
    ######### Getting TRs

    TRs_pkl = os.path.join(working_dir, "TRs_raw", seq_name + "_TRs.pkl")

    if os.path.exists(TRs_pkl):
        # Getting back the sequences:
        with open(TRs_pkl,'rb') as f: # files saved before  
            denovo_list = pickle.load(f)
    else:
        denovo_list = seq.detect(denovo=True)
        # Saving this sequences as binary files:
        with open(TRs_pkl, 'wb') as f: 
            pickle.dump(denovo_list, f)

    # print("Found", len(denovo_list.repeats), "denovo repeats in", seq_name)

    ##########################################################################
    ######### Filtering TRs


    output_pickle_file = os.path.join(result_dir, seq_name + ".pkl")
    output_tsv_file = os.path.join(result_dir, seq_name + ".tsv")

    if os.path.exists(output_pickle_file) and os.path.exists(output_tsv_file):
        with open(output_pickle_file,'rb') as f: 
            denovo_list_remastered = pickle.load(f)
    else:
        # filtering for pvalue
        denovo_list = denovo_list.filter(
            "pvalue",
            score,
            pvalue_threshold)

        ## filtering for divergence
        denovo_list = denovo_list.filter(
            "divergence",
            score,
            divergence_threshold)

        ## filtering for number of repeat units
        denovo_list = denovo_list.filter(
            "attribute", 
            "n_effective", 
            "min", 
            n_threshold)

        ## filtering for length of repeat units
        denovo_list = denovo_list.filter(
            "attribute", 
            "l_effective", 
            "max", 
            l_threshold)

        ##########################################################################
        #########  Building HMM with hmmbuild

        # # De novo TRs were remastered with HMM
        denovo_hmm = [hmm.HMM.create(input_format = 'repeat', repeat=iTR) for iTR in denovo_list.repeats] # only possible with hmmbuild
        denovo_list_remastered = seq.detect(lHMM=denovo_hmm)

        ##########################################################################
        ######### Clustering

        # De novo TRs were clustered for overlap (common ancestry). Only best =
        # lowest p-Value and lowest divergence were retained.
        denovo_list_remastered = denovo_list.filter(
            "none_overlapping", ["common_ancestry"], {
                "pvalue": score, "divergence": score})

        ##########################################################################
        ######### Save Tandem Repeats

        denovo_list_remastered.write(output_format = "pickle", file = output_pickle_file)
        denovo_list_remastered.write(output_format = "tsv", file = output_tsv_file)

        # function to save as fasta has to be integrated

    print("\n***", seq_name, "***")
    print("denovo repeats:", len(denovo_list.repeats))
    print("repeats after filtering and clustering:", len(denovo_list_remastered.repeats))

    for i in range(len(denovo_list_remastered.repeats)):
        print(denovo_list_remastered.repeats[i])
