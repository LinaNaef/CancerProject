import sys

##########################################################################
### Importing required modules
##########################################################################

import os
import pickle

# # modules in this folder
import get_sequences
import find_TRs_in_genes
import pyfaidx


import logging
import logging.config
import os

from tral.paths import config_file, PACKAGE_DIRECTORY
from tral import configuration

from tral.sequence import sequence
from tral.hmm import hmm

logging.config.fileConfig(config_file("logging.ini"))
log = logging.getLogger('root')

seq_type = "AA" # change sequence type in config file as well!!

if seq_type == "AA":
    ## AA reference
    working_dir = "/home/lina/Desktop/TRAL_Masterthesis/references/Uniprot_data/pickles"
    sequences_path = "/home/lina/Desktop/TRAL_Masterthesis/references/Uniprot_data"
    output_path = "/home/lina/SynologyDrive/TRAL_Masterthesis/TRAL_Pipeline_Analytics/test_output/AA"
    detectors = ["XSTREAM","HHrepID","T-REKS","TRUST"] # AA compatible detectors

    # Thresholds for filtering
    pvalue_threshold = 0.05
    divergence_threshold = 0.8
    n_threshold = 2.5 # minimun repeat unit count
    l_threshold = 3 # maximum repeat unit length

elif seq_type == "DNA":
    # Since TRAL is not ready to calculate the pvalue for DNA, I cannot longer focus on this part

    ## DNA reference
    working_dir = "/home/lina/Desktop/TRAL_Masterthesis/references/NCBI/pickles"
    sequences_path = "/home/lina/Desktop/TRAL_Masterthesis/references/NCBI"
    output_path = "/home/lina/SynologyDrive/TRAL_Masterthesis/TRAL_Pipeline_Analytics/test_output/DNA"
    detectors = ["T-REKS", "TRF", "XSTREAM"] # DNA compatible detectors (without TRED and Phobos)

    # Thresholds for filtering
    pvalue_threshold = 0.01
    divergence_threshold = 0.8
    n_threshold = 2.5 # minimun repeat unit count
    l_threshold = 5 # maximum repeat unit length

else:
    print("This sequence type does not exist.")

genes = ["APC", "CDKL1", "TGFBR2", "TP53I3", "TP53I11", "HTT"]
# genes = ["TGFBR2"]

# define this as in the config
CONFIG_GENERAL = configuration.Configuration.instance().config
CONFIG = CONFIG_GENERAL["repeat_list"]
score = CONFIG["model"]

##########################################################################
######### Getting sequences

####### Set of proteins (maybe better version)
# to get the sequences from one file (different proteins with one sequence per file)
# l_sequence = Fasta(sequences_file)
# for iS_pyfaidx in l_sequence:
#     iS = sequence.Sequence(seq=str(iS_pyfaidx), name=iS_pyfaidx.name.split("|")[1]) # name is protein identifier

# with this method we obtain the sequences differently

##########################################################################
######### used for only single sequences per fasta file (one protein per file!)
for gene in genes:
    ####### Some single proteins (own function to get sequence used)
    # to get the sequences from different files (file per protein)

    sequence_pkl = os.path.join(working_dir, gene + "_sequences.pkl")

    if os.path.exists(sequence_pkl):
        # Getting back the sequences:
        with open(sequence_pkl,'rb') as f: # files saved before  
            sequence = pickle.load(f)
    else:
        sequence = get_sequences.get_sequences(sequences_path, gene)[0]
        # Saving this sequences as binary files:
        with open(sequence_pkl, 'wb') as f: 
            pickle.dump(sequence, f)

    ##########################################################################
    ######### Getting TRs

    TRs_pkl = os.path.join(working_dir, "TRs_raw", gene + "_TRs.pkl")

    if os.path.exists(TRs_pkl):
        # Getting back the sequences:
        with open(TRs_pkl,'rb') as f: # files saved before  
            denovo_list = pickle.load(f)
    else:
        denovo_list = sequence.detect(denovo=True)
        # Saving this sequences as binary files:
        with open(TRs_pkl, 'wb') as f: 
            pickle.dump(denovo_list, f)

    # print("Found", len(denovo_list.repeats), "denovo repeats in", gene)

    ##########################################################################
    ######### Filtering TRs

    output_pickle_file = os.path.join(output_path, gene + ".pkl")
    output_tsv_file = os.path.join(output_path, gene + ".tsv")

    if os.path.exists(output_pickle_file) and os.path.exists(output_tsv_file):
        with open(output_pickle_file,'rb') as f: 
            denovo_list_remastered = pickle.load(f)
    else:
        # filtering for pvalue
        denovo_list = denovo_list.filter(
            "pvalue",
            score,
            pvalue_threshold)
        # print("Repeats after filtering for pvalue 0.01:", len(denovo_list.repeats))

        ## filtering for divergence
        denovo_list = denovo_list.filter(
            "divergence",
            score,
            divergence_threshold)
        # print("Repeats after filtering for divergence 0.8:", len(denovo_list.repeats))

        ## filtering for number of repeat units
        denovo_list = denovo_list.filter(
            "attribute", 
            "n_effective", 
            "min", 
            n_threshold)
        # print("Repeats after filtering for min 2.5 repeat units:", len(denovo_list.repeats))

        ## filtering for length of repeat units
        denovo_list = denovo_list.filter(
            "attribute", 
            "l_effective", 
            "max", 
            l_threshold)
        # print("Repeats after filtering for a length of at least 10:", len(denovo_list.repeats))

        ##########################################################################
        #########  Building HMM with hmmbuild

        # # De novo TRs were remastered with HMM
        denovo_hmm = [hmm.HMM.create(input_format = 'repeat', repeat=iTR) for iTR in denovo_list.repeats] # only possible with hmmbuild
        denovo_list_remastered = sequence.detect(lHMM=denovo_hmm)

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

    print("\n***", gene, "***")
    print("denovo repeats:", len(denovo_list.repeats))
    print("repeats after filtering and clustering:", len(denovo_list_remastered.repeats))

    for i in range(len(denovo_list_remastered.repeats)):
        print(denovo_list_remastered.repeats[i])
