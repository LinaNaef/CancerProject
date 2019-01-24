#!/home/lina/envs/tral_env/bin/python3.7

##########################################################################
### Importing required modules
##########################################################################

import sys
import os
import pickle
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
from collections import Counter, OrderedDict

from pyfaidx import Fasta

import logging
import logging.config

from tral.paths import config_file, PACKAGE_DIRECTORY
from tral import configuration
from tral.sequence import sequence
from tral.hmm import hmm

import analyzing_functions

logging.config.fileConfig(config_file("logging.ini"))
log = logging.getLogger('root')

# define this as in the config
CONFIG_GENERAL = configuration.Configuration.instance().config
CONFIG = CONFIG_GENERAL["repeat_list"]
score = CONFIG["model"]

# Create function to count AAs
def AACount(text):
    return Counter(c for c in text if c.isalpha())

##########################################################################
######### Defining Paths and Parameters

## AA reference
working_dir = "/home/lina/Desktop/TRAL_Masterthesis/references/Uniprot_data/proteom_reference/pickles"
sequences_path = "/home/lina/Desktop/TRAL_Masterthesis/references/Uniprot_data/proteom_reference"
# https://www.uniprot.org/proteomes/UP000005640
output_path = "/home/lina/SynologyDrive/TRAL_Masterthesis/TRAL_Pipeline_Analytics/test_output/proteom"
output_figures = "/home/lina/SynologyDrive/TRAL_Masterthesis/TRAL_Pipeline_Analytics/output_figures"

list_proteins_with_TR = [] # initialize list to get all proteins containing TRs
all_AA = "" # initialize string for all AA in TRs
count_dist = {1:{},2:{},3:{}} # initialize dictionaries for length count
max_unit_count = 20 # expected max count for visualization of unit length, unit count distribution

number_proteins = 0
number_TR_proteins = 0
number_TRs = 0

set_name = "Proteome"

all_chr = [str(i) for i in range(1,23)] + ["X","Y"]
file_protein_list = os.path.join(output_path, set_name + '.txt')

for chr_nr in all_chr: # analyze each chromosome
    set_file = "AUP000005640_Chromosome" + chr_nr + ".fasta"
    set_name = set_file.split(".fasta")[0]
    sequences_file = os.path.join(sequences_path, set_file)
    result_dir = os.path.join(output_path, set_name)

    ##########################################################################
    ######### Getting sequences

    ## obtaining sequences from fasta format
    proteins = Fasta(sequences_file)
    # print(proteins.keys())
    # proteins['sp|P03886|NU1M_HUMAN'][:].seq

    for pyfaidx in proteins:
        number_proteins += 1
        UniqueIdentifier = pyfaidx.name.split("|")[1]
        EntryName = pyfaidx.name.split("|")[2]
        output_pickle_file = os.path.join(result_dir, UniqueIdentifier + ".pkl")
        output_tsv_file = os.path.join(result_dir, UniqueIdentifier + ".tsv")

        # open remastered TR file
        if os.path.exists(output_pickle_file) and os.path.exists(output_tsv_file):
            with open(output_pickle_file,'rb') as f: 
                denovo_list_remastered = pickle.load(f)

        # analyze all proteins that contain TRs
        if len(denovo_list_remastered.repeats) > 0:
            protein_tuple = UniqueIdentifier, EntryName, len(denovo_list_remastered.repeats)
            list_proteins_with_TR.append(protein_tuple)
            number_TR_proteins += 1
            number_TRs += len(denovo_list_remastered.repeats)

            for TR in range(len(denovo_list_remastered.repeats)):

                l = denovo_list_remastered.repeats[TR].l_effective
                n = denovo_list_remastered.repeats[TR].n

                if n in count_dist[l]:
                    count_dist[l][n] += 1
                else:
                    count_dist[l][n] = 1

                all_AA += denovo_list_remastered.repeats[TR].textD_standard_aa # string with all AA

chr_nr = "whole proteome"
chr_name = "whole proteome"

# Calculate Figures and 
analyzing_functions.TR_list_txt(list_proteins_with_TR, file_protein_list)
analyzing_functions.AA_frequency(all_AA, chr_name, output_figures)
analyzing_functions.l_n_distribution(count_dist, chr_name, output_figures)
analyzing_functions.overview(chr_name, number_proteins, number_TR_proteins, number_TRs)
