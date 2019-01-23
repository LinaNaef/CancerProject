##########################################################################
### Importing required modules
##########################################################################

import sys
import os
import pickle

from pyfaidx import Fasta

import logging
import logging.config

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
working_dir = "/home/lina/Desktop/TRAL_Masterthesis/references/Uniprot_data/proteom_reference/pickles"
sequences_path = "/home/lina/Desktop/TRAL_Masterthesis/references/Uniprot_data/proteom_reference"
output_path = "/home/lina/SynologyDrive/TRAL_Masterthesis/TRAL_Pipeline_Analytics/test_output/proteom"

# set_file = "AUP000005640_Mitochondrion.fasta"
set_file = "AUP000005640_Chromosome21.fasta"
set_name = set_file.split(".fasta")[0]
sequences_file = os.path.join(sequences_path, set_file)
result_dir = os.path.join(output_path, set_name)

##########################################################################
######### Getting sequences

## obtaining sequences from fasta format
proteins = Fasta(sequences_file)
# print(proteins.keys())
# proteins['sp|P03886|NU1M_HUMAN'][:].seq

number_proteins = 0
number_TR_proteins = 0
number_TRs = 0

for pyfaidx in proteins:
    number_proteins += 1
    seq_name = pyfaidx.name.split("|")[1]

    output_pickle_file = os.path.join(result_dir, seq_name + ".pkl")
    output_tsv_file = os.path.join(result_dir, seq_name + ".tsv")

    if os.path.exists(output_pickle_file) and os.path.exists(output_tsv_file):
        with open(output_pickle_file,'rb') as f: 
            denovo_list_remastered = pickle.load(f)

    if len(denovo_list_remastered.repeats) > 0:
        number_TR_proteins += 1
        number_TRs += len(denovo_list_remastered.repeats)


print("\nOf {} proteins, {} contain a total of {} repeats.".format(number_proteins, number_TR_proteins, number_TRs))
print("That is {} % of all proteins examined in this run.".format(round((number_TR_proteins / number_proteins * 100), 2)))
print("\nFiltering criteria has been: \
        \n --> pvalue < 0.05 \
        \n --> divergence < 0.1 \
        \n --> minimun repeat unit count: 2.5 \
        \n --> maximum repeat unit length: 3\n")
