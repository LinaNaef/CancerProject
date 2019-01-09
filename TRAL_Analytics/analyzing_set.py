import sys

##########################################################################
### Importing required modules
##########################################################################

import os
import pickle

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

    # print("\n***", seq_name, "***")
    # print("repeats after filtering and clustering:", len(denovo_list_remastered.repeats))

print("Of {} proteins, {} contain a total of {} repeats.".format(number_proteins, number_TR_proteins, number_TRs))