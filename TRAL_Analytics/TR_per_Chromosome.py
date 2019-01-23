#!/home/lina/envs/tral_env/bin/python3.7

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

# working_dir = "/home/naefpau1/proteome_reference/pickles"
# sequences_path = "/home/naefpau1/proteome_reference"
# output_path = "/home/naefpau1/output

# Thresholds for filtering
pvalue_threshold = 0.05
divergence_threshold = 0.1
n_threshold = 2.5 # minimun repeat unit count
l_threshold = 3 # maximum repeat unit length

# set_file = "AUP000005640_Mitochondrion.fasta"
# set_file = "AUP000005640_Chromosome22.fasta"

if len(sys.argv) > 1:
    chr_nr = sys.argv[1] # check commmand line argument
else:
    sys.exit("Please provide the chromosome number as input.")

set_file = "AUP000005640_Chromosome" + chr_nr + ".fasta"
print(set_file)

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

all_denovo_repeats = 0
all_filtered_repeats = 0

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
        for TR in denovo_list.repeats:
            TR.calculate_pvalues()
        # Saving this sequences as binary files:
        with open(TRs_pkl, 'wb') as f: 
            pickle.dump(denovo_list, f)

    ##########################################################################
    ######### Filtering TRs

    all_denovo_repeats += len(denovo_list.repeats) # add number of denovo found repeats

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

    all_filtered_repeats += len(denovo_list_remastered.repeats) # add number of clustered repeats
    print("\n***", seq_name, "***")
    print("denovo repeats:", len(denovo_list.repeats))
    print("repeats after filtering and clustering:", len(denovo_list_remastered.repeats))

    for i in range(len(denovo_list_remastered.repeats)):
        print(denovo_list_remastered.repeats[i])

print("\nThere where {} repeats found de novo.".format(all_denovo_repeats))
print("After filtering and clustering there where only {} repeats left.\n".format(all_filtered_repeats))
