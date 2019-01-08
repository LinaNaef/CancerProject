import sys
# sys.path.append('/home/lina/SynologyDrive/TRAL_Masterthesis/TRAL_Pipeline_Analytics/TRAL_Analytics') 

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

# ## AA reference
working_dir = "/home/lina/Desktop/TRAL_Masterthesis/references/Uniprot_data/pickles"
sequences_path = "/home/lina/Desktop/TRAL_Masterthesis/references/Uniprot_data"
output_path = "/home/lina/SynologyDrive/TRAL_Masterthesis/TRAL_Pipeline_Analytics/test_output"

## DNA reference
#working_directory = "/home/lina/Desktop/TRAL_Masterthesis/references/NCBI/pickles"
#sequences_path = "/home/lina/Desktop/TRAL_Masterthesis/references/NCBI"

# genes = ["APC", "CDKL1", "TGFBR2","TP53I3","TP53I11"]
genes = ["APC"]

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

    # sequences = get_sequences.sequences_per_gene(genes, sequences_path, show_time=False, as_pickle=False, patient=False)

    ##########################################################################
    ######### Getting TRs

    detectors_AA = ["XSTREAM","HHrepID","T-REKS","TRUST"] # AA compatible detectors
    detectors_DNA = ["T-REKS", "TRF", "XSTREAM"] # DNA compatible detectors (without TRED and Phobos)

    # de novo detection methods (Trust, T-reks, Xstream, HHrepID) are used to search the
    # INSERT OWN PARAMTERS USING: test_denovo_list = test_seq.detect(denovo =
    # True, **TEST_DENOVO_PARAMETERS)

    # test_denovo_list = test_seq.detect(denovo=True)
    # # Saving this TRs as binary files:
    # with open(working_dir + '/denovo_list.pkl', 'wb') as f: 
    #     pickle.dump(test_denovo_list, f)

    # Getting back the TRs:
    with open(working_dir + '/denovo_list.pkl','rb') as f: # files saved before  
        test_denovo_list = pickle.load(f)

    # When Trust is part of the detectors, the number of found repeats may
    # differ between runs...
    print("Found repeats:", len(test_denovo_list.repeats))

    ##########################################################################

    test_denovo_list = test_denovo_list.filter(
        "pvalue",
        score,
        0.01)
    print("Repeats after filtering for pvalue 0.01:", len(test_denovo_list.repeats))
    test_denovo_list = test_denovo_list.filter(
        "divergence",
        score,
        0.8)
    print("Repeats after filtering for divergence 0.8:", len(test_denovo_list.repeats))
    test_denovo_list = test_denovo_list.filter("attribute", "n_effective", "min", 2.5)
    print("Repeats after filtering for min 2.5 repeat units:", len(test_denovo_list.repeats))
    test_denovo_list = test_denovo_list.filter("attribute", "l_effective", "max", 3) # here an exception would be raise if the attribute is l
    print("Repeats after filtering for a length of at least 10:", len(test_denovo_list.repeats))
    # for i in range(len(test_denovo_list.repeats)):
    # 	print(test_denovo_list.repeats[i])


    ##########################################################################
    #########  Building HMM with hmmbuild

    # # De novo TRs were remastered with HMM
    test_denovo_hmm = [hmm.HMM.create(input_format = 'repeat', repeat=iTR) for iTR in test_denovo_list.repeats] # only possible with hmmbuild
    test_denovo_list_remastered = test_seq.detect(lHMM=test_denovo_hmm)

    ##########################################################################
    ######### Clustering


    # De novo TRs were clustered for overlap (common ancestry). Only best =
    # lowest p-Value and lowest divergence were retained.
    test_denovo_list = test_denovo_list.filter(
        "none_overlapping", ["common_ancestry"], {
            "pvalue": score, "divergence": score})
    print("Repeats after clustering:", len(test_denovo_list.repeats))


    ##########################################################################
    ######### Save Tandem Repeats


    output_pickle_file = os.path.join(output_path, gene + ".pkl")
    test_denovo_list.write(output_format = "pickle", file = output_pickle_file)


    output_tsv_file = os.path.join(output_path, gene + ".tsv")
    test_denovo_list.write(output_format = "tsv", file = output_tsv_file)

    # function to save as fasta has to be integrated