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

from tral.sequence import sequence
from tral.hmm import hmm

logging.config.fileConfig(config_file("logging.ini"))
log = logging.getLogger('root')

# ## AA reference
working_dir = "/home/lina/Desktop/TRAL_Masterthesis/references/Uniprot_data/pickles"
sequences_path = "/home/lina/Desktop/TRAL_Masterthesis/references/Uniprot_data"

## DNA reference
#working_directory = "/home/lina/Desktop/TRAL_Masterthesis/references/NCBI/pickles"
#sequences_path = "/home/lina/Desktop/TRAL_Masterthesis/references/NCBI"

genes = ["APC"]
gene = "APC"

##########################################################################
######### Getting sequences

####### Set of proteins
# to get the sequences from one file (different proteins with one sequence per file)
# l_sequence = Fasta(sequences_file)
# for iS_pyfaidx in l_sequence:
#     iS = sequence.Sequence(seq=str(iS_pyfaidx), name=iS_pyfaidx.name.split("|")[1]) # name is protein identifier

# with this method we obtain the sequences differently

####### Some single proteins
# to get the sequences from different files (file per protein)
# sequences = get_sequences.sequences_per_gene(genes, sequences_path, show_time=False, as_pickle=False, patient=False)
# Saving this sequences as binary files:
# with open(working_dir + '/sequences.pkl', 'wb') as f: 
#     pickle.dump(sequences, f)

# Getting back the sequences:
with open(working_dir + '/sequences.pkl','rb') as f: # files saved before  
    sequences = pickle.load(f)

##########################################################################
######### Getting TRs

# detectors_AA = ["XSTREAM","HHrepID","T-REKS","TRUST"] # AA compatible detectors
# detectors_DNA = ["T-REKS", "TRF", "XSTREAM"] # DNA compatible detectors (without TRED and Phobos)

test_seq = sequences["APC"][0]
score = "phylo_gap001"

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

# De novo TRs with dTR_units (divergence) > 0.8; n_effective < 2.5; l < 10 or
# pvalue "phylo_gap01_ignore_trailing_gaps_and_coherent_deletions" > 0.01
# are discarded.
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


# Write result set of Pfam TRs
#test_entire_set.write(format = "tsv,...")





# # for TR in sequence_test.detect(denovo=True).repeats:
# #     TR.calculate_pvalues()

# # for i in range(len(sequence_test.detect(denovo=True).repeats)):
# # 	sequence_test.detect(denovo=True).repeats[i].calculate_pvalues()

# # print(sequence.detect(denovo=True).repeats[0])

# def get_TRs_per_detector(detectors, sequences, gene):
#     TRs_per_detector = {}
#     for detector in detectors:
#         TRs = find_TRs_in_genes.get_TRs(detector, sequences[gene])
#         TRs_per_detector.update({detector:TRs[0]}) # because there is only one sequence i just put [0]
#     return TRs_per_detector

# TRs_per_detector = get_TRs_per_detector(detectors_DNA, sequences, gene)

# # for detector, TRs in TRs_per_detector.items():
# #     print("\n",detector,"\n")
# #     find_TRs_in_genes.loop_list_TR(TRs)

# # for gene, sequence_list in sequences.items():
# #     for sequence in sequence_list:
# #         sequence.get_repeatlist("denovo").cluster("shared_char")
    

# # for detector, TR_list_list in TRs_per_detector.items():
# #     print("\n",detector,"\n")
# #     for TR_list in TR_list_list:
# #         # print(TR_list)
# #         TR_list.cluster("shared_char")


# # TRs_per_detector["T-REKS"].cluster("shared_char")

# type(TRs_per_detector["T-REKS"].repeats[1])



# for i in range(len(TRs_per_detector["T-REKS"].repeats)):
# 	TRs_per_detector["T-REKS"].repeats[i].calculate_pvalues()

# for i in range(len(TRs_per_detector["T-REKS"].repeats)):
# 	print(TRs_per_detector["T-REKS"].repeats[i])


# # # filter for n_effective
# # filtered_TRs = TRs_per_detector["T-REKS"].filter("attribute", "n_effective", "min", 3.5)

# # print(len(TRs_per_detector["T-REKS"].repeats))
# # print(len(filtered_TRs.repeats))

# # for i in range(len(filtered_TRs.repeats)):
# # 	print(filtered_TRs.repeats[i])


# # # filter for l_effective (max 3.5)
# # filtered_TRs = TRs_per_detector["T-REKS"].filter("attribute", "l_effective", "max", 3.5)

# # print(len(TRs_per_detector["T-REKS"].repeats))
# # print(len(filtered_TRs.repeats))

# # for i in range(len(filtered_TRs.repeats)):
# # 	print(filtered_TRs.repeats[i])

# # # filter for l_msa (max 3.5)
# # filtered_TRs = TRs_per_detector["T-REKS"].filter("attribute", "l_effective", "max", 3.5)

# # print(len(TRs_per_detector["T-REKS"].repeats))
# # print(len(filtered_TRs.repeats))

# # for i in range(len(filtered_TRs.repeats)):
# # 	print(filtered_TRs.repeats[i])


# # filter for pvalue (min 0.05)
# filtered_TRs = TRs_per_detector["T-REKS"].filter(func_name = "pvalue", score = "phylo_gap001", threshold = 0.05)

# print(len(TRs_per_detector["T-REKS"].repeats))
# print(len(filtered_TRs.repeats))

# for i in range(len(filtered_TRs.repeats)):
# 	print(filtered_TRs.repeats[i])