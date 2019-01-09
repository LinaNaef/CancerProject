##########################################################################
### Importing required modules
##########################################################################

import os

# modules in this folder
import get_sequences
import find_TRs_in_genes

## AA reference
working_directory = "/home/lina/Desktop/TRAL_Masterthesis/references/Uniprot_data/pickles"
sequences_path = "/home/lina/Desktop/TRAL_Masterthesis/references/Uniprot_data"

# ## DNA reference
# working_directory = "/home/lina/Desktop/TRAL_Masterthesis/references/NCBI/pickles"
# sequences_path = "/home/lina/Desktop/TRAL_Masterthesis/references/NCBI"

##########################################################################
### Defining paths and variables
##########################################################################

# genes = ["TGFBR2"]
# gene = "TGFBR2"

genes = ["APC"]
gene = "APC"

sequences = get_sequences.sequences_per_gene(genes, sequences_path, show_time=False, as_pickle=False, patient=False)
detectors_AA = ["XSTREAM","HHrepID","T-REKS","TRUST"] # AA compatible detectors
detectors_DNA = ["T-REKS", "TRF", "XSTREAM"] # DNA compatible detectors (without TRED and Phobos)


def get_TRs_per_detector(detectors, sequences, gene):
    TRs_per_detector = {}
    for detector in detectors:
        TRs = find_TRs_in_genes.get_TRs(detector, sequences[gene])
        TRs_per_detector.update({detector:TRs})
    return TRs_per_detector

TRs_per_detector = get_TRs_per_detector(detectors_DNA, sequences, gene)

# for detector, TRs in TRs_per_detector.items():
#     print("\n",detector,"\n")
#     find_TRs_in_genes.loop_list_TR(TRs)

# for gene, sequence_list in sequences.items():
#     for sequence in sequence_list:
#         sequence.get_repeatlist("denovo").cluster("shared_char")
    

for detector, TR_list_list in TRs_per_detector.items():
    print("\n",detector,"\n")
    for TR_list in TR_list_list:
        # print(TR_list)
        TR_list.cluster("shared_char")


type(sequences)

len(TRs_per_detector["T-REKS"]) # list
len(TRs_per_detector["T-REKS"][0].repeats) # repeatlist

TRs_per_detector["T-REKS"][0].cluster("shared_char")

type(TRs_per_detector["T-REKS"][0].repeats[1])


# for i range(len(TRs_per_detector["T-REKS"])):
#     TRs_per_detector["T-REKS"][i].repeats      

# type(TRs_per_detector["T-REKS"])
# for i in TRs_per_detector["T-REKS"]:
# 	i.repeats

for i in range(len(TRs_per_detector["T-REKS"][0].repeats)):
	TRs_per_detector["T-REKS"][0].repeats[i].calculate_pvalues()
        