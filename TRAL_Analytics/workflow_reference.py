##########################################################################
### Importing required modules
##########################################################################

import os

# modules in this folder
import get_sequences
import find_TRs_in_genes

working_directory = "/home/lina/SynologyDrive/TRAL_Masterthesis/Uniprot_data/pickles"
sequences_path = "/home/lina/SynologyDrive/TRAL_Masterthesis/Uniprot_data"

##########################################################################
### Defining paths and variables
##########################################################################

genes = ["TGFBR2"]
data_path = "/home/lina/SynologyDrive/TRAL_Masterthesis/Uniprot_data"
sequences = get_sequences.sequences_per_gene(genes, sequences_path, show_time=False, as_pickle=False, patient=False)
detectors = ["XSTREAM","HHrepID","T-REKS","TRUST"] # AA compatible detectors
gene = "TGFBR2"



def get_TRs_per_detector(detectors, sequences, gene):
    TRs_per_detector = {}
    for detector in detectors:
        TRs = find_TRs_in_genes.get_TRs(detector, sequences[gene])
        TRs_per_detector.update({detector:TRs})
    return TRs_per_detector

TRs_per_detector = get_TRs_per_detector(detectors, sequences, gene)

for detector, TRs in TRs_per_detector.items():
    print("\n",detector,"\n")
    find_TRs_in_genes.loop_list_TR(TRs)


# TRs_XSTREAM = find_TRs_in_genes.get_TRs("XSTREAM", sequences["TGFBR2"])
# TRs_HHrepID = find_TRs_in_genes.get_TRs("HHrepID", sequences["TGFBR2"])
# TRs_TREKS = find_TRs_in_genes.get_TRs("T-REKS", sequences["TGFBR2"])
# TRs_TRUST = find_TRs_in_genes.get_TRs("TRUST", sequences["TGFBR2"])
# find_TRs_in_genes.loop_list_TR(TRs_XSTREAM)
# find_TRs_in_genes.loop_list_TR(TRs_HHrepID)
# find_TRs_in_genes.loop_list_TR(TRs_TREKS)
# find_TRs_in_genes.loop_list_TR(TRs_TRUST)