##########################################################################
### Importing required modules
##########################################################################

import os

# modules in this folder
import get_sequences
import find_TRs_in_genes

# ## AA reference
# working_directory = "/home/lina/Desktop/TRAL_Masterthesis/references/Uniprot_data/pickles"
# sequences_path = "/home/lina/Desktop/TRAL_Masterthesis/references/Uniprot_data"

## DNA reference
working_directory = "/home/lina/Desktop/TRAL_Masterthesis/references/NCBI/pickles"
sequences_path = "/home/lina/Desktop/TRAL_Masterthesis/references/NCBI"

##########################################################################
### Defining paths and variables
##########################################################################

genes = ["TGFBR2"]
gene = "TGFBR2"

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


## does not work till now: will try to perform statistical significance test

for detector, TRs in TRs_per_detector.items():
    print("\n",detector,"\n")

    for TR in TRs.repeats:
        TR.calculate_pvalues()


