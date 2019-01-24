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
output_path = "/home/lina/SynologyDrive/TRAL_Masterthesis/TRAL_Pipeline_Analytics/test_output/proteom"

if len(sys.argv) > 1:
    chr_nr = sys.argv[1] # check commmand line argument
else:
    chr_nr = "2"
    # sys.exit("Please provide the chromosome number as input.")

set_file = "AUP000005640_Chromosome" + chr_nr + ".fasta"
chr_name = "Chromosome " + chr_nr
print(set_file)

set_name = set_file.split(".fasta")[0]
sequences_file = os.path.join(sequences_path, set_file)
result_dir = os.path.join(output_path, set_name)

all_AA = ""

##########################################################################
######### Getting sequences

## obtaining sequences from fasta format
proteins = Fasta(sequences_file)
# print(proteins.keys())
# proteins['sp|P03886|NU1M_HUMAN'][:].seq

number_proteins = 0
number_TR_proteins = 0
number_TRs = 0

max_length = 41
l_n = np.zeros((4,max_length))
print(l_n)

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
        for TR in range(len(denovo_list_remastered.repeats)):
            # print(denovo_list_remastered.repeats[TR])
            l = denovo_list_remastered.repeats[TR].l_effective
            n = denovo_list_remastered.repeats[TR].n
            all_AA += denovo_list_remastered.repeats[TR].textD_standard_aa # string with all AA
            # print(l,n)
            l_n[l,n] += 1  

############ Create Dictionary with all AA appended ###############

# Create dict which counts AAs
AA_dict = (Counter(all_AA))
# remove whitespace and gap
AA_dict.pop(' ', None)
AA_dict.pop('-', None)

# add missing AAs
AAs = "ACDEFGHIKLMNPQRSTVWY"
for AA in AAs:
    if AA not in AA_dict:
        AA_dict[AA] = 0

# sort by count
sorted_AA_count = sorted(AA_dict.items(), key=lambda x: x[1], reverse=True)
AA_names = list(zip(*sorted_AA_count))[0]
AA_count = list(zip(*sorted_AA_count))[1]
normed_AA_count = [float(i)/sum(AA_count) for i in AA_count]

# Plot
plt.bar(AA_names, normed_AA_count, color='g')
plt.title("Relative apearance of amino acids in tandem repeats in " + chr_name)
plt.xlabel('Amino Acid')
plt.ylabel('Relative appearance in TRs')
plt.show()
plt.close()

############ Lenght/Unit Count Distribution ###############

# Three subplots sharing both x/y axes
# f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True)
# plt.imshow(l_n[1,], cmap='hot')


# print(l_n[1,])
# # shiftig to have a better visualisation
# l_n_roll = np.roll(l_n,1) # shift to the right 
# l_n_roll2 = np.roll(l_n_roll,1,0) # shift down
# print(l_n_roll2)

# plt.xlabel("Unit Number")
# plt.ylabel("Unit Lenght")
# plt.xlim(1,max_length-1)
# # plt.ylim(1,3)
# plt.imshow(l_n_roll2, cmap='hot')
# plt.show()
# plt.close()

########### Overview ###############

print("\nOf {} proteins, {} contain a total of {} repeats.".format(number_proteins, number_TR_proteins, number_TRs))
print("That is {} % of all proteins examined in this run.".format(round((number_TR_proteins / number_proteins * 100), 2)))
print("\nFiltering criteria has been: \
        \n --> pvalue < 0.05 \
        \n --> divergence < 0.1 \
        \n --> minimun repeat unit count: 2.5 \
        \n --> maximum repeat unit length: 3\n")
