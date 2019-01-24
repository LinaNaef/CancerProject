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
set_name = set_file.split(".fasta")[0]
sequences_file = os.path.join(sequences_path, set_file)
result_dir = os.path.join(output_path, set_name)
file_protein_list = os.path.join(output_path, set_name + '.txt')


list_proteins_with_TR = []
all_AA = ""
count_dist = {1:{},2:{},3:{}}
max_unit_count = 20 # expected max count for visualization of unit length, unit count distribution


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

for pyfaidx in proteins:
    number_proteins += 1
    UniqueIdentifier = pyfaidx.name.split("|")[1]
    EntryName = pyfaidx.name.split("|")[2]
    output_pickle_file = os.path.join(result_dir, UniqueIdentifier + ".pkl")
    output_tsv_file = os.path.join(result_dir, UniqueIdentifier + ".tsv")

    if os.path.exists(output_pickle_file) and os.path.exists(output_tsv_file):
        with open(output_pickle_file,'rb') as f: 
            denovo_list_remastered = pickle.load(f)

    # all proteins that contain TRs
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

############ Create txt with Proteins that contain TRs ###############

f = open(file_protein_list, 'w')
f.write('UniqueIdentifier EntryName NumberTRs\n')
for protein in list_proteins_with_TR:
    f.write(' '.join(str(s) for s in protein) + '\n')
f.close()

############ Create Dictionary with all AA appended ###############

# # Create dict which counts AAs
# AA_dict = (Counter(all_AA))
# # remove whitespace and gap
# AA_dict.pop(' ', None)
# AA_dict.pop('-', None)

# # add missing AAs
# AAs = "ACDEFGHIKLMNPQRSTVWY"
# for AA in AAs:
#     if AA not in AA_dict:
#         AA_dict[AA] = 0

# # sort by count
# sorted_AA_count = sorted(AA_dict.items(), key=lambda x: x[1], reverse=True)
# AA_names = list(zip(*sorted_AA_count))[0]
# AA_count = list(zip(*sorted_AA_count))[1]
# normed_AA_count = [float(i)/sum(AA_count) for i in AA_count]

# # Plot
# plt.bar(AA_names, normed_AA_count, color='g')
# # plt.title("Relative apearance of amino acids in tandem repeats in " + chr_name)
# plt.xlabel('Amino Acid')
# plt.ylabel('Relative appearance in TRs')
# plt.show()
# plt.close()

############ Lenght/Unit Count Distribution ###############

# # get max count to define x-axis
# max_count = 0
# for n in range(1,4):
#     if max(count_dist[n]) > max_count:
#         max_count = max(count_dist[n])

# # add missing lengts
# for n in range(1,4):
#     for l in range(1,max_count+1):
#         if l not in count_dist[n]:
#             count_dist[n][l] = 0

# # sort values
# one_AA = list(sorted(count_dist[1].items()))
# two_AA = list(sorted(count_dist[2].items()))
# three_AA = list(sorted(count_dist[3].items()))

# length = list(zip(*one_AA))[0]
# one_AA_count = list(zip(*one_AA))[1]
# two_AA_count = list(zip(*two_AA))[1]
# three_AA_count = list(zip(*three_AA))[1]

# # Three subplots sharing both x/y axes
# f, (ax) = plt.subplots(3, sharey=True)

# ax[0].set(ylabel='Count', xlabel='Number of repetitive units with 1 Amino Acid.')
# ax[1].set(ylabel='Count', xlabel='Number of repetitive units with 2 Amino Acids.')
# ax[2].set(ylabel='Count', xlabel='Number of repetitive units with 3 Amino Acids.')

# ax[0].bar(length, one_AA_count, color='g')
# # ax[0].set_title('1 Amino Acid')
# ax[1].bar(length, two_AA_count, color='g')
# # ax[1].set_title('2 Amino Acids')
# ax[2].bar(length, three_AA_count, color='g')
# # ax[2].set_title('3 Amino Acids')
# f.subplots_adjust(hspace=0.5)

# plt.show()

########### Overview ###############

print('******************',chr_name,'******************')
print("\nOf {} proteins, {} contain a total of {} repeats.".format(number_proteins, number_TR_proteins, number_TRs))
print("That is {} % of all proteins in Chromosome {}.".format(round((number_TR_proteins / number_proteins * 100), 2), chr_nr))
print("\nFiltering criteria has been: \
        \n --> pvalue < 0.05 \
        \n --> divergence < 0.1 \
        \n --> minimun repeat unit count: 2.5 \
        \n --> maximum repeat unit length: 3\n")
print('************************************************')
