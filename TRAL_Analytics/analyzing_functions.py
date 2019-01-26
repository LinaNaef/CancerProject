#########################################################################
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

##########################################################################
######### Create txt with Proteins that contain TRs

def TR_list_txt(list_proteins_with_TR,file_protein_list):
    f = open(file_protein_list, 'w')
    f.write('UniqueIdentifier EntryName GeneName NumberTRs\n')
    for protein in list_proteins_with_TR:
        f.write(' '.join(str(s) for s in protein) + '\n')
    f.close()

##########################################################################
######### Create Figure with AA frequency

def AA_frequency(all_AA, chr_name, output_statistics):
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
    # plt.title("Relative apearance of amino acids in tandem repeats in " + chr_name)
    plt.xlabel('Amino acid')
    plt.ylabel('Relative frequency in TRs')

    # Save Figure
    output_AA_frequency = os.path.join(output_statistics, "aa_frequency" , chr_name + ".png")

    plt.savefig(output_AA_frequency)
    plt.close()

##########################################################################
######### Lenght/Unit Count Distribution Figure

def l_n_distribution(count_dist, chr_name, output_statistics):
    # get max count to define x-axis
    max_count = 0
    for n in range(1,4):
        if max(count_dist[n]) > max_count:
            max_count = max(count_dist[n])

    # add missing lengts
    for n in range(1,4):
        for l in range(1,max_count+1):
            if l not in count_dist[n]:
                count_dist[n][l] = 0

    # sort values
    one_AA = list(sorted(count_dist[1].items()))
    two_AA = list(sorted(count_dist[2].items()))
    three_AA = list(sorted(count_dist[3].items()))

    length = list(zip(*one_AA))[0]
    one_AA_count = list(zip(*one_AA))[1]
    two_AA_count = list(zip(*two_AA))[1]
    three_AA_count = list(zip(*three_AA))[1]

    # Three subplots sharing both x/y axes
    f, (ax) = plt.subplots(3, sharey=True)

    # only include TRs that have a max n of 25, others are relatively few
    ax[0].set(ylabel='Count', xlabel='Number of repetitive units with 1 amino acid.', xlim = (1,25))
    ax[1].set(ylabel='Count', xlabel='Number of repetitive units with 2 amino acids.', xlim = (1,25))
    ax[2].set(ylabel='Count', xlabel='Number of repetitive units with 3 amino acids.', xlim = (1,25))

    ax[0].bar(length, one_AA_count, color='g')
    # ax[0].set_title('1 Amino Acid')
    ax[1].bar(length, two_AA_count, color='g')
    # ax[1].set_title('2 Amino Acids')
    ax[2].bar(length, three_AA_count, color='g')
    # ax[2].set_title('3 Amino Acids')
    f.subplots_adjust(hspace=0.5)

    # Save Figure
    output_l_n_distribution = os.path.join(output_statistics, "length_unit_distribution" , chr_name + ".png")
    
    plt.savefig(output_l_n_distribution)
    plt.close()

##########################################################################
######### Overview

def overview(chr_name, number_proteins, number_TR_proteins, number_TRs):
    print('\n******************',chr_name,'******************')
    print("\nOf {} proteins, {} contain a total of {} repeats.".format(number_proteins, number_TR_proteins, number_TRs))
    print("That is {} % of all proteins in {}.".format(round((number_TR_proteins / number_proteins * 100), 2), chr_name))
    print("\nFiltering criteria has been: \
            \n --> pvalue < 0.05 \
            \n --> divergence < 0.1 \
            \n --> minimun repeat unit count: 2.5 \
            \n --> maximum repeat unit length: 3\n")
    print('************************************************')
