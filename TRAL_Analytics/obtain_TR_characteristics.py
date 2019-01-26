#!/home/lina/envs/tral_env/bin/python3.7

##########################################################################
### Importing required modules
##########################################################################

import sys
import os
import pickle
import numpy as np
from numpy.core import multiarray

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

import analyzing_functions

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
######### Defining Paths and Parameters, Initialize lists etc

## AA reference
working_dir = "/home/lina/Desktop/TRAL_Masterthesis/references/Uniprot_data/proteom_reference/pickles"
sequences_path = "/home/lina/Desktop/TRAL_Masterthesis/references/Uniprot_data/proteom_reference"
# https://www.uniprot.org/proteomes/UP000005640
# output_path = "/home/lina/SynologyDrive/TRAL_Masterthesis/TRAL_Pipeline_Analytics/test_output/proteom"
output_path = "/home/lina/Desktop/TRAL_Masterthesis/TRAL_Pipeline_Analytics/output_all_TRs_Cluster"
output_statistics = "/home/lina/SynologyDrive/TRAL_Masterthesis/TRAL_Pipeline_Analytics/output_statistics"

list_proteins_with_TR = [] # initialize list to get all proteins containing TRs
all_AA = "" # initialize string for all AA in TRs
count_dist = {1:{},2:{},3:{}} # initialize dictionaries for length count

number_proteins = 0
number_TR_proteins = 0
number_TRs = 0

set_name = "proteome"
chr_nr = "whole proteome"
chr_name = "whole_proteome"

all_chr = [str(i) for i in range(1,23)] + ["X","Y"] # defining each chromosome

pvalue_above_zero = 0
divergence_above_zero = 0
max_TR = 0
list_long_TR_25 = []
list_long_TR_50 = []
EOFError_Count = 0

# from https://www.proteinatlas.org/search/protein_class:COSMIC+somatic+mutations+in+cancer+genes
from cancer_analysis.import_tsv import cancer_proteome
cancer_TRs = {}
list_colorect_cancer_proteins_with_TRs = []

# from https://www.proteinatlas.org/search/cancer_specificity_rna%3Acolorectal+cancer%3Belevated+AND+sort_by%3Atissue+specific+score+AND+show_columns%3Acancerrnacategory
# elevated mRNA in colorectal cancer tissue
from cancer_analysis.import_tsv import colorec_canc_proteome
colorec_cancer_TRs = {}
list_cancer_proteins_with_TRs = []

for chr_nr in all_chr: # analyze each chromosome
    set_file = "AUP000005640_Chromosome" + chr_nr + ".fasta"
    set_name = set_file.split(".fasta")[0]
    sequences_file = os.path.join(sequences_path, set_file)
    result_dir = os.path.join(output_path, set_name)

    ##########################################################################
    ######### Getting sequences

    ## obtaining sequences from fasta format
    proteins = Fasta(sequences_file)
    # print(proteins.keys())
    # proteins['sp|P03886|NU1M_HUMAN'][:].seq

    for pyfaidx in proteins:
        number_proteins += 1
        UniqueIdentifier = pyfaidx.name.split("|")[1]
        EntryName = pyfaidx.name.split("|")[2]
        try:
            GeneName = pyfaidx.long_name.split("GN=")[1].split(" ")[0]
        except IndexError:
            GeneName = "NoGeneName"
        output_pickle_file = os.path.join(result_dir, UniqueIdentifier + ".pkl")
        output_tsv_file = os.path.join(result_dir, UniqueIdentifier + ".tsv")

        # open remastered TR file
        if os.path.exists(output_pickle_file) and os.path.exists(output_tsv_file):
            try:
                with open(output_pickle_file,'rb') as f: 
                    denovo_list_remastered = pickle.load(f)
            except EOFError:
                EOFError_Count += 1
                pass            

        # analyze all proteins that contain TRs
        if len(denovo_list_remastered.repeats) > 0:
            protein_tuple = UniqueIdentifier, EntryName, GeneName, len(denovo_list_remastered.repeats)
            list_proteins_with_TR.append(protein_tuple)
            number_TR_proteins += 1
            number_TRs += len(denovo_list_remastered.repeats)

            # TR in cancer related gene?
            for gene in cancer_proteome:
                if GeneName == gene:
                    if gene not in cancer_TRs:
                        cancer_TRs[gene] = denovo_list_remastered.repeats
                        list_cancer_proteins_with_TRs.append(protein_tuple)
                        break

            # TR in colorectal cancer related gene?
            for gene in colorec_canc_proteome:
                if GeneName == gene:
                    if gene not in colorec_cancer_TRs:
                        colorec_cancer_TRs[gene] = denovo_list_remastered.repeats
                        list_colorect_cancer_proteins_with_TRs.append(protein_tuple)
                        break

            for TR in range(len(denovo_list_remastered.repeats)):
                
                # some length statistics
                l = denovo_list_remastered.repeats[TR].l_effective
                n = denovo_list_remastered.repeats[TR].n
                if n in count_dist[l]:
                    count_dist[l][n] += 1
                else:
                    count_dist[l][n] = 1
                all_AA += denovo_list_remastered.repeats[TR].textD_standard_aa # string with all AA

                # just have a look if not all TRs have a pvalue of 0 and a very small divergence
                if denovo_list_remastered.repeats[TR].d_pvalue['phylo_gap01'] > 0:
                    pvalue_above_zero += 1
                if denovo_list_remastered.repeats[TR].d_divergence['phylo_gap01'] > 1e-10:
                    divergence_above_zero += 1

                # find long repeats
                TR_length = denovo_list_remastered.repeats[TR].n * denovo_list_remastered.repeats[TR].l_effective
                if max_TR < TR_length:
                    max_TR = TR_length
                    longest_TR, longest_TR_protein, max_TR = denovo_list_remastered.repeats[TR], pyfaidx, TR_length

                if TR_length > 50:
                    list_long_TR_50.append((denovo_list_remastered.repeats[TR], pyfaidx))

                if TR_length > 25:
                    list_long_TR_25.append((denovo_list_remastered.repeats[TR], pyfaidx))

