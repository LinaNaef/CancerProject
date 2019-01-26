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
output_path = "/home/lina/SynologyDrive/TRAL_Masterthesis/TRAL_Pipeline_Analytics/test_output/proteom"
# output_path = "/home/lina/Desktop/TRAL_Masterthesis/TRAL_Pipeline_Analytics/output_all_TRs_Cluster"
output_figures = "/home/lina/SynologyDrive/TRAL_Masterthesis/TRAL_Pipeline_Analytics/output_figures"

list_proteins_with_TR = [] # initialize list to get all proteins containing TRs
all_AA = "" # initialize string for all AA in TRs
count_dist = {1:{},2:{},3:{}} # initialize dictionaries for length count

number_proteins = 0
number_TR_proteins = 0
number_TRs = 0

set_name = "Proteome"

all_chr = [str(i) for i in range(1,23)] + ["X","Y"] # defining each chromosome
file_protein_list = os.path.join(output_path, set_name + '.txt')

pvalue_above_zero = 0
divergence_above_zero = 0
max_TR = 0
list_long_TR = []

# from https://www.proteinatlas.org/search/protein_class:COSMIC+somatic+mutations+in+cancer+genes
from cancer_analysis.import_tsv import cancer_proteome
cancer_TRs = {}

# from https://www.proteinatlas.org/search/cancer_specificity_rna%3Acolorectal+cancer%3Belevated+AND+sort_by%3Atissue+specific+score+AND+show_columns%3Acancerrnacategory
# elevated mRNA in colorectal cancer tissue
from cancer_analysis.import_tsv import colorec_canc_proteome
colorec_cancer_TRs = {}

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
            with open(output_pickle_file,'rb') as f: 
                denovo_list_remastered = pickle.load(f)

        # analyze all proteins that contain TRs
        if len(denovo_list_remastered.repeats) > 0:
            protein_tuple = UniqueIdentifier, EntryName, len(denovo_list_remastered.repeats)
            list_proteins_with_TR.append(protein_tuple)
            number_TR_proteins += 1
            number_TRs += len(denovo_list_remastered.repeats)

            # TR in cancer related gene?
            for gene in cancer_proteome:
                if GeneName == gene:
                    if gene not in cancer_TRs:
                        cancer_TRs[gene] = denovo_list_remastered.repeats
                        break

            # TR in colorectal cancer related gene?
            for gene in colorec_canc_proteome:
                if GeneName == gene:
                    if gene not in colorec_cancer_TRs:
                        colorec_cancer_TRs[gene] = denovo_list_remastered.repeats
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
                    longest_TR, protein, max_TR = denovo_list_remastered.repeats[TR], pyfaidx, TR_length

                if TR_length > 50:
                    list_long_TR.append((denovo_list_remastered.repeats[TR], pyfaidx))

                

#print("The longest TR is:", protein.name)
#print(protein)
#print(longest_TR)
print("length of longest TR:", max_TR)
# sp|Q9NZW4|DSPP_HUMAN longest TR in a first run, why does it find something like this?

print("There are {} TRs longer than 50.".format(len(list_long_TR)))

# for TR, protein in list_long_TR:
#     print(protein.name,protein)

print("pvalue not zero:", pvalue_above_zero)
print("divergence not zero",divergence_above_zero)

print("Cancer Genes with STRs:",len(cancer_TRs))
print("Colorectal Cancer Genes with STRs:",len(colorec_cancer_TRs))

for gene, TRs in colorec_cancer_TRs.items():
    for TR in TRs:
        print("************************")
        print(gene)
        print(TR)

colorec_cancer_TRs
chr_nr = "whole proteome"
chr_name = "whole proteome"



print("From {} cancer related genes contain {} STRs.".format(len(cancer_proteome),len(cancer_TRs)))
print("From {} colorectal cancer related genes contain {} STRs.".format(len(colorec_canc_proteome),len(colorec_cancer_TRs)))

# Calculate Figures and Overview
# analyzing_functions.TR_list_txt(list_proteins_with_TR, file_protein_list)
# analyzing_functions.AA_frequency(all_AA, chr_name, output_figures)
# analyzing_functions.l_n_distribution(count_dist, chr_name, output_figures)
analyzing_functions.overview(chr_name, number_proteins, number_TR_proteins, number_TRs)
