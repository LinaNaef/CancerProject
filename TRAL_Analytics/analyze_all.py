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
######### Defining Paths and Parameters

## AA reference
working_dir = "/home/lina/Desktop/TRAL_Masterthesis/references/Uniprot_data/proteom_reference/pickles"
sequences_path = "/home/lina/Desktop/TRAL_Masterthesis/references/Uniprot_data/proteom_reference"
# https://www.uniprot.org/proteomes/UP000005640
output_path = "/home/lina/SynologyDrive/TRAL_Masterthesis/TRAL_Pipeline_Analytics/test_output/proteom"
output_figures = "/home/lina/SynologyDrive/TRAL_Masterthesis/TRAL_Pipeline_Analytics/output_figures"

list_proteins_with_TR = [] # initialize list to get all proteins containing TRs
all_AA = "" # initialize string for all AA in TRs
count_dist = {1:{},2:{},3:{}} # initialize dictionaries for length count
max_unit_count = 20 # expected max count for visualization of unit length, unit count distribution

number_proteins = 0
number_TR_proteins = 0
number_TRs = 0

set_name = "Proteome"

all_chr = [str(i) for i in range(1,23)] + ["X","Y"]
file_protein_list = os.path.join(output_path, set_name + '.txt')

pvalue_above_zero = 0
divergence_above_zero = 0
max_TR = 0
list_long_TR = []

# from https://www.proteinatlas.org/search/protein_class:COSMIC+somatic+mutations+in+cancer+genes
cancer_proteome = ['ABI1', 'ABL1', 'ABL2', 'ACKR3', 'ACSL3', 'ACSL6', 'ACVR1', 'ACVR2A', 'AFDN', 'AFF1', 'AFF3', 'AFF4', 'AKAP9', 'AKT1', 'AKT2', 'ALDH2', 'ALK', 'AMER1', 'APC', 'APOBEC3B', 'AR', 'ARHGAP26', 'ARHGEF12', 'ARID1A', 'ARID1B', 'ARID2', 'ARNT', 'ASPSCR1', 'ASXL1', 'ATF1', 'ATIC', 'ATM', 'ATP1A1', 'ATP2B3', 'ATR', 'ATRX', 'AXIN1', 'AXIN2', 'BAP1', 'BCL10', 'BCL11A', 'BCL11B', 'BCL2', 'BCL3', 'BCL6', 'BCL7A', 'BCL9', 'BCL9L', 'BCOR', 'BCORL1', 'BCR', 'BIRC3', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD3', 'BRD4', 'BRIP1', 'BTG1', 'BTK', 'BUB1B', 'C15orf65', 'CACNA1D', 'CALR', 'CAMTA1', 'CANT1', 'CARD11', 'CARS', 'CASP8', 'CBFA2T3', 'CBFB', 'CBL', 'CBLB', 'CBLC', 'CCDC6', 'CCNB1IP1', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD74', 'CD79A', 'CD79B', 'CDC73', 'CDH1', 'CDH11', 'CDK12', 'CDK4', 'CDK6', 'CDKN1B', 'CDKN2A', 'CDKN2C', 'CDX2', 'CEBPA', 'CEP89', 'CHCHD7', 'CHD4', 'CHEK2', 'CHIC2', 'CIC', 'CIITA', 'CLIP1', 'CLP1', 'CLTC', 'CLTCL1', 'CNBP', 'CNOT3', 'CNTRL', 'COL1A1', 'COL2A1', 'COX6C', 'CREB1', 'CREB3L1', 'CREB3L2', 'CREBBP', 'CRLF2', 'CRTC1', 'CRTC3', 'CSF3R', 'CTCF', 'CTNNB1', 'CUX1', 'CXCR4', 'CYLD', 'DAXX', 'DCTN1', 'DDB2', 'DDIT3', 'DDR2', 'DDX10', 'DDX3X', 'DDX5', 'DDX6', 'DEK', 'DICER1', 'DNAJB1', 'DNM2', 'DNMT3A', 'DROSHA', 'EBF1', 'ECT2L', 'EGFR', 'EIF3E', 'EIF4A2', 'ELF4', 'ELK4', 'ELL', 'ELN', 'EML4', 'EP300', 'EPAS1', 'EPS15', 'ERBB2', 'ERBB3', 'ERBB4', 'ERC1', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERG', 'ESR1', 'ETNK1', 'ETV1', 'ETV4', 'ETV5', 'ETV6', 'EWSR1', 'EXT1', 'EXT2', 'EZH2', 'EZR', 'FAM131B', 'FAM46C', 'FANCA', 'FANCC', 'FANCD2', 'FANCE', 'FANCF', 'FANCG', 'FAS', 'FAT1', 'FAT4', 'FBXO11', 'FBXW7', 'FCGR2B', 'FCRL4', 'FES', 'FEV', 'FGFR1', 'FGFR1OP', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FHIT', 'FIP1L1', 'FLCN', 'FLI1', 'FLT3', 'FLT4', 'FNBP1', 'FOXA1', 'FOXL2', 'FOXO1', 'FOXO3', 'FOXO4', 'FOXP1', 'FSTL3', 'FUBP1', 'FUS', 'GAS7', 'GATA1', 'GATA2','GATA3', 'GMPS', 'GNA11', 'GNAQ', 'GNAS', 'GOLGA5', 'GOPC', 'GPC3', 'GPHN', 'GRIN2A', 'H3F3A', 'H3F3B', 'HERPUD1', 'HEY1', 'HIF1A', 'HIP1', 'HIST1H3B', 'HIST1H4I', 'HLA-A', 'HLF', 'HMGA1', 'HMGA2', 'HNF1A', 'HNRNPA2B1', 'HOOK3', 'HOXA11', 'HOXA13', 'HOXA9', 'HOXC13', 'HOXD11', 'HOXD13', 'HRAS', 'HSP90AA1', 'HSP90AB1', 'IDH1', 'IDH2', 'IKBKB', 'IKZF1', 'IL2', 'IL21R', 'IL6ST', 'IL7R', 'IRF4', 'ITK', 'JAK1', 'JAK2', 'JAK3', 'JAZF1', 'JUN', 'KAT6A', 'KAT6B', 'KCNJ5', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KDSR', 'KEAP1', 'KIAA1549', 'KIF5B', 'KIT', 'KLF4', 'KLF6', 'KLK2', 'KMT2A', 'KMT2C', 'KMT2D', 'KNL1', 'KRAS', 'KTN1', 'LASP1', 'LCK', 'LCP1', 'LEF1', 'LHFP', 'LIFR', 'LMNA', 'LMO1', 'LMO2', 'LPP', 'LRIG3', 'LRP1B', 'LSM14A', 'LYL1', 'LZTR1', 'MAF', 'MAFB', 'MALT1', 'MAML2', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAPK1', 'MAX', 'MDM2', 'MDM4', 'MDS2', 'MECOM', 'MED12', 'MEN1', 'MET', 'MITF', 'MKL1', 'MLF1', 'MLH1', 'MLLT1', 'MLLT10', 'MLLT11', 'MLLT3', 'MLLT6', 'MN1', 'MNX1', 'MPL', 'MSH2', 'MSH6', 'MSI2','MSN', 'MTCP1', 'MTOR', 'MUC1', 'MUTYH', 'MYB', 'MYC', 'MYCL', 'MYCN', 'MYD88', 'MYH11', 'MYH9', 'MYO5A', 'MYOD1', 'NAB2', 'NACA', 'NBN', 'NCKIPSD', 'NCOA1', 'NCOA2', 'NCOA4', 'NCOR1', 'NCOR2', 'NDRG1', 'NF1', 'NF2', 'NFATC2', 'NFE2L2', 'NFIB', 'NFKB2', 'NFKBIE', 'NIN', 'NKX2-1', 'NONO', 'NOTCH1', 'NOTCH2', 'NPM1', 'NR4A3', 'NRAS', 'NRG1', 'NSD1', 'NSD2', 'NSD3', 'NT5C2', 'NTRK1', 'NTRK3', 'NUMA1', 'NUP214', 'NUP98', 'NUTM1', 'NUTM2A', 'NUTM2B', 'OLIG2', 'OMD', 'P2RY8', 'PAFAH1B2', 'PALB2', 'PATZ1', 'PAX3', 'PAX5', 'PAX7', 'PAX8', 'PBRM1', 'PBX1', 'PCM1', 'PCSK7', 'PDCD1LG2', 'PDE4DIP', 'PDGFB', 'PDGFRA', 'PDGFRB', 'PER1', 'PHF6', 'PHOX2B', 'PICALM', 'PIK3CA', 'PIK3R1', 'PIM1', 'PLAG1', 'PLCG1', 'PML', 'PMS1', 'PMS2', 'POLE', 'POT1', 'POU2AF1', 'POU5F1', 'PPARG', 'PPFIBP1', 'PPM1D', 'PPP2R1A', 'PPP6C', 'PRCC', 'PRDM1', 'PRDM16', 'PREX2', 'PRF1', 'PRKACA', 'PRKAR1A', 'PRRX1', 'PSIP1', 'PTCH1', 'PTEN', 'PTK6', 'PTPN11', 'PTPN13', 'PTPRB', 'PTPRC', 'PTPRK', 'PTPRT', 'PWWP2A', 'QKI', 'RABEP1', 'RAC1', 'RAD21', 'RAD51B', 'RAF1', 'RALGDS', 'RANBP17', 'RANBP2', 'RAP1GDS1', 'RARA', 'RB1', 'RBM10', 'RBM15', 'RECQL4', 'REL', 'RET', 'RHOA', 'RHOH', 'RMI2', 'RNF213', 'RNF43', 'ROS1', 'RPL10', 'RPL22', 'RPL5', 'RPN1', 'RSPO2', 'RSPO3', 'RUNX1', 'RUNX1T1', 'SALL4', 'SBDS', 'SDC4', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SEPT5', 'SEPT6', 'SEPT9', 'SET', 'SETBP1', 'SETD2', 'SF3B1', 'SFPQ', 'SFRP4', 'SH2B3', 'SH3GL1', 'SHTN1', 'SLC34A2', 'SLC45A3', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMARCE1', 'SMO', 'SND1', 'SOCS1', 'SOX2', 'SPECC1', 'SPEN', 'SPOP', 'SRC', 'SRGAP3', 'SRSF2', 'SRSF3', 'SS18', 'SS18L1', 'SSX1', 'SSX2', 'SSX4', 'STAG2', 'STAT3', 'STAT5B', 'STAT6', 'STIL', 'STK11', 'STRN', 'SUFU', 'SUZ12', 'SYK', 'TAF15', 'TAL1', 'TAL2', 'TBL1XR1', 'TBX3', 'TCEA1', 'TCF12', 'TCF3', 'TCF7L2', 'TCL1A', 'TERT', 'TET1', 'TET2', 'TFE3', 'TFEB', 'TFG', 'TFPT', 'TFRC', 'TGFBR2', 'THRAP3', 'TLX1', 'TLX3', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TNFRSF17', 'TOP1', 'TP53', 'TP63', 'TPM3', 'TPM4', 'TPR', 'TRAF7', 'TRIM24', 'TRIM27', 'TRIM33', 'TRIP11', 'TRRAP', 'TSC1', 'TSC2', 'TSHR', 'TTL', 'U2AF1', 'UBR5', 'USP6', 'USP8', 'VHL', 'VTI1A', 'WAS', 'WDCP', 'WIF1', 'WRN', 'WT1', 'WWTR1', 'XPA', 'XPC', 'XPO1', 'YWHAE', 'ZBTB16', 'ZCCHC8', 'ZFHX3', 'ZMYM2', 'ZNF331', 'ZNF384', 'ZNF521', 'ZRSR2']
cancer_TR_list = {}

# from https://www.proteinatlas.org/search/cancer_specificity_rna%3Acolorectal+cancer%3Belevated+AND+sort_by%3Atissue+specific+score+AND+show_columns%3Acancerrnacategory
# elevated mRNA in colorectal cancer tissue
colorec_canc_proteome = ['DEFA6', 'DEFA5', 'INSL5', 'AQP8', 'C17orf78', 'PRSS3', 'TEX43', 'MEP1A', 'ISX', 'KLK1', 'TM4SF20', 'TFF1', 'CA1', 'TBX10', 'BTNL3', 'MUC12', 'TMIGD1', 'APOBEC1', 'C5orf52', 'MS4A12', 'IGF2', 'MYADML2', 'NOX1', 'SPACA3', 'DUSP27', 'GUCA2A', 'CHST5', 'GUCY2C', 'MYO1A', 'NAT2', 'NXPE1', 'CDHR2', 'FABP2', 'LYPD8', 'XPNPEP2', 'SLC18A1', 'CHP2', 'HHLA2', 'NR1I2', 'ATOH1', 'CDHR5', 'CDX1', 'LRRC19', 'NXPE4', 'REG1A', 'REG4', 'TMEM236', 'CDH17', 'CLCA1', 'CLRN3', 'GAL3ST2', 'GPA33', 'MUC17', 'SOX1', 'CEACAM7', 'EPS8L3', 'FABP1', 'PDX1', 'PHGR1', 'PPP1R14D', 'VIL1', 'ACSL6', 'BTNL8', 'CASP5', 'CDX2', 'EFNA2', 'EREG', 'GJB4', 'HNF4A', 'IHH', 'MISP', 'MOGAT3', 'MUC13', 'NEU4', 'PRAP1', 'SPINK4', 'TRIM31', 'URAD', 'ZG16', 'AKR7A3', 'ANXA13', 'CDHR1', 'CLC', 'MGAM2', 'NOS2', 'PF4', 'SATB2', 'SULT1B1', 'VIP', 'AGMAT', 'AIFM3', 'ANKS4B', 'ANO9', 'AOC1', 'ARL14', 'ASCL2', 'ATP10B', 'B3GNT6', 'BCL2L14', 'BCL2L15', 'BEST2', 'C6orf222', 'CA7', 'CCL24', 'CCNI2', 'CEACAM5', 'CFTR', 'CLDN2', 'CXCL3', 'CYP2W1', 'CYP4F12', 'DMBT1', 'DPEP1', 'DUOXA2', 'EPHB2', 'FAM3D', 'FFAR4', 'FLJ22763', 'FRMD1', 'FUT3', 'GPR35', 'GRM8', 'HAND1', 'HIST1H1B', 'HIST1H4L', 'HOXB9', 'ITLN1', 'KRT20', 'KRT40', 'LCN15', 'LEFTY1', 'LGALS4', 'LGR5', 'LINC00675', 'LY6G6D', 'MAB21L2', 'MYO7B', 'NKD1', 'NMUR2', 'NXPE2', 'OLFM4', 'PIGR', 'PLCB4', 'POF1B', 'POU5F1B', 'PPBP', 'PTPRO', 'R3HDML', 'RETNLB', 'RNF186', 'RNF43', 'RUBCNL', 'RXFP4', 'SGK2', 'SLC26A3', 'SLC35D3', 'SLC6A7', 'TDGF1', 'TLDC2', 'TMEM150B', 'TMEM211', 'TPSG1', 'TRABD2A', 'TRIM15', 'TRIM40', 'TSPAN8', 'TUBAL3', 'UGT2B17', 'USH1C']
colorec_canc_cancer_TR_list = {}

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
                    if gene not in cancer_TR_list:
                        cancer_TR_list[gene] = denovo_list_remastered.repeats
                        break

            # TR in colorectal cancer related gene?
            for gene in colorec_canc_proteome:
                if GeneName == gene:
                    if gene not in colorec_canc_cancer_TR_list:
                        colorec_canc_cancer_TR_list[gene] = denovo_list_remastered.repeats
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

                # print(denovo_list_remastered.repeats[TR])

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

print("Cancer Genes with STRs:",len(cancer_TR_list))


print("Colorectal Cancer Genes with STRs:",len(colorec_canc_cancer_TR_list))
for gene, TRs in colorec_canc_cancer_TR_list.items():
    for TR in TRs:
        print("************************")
        print(gene)
        print(TR)

colorec_canc_cancer_TR_list
chr_nr = "whole proteome"
chr_name = "whole proteome"



print("From {} cancer related genes contain {} STRs.".format(len(cancer_proteome),len(cancer_TR_list)))
print("From {} colorectal cancer related genes contain {} STRs.".format(len(colorec_canc_proteome),len(colorec_canc_cancer_TR_list)))

# Calculate Figures and Overview
# analyzing_functions.TR_list_txt(list_proteins_with_TR, file_protein_list)
# analyzing_functions.AA_frequency(all_AA, chr_name, output_figures)
# analyzing_functions.l_n_distribution(count_dist, chr_name, output_figures)
analyzing_functions.overview(chr_name, number_proteins, number_TR_proteins, number_TRs)
