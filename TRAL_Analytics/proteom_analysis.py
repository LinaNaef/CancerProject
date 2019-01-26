from obtain_TR_characteristics import *

#################################################
######### Analysis long TRs

#print("The longest TR is:", longest_TR_protein.name)
#print(longest_TR_protein)
#print(longest_TR)
print("The length of the longest TR is:", max_TR)
# sp|Q9NZW4|DSPP_HUMAN longest TR in a first run, why does it find something like this?

print("There are {} TRs longer than 50.".format(len(list_long_TR)))

# for TR, protein in list_long_TR: # print name of protein + whole protein
#     print(protein.name,protein)

# for TR, protein in list_long_TR: # print TRs
#     print("**********************************")
#     print(protein.long_name)
#     print(TR)

#################################################
######### Analyse Cancer Proteins of TRs

# for gene, TRs in colorec_cancer_TRs.items():
#     for TR in TRs:
#         print("************************")
#         print(gene)
#         print(TR)


print("From {} cancer related genes contain {} STRs.".format(len(cancer_proteome),len(cancer_TRs)))
print("From {} colorectal cancer related genes contain {} STRs.".format(len(colorec_canc_proteome),len(colorec_cancer_TRs)))

#################################################
######### Output txt files with genes conataining TRs

# all proteins containing TRs
file_protein_list = os.path.join(output_statistics, 'proteome_with_TRs.txt')
analyzing_functions.TR_list_txt(list_proteins_with_TR, file_protein_list)

# all colorectal cancer related proteins containing TRs
file_colorect_cancer_protein_list = os.path.join(output_statistics, 'colorect_cancer_proteome_with_TRs.txt')
analyzing_functions.TR_list_txt(list_colorect_cancer_proteins_with_TRs, file_colorect_cancer_protein_list)

# all cancer related proteins containing TRs
file_cancer_protein_list = os.path.join(output_statistics, 'cancer_proteome_with_TRs.txt')
analyzing_functions.TR_list_txt(list_cancer_proteins_with_TRs, file_cancer_protein_list)

#################################################
######### Calculate Figures and Overview

analyzing_functions.AA_frequency(all_AA, chr_name, output_statistics)
analyzing_functions.l_n_distribution(count_dist, chr_name, output_statistics)
analyzing_functions.overview(chr_name, number_proteins, number_TR_proteins, number_TRs)

#################################################
######### Analysis Statistics

print("pvalue not zero:", pvalue_above_zero)
print("divergence not zero",divergence_above_zero)



