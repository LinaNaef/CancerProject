# cd /home/lina/Desktop/tral/easy_setup
# . activateTRAL.sh
# cd /home/lina/SynologyDrive/TRAL_Masterthesis/TRAL_Pipeline_Analytics/TRAL_Analytics
# python3


import sys
import pickle


# manipulating system path for this session
# only used if module handlingTRALsequencesTRs.py is not in the sytem path

working_dir = "/home/lina/SynologyDrive/TRAL_Masterthesis/test_TRAL" # adapt to your own path
sys.path.append(working_dir)
#print('\n'.join(sys.path))

import handlingTRALsequencesTRs



""" 12.06.2018: Delucchi, Matteo; Näf, Paulina
    This skript provides an example workflow, how to use handlingTRALsequencesTRs
    To be used with TRAL (Schaper et. al. 2015)
    """


##########################################################################
### Getting sequences from file system (sequences in FASTA-format)
##########################################################################

### Defining paths for sequences which should be obtained
# all paths should be adapted to your own paths

main_directory = "/home/lina/SynologyDrive/Maria_Files_20180604/TandemRepeat_Detection_Benchmark/"
subdirectory_DNA_40 = "/DNA/Simulated_TandemRepeats/40/0.0"
subdirectory_DNA_80 = "/DNA/Simulated_TandemRepeats/80/0.0"
subdirectory_DNA_120_0 = "/DNA/Simulated_TandemRepeats/120/0.0"

# does not work for simulated TRs with indels (problem within TRAL)
#subdirectory_DNA_120_01 = "/DNA/Simulated_TandemRepeats/120/0.01" 
#subdirectory_DNA_120_001 = "/DNA/Simulated_TandemRepeats/120/0.001"

subdirectory_AA_40 = "/AA/Simulated_TandemRepeats/40/0.0"
subdirectory_AA_80 = "/AA/Simulated_TandemRepeats/80/0.0"
subdirectory_AA_120_0 = "/AA/Simulated_TandemRepeats/120/0.0"
# does not work for simulated TRs with indels (problem within TRAL)
#subdirectory_AA_120_01 = "/AA/Simulated_TandemRepeats/120/0.01"
#subdirectory_AA_120_001 = "/AA/Simulated_TandemRepeats/120/0.001"

subdirectories_DNA = [subdirectory_DNA_40, subdirectory_DNA_80, subdirectory_DNA_120_0]
subdirectories_AA = [subdirectory_AA_40, subdirectory_AA_80, subdirectory_AA_120_0]


list_l = list(range(1,(25+1),4)) # length of TR units (every fourth)
#list_l = list(range(1,(25+1))) # all
list_n = list(range(2,(15+1),2)) # number of repeats

## Getting Sequences from DNA and from AA
        
list_DNA = handlingTRALsequencesTRs.parse_sequence_files(main_directory,subdirectories_DNA,list_l, list_n)
list_AA = handlingTRALsequencesTRs.parse_sequence_files(main_directory,subdirectories_AA,list_l, list_n)


# Saving this lists as binary files:
with open(working_dir + 'list_sequences.pkl', 'wb') as f: 
    pickle.dump([list_DNA,list_AA], f)

##########################################################################
### de novo TR detection
##########################################################################

# Getting back the objects (later):
with open(working_dir + '/analysislist_sequences.pkl','rb') as f: # files saved before  
    list_DNA, list_AA = pickle.load(f)


# Creating TR-dicts for detectors detectors

#detectors = ["TRF","T-REKS","TRUST","Phobos","HHrepID"]
detectors = ["Phobos"]

TR_DNA = handlingTRALsequencesTRs.TR_dict_per_detector_per_PAM(list_DNA,3, detectors, show_time= True)
TR_AA = handlingTRALsequencesTRs.TR_dict_per_detector_per_PAM(list_AA,3, detectors, show_time= True) # does work for TRUST with 40 but not 80 sequences

# Saving files
# ATTENTION: change file name if you don't want to overwrite a previous result!!!
with open(working_dir + 'dict_TR_DNA_3seq.pkl', 'wb') as f: 
    pickle.dump(TR_DNA, f)

    
with open(working_dir + 'dict_TR_AA_3seq.pkl', 'wb') as f: 
    pickle.dump(TR_AA, f)


# small example how to look at this data:

# shows all detected TRs for the defined properties:  
handlingTRALsequencesTRs.loop_list_TR(TR_DNA["T-REKS"]["120"]["T-REKS_25_4"])

handlingTRALsequencesTRs.loop_list_TR(TR_DNA["Phobos"]["120"]["Phobos_25_4"])
"""

Sample Output:

... 


 Sequence number 48 :

 1 :
> begin:59 l_effective:3 n:8
-AG--A-
-AG--AG
-GGTCA-
-AG--G-
-AG--GT
-TC-CA-
-TA--T-
TAA--A-

 2 :
> begin:8 l_effective:6 n:8
GCG-AACA-
--G-CTAAC
-C--CAATT
-CTAC-GG-
-CT-TTGG-
-GA-CAGC-
-CT--TGA-
-CT-CTGA-

 Sequence number 49 :

 1 :
> begin:3 l_effective:12 n:2
TGCGTACACGAAAG-G
TG--TA-ACAATTGAG

 2 :
> begin:35 l_effective:8 n:3
CCGA-TGA
CCGTGTGA
CCTCG-GA

 3 :
> begin:59 l_effective:7 n:4
ATGA-AG-
ATTTTAGT
ATAGTCG-
CTTATA-T

 Sequence number 50 :

 1 :
> begin:12 l_effective:5 n:3
-CTACA--
-CT-CGTC
ACA-CG-C

 2 :
> begin:42 l_effective:16 n:3
GCC-GTCATCAGAATGT-
GCTACAGGTGAA-ATACC
GCTTATAAT-AA-ATAA-


...

"""



# shows all lenghts and numbers of the detected TRs for the defined properties
handlingTRALsequencesTRs.loop_l_n_list_TR(TR_DNA["TRUST"]["40"]["TRUST_21_10"])  

""" Sample Output:

 Sequence number 1 :
1 :
l:,  42 n:  3 


 Sequence number 2 :
No Tandem Repeats detected.

 Sequence number 3 :
1 :
l:,  42 n:  4 


 Sequence number 4 :
1 :
l:,  21 n:  3 




 ....
 
 
 
 
 Sequence number 78 :
1 :
l:,  42 n:  5 


 Sequence number 79 :
1 :
l:,  42 n:  3 


 Sequence number 80 :
1 :
l:,  42 n:  5 

"""  


##########################################################################


# Getting back the objects (later):with open(working_dir + 'dict_TR_40seq_DNA.pkl','rb') as f:  
with open(working_dir + 'dict_TR_DNA_80seq.pkl','rb') as f:  
    TR_DNA = pickle.load(f)

with open(working_dir + 'dict_TR_AA_40seq.pkl','rb') as f:  
    TR_AA = pickle.load(f)
    
