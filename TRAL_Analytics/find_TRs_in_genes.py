##########################################################################
### Importing required modules
##########################################################################

import os
import pickle
import time
import importlib

# to use this modules, TRAL has to be installed properly
from tral.sequence import sequence
from tral.repeat import repeat
from tral.repeat_list import repeat_list

# modules in this folder
import get_sequences

##########################################################################
### Defining paths and variables
##########################################################################

working_directory = "/home/lina/SynologyDrive/TRAL_Masterthesis/TRAL_Pipeline_Analytics/TRAL_Analytics/working_files"
sequences_path = "/home/lina/SynologyDrive/TRAL_Masterthesis/IBM_files/Assembled_genes"

## tissues to analyze
tissues = ["blood_derived_normal","primary_tumor","solid_tissue_normal"]

## genes to analyze
# genes_list = "/home/lina/SynologyDrive/TRAL_Masterthesis/TRAL_Pipeline_Analytics/Prepare_PASS/colorectal_msi_genes.txt"
# with open(genes_list) as g:
#     genes = g.readlines()
# genes = [x.strip() for x in genes]
# or define list with only few genes:

genes = ["APC","MLH1"]
# genes = ["APC"]
patients = [ patient for patient in os.listdir(sequences_path)]
# patients = ["Pat1","Pat10","Pat34"]


##########################################################################       
### use detector for sequences and get tandem repeats
##########################################################################


def loop_list_TR(list_TR, ignore_empty=True, **function):    
    """ returns length and number of TR units of a List of Repeatlists 
        Example of output:
            ...
             Sequence number 17 :
            1 :
            l:,  11 n:  5 
            
            
             Sequence number 18 :
            1 :
            l:,  5 n:  16 
            ...
Args:
        TR_list (list)
Keywords:
    ''    
    """     
    count = 0
    for i in list_TR:
        count +=1
        if ignore_empty:
            print_TR(i,count)
        else:
            print_TR(i,count,ignore_empty=False)



def print_TR_l_n(TR):
    """ returns length and number of TR units of a RepeatList
    Args:
        TR (class RepeatList)   
    """
    
    if len(TR.repeats) == 0:
        print("No Tandem Repeats detected.")
        
    else:
        for i in range(len(TR.repeats)):
            print(i+1,":")
            print( "l:, ",TR.repeats[i].l_effective,"n: ",TR.repeats[i].n,"\n")

def print_TR(TR, count, ignore_empty=True):
    """ returns a TRs from a RepeatList
    Args:
        TR (class RepeatList)
    
    """
    
    if len(TR.repeats) == 0:
        if not ignore_empty:
            print("\n Sequence number",count,":")
            print("No Tandem Repeats detected.")
        else:
            pass
        
    else:
        print("\n Sequence number",count,":")
        for i in range(len(TR.repeats)):
            print("\n",i+1,":")
            print(TR.repeats[i])


def get_TRs(detector,sequences,show_time=False):
    """ returns a list with lists of tandem repeats per sequence (class RepeatList)
    
    Args:
        detector (str): detector to be used
        sequence (TRAL class Sequence)
           
    """
    
    list_TR = []
    start = time.time()
    
    for i in range(len(sequences)):
        TR = sequences[i].detect(denovo = True, detection = {"detectors": [detector]})
        list_TR.append(TR)
        
    end = time.time()
    if show_time: 
        print(end - start)   

    return list_TR
