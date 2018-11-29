""" 12.06.2018: Näf, Paulina
    
    This skript provides functions to parse directories with fasta files
    containing DNA or AA TR sequences and to apply TRDs on these.
    Finally, one can create a dictionary structure which can be used easily  to apply statistics on.

    To be used with TRAL (Schaper et. al. 2015)
    
    TODO:
        Error Handling
        
        """

##########################################################################
### Importing required modules
##########################################################################

import os
import pickle
import time

# to use this modules, TRAL has to be installed properly
from tral.sequence import sequence
from tral.repeat import repeat
from tral.repeat_list import repeat_list

##########################################################################
### Create dictionairy with TRs per detector, per divergence 
### and per TR length and number of repeat units
##########################################################################

def TR_dict_per_detector_per_PAM(list_sequences,number_sequences, list_detectors, show_time = False):
    """ returns dict with detectors as keys, 
        within a dict with PAM as keys,
        within a dict with TR_types as keys,
        within a list with RepeatLists of these TR_types
    
     
        Args:
            list_sequences (list): output from parse_sequence_files().
            number_sequences (int): how many sequences should be taken (at most as much as can be find in the input list).
            list_detectors (list): list of detectors which should be used to detect TRs.
       
       
        Example to look at generated data:
            loop_l_n_list_TR(TR_DNA["T-REKS"]["120"]["T-REKS_21_4"])
            loop_list_TR(TR_DNA["T-REKS"]["80"]["T-REKS_25_4"])
    
       """
       
    start = time.time()
    detectors = list_detectors
    dict_detectors_dict_per_PAM = dict()# new function
   


    for detector in detectors:
        
        PAM = ["40","80","120"]
        dict_per_PAM = dict()
        
        for i in range(len(list_sequences)):  
            
            dict_per_PAM[PAM[i]] = TR_per_sequence(list_sequences[i],detector,number_sequences)
            
        dict_detectors_dict_per_PAM[detector] = dict_per_PAM    
        print("Detector",detector,"done")
        
    end = time.time()
    
    if show_time:
        print(end - start)        
        
    return(dict_detectors_dict_per_PAM)   



##########################################################################
### Get TRs from sequences
##########################################################################

def TR_per_sequence(sequences_per_TR_type,detector,number_sequences,show_time=False):
    
    """ takes a dictionary with keys in format l_n (definition of TR_type) and corresponding sequences
    and a certain detector and how many sequences should be analysed 
    and returns a dictionary with detected TRs
    
    
    Args:
        sequences_per_TR_type (dict):  keys in format l_n (definition of TR_type) and corresponding sequences
        number_sequences (int): how many sequences should be taken (at most as much as can be find in the input list). 
        
        """
    
    start = time.time()
    predicted_TR_per_TR_type = dict()
    
    for TR_type in sorted(sequences_per_TR_type.keys()):
        TRs = get_TR(detector,sequences_per_TR_type[TR_type],number_sequences)
        predicted_TR_per_TR_type[detector + "_" + TR_type] = TRs
        print("Calculated TRs from sequences with TR_type",TR_type)
        
    end = time.time()   
    
    if show_time:
        print("Function TR_per_sequence() took",end - start,"seconds")  
        
    return predicted_TR_per_TR_type
    

##########################################################################
### Iteration through files
##########################################################################

def get_sequences_in_directory_per_TR_type(directory_in_str,list_l,list_n,show_time=False):
    
    """ returns a dictionary with detected TRs
            key in form: [l_n]
    
    
    Args:
        directory_in_str (str): path with tandem repeat files
        list_l (list): length of simulated TR unit sequences which should be taken
        list_n (list): number of simulated TR unit sequences which should be taken
        
        """
     
    start = time.time()
    sequences_per_TR_type = dict()
    
    for l in list_l:
        
        for n in list_n:
            
            TR_type = "_".join([str(l),str(n)])
            filename = TR_type + ".faa"
            path = os.path.join(directory_in_str, filename)
            sequences = sequence.Sequence.create(file = path, input_format = 'fasta')
            sequences_per_TR_type[TR_type] = sequences
            print("Got sequences from TR_type",TR_type)
            
    end = time.time()   
    
    if show_time:
        print("Function get_sequences_in_directory_per_TR_type() took",end - start,"seconds")  
        
    return sequences_per_TR_type


def parse_sequence_files(main_directory,subdirectories,list_l, list_n):
    
    """ returns a list with dictonaries with detected TRs
            key in form: [l_n]
    
    
    Args:
        main_directory (str): path to subdirectories
        subdirectories (list): list with subdirectories as str (part behind main directory)
        list_l (list): length of simulated TR unit sequences which should be taken
        list_n (list): number of simulated TR unit sequences which should be taken
        
        """
        
    all_sequences = []
    
    for subdirectory in subdirectories:
        
        print(subdirectory)
        path = main_directory + subdirectory
        sequences = get_sequences_in_directory_per_TR_type(path,list_l,list_n,show_time=True)
        all_sequences.append(sequences)
        
    return(all_sequences)    


##########################################################################       
### use detector for sequences and get tandem repeats as a list
##########################################################################

def get_all_TR(detector,sequences,show_time=False):
    
    """ returns a list with lists of tandem repeats per sequence (class RepeatList)
    
    Args:
        detector (str): detector to be used
        sequence (class Sequence)
           
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
    
def get_TR(detector,sequences,number_sequences,show_time=False):
    
    """ returns a list with lists of tandem repeats per sequence (class RepeatList)
    
    Args:
        detector (str): detector to be used
        sequence (class Sequence)
        number_sequences (int): how many sequences should be taken (at most as much as can be find in the input list). 
           
    """
    list_TR = []
    start = time.time()
    
    for i in range(number_sequences):
        
        TR = sequences[i].detect(denovo = True, detection = {"detectors": [detector]})
        list_TR.append(TR)
        
    end = time.time()
    
    if show_time:
        
        print(end - start)       
        
    return list_TR
    
### printing of tandem repeats (Class: Repeat)    
    
def print_TR(TR):
    
    """ returns a TRs from a RepeatList
    
    Args:
        TR (class RepeatList)
    
    """
    
    
    if len(TR.repeats) == 0:
        
        print("No Tandem Repeats detected.")
        
    else:
        
        for i in range(len(TR.repeats)):
            print("\n",i+1,":")
            print(TR.repeats[i]) 
            
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
            
### printing of tandem repeats (Class: RepeatList)
# iterates through Lists and then uses functions for Repeat            

def loop_l_n_list_TR(list_TR): 
    
    """ returns detected TRs of a List of Repeatlists
    
    
        Example of output:
            
            ...
    
            Sequence number 48 :
            
             1 :
            > begin:1 l_effective:3 n:5
            -GT-
            -GTC
            CATT
            -GGC
            -GGG
            
             2 :
            > begin:27 l_effective:18 n:4
            TGT-C--TAATCTGAGGTAAT
            -GAG---GAATT-G-GCCT-T
            TGCGCGGTTATC--AGG-AAG
            TGCGC--TTATA-GAGACAAT
            
             Sequence number 49 :
            
             1 :
            > begin:51 l_effective:6 n:6
            TG---TA-TC--
            TG--GTA---A-
            -GATCGA-C---
            ----CGAGC-AG
            CG--CT--C-T-
            -G--CTA-C-A-
            
            ...
            
                
    Args:
        TR_list (list)
    
    """
        
    count = 0
    
    for i in list_TR:
        
        count +=1
        print("\n Sequence number",count,":")
        print_TR_l_n(i)

def loop_list_TR(list_TR):
    
    """ returns length and number of TR units of a List of Repeatlists
    
    
        Example of output:
            
            ...
    
             Sequence number 15 :
            1 :
            l:,  3 n:  8 
            
            
             Sequence number 16 :
            1 :
            l:,  6 n:  4 
            
            2 :
            l:,  8 n:  4 
            
            
             Sequence number 17 :
            1 :
            l:,  11 n:  5 
            
            
             Sequence number 18 :
            1 :
            l:,  5 n:  16 
            
            ...
            
                
    Args:
        TR_list (list)
    
    """    
    
    count = 0
    
    for i in list_TR:
        
        count +=1
        print("\n Sequence number",count,":")
        print_TR(i)
            

            

    
