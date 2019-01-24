from collections import Counter

def AACount(text):
    return Counter(c for c in text if c.isalpha())


print(AACount("AAAAAAAAAAATTTTTTDCFSDFSDFE"))


# import os
# from tral.sequence import sequence
# from tral.paths import PACKAGE_DIRECTORY

# proteome_HIV = os.path.join(PACKAGE_DIRECTORY, "examples", "data", "HIV-1_388796.faa")
# sequences_HIV = sequence.Sequence.create(file = proteome_HIV, input_format = 'fasta')

# tandem_repeats = sequences_HIV[0].detect(denovo = True, detection = {"detectors": ["XSTREAM"]})

# for iSequence in sequences_HIV:
#     iTandem_repeats = iSequence.detect(denovo = True)
#     iSequence.set_repeatlist(iTandem_repeats, "denovo")
    

# print([i for i in sequences_HIV[1].get_repeatlist('denovo').repeats if i.TRD == "T-REKS"][0])

