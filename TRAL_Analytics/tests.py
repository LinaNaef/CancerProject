import os
from tral.sequence import sequence
from tral.paths import PACKAGE_DIRECTORY

proteome_HIV = os.path.join(PACKAGE_DIRECTORY, "examples", "data", "HIV-1_388796.faa")
sequences_HIV = sequence.Sequence.create(file = proteome_HIV, input_format = 'fasta')

tandem_repeats = sequences_HIV[0].detect(denovo = True, detection = {"detectors": ["XSTREAM"]})

for iSequence in sequences_HIV:
    # iTandem_repeats = iSequence.detect(denovo = True)
    iSequence.set_repeatlist(iTandem_repeats, "denovo")
    for iTandemRepeat in iSequence.get_repeatlist('denovo').repeats:
        iTandemRepeat.calculate_pvalues()
