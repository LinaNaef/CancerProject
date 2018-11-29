# cd home/lina/Desktop/tral/easy_setup
# . activateTRAL.sh
# cd ..
# python3



# Test run and parse TRAL

import os
print(os.getcwd())

from tral.sequence import sequence
from tral.paths import PACKAGE_DIRECTORY

proteome_HIV = os.path.join(PACKAGE_DIRECTORY, "examples", "data", "HIV-1_388796.faa")
sequences_HIV = sequence.Sequence.create(file = proteome_HIV, input_format = 'fasta')

tandem_repeats = sequences_HIV[0].detect(denovo = True, detection = {"detectors": ["T-REKS"]})
print(len(tandem_repeats.repeats))
print(tandem_repeats.repeats[0])
# 1 repeat

tandem_repeats = sequences_HIV[0].detect(denovo = True, detection = {"detectors": ["Phobos"]})
print(len(tandem_repeats.repeats))
# 0 repeats

tandem_repeats = sequences_HIV[0].detect(denovo = True, detection = {"detectors": ["HHrepID"]})
print(len(tandem_repeats.repeats))
# 0 repeats

tandem_repeats = sequences_HIV[0].detect(denovo = True, detection = {"detectors": ["TRF"]})
print(len(tandem_repeats.repeats))
# 0 repeats

tandem_repeats = sequences_HIV[0].detect(denovo = True, detection = {"detectors": ["TRUST"]})
print(len(tandem_repeats.repeats))
# 0 repeats



# tandem_repeats = sequences_HIV[0].detect(denovo = True, detection = {"detectors": ["TRED"]})
# Did not find TRED result file in /tmp/tmpi_9sdpkc/001/TRED/tred.o
# Problem: TRED only works with python2! you have to run this with python2!


# Test for tred in python2
# Ctrl D
# python2

import os
print(os.getcwd())

from tral.sequence import sequence
from tral.paths import PACKAGE_DIRECTORY

proteome_HIV = os.path.join(PACKAGE_DIRECTORY, "examples", "data", "HIV-1_388796.faa")
sequences_HIV = sequence.Sequence.create(file = proteome_HIV, input_format = 'fasta')

tandem_repeats = sequences_HIV[0].detect(denovo = True, detection = {"detectors": ["TRED"]})