import os
import sys
from itertools import islice
import argparse

'''This script can parallelize a script -i into -n jobs (as long as the script supports --num_jobs and --index arguments'''

#Prepare input arguments
#num_jobs: split our input files into n different processes
#index: the index of the current running process, defines which chunk of input files to fetch
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", required = True)
parser.add_argument("-n", "--num_jobs", required = True)
args = parser.parse_args()
macro_name = args.input
num_jobs = args.num_jobs



#Configure submission file argument template
submit = 'universe = vanilla\n' ## writing .sub file
submit += 'arguments = "$(argument)"\n'
submit += 'output = submit01.out\n'
submit += 'error = submit01.err\n'
submit += 'log = submit01.log\n'
submit += '+JobFlavour = "microcentury"\n' 
submit += 'queue\n'
submitName = "plots/tmp/" + args.input + '.sub'
sub1 = open(submitName,'w')
sub1.write(submit+'\n')
sub1.close() ##finish writing .sub file
nfile = 1 ##number of sample file processed per each condor job


#Write executable command to sub-file for every job
for i in range(int(num_jobs)):

    #Set CMSSW environment
    create = '#!/bin/bash\n' ##writng .sh file
    create += 'export CMSSW_PROJECT_SRC=/afs/cern.ch/user/n/nhurley/CMSSW_12_4_6_emtfReEmul/src\n' 
    create += 'cd /afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/ \n' ## go to the src directly 

    #Run python script
    create += 'python ' + args.input + ' -n ' + args.num_jobs + ' -i ' + str(i) + ' \n'
    createName = 'plots/tmp/submit'+str(i)+'.sh' #Job's submission file output

    sub2 = open(createName,'w')
    sub2.write(create+'\n')
    sub2.close() ## finish writing .sh file
    os.system('chmod 755 '+createName + " &") ## make .sh file executable

for i in range(int(num_jobs)):
    if not i == int(num_jobs) - 1:
        os.system('condor_submit '+ submitName +' executable='+createName + " &") ## submit the job using condor_submit command
    else:
        os.system('condor_submit '+ submitName +' executable='+createName) ## submit the job using condor_submit command