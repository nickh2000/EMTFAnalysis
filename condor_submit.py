import os
import sys
from itertools import islice
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", required = True)
parser.add_argument("-n", "--num_jobs", required = True)

args = parser.parse_args()
macro_name = args.input
num_jobs = args.num_jobs


submit = 'universe = vanilla\n' ## writing .sub file
submit += 'arguments = "$(argument)"\n'
submit += 'output = submit01.out\n'
submit += 'error = submit01.err\n'
submit += 'log = submit01.log\n'
submit += '+JobFlavour = "microcentury"\n' 
submit += 'queue\n'
submitName = args.input + '.sub'
sub1 = open(submitName,'w')
sub1.write(submit+'\n')
sub1.close() ##finish writing .sub file
nfile = 1 ##number of sample file processed per each condor job

for counter1 in range(int(num_jobs)):
    create = '#!/bin/bash\n' ##writng .sh file
    create += 'export CMSSW_PROJECT_SRC=/afs/cern.ch/user/n/nhurley/CMSSW_12_4_6_emtfReEmul/src\n' 
    create += 'cd /afs/cern.ch/user/n/nhurley/EMTFAnalyzer/AWBTools/macros/ \n' ## go to the src directly 
    #create += 'eval `scramv1 runtime -sh`\n' ## this is effectively "cmsenv"
    #create += 'export X509_USER_PROXY=/afs/cern.ch/user/r/rmccarth/x509up_u120948\n' ## exporting the proxy
    #create += 'cd /afs/cern.ch/user/r/rmccarth/private/dispVert/CMSSW_12_5_0_pre2/src/L1Trigger/TrackFindingTracklet/test\n'## go the directly where you have the producer and python config file
    create += 'python ' + args.input + ' -n ' + args.num_jobs + ' -i ' + str(counter1) + ' \n'
    createName = 'plots/tmp/submit'+str(counter1)+'.sh'
    sub2 = open(createName,'w')
    sub2.write(create+'\n')
    sub2.close() ## finish writing .sh file
    counter1+=1
    os.system('chmod 755 '+createName) ## make .sh file executable
    os.system('condor_submit '+ submitName +' executable='+createName) ## submit the job using condor_submit command