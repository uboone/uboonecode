import os, commands

jobs=[j.split()[0] for j in commands.getoutput('jobsub_q --user=kterao').split('\n') if j.find('fifebatch')>=0]
for j in jobs:
    os.system('jobsub_rm --jobid=%s' % j)

