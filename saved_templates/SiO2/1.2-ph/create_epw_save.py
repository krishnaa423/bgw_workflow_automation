
#
# Post-processing script QE --> EPW
# 14/07/2015 - Samuel Ponce
#

from builtins import input
import numpy as np
import os

# Enter the number of irr. q-points
# prefix = input('Enter the prefix used for PH calculations (e.g. diam)\n')
prefix = 'SiO2'
outdir = './tmp'

# Enter the number of irr. q-points
nqpt = 8

# try:
#   nqpt = int(nqpt)
# except ValueError:
#   raise Exception('The value you enter is not an integer!')

if not os.path.exists('save'):
        os.mkdir('save')

for iqpt in np.arange(1,nqpt+1):
  label = str(iqpt)

  os.system('cp '+prefix+'.dyn'+str(iqpt)+' save/'+prefix+'.dyn_q'+label)
  if (iqpt == 1):
    os.system('cp ' + outdir + '/_ph0/'+prefix+'.dvscf1 save/'+prefix+'.dvscf_q'+label)
    os.system('cp -r ' + outdir + '/_ph0/'+prefix+'.phsave save/')

    # # KV: Maybe. 
    # os.system('cp '+prefix+'.dyn'+' save/'+prefix+'.dyn_q'+label)
  else:
    os.system('cp ' + outdir + '/_ph0/'+prefix+'.q_'+str(iqpt)+'/'+prefix+'.dvscf1 save/'+prefix+'.dvscf_q'+label)
    os.system('rm ' + outdir + '/_ph0/'+prefix+'.q_'+str(iqpt)+'/*wfc*' )

