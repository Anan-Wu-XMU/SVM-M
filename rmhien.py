#!/share/apps/python3.7.1/bin/python3.7
# remove the conformers with high free energies 

import sys
import os
import math
import string
import numpy as np
import pandas as pd

def main():
   fname="freeG.txt"
   ninput1=sys.argv[1]
   ien=int(ninput1)
   dataen=pd.read_table(fname,sep=",",header=None)
   nrow=dataen[0].shape[0]
   if ien==1:
     ecut=50.0
   if ien==2:
     ecut=20.0
   if ien==3:
     ecut=10.0
   enr=np.zeros(nrow)
# factor from atomic unit to kj/mol
   factor=2625.502
#evaluate energy difference in kj/mol
   enmin=min(dataen[1])
   for i in range(0,nrow):
     edif=(dataen[1][i]-enmin)*factor
     enr[i]=edif
     if edif >= ecut:
#       print(data[i,0], "with high energy!")   
       os.system("rm "+dataen[0][i])

if __name__ == "__main__":
  main()

   
