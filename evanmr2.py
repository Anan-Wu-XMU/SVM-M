#!/share/apps/python3.7.1/bin/python3.7
# evaluate the averaged NMR chemical shifts
import math
import sys
import os
import numpy as np
import pandas as pd
from rdkit import Chem

def to_ch3(sdf):
   suppl=Chem.SDMolSupplier(sdf, removeHs=False) 
   mol=suppl[0]
   smartsrule = '[CH3]'
   patt = Chem.MolFromSmarts(smartsrule)
   tC_list = mol.GetSubstructMatches(patt)
   atomC_list = []
   atomH_list = []
   for a in tC_list:
      atomC_list.append(a[0])

   for c in atomC_list:
      iatom = mol.GetAtomWithIdx(c)
      tH_list = []
      for i in iatom.GetNeighbors():
         if i.GetAtomicNum() == 1:
            tH_list.append(i.GetIdx()+1)
      atomH_list.append(tH_list)
   return atomH_list

def to_oh(sdf):
   suppl=Chem.SDMolSupplier(sdf, removeHs=False)
   mol=suppl[0]
   smartsrule = '[OH]'
   patt = Chem.MolFromSmarts(smartsrule)
   O_list = mol.GetSubstructMatches(patt)
   atomO_list = []
   atomH_list = []
   for a in O_list:
      atomO_list.append(a[0])

   for c in atomO_list:
      iatom = mol.GetAtomWithIdx(c)
      for i in iatom.GetNeighbors():
         if i.GetAtomicNum() == 1:
            atomH_list.append(i.GetIdx()+1)
   return atomH_list

def main():
   fname1="freeG.txt"
   n_C=int(sys.argv[1])
   n_H=int(sys.argv[2])
   fn=open("avNMR.txt","w")
   fHn=open("avHNMR.txt","w")
   dataen=pd.read_table(fname1,sep=",",header=None)
   sdfName=dataen[0][0].split("geom")[0]+".sdf"

   emin=min(dataen[1])
   factor=627.51
   R=0.00198718
   T=298.15
   TMS=191.6379
   TMS_H=31.6474
   en=(dataen[1]-emin)*factor
   nmol=en.shape[0]
   prob=np.zeros(nmol)
   for i in range(nmol):
     prob[i]=math.exp(-en[i]/R/T)
   total=prob.sum()
# prob contsist of probability for each conformer 
   prob=prob/total
   nmol=dataen[1].shape[0]
   shielding=np.zeros(n_C)
   Hshielding=np.zeros(n_H)
   shift=np.zeros(n_C)
   Hshift=np.zeros(n_H)

# Boltzmann average for each conformer
   CH3_l=to_ch3(sdfName)
   n_CH3=len(CH3_l)
   for i in range(0,nmol,1):
     str=dataen[0][i].split(".")
     iprob=prob[i]
     # C shielding
     fnmr=str[0]+"-NMR.NMR"
     dorig=pd.read_table(fnmr,sep=" ",header=None)
     t1=dorig[1]*iprob
     shielding=shielding+t1
     # H shielding
     fHnmr=str[0]+"-NMR.HNMR"
     if n_CH3 > 0:
        with open(fHnmr, 'r') as f:
           origH_l=f.readlines()
        shield_av=[]
        Hshdict={}
        for line in origH_l:
           Hshdict[int(line.split()[0])] = line.split()[1]
       # print("orig_dic:", Hshdict)
        for l in CH3_l:
           sumShield=0
           for Hs in l:
              sumShield += float(Hshdict[Hs])
           meanShield=sumShield/3
           shield_av.append(format(meanShield, '0.4f'))

        for a in range(len(shield_av)):
           num=CH3_l[a]
           for b in num:
              Hshdict[b] = shield_av[a]
#        print("new_dic:", Hshdict)
 
        dHorig=pd.DataFrame(list(Hshdict.items()))
        dH_tmp=pd.to_numeric(dHorig[1],errors='coerce')
        t2=dH_tmp*iprob
       # dHorig=pd.DataFrame(pd.Series(Hshdict))
     else:
        dHorig=pd.read_table(fHnmr,sep=" ",header=None)    
        t2=dHorig[1]*prob[i]
     Hshielding=Hshielding+t2
   shift=TMS-shielding
   Hshift=TMS_H-Hshielding

   print("=====================================================")
   print("      Conformer   Relative En(kcal/mol)   Probability")
   for i in range(0,nmol,1): 
      print(f'  {dataen[0][i]:13s}         {en[i]:6.2f}             {prob[i]:8.4f}')
   print("=====================================================")

   print("Boltzmann averaged NMR shieldings and Chemical shifts!")
   print("=====================================================")
   print(" Atom    Shielding   Shifts")
   fn.write("Atom  Shieldings  Shifts\n")
   for i in range(0,n_C,1):
     print(f'C{dorig[0][i]:3d}   {shielding[i]:8.4f}  {shift[i]:8.2f}')
     fn.write("C%3d  %8.4f  %8.2f\n" % (dorig[0][i],shielding[i],shift[i]))
   print("=====================================================")
   print("      ")
   print("=====================================================")
   print(" Atom    Shielding   Shifts  Note")
   fHn.write("Atom  Shieldings  Shifts  Note\n")

   OH_l=to_oh(sdfName) 
   for i in range(0,n_H,1):
     flag = 0
     for j in CH3_l:
        if dHorig[0][i] in j :
           print(f'H{dHorig[0][i]:3d}    {Hshielding[i]:8.4f}  {Hshift[i]:8.2f}   CH3_{CH3_l.index(j):1d}')
           fHn.write("H%3d  %8.4f  %7.2f   CH3_%1d\n" % (dHorig[0][i],Hshielding[i],Hshift[i],CH3_l.index(j)))
           flag = 1
     if flag==0:
        if dHorig[0][i] in OH_l :
           print(f'H{dHorig[0][i]:3d}    {Hshielding[i]:8.4f}  {Hshift[i]:8.2f}   OH')
           fHn.write("H%3d  %8.4f  %7.2f   OH\n" % (dHorig[0][i],Hshielding[i],Hshift[i]))
        else:   
           print(f'H{dHorig[0][i]:3d}    {Hshielding[i]:8.4f}  {Hshift[i]:8.2f}')
           fHn.write("H%3d  %8.4f  %7.2f\n" % (dHorig[0][i],Hshielding[i],Hshift[i]))

   print("=====================================================")


if __name__ == "__main__":
  main()

   
