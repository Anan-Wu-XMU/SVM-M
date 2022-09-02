* SVM-M
The SVM-M repository incldues all scripts for stucture elucidation of natural products. It contains 
two parts: 

   1. Given all intial conformers of a candidate strucure, automated calculation of 13C/1H NMR 
      cheimical shifts at the xOPBE/6-311+G(2D,P) level of theory. Geometry optimzations, conformational
      analysis, NMR calculations and Boltzman average are performed automatically.

   2. Structure elucidation with the SVM-M protocol. 

   autonmr.sh/autopcm.sh: perform automatical geometry optimization, conformation analysis, NMR calculations 
                          and Boltzman average in gas/solution phase. 

   rmifreq.sh           : Remove failed geometry optimization and transition states. 
   rmhien.py            : Remove the conformers with high free energies (50, 20 and 10 kJ/mol)
   retnmr.sh            : Retrieve all NMR data from gaussian log files 
   evanmr2.py           : Evaluate the Boltzman averaged 13C/1H NMR cheimical shifts
   Checkgeom            : Perform the conformational analysis(Fortran code), in return 
                                 0: non-identical conformers
                                 1: identical conformers                                       
   nmrstat.py           : Perform the structure elucidation. 
      
* Requirment 

Gaussian 16 package, OpenPBS compatible queue system, Linux (CentOS compatible),  Python3 and rdkit

If we are running queue sysstem other than OpenPBS, we may need to modify the script q16 and autonmr.sh/autopcm.sh

Since redundant coordinates are used for conformational analysis, this repository requires a CentOS compatible system
to run the precomplied Fortran code (Checkgeom)

* Data availability

  The dataset are stored in SVM-correct-380.pickle for the correct molecules and SVM-inco-380.pickle for the incorrect molecules.

* Download the model and run directly

   autonmr.sh xxx.sdf queue  (Geometry optimization in the gas phase)
   autopcm.sh xxx.sdf queue  (Geometry optimization in the solution phase)

   We note that all NMR calculations are performed in the solution phase.

* Contact
  
  Anan Wu: ananwu@xmu.edu.cn
