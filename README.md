SVM-M
===

The SVM-M repository includes all scripts for the SVM-M structure elucidation of natural products. It contains two parts: 

   1. Given all intial conformers of a candidate strucure, automated evaluation of 13C/1H NMR cheimical shifts at the xOPBE/6-311+G(2D,P) level of theory. Geometry optimzations, conformational analysis, NMR calculations and Boltzman average are performed automatically

   2. Structure elucidation with the SVM-M protocol

   * q16                  : Linux script for running g16 (g16_model for the g16 environment setting)

   * autonmr.sh/autopcm.sh: perform automatical geometry optimization, conformation analysis, NMR calculations and Boltzman average in gas/solution phase

   * rmifreq.sh           : Remove failed geometry optimization and transition states. 
   
   * rmhien.py            : Remove the conformers with high free energies (50, 20 and 10 kJ/mol)
   
   * retnmr.sh            : Retrieve all NMR data from gaussian log files 
   
   * evanmr2.py           : Evaluate the Boltzman averaged 13C/1H NMR cheimical shifts
   
   * Checkgeom            : Perform the conformational analysis(Fortran code), in return 
   
                                 0: non-identical conformers
                                 
                                 1: identical conformers                    
                                 
   * nmrstat.py           : Perform the structure elucidation. 
      
Requirment 
===

Gaussian 16 package, OpenPBS compatible queue system, Linux (CentOS compatible),  Python3 and RDkit

If you are running queue sysstem other than OpenPBS, you may need to modify the script q16 and the corresponding environments.

Since redundant coordinates are used for conformational analysis, this repository requires a CentOS compatible system to run the precomplied Fortran code (Checkgeom)

Data availability
===

  The dataset are stored in SVM-correct-380.pickle for the correct molecules and SVM-inco-380.pickle for the incorrect molecules.

Download the model and run directly
===

   * autonmr.sh xxx.sdf queue  (Geometry optimization in the gas phase)
   
   * autopcm.sh xxx.sdf queue  (Geometry optimization in the solution phase)

   We note that all NMR calculations are performed in the solution phase.
   
   Once all calculations are done, you will find two files, named as allNMR.txt and allHNMR.txt, in which contain the Boltzman averaged 13C/1H chemical shieldings and chemical shifts. To perform the SVM-M structure elucidation, you need to put the experimental 13C chemical shifts into a file "exp.txt" and run the python script
   
   * nmrstat.py

Contact
===
  
  Anan Wu: ananwu@xmu.edu.cn
