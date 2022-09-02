#!/bin/bash 
# retrive NMR data
for a in `ls *NMR.log`
  do
  pref=`echo $a | cut -d '.' -f1`
  grep "C    Isotropic =" $a | awk -e '{print $1 " " $5}' > $pref.NMR
  grep "H    Isotropic =" $a | awk -e '{print $1 " " $5}' > $pref.HNMR
done

ls *.NMR > dir
sed -i 's/\.NMR//g' dir
# put all NMR data into one file
nmol=`ls -la *.NMR | wc -l`
if [ $nmol -gt 1 ]
  then
  ncount=1
  for a in `cat dir`
    do
    if [ $ncount -eq 1 ]
      then 
        cp $a".NMR" allNMR
        cp $a".HNMR" allHNMR
        ncount=$(($ncount+1))
    else
      join allNMR $a".NMR" > t1
      join allHNMR $a".HNMR" > t2
      mv t1 allNMR
      mv t2 allHNMR
    fi
  done
else
   cp *.NMR allNMR
   cp *.HNMR allHNMR
fi

rm dir
#rm *.NMR
