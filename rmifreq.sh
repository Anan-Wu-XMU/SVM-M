#!/bin/bash 
#remove failed jobs
for a in `ls *.log`
  do
  normal=`grep "Normal termi" $a | wc -l`
  if [ $normal -lt 2 ]
    then
      pre1=`echo $a | cut -d "." -f1`
      rm $a
      rm $pre1.com
      rm $pre1.chk
      echo "geometry optimization failed for $a!!!"
  fi
done

#remove transition states
for a in `ls *.log`
  do
  freq1=`grep Frequencies $a | head -n 1 | awk -e '{print $3}'`
#freq is a real number, need to round to a integer number.
  if [ ${freq1//.*/} -lt 0 ]
    then
      pre1=`echo $a | cut -d "." -f1`
      rm $a
      rm $pre1.com
      rm $pre1.chk
      echo "$a is a TS!"
  fi
done
#

if [ -e "freeG.txt" ] 
  then
    rm freeG.txt
fi
for a in `ls *.log`
  do 
  freeG=`grep "thermal Free Energies=" $a | awk -e '{print $8}'`
  echo "$a,$freeG" >> freeG.txt
done
nrow=`wc -l freeG.txt | awk -e '{print $1}'`
#remove conformers with ennergy 50 kj/mol higher than the most stable conformer
rmhien.py $1

#extract the new geometries
for a in `ls *.log`
  do 
    fname=`echo $a | cut -d "." -f1`
    sed -i '/Redundant internal coordinates found in file/,/Recover connectivity data from disk./w test.temp' $a
    mv test.temp $fname.temp
    sed -i '/Redundant internal coordinates found in file/d' $fname.temp
    sed -i '/Recover connectivity data from disk./d' $fname.temp
    sed -i 's/,/ /g' $fname.temp
    cat $fname.temp | awk -e '{print $1 "  " $3 "  " $4 "  " $5}' > $fname.new
    rm $fname.temp
done

#remove duplicated geometries
fname=`ls *.new | head -n 1`
natom=`cat $fname | wc -l`
#echo $natom
for a in `ls *.new`
  do
  if [ -e $a ]
    then 
    for b in `ls *.new`
      do
      if [ -e $b ]
        then 
        if [ "$a" != "$b" ]
          then
            idup=`Checkgeom $natom $a $b | awk -e '{print $1}'`
            if [[ $idup = 1 || $idup = 2 ]] 
              then 
                pre1=`echo $b | cut -d "." -f1`
#                echo "$b is duplicated with $a!"
                rm $b
                rm $pre1.log
                rm $pre1.com
                rm $pre1.chk
            fi
        fi
      fi
    done
  fi
done 
#rm *.err; rm *.out; rm head
