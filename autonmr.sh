#!/bin/bash 
# script to perform NMR calculations over different conformers and 
# give the statistic data of calculated NMR. 
module load python

if [ $# -lt 2 ]
  then 
    echo " Usage  : autonmr.sh fname queue" 
    echo " Usage  : autonmr.sh fname queue charge"
    echo "     fname: output from spartan"
    echo "     queue: gold1,E2680,gold2,E5 and snode"
    echo "====================================================="
    echo " recommendation Usage: nohup autonmr.sh fname queue &"
    echo "====================================================="
    echo "       Final results in avNMR.txt and nohup.out"
    echo "       Any problem, pls conact ananwu@xmu.edu.cn"
    exit 1
fi

queue=$2

#get prefix of $1 and create the directory $pref
pref=`echo $1 | cut -d "." -f1`
suffix=`echo $1 | cut -d "." -f2`
if [ ! -d $pref ]
  then
    mkdir $pref
    cd $pref
    cp ../$1 .
else
  echo "Dir $pref exist!"
  exit 1
fi

if [ $# -eq 2 ]
  then
    charge='0'
fi

if [ $# -eq 3 ]
  then
    charge=$3
fi

if [ $suffix = 'mol2' ]
  then 
#number of confermors in $1
    nmol=`grep MOLECULE $1 | wc -l`
#number of lines in $1
    nlin=`wc -l $1 | awk -e '{print $1}'`
#number of lines for split $1
    nsplit=$(($nlin/$nmol))
    split -l $nsplit $1 mmg
    ncount=1
    for a in `ls mmg*`
      do
        sed -i '/ATOM/,/BOND/w geomt' $a
        sed -i '/ATOM/d' geomt
        sed -i '/BOND/d' geomt
        awk -e '{print $2, " ",$3," ", $4," ",$5}' geomt > $pref"geom"$ncount
        echo " " >> $pref"geom"$ncount
        ncount=$(($ncount+1))
    done
    rm geomt
    rm mmg*
fi

if [ $suffix = 'sdf' ]
  then 
    nmol=`grep "M  END" $1 | wc -l`
    nlin=`wc -l $1 | awk -e '{print $1}'`
    nsplit=$(($nlin/$nmol))
    split -l $nsplit $1 mmg
    sed -i '4,4w t1' $1
    cut -c1-3 t1 > t2; mv t2 t1
    na=`awk -e '{print $1}' t1`
    ncount=1
    for a in `ls mmg*`
      do
        naend=$(($na+4))
        sed -i "5,$naend w geomt" $a
        awk -e '{print $4, " ",$1," ",$2," ",$3}' geomt > $pref"geom"$ncount
#        echo " " >> $pref"geom"$ncount
        ncount=$(($ncount+1))
    done
    rm t1
    rm geomt 
    rm mmg*
fi

if [ $nmol = 1 ]
  then 
cat > head << EOF
%mem=10GB
#P wB97XD/6-31G* OPT=(CalcFC,,MaxCyc=280) Freq nosym

wB97XD geometry optimization

$charge 1
EOF
echo " " > tail
for a in `ls *geom*`
  do
    cat head $a tail > $a.com
    rm $a 
    q16 $queue $a > /dev/null 2>&1
done

IFin=1
while [ $IFin != 0 ]
do
  sleep 60
  IFin=`qstat -u $USER | grep " "$pref"g" | wc -l`
done

echo "Number of conformers in wB97XD/6-31G*: $nmol"
rmifreq.sh 3
#need to run rmifreq.sh 3 twice to get the corrent free energies
rm *.new
nmol=`ls -la *.log | wc -l`
echo "Number of conformers in NMR calculations: $nmol"
rmifreq.sh 3
rm *.err; rm *.out; rm head

# evaluate NMR 
cat > head << EOF
%mem=10GB
#P OPBE/6-311+G(2D,P) NMR SCRF=(PCM,Solvent=Chloroform) IOP(3/76=0770002300,3/77=1000010000,3/78=1000010000)

xOPBE NMR

$charge 1
EOF

for a in `ls *geom*.new`
  do
    fname=`echo $a | cut -d "." -f1`
    cat head $fname.new > $fname-NMR.com
    echo "  " >> $fname-NMR.com
#    rm $a 
    q16 $queue $fname-NMR > /dev/null 2>&1
done
IFin=1
while [ $IFin != 0 ]
do
  sleep 60
  IFin=`qstat -u $USER | grep " "$pref"g" | wc -l`
done

#echo "Finished NMR calculations in Chloroform!"
rm *.err; rm *.out; rm *.chk; rm head; rm tail
# retrive NMR data into NMRall (shieldings)
retnmr.sh
natm=`wc -l *.NMR | head -n 1 | awk -e '{print $1}'`
evanmr2.py $natm
exit 1
fi

#generator Gaussian input, first run with PM7 geometry optimization!
cat > head << EOF
%mem=10GB
#P PM7 OPT=(CalcFC,MaxCyc=280) Freq

PM7 geometry optimization

$charge 1
EOF
echo " " > tail
for a in `ls *geom*`
  do 
    echo "%chk=$a.chk" > tchk
    cat tchk head $a tail > $a.com
    rm $a tchk
    q16 $queue $a > /dev/null 2>&1
done

IFin=1
while [ $IFin != 0 ]
do
  sleep 60
  IFin=`qstat -u $USER | grep " "$pref"g" | wc -l`  
done

#echo "Finshed geometry optimization at PM7!"
nmol=`ls -la *.log | wc -l`
echo "Number of conformers in PM7 calculations: $nmol"
rmifreq.sh 1
rm *.err; rm *.out; rm head
rm *.com; rm *.log; rm freeG.txt

cat > head << EOF
%mem=10GB
#P wB97XD/3-21G OPT=(ReadFC,MaxCyc=280) Freq

wB97XS geometry optimization

$charge 1
EOF
for a in `ls *geom*.new`
  do
    fname=`echo $a | cut -d "." -f1`
    echo "%chk=$fname.chk" > tchk
    cat tchk head $fname.new > $fname.com
    echo "  " >> $fname.com
    rm $a tchk
    q16 $queue $fname > /dev/null 2>&1
done

IFin=1
while [ $IFin != 0 ]
do
  sleep 60
  IFin=`qstat -u $USER | grep " "$pref"g" | wc -l`
done

#echo "Finished geometry optimization at wB97XD/3-21G"
nmol=`ls -la *.log | wc -l`
echo "Number of conformers in wB97XD/3-21G: $nmol"
rmifreq.sh 2
rm *.err; rm *.out; rm head
rm *.com; rm *.log; rm freeG.txt

cat > head << EOF
%mem=10GB
#P wB97XD/6-31G* OPT=(ReadFC,,MaxCyc=280) Freq

wB97XS geometry optimization

$charge 1
EOF
for a in `ls *geom*.new`
  do
    fname=`echo $a | cut -d "." -f1`
    echo "%chk=$fname.chk" > tchk
    cat tchk head $fname.new > $fname.com
    echo "  " >> $fname.com
    rm $a tchk
    q16 $queue $fname > /dev/null 2>&1
done

IFin=1
while [ $IFin != 0 ]
do
  sleep 60
  IFin=`qstat -u $USER | grep " "$pref"g" | wc -l`
done

#echo "Finished final geometry optimization at wB97XD/6-31G*"
nmol=`ls -la *.log | wc -l`
echo "Number of conformers in wB97XD/6-31G*: $nmol"
rmifreq.sh 3
#need to run rmifreq.sh 3 twice to get the corrent free energies
rm *.new
nmol=`ls -la *.log | wc -l`
echo "Number of conformers in NMR calculations: $nmol"
rmifreq.sh 3
rm *.err; rm *.out; rm head

# evaluate NMR 
cat > head << EOF
%mem=10GB
#P OPBE/6-311+G(2D,P) NMR SCRF=(PCM,Solvent=Chloroform) IOP(3/76=0770002300,3/77=1000010000,3/78=1000010000)

xOPBE NMR

$charge 1
EOF

for a in `ls *geom*.new`
  do
    fname=`echo $a | cut -d "." -f1`
    cat head $fname.new > $fname-NMR.com
    echo "  " >> $fname-NMR.com
#    rm $a 
    q16 $queue $fname-NMR > /dev/null 2>&1
done
IFin=1
while [ $IFin != 0 ]
do
  sleep 60
  IFin=`qstat -u $USER | grep " "$pref"g" | wc -l`
done

#echo "Finished NMR calculations in Chloroform!"
rm *.err; rm *.out; rm *.chk; rm head; rm tail
# retrive NMR data into NMRall (shieldings)
retnmr.sh
# modified by Anan Wu, 20220607, implementation of H_shifts
#natm=`wc -l *.NMR | head -n 1 | awk -e '{print $1}'`
#evanmr2.py $natm
n_C=`wc -l *.NMR | head -n 1 | awk -e '{print $1}'`
n_H=`wc -l *.HNMR | head -n 1 | awk -e '{print $1}'`
evanmr2.py $n_C $n_H
