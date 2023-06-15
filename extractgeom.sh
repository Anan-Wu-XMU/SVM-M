#!/bin/bash
fname=`echo $1 | cut -d "." -f1`
sed -i '/Redundant internal coordinates found in file/,/Recover connectivity data from disk./w test.temp' $1
mv test.temp $fname.temp
sed -i '/Redundant internal coordinates found in file/d' $fname.temp
sed -i '/Recover connectivity data from disk./d' $fname.temp
sed -i 's/,/ /g' $fname.temp
cat $fname.temp | awk -e '{print $1 "  " $3 "  " $4 "  " $5}' > $fname.new
rm *.temp
