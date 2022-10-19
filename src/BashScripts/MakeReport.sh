#!/bin/bash
if [ -f Report.txt ] ; then rm -rf Report.txt ; touch Report.txt ; fi
if [ -f tmp ] ; then rm -rf tmp ; touch tmp ; fi
for file in *_Staggering_2nd_histogram.dat ; do 
 name=$(echo $file | sed 's/_Staggering_2nd_histogram.dat//g')
 n=$(sed '/#/d' $file | awk '{sum+= $2*$1} END {print sum}')
 m=$(sed '/#/d' ${name}_Staggering_1st_histogram.dat | awk '{sum+= $2*$1} END {print sum}')
 s_tot=$(grep "moments" ${name}_siosi_histogram.dat | awk '{print $4}')
 #s_ttt=$(grep "moments" ${name}_sisisi_histogram.dat | awk '{print $4}')
 echo $name $n $m $s_tot >> tmp
done
max_2=$(sort -gk2 tmp | awk '{print $2}' | tail -n1) ; max_4=$(sort -gk3 tmp | awk '{print $4}' | tail -n1) ; max_3=$(sort -gk3 tmp | awk '{print $2}' | tail -n1)
awk -v max_2=$max_2 -v max_3=$max_3 -v max_4=$max_4 '{print $1,$2/max_2 + $3/max_3 + $4/max_4}' tmp | sort -gk2 > Report.txt
