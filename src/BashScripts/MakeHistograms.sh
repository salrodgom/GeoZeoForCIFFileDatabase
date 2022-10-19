#!/bin/bash 
# By Salvador R.G. Balestra, 2020
# This code calculates angle distributions in zeolitic structures
gfortran histograms.f90 -o histograms.exe -O0
# all_zeo_osio.dat
for code in ooo osio siosi sisisi sisi Staggering ; do
 if [ -f all_zeo_${code}.dat ] ; then rm -rf all_zeo_${code}.dat ; touch all_zeo_${code}.dat ; fi
 if [ $code == "ooo" ] ; then
   col=4; minmax="64 55.0 185"
 elif [ $code == "osio" ] ; then
   col=4; minmax="64 100.0 120.0"
 elif [ $code == "siosi" ] ; then
   col=4; minmax="64 115.0 185.0"
 elif [ $code == "sisisi" ] ; then
   col=4; minmax="64 55.0 165.0"
 elif [ $code == "Staggering" ] ; then
   col=4; minmax="128 0.0 60.01"
 elif [ $code == "sisi" ] ; then
   col=3; minmax="64 2.9 3.3" 
 fi
#
 ls *.cif | sed "s/\.cif//g" | while read abc ; do
 #for abc in "MTN_111" "MTN_222" ; do # <- debugging: check that both 111 and 222 supercells have the same histograms.
  if [ -f ${abc}_$code.dat ] ; then
  echo $abc $code
  if [ $code == "Staggering" ] ; then
   awk -v col=3 '{print $col}' ${abc}_$code.dat > input ; echo ${minmax} > bins ; ./histograms.exe < bins > ${abc}_${code}_1st_histogram.dat
   awk -v col=4 '{print $col}' ${abc}_$code.dat > input ; echo ${minmax} > bins ; ./histograms.exe < bins > ${abc}_${code}_2nd_histogram.dat
  else
   awk -v col=$col '{print $col}' ${abc}_${code}.dat > input ; echo ${minmax} > bins ; ./histograms.exe < bins > ${abc}_${code}_histogram.dat
  fi
  cat ${abc}_${code}.dat >> all_zeo_${code}.dat
  rm -rf input bins
  fi
 done
 if [ $code == "Staggering" ] ; then
  awk '{print $3}' all_zeo_${code}.dat > input ; echo ${minmax} > bins ; ./histograms.exe < bins > all_zeo_${code}_1st_histogram.dat ; rm -rf bins input
  awk '{print $4}' all_zeo_${code}.dat > input ; echo ${minmax} > bins ; ./histograms.exe < bins > all_zeo_${code}_2nd_histogram.dat ; rm -rf bins input
 else
  awk -v col=$col '{print $col}' all_zeo_${code}.dat > input ; echo ${minmax} > bins ; ./histograms.exe < bins > all_zeo_${code}_histogram.dat ; rm -rf bins input
 fi
done
if [ ! -d AllZeo ] ; then mkdir AllZeo ; fi
mv all_zeo_* AllZeo/.
rm -rf histograms.exe
for file in  *_Staggering_1st_histogram.dat ; do echo $file ; name=$(echo $file | sed 's/_1st_histogram\.dat//g') ; sed "s/NAME/${name}/g" gp_Staggering > gp_tmp ; gnuplot gp_tmp ; rm -rf gp_tmp ; done
for file in  *_siosi_histogram.dat ; do echo $file ; name=$(echo $file | sed 's/\.dat//g') ; sed "s/NAME/$name/g" gp_SiOSi > gp_tmp ; sed -i "s/FILE/$file/g" gp_tmp ; gnuplot gp_tmp ; rm -rf gp_tmp ; done
for file in  *_ooo_histogram.dat ; do echo $file ; name=$(echo $file | sed 's/\.dat//g') ; sed "s/NAME/$name/g" gp_OOO > gp_tmp ; sed -i "s/FILE/$file/g" gp_tmp ; gnuplot gp_tmp ; rm -rf gp_tmp ; done
for file in  *_sisisi_histogram.dat ; do echo $file ; name=$(echo $file | sed 's/\.dat//g') ; sed "s/NAME/$name/g" gp_SiSiSi > gp_tmp ; sed -i "s/FILE/$file/g" gp_tmp ; gnuplot gp_tmp ; rm -rf gp_tmp ; done
for file in  *_osio_histogram.dat ; do echo $file ; name=$(echo $file | sed 's/\.dat//g') ; sed "s/NAME/$name/g" gp_OSiO > gp_tmp ; sed -i "s/FILE/$file/g" gp_tmp ; gnuplot gp_tmp ; rm -rf gp_tmp ; done
for file in  *_sisi_histogram.dat ; do echo $file ; name=$(echo $file | sed 's/\.dat//g') ; sed "s/NAME/$name/g" gp_SiSi > gp_tmp ; sed -i "s/FILE/$file/g" gp_tmp ; gnuplot gp_tmp ; rm -rf gp_tmp ; done
if [ ! -d Figures ] ; then mkdir Figures ; fi  
mv *.eps Figures/.
cd Figures
 ls *.eps | sed 's/_/ /g' | awk '{print $1}' | uniq | while read abc ; do
  if [ ! -d $abc ] ; then mkdir $abc ; fi
  mv ${abc}*.eps $abc/.
 done
cd ..
exit 0
