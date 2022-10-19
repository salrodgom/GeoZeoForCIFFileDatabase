#!/bin/bash
# make the header of the CSV file:
if [ -f data.csv ] ; then rm -rf data.csv ; touch data.csv ; fi
#if [ ! -d Download_Database ] ; then mkdir Download_Database ; fi
echo "NameCode;TSites;RingSizes;TopologicalDensity;StaggeringAngle;Population_StaggeringAngle;TOTAngle;Population_TOTAngle;TTDistance;Population_TTDistance;TTTAngle;Population_TTTAngle;OOO_Angle;Population_OOOAngle" > data.csv
for CIFFile in *.cif ; do
 Name=$(echo $CIFFile | sed 's/\.cif//g' | sed 's/_/ /g' | awk '{print $1}')
 SuperCell=$(echo $CIFFile | sed 's/\.cif//g' | sed 's/_/ /g' | awk '{print $2}')
 #n_T=$(cat Download_Database/Original_IZA_Structures/${Name}.cif | sed 's/\t/    /g' | grep " Si " | wc -l | awk '{print $1}')
 # Ring Sizes (from IZA database)
 #if [ ! -f Download_Database/Original_IZA_Structures/${Name}.html ] ; then
  ## curl https://europe.iza-structure.org/IZA-SC/framework.php?STC=${Name} --output Download_Database/Original_IZA_Structures/${Name}.html
  #echo "Dummy data" > 
 #fi
 #Rings=$(grep -A3 "Ring sizes" Download_Database/Original_IZA_Structures/${Name}.html | grep "Normal" | tail -n1 | sed -e 's/>/ /g' -e 's/</ /g' | awk '{print $3}' | sed 's/&nbsp;&nbsp;/ /g')
 #TD=$(grep -A3 "Topological density:" Download_Database/Original_IZA_Structures/${Name}.html | tail -n1 | awk '{print $5}' | sed -e 's/<\/td>//g')
 n_geom=0
 for geom in "Staggering_1st" "siosi" "sisi" "sisisi" "ooo" ; do
  let n_geom++
  file=${Name}_${SuperCell}_${geom}_histogram.dat
  sed '/#/d' $file | awk '{if ($2>0)  print $0}' > tmp
  n_lines=0
  while read line ; do
   let n_lines++
   Staggering_angle[${n_lines}]=$(echo $line | awk '{print $1}')
   Staggering_angle_population[${n_lines}]=$(echo $line | awk '{print $2}')
  done < tmp
  meassure[${n_geom}]=$(awk '{S=S?S OFS s1 $1 s1:s1 $1 s1} END{print S}' OFS=" " tmp)
  population[${n_geom}]=$(awk '{S=S?S OFS s1 $2 s1:s1 $2 s1} END{print S}' OFS=" " tmp)
 done
 echo "${Name};${n_T};${Rings};${TD};[${meassure[1]}];[${population[1]}];[${meassure[2]}];[${population[2]}];[${meassure[3]}];[${population[3]}];[${meassure[4]}];[${population[4]}];[${meassure[5]}];[${population[5]}]" >> data.csv
 rm -rf tmp
done
exit 0
