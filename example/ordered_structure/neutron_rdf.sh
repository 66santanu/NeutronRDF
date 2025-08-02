#!/bin/bash
## Scattering legths
b_F=5.654 
b_Be=7.79
## populations of elements
cF=0.6667
cBe=0.3333
## for getting the last frame
let nline=2306
let natoms=2304
echo "$natoms" > natoms
## box lenths for extended XYZ format
echo "33.14484 33.14484 33.14484" > lattice

##########################################################################
for((i=0; i<1; i++))
do
 tail -$natoms sample${i}.xyz > coordinates
 cat natoms lattice coordinates > Structure.xyz
 #####################################
 ../../bin/rdf <input
 awk '{printf "%6.3f %6.3f\n", $1,$2}' gofr-AA.txt >gofrBeBe${i}.dat
 awk '{print $2}' gofr-AA.txt > BeBe

 #####################################
 awk '{printf "%6.3f %6.3f\n", $1,$2}' gofr-AB.txt >gofrBeF${i}.dat
 awk '{print $2}' gofr-AB.txt > BeF

 #####################################
 awk '{printf "%6.3f %6.3f\n", $1,$2}' gofr-BB.txt >gofrFF${i}.dat
 awk '{print $2}' gofr-BB.txt > FF

 #####################################
 awk '{printf "%6.3f\n",$1}' gofr-AB.txt > dist
 cat BeBe | awk '{print $1*('$cBe'^2)*('$b_Be'^2)}'>gBeBe
 cat FF | awk '{print $1*('$cF'^2)*('$b_F'^2)}'>gFF
 cat BeF | awk '{print $1*2*('$cBe')*('$cF')*('$b_Be')*('$b_F')}'>gBeF

 NormConst=$(awk "BEGIN {print (($cBe*$b_Be)+($cF*$b_F))^2}") 
 echo "Normalization constant is : $NormConst"
 
 paste gBeBe gFF gBeF | awk '{print ($1+$2+$3)/('$NormConst')}' > weighted_sum
 paste dist  weighted_sum  | awk '{print $1,$2}' > gofr_total${i}.dat
 rm dist weighted_sum gBeBe gFF gBeF BeBe FF BeF gofr-AA.txt gofr-AB.txt gofr-BB.txt lattice  natoms 
 rm Structure.xyz coordinates  
done
