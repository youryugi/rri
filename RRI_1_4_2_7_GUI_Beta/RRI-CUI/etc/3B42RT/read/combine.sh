#!/bin/sh

# variables
jleft=269
ibottom=106
jright=281
itop=95

outfile="./rain_ind.txt"
rm -f $outfile

t=0

for year in 2015
do

 for month in 12
 do

  #for day in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31
  for day in 27 28 29 30
  do

   for hour in 00 03 06 09 12 15 18 21
   do

    ofname=./outfile/3B42RT.$year$month$day$hour.txt

    if [ -f $ofname ]
    then

      echo "combine $ofname"

      col=`expr $jright - $jleft + 1`
      row=`expr $ibottom - $itop + 1`

      echo "$t $year $month $day $hour"
      echo "$t $col $row" >> $outfile
      cat < $ofname >> $outfile
      t=`expr $t + 10800`

    fi
   done
  done
 done
done
