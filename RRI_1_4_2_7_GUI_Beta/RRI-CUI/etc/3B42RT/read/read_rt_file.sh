#!/bin/sh

# variables
jleft=269
ibottom=106
jright=281
itop=95

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

    ifname=./infile/3B42RT.$year$month$day$hour.7.bin
    ofname=./outfile/3B42RT.$year$month$day$hour.txt

    if [ -f $ifname ]
    then

     if [ -f $ofname ]
     then
      echo "skip $ifname"
     else

      echo "process $ifname"

      echo "$ifname" >  infile.txt
      echo "3B42RT"  >> infile.txt
      echo "$jleft"   >> infile.txt
      echo "$ibottom" >> infile.txt
      echo "$jright"  >> infile.txt
      echo "$itop"    >> infile.txt

      #col=`expr $jright - $jleft + 1`
      #row=`expr $ibottom - $itop + 1`

      ./read_rt_file.exe < infile.txt > $ofname

     fi
    fi
   done
  done
 done
done
